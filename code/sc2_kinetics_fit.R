#!/usr/bin/env Rscript
# =============================================================================
# sc2_kinetics_fit.R
#
# Fits the piecewise-linear ("tent") viral-kinetics model to EVERY infection in
# ct_dat_refined.csv under a SINGLE analysis category (all infections shrunk
# toward one common population mean — no first/second/variant/paired structure),
# with STANDARD observation noise, then runs the time-below-Ct-threshold
# ("infectious window") analysis. The Stan model lives in a separate file
# (code/sc2_kinetics_fit.stan); the data-prep helpers (from utils.R) and priors
# (set_run_pars.R index 1) are inlined below.
#
# Usage:
#   Rscript refit_single_category.R
# Config via env vars (all optional):
#   CT_DATA      path to ct_dat_refined.csv     (default: data/ct_dat_refined.csv)
#   CT_OUT       output dir                      (default: output/single_category)
#   CT_STAN      path to the Stan model file     (default: code/sc2_kinetics_fit.stan)
#   CT_ITER      Stan iterations per chain       (default: 2000)
#   CT_CHAINS    Stan chains                      (default: 4)
#   CT_FIRST_ONLY  "TRUE" to fit first infections only (default: FALSE = all)
#   CT_MAXINF    cap #infections for a quick smoke test (default: 0 = all)
#
# NOTE: a full fit of ~2,000 infections takes tens of minutes to hours and
# several GB RAM. Use CT_MAXINF + small CT_ITER to smoke-test first.
# =============================================================================

suppressMessages({library(tidyverse); library(rstan)})
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("code/sc2_kinetics_utils.R")   # window computations + figure functions

## ---- config -----------------------------------------------------------------
DATA       <- Sys.getenv("CT_DATA", "data/ct_dat_refined.csv")
OUTDIR     <- Sys.getenv("CT_OUT",  "output/")
STAN_FILE  <- Sys.getenv("CT_STAN", "code/sc2_kinetics_fit.stan")
ITER       <- as.integer(Sys.getenv("CT_ITER",   "2000"))
CHAINS     <- as.integer(Sys.getenv("CT_CHAINS", "4"))
FIRST_ONLY <- toupper(Sys.getenv("CT_FIRST_ONLY","FALSE")) %in% c("TRUE","1","YES")
MAXINF     <- as.integer(Sys.getenv("CT_MAXINF", "0"))
CT_THRESHOLDS <- c(35, 32, 26)   # 32 = "bare" infectiousness, 26 = ~50% culture positivity
GI_SD      <- 2.0            # SARS-CoV-2 intrinsic generation-interval SD (days), reference
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

convert_Ct_logGEML <- function(Ct, m_conv=-3.609714286, b_conv=40.93733333){
	out <- (Ct-b_conv)/m_conv * log10(10) + log10(250)
	return(out) 
}

global_pars <- c(lod = 40)  # limit of detection (Ct)

## priors — standard-noise config (matches set_run_pars.R index 1) -------------
prior_pars_base <- list(
	tp_prior    = c(0, 2),
	dp_midpoint = 20,
	wp_midpoint = 5,
	wr_midpoint = 12,
	sigma_prior = c(0, 0.5),     # standard observation noise (NOT the "lowsig" variant)
	lambda      = 0.01,
	fpmean      = 1/log(10),
	priorsd     = 0.25
)

## ---- data-prep helpers (verbatim logic from utils.R) ------------------------
clean_infection_events <- function(ct_dat){
	ie_list <- sort(unique(ct_dat$InfectionEvent))
	ie_df <- data.frame(InfectionEvent = ie_list, InfectionEventClean = seq_along(ie_list))
	left_join(ct_dat, ie_df, by = "InfectionEvent")
}
trim_negatives <- function(indiv_data, global_pars){
	indiv_data %>% split(.$id) %>%
		map(~ arrange(., t)) %>%
		map(~ mutate(., ispositive = if_else(y < global_pars[["lod"]], 1, 0),
										 il  = lag(ispositive), il2 = lag(ispositive, 2),
										 ld  = lead(ispositive), ld2 = lead(ispositive, 2))) %>%
		map(~ filter(., ispositive == 1 | il == 1 | il2 == 1 | ld == 1 | ld2 == 1)) %>%
		map(~ select(., -ispositive, -il, -il2, -ld, -ld2)) %>%
		bind_rows()
}

## ---- build single-category analysis data ------------------------------------
ct <- read_csv(DATA, show_col_types = FALSE)
if (nrow(ct) == 0) stop("ct_dat_refined.csv is empty at: ", DATA,
												" (the in-repo copy is a gitignored stub; point CT_DATA at the real file).")
if (FIRST_ONLY) ct <- filter(ct, InfNum == 1)
if (MAXINF > 0) {
	keep <- ct %>% distinct(InfectionEvent) %>% slice_head(n = MAXINF) %>% pull(InfectionEvent)
	ct <- filter(ct, InfectionEvent %in% keep)
}

indiv_data <- ct %>%
	mutate(id = InfectionEvent, t = TestDateIndex, y = CtT1, category = 1L) %>%
	trim_negatives(global_pars) %>%
	clean_infection_events() %>%
	mutate(id_clean = InfectionEventClean) %>%
	ungroup() %>%
	select(id, id_clean, t, y, category) %>%
	arrange(id_clean, t)

n_indiv <- n_distinct(indiv_data$id_clean)
catlist <- rep(1L, n_indiv)          # single analysis category
adjmat  <- matrix(1L, n_indiv, 1)    # single (null) adjustment category
cat(sprintf("Fitting %d infections | %d measurements | single category | %s noise | iter=%d chains=%d\n",
						n_indiv, nrow(indiv_data), "standard", ITER, CHAINS))


## ---- fit --------------------------------------------------------------------
model <- stan_model(file = STAN_FILE)
t0 <- Sys.time()
fit <- sampling(model, data = list(
	N = nrow(indiv_data), n_id = n_indiv, lod = global_pars[["lod"]],
	id = indiv_data$id_clean, category = catlist, max_category = max(catlist),
	n_adj = ncol(adjmat), adjust = adjmat, max_adjust = max(adjmat),
	t = indiv_data$t, y = indiv_data$y,
	tp_prior = prior_pars_base$tp_prior, dp_midpoint = prior_pars_base$dp_midpoint,
	wp_midpoint = prior_pars_base$wp_midpoint, wr_midpoint = prior_pars_base$wr_midpoint,
	sigma_prior = prior_pars_base$sigma_prior, lambda = prior_pars_base$lambda,
	fpmean = prior_pars_base$fpmean, priorsd = prior_pars_base$priorsd),
	iter = ITER, chains = CHAINS)
cat(sprintf("Fit time: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
# quick convergence check on the REAL parameters (single-category empties give NA R-hat -> ignore)
print(summary(fit, pars=c("log_dp_mean","log_dp_sd","log_wp_mean","log_wp_sd",
													"log_wr_mean","log_wr_sd","sigma"))$summary[,c("mean","n_eff","Rhat")])

## ---- reconstruct per-infection posteriors (single category => no adjustment)-
## Extract only the small params we need, save the SMALL features first, then
## release the (large) fit. The full fit is NOT saved by default (set CT_SAVEFIT=TRUE).
fl <- rstan::extract(fit, pars = c("log_dp_mean","log_dp_sd","dp_raw",
																	 "log_wp_mean","log_wp_sd","wp_raw",
																	 "log_wr_mean","log_wr_sd","wr_raw","tp"))
# dp[i] = exp(log_dp_mean + log_dp_sd*dp_raw[i]) * midpoint   (adj/category terms are 0)
recon <- function(mean_v, sd_v, raw_m, mid){                                    # [ndraws x n_id]
	mean_v <- as.numeric(mean_v); sd_v <- as.numeric(sd_v); raw_m <- as.matrix(raw_m)
	exp(matrix(mean_v, nrow(raw_m), ncol(raw_m)) +
			matrix(sd_v,  nrow(raw_m), ncol(raw_m)) * raw_m) * mid
}
dp <- recon(fl$log_dp_mean, fl$log_dp_sd, fl$dp_raw, prior_pars_base$dp_midpoint)
wp <- recon(fl$log_wp_mean, fl$log_wp_sd, fl$wp_raw, prior_pars_base$wp_midpoint)
wr <- recon(fl$log_wr_mean, fl$log_wr_sd, fl$wr_raw, prior_pars_base$wr_midpoint)

idmap <- indiv_data %>% distinct(id_clean, id) %>% arrange(id_clean)   # id = InfectionEvent
tp <- fl$tp                                                            # [ndraws x n_id] peak timing
colnames(dp) <- colnames(wp) <- colnames(wr) <- colnames(tp) <- idmap$id
key <- ct %>% group_by(InfectionEvent) %>% slice(1) %>% ungroup() %>%
	select(InfectionEvent, PersonID, InfNum, AgeGrp, VaccinationStatus, BoosterStatus, LineageBroad)
# tp + raw points + lod are stashed so the window analysis and every figure can
# be regenerated from the rds alone (see sc2_kinetics_figures.R) without re-fitting.
feats <- list(dp = dp, wp = wp, wr = wr, tp = tp, key = key,
							raw = indiv_data %>% transmute(InfectionEvent = id, id_clean, t, y),
							lod = global_pars[["lod"]], run = "single_category")
saveRDS(feats, file.path(OUTDIR, "posterior_features_single_category.rds"))
if (toupper(Sys.getenv("CT_SAVEFIT","FALSE")) %in% c("TRUE","1","YES")) {
	ok <- tryCatch({ save(fit, indiv_data, prior_pars_base, file=file.path(OUTDIR,"fit.RData")); TRUE },
								 error=function(e){ message("WARNING: could not save full fit (", conditionMessage(e),
																						 ") -- features already saved, continuing."); FALSE })
}
rm(fit); gc(verbose=FALSE)   # free memory before the window analysis

## ---- time below Ct threshold (infectious window) ----------------------------
## All window computations + figures are factored into sc2_kinetics_utils.R and
## operate on `feats`, so they can be reproduced from the saved rds without a
## re-fit (see sc2_kinetics_figures.R).
report_ct_windows(feats, thresholds = CT_THRESHOLDS, gi_sd = GI_SD)

## window-width summary for every integer Ct 40..20 (thresholds can be changed
## later without re-fitting)
window_summary <- window_summary_by_ct(feats, cts = 40:20)
write_csv(window_summary, file.path(OUTDIR, "window_summary_by_Ct.csv"))
cat(sprintf("Saved: %s/window_summary_by_Ct.csv (window width + #trajectories exceeding, Ct 40..20)\n", OUTDIR))

## per-infection feature table (window at each threshold + peak params + covariates)
write_csv(infection_window_features(feats, thresholds = c(32, 26)),
					file.path(OUTDIR, "infection_windows_single_category.csv"))

## ---- fitted tent trajectories vs raw data (25 infections) -------------------
ptraj <- plot_trajectories(feats, n = 25, seed = 1) +
	labs(title = "Fitted tent trajectories vs raw Ct (25 randomly-drawn infections, >=4 obs)",
			 subtitle = "line = posterior-mean tent; points = observed Ct; dashed = Ct 26 (~50% culturable)")
ggsave(file.path(OUTDIR, "trajectories_ct26.png"), ptraj, width = 10, height = 8, dpi = 150)
cat(sprintf("Saved: %s/trajectories_ct26.png (25 trajectories + raw data, Ct=26 line)\n", OUTDIR))

cat(sprintf("\nSaved: %s/{posterior_features_single_category.rds, infection_windows_single_category.csv}%s\n",
		OUTDIR, if(toupper(Sys.getenv("CT_SAVEFIT","FALSE")) %in% c("TRUE","1","YES")) " (+ fit.RData if space allowed)" else ""))
