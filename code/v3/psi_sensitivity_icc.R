# ==============================================================================
# ICC-based psi estimation: bootstrap uncertainty and parameter sensitivity
# ==============================================================================
#
# A transparent complement to the full Bayesian likelihood in
# psi_empirical_litgi.R. The core idea: within a cluster, all siblings share
# the same infector timing, so their serial intervals are correlated. The
# intraclass correlation (ICC) of serial intervals grouped by infector measures
# this clustering. After correcting for incubation period variance (which adds
# within-cluster noise regardless of psi), ICC maps directly onto psi.
#
# Three outputs:
#   1. Bootstrap CI on psi for each pathogen (resampling clusters)
#   2. Comparison figure: ICC-based vs. likelihood-based psi estimates
#   3. Sensitivity surface: how psi changes as GI and incubation SDs vary
#
# Depends on: psi_empirical_litgi.R (saves output/psi_cluster_data.RDS)
# ==============================================================================

cat("=== ICC-based psi sensitivity analysis ===\n")

library(tidyverse)
source("code/utils.R")

# Load cluster data saved by psi_empirical_litgi.R
cluster_file <- file.path("output", "psi_cluster_data.RDS")
if (!file.exists(cluster_file))
	stop("Run psi_empirical_litgi.R first to generate output/psi_cluster_data.RDS")

dat                  <- readRDS(cluster_file)
all_disease_clusters <- dat$all_disease_clusters
pathogen_params      <- dat$pathogen_params

# ==============================================================================
# Core functions
# ==============================================================================

#' Intraclass correlation of serial intervals, grouped by infector.
#'
#' Uses a one-way ANOVA. Only clusters with >= 2 siblings are informative;
#' singletons are excluded. Returns NA if fewer than 2 multi-offspring clusters.
#'
#' @param clusters List of numeric vectors (serial intervals per infector).
#' @return Scalar ICC in (-1, 1).
compute_icc <- function(clusters) {
	multi <- Filter(function(cl) length(cl) >= 2, clusters)
	if (length(multi) < 2) return(NA_real_)

	df <- map_dfr(seq_along(multi), ~tibble(cluster = .x, si = multi[[.x]]))

	ns <- lengths(multi)
	n0 <- length(ns) / sum(1 / ns)          # harmonic mean cluster size
	ms <- summary(aov(si ~ factor(cluster), data = df))[[1]][, "Mean Sq"]

	(ms[1] - ms[2]) / (ms[1] + (n0 - 1) * ms[2])
}

#' Convert ICC to psi given GI variance and incubation period variance.
#'
#' Derivation: within-cluster SI variance = psi * var_gi + var_inc;
#' total SI variance = var_gi + 2 * var_inc; ICC = (1-psi)*var_gi + var_inc)
#' / (var_gi + 2*var_inc). Solving for psi gives this formula.
#' Result is clamped to [0, 1] to handle finite-sample noise.
#'
#' @param icc    Scalar or vector of ICC values.
#' @param var_gi Generation interval variance (days^2).
#' @param var_inc Incubation period variance (days^2).
#' @return Scalar or vector of psi values in [0, 1].
icc_to_psi <- function(icc, var_gi, var_inc) {
	psi <- 1 - (icc * (var_gi + 2 * var_inc) - var_inc) / var_gi
	pmax(0, pmin(1, psi))
}

# ==============================================================================
# Bootstrap psi estimates
# ==============================================================================

cat("\n--- Bootstrap psi estimates ---\n")

set.seed(42)
n_boot            <- 2000
min_multi_clust   <- 5

boot_results <- list()

for (disease in names(all_disease_clusters)) {
	pp       <- pathogen_params[[disease]]
	clusters <- all_disease_clusters[[disease]]$clusters
	multi    <- Filter(function(cl) length(cl) >= 2, clusters)

	if (length(multi) < min_multi_clust) {
		cat(sprintf("  Skipping %-20s (only %d multi-offspring clusters)\n",
		    disease, length(multi)))
		next
	}

	var_gi  <- pp$alpha_gi / pp$beta_gi^2
	var_inc <- pp$a_obs    / pp$b_obs^2

	icc_obs <- compute_icc(clusters)
	psi_obs <- icc_to_psi(icc_obs, var_gi, var_inc)

	# Resample clusters with replacement; compute psi for each bootstrap replicate
	boot_psis <- replicate(n_boot, {
		boot_cls <- sample(clusters, length(clusters), replace = TRUE)
		icc_to_psi(compute_icc(boot_cls), var_gi, var_inc)
	})

	ci_lo <- quantile(boot_psis, 0.025, na.rm = TRUE)
	ci_hi <- quantile(boot_psis, 0.975, na.rm = TRUE)

	cat(sprintf("  %-20s  psi = %.2f  95%% CI [%.2f, %.2f]  (ICC = %.2f, n_clusters = %d)\n",
	    disease, psi_obs, ci_lo, ci_hi, icc_obs, length(multi)))

	boot_results[[disease]] <- tibble(
		disease  = disease,
		psi_obs  = psi_obs,
		icc_obs  = icc_obs,
		ci_lo    = ci_lo,
		ci_hi    = ci_hi,
		var_gi   = var_gi,
		var_inc  = var_inc,
		n_multi  = length(multi)
	)
}

boot_df <- bind_rows(boot_results)

# ==============================================================================
# Sensitivity surface: vary GI SD and incubation SD
# ==============================================================================

cat("\n--- Computing sensitivity surfaces ---\n")

# Each axis spans 70%–130% of the published SD value (variance varies as square)
frac_grid    <- seq(0.7, 1.3, length.out = 30)
sens_results <- list()

for (disease in names(boot_results)) {
	pp      <- pathogen_params[[disease]]
	icc_obs <- boot_results[[disease]]$icc_obs

	gi_sd_pub  <- sqrt(pp$alpha_gi / pp$beta_gi^2)
	inc_sd_pub <- sqrt(pp$a_obs    / pp$b_obs^2)

	sens_results[[disease]] <- expand_grid(
		gi_frac  = frac_grid,
		inc_frac = frac_grid
	) %>%
		mutate(
			psi     = icc_to_psi(icc_obs,
			                     var_gi  = (gi_sd_pub  * gi_frac)^2,
			                     var_inc = (inc_sd_pub * inc_frac)^2),
			disease = disease
		)
}

sens_df <- bind_rows(sens_results)

# ==============================================================================
# Figures
# ==============================================================================

cat("\n--- Generating figures ---\n")

# Consistent disease ordering (by ICC-based psi, low to high)
disease_order <- boot_df %>% arrange(psi_obs) %>% pull(disease)
boot_df$disease <- factor(boot_df$disease, levels = disease_order)

# --------------------------------------------------------------------------
# Figure 1: ICC-based vs. likelihood-based psi estimates
# --------------------------------------------------------------------------

# Load likelihood posteriors for comparison
lik_file <- file.path("output", "psi_empirical_litgi_results.csv")
if (file.exists(lik_file)) {
	lik_df <- read_csv(lik_file, show_col_types = FALSE) %>%
		filter(disease %in% boot_df$disease) %>%
		group_by(disease) %>%
		summarise(
			psi_est = first(post_mean),
			ci_lo   = first(ci_lo),
			ci_hi   = first(ci_hi),
			.groups = "drop"
		) %>%
		mutate(method = "Likelihood (posterior mean)",
		       disease = factor(disease, levels = disease_order))

	icc_df <- boot_df %>%
		transmute(disease, psi_est = psi_obs, ci_lo, ci_hi,
		          method = "ICC-based (bootstrap CI)")

	compare_df <- bind_rows(icc_df, lik_df)

	fig_compare <- ggplot(compare_df,
	    aes(x = psi_est, y = disease, color = method, shape = method)) +
		geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
		    height = 0.25, linewidth = 0.6,
		    position = position_dodge(width = 0.55)) +
		geom_point(size = 3,
		    position = position_dodge(width = 0.55)) +
		scale_color_manual(
		    values = c("ICC-based (bootstrap CI)"    = "darkorange",
		               "Likelihood (posterior mean)" = "steelblue"),
		    name = NULL) +
		scale_shape_manual(
		    values = c("ICC-based (bootstrap CI)"    = 17,
		               "Likelihood (posterior mean)" = 16),
		    name = NULL) +
		xlim(0, 1) +
		labs(x = expression(hat(psi)),
		     y = NULL,
		     title = expression("ICC-based vs. likelihood-based estimates of" ~ psi),
		     subtitle = "Points = posterior mean (likelihood) or ICC point estimate; bars = 95% intervals") +
		theme_minimal(base_size = 13) +
		theme(panel.grid.major.y = element_blank(),
		      legend.position    = "bottom")

	save_fig(fig_compare, "psi_sensitivity_comparison",
	         width = 9, height = max(4, nrow(boot_df) * 0.65))
	cat("  Saved comparison figure\n")
} else {
	cat("  Skipping comparison figure (likelihood results not found)\n")
}

# --------------------------------------------------------------------------
# Figure 2: Sensitivity surfaces (one facet per disease)
# --------------------------------------------------------------------------

sens_df$disease <- factor(sens_df$disease, levels = disease_order)
n_dis <- length(unique(sens_df$disease))

fig_sens <- ggplot(sens_df, aes(x = gi_frac, y = inc_frac, fill = psi)) +
	geom_tile() +
	geom_contour(aes(z = psi), color = "white", alpha = 0.6,
	             linewidth = 0.3, breaks = seq(0, 1, by = 0.2)) +
	# Mark the published parameter values
	geom_point(data = tibble(gi_frac = 1, inc_frac = 1),
	           inherit.aes = FALSE, shape = 3, size = 4,
	           color = "white", stroke = 1.5) +
	scale_fill_viridis_c(name = expression(hat(psi)),
	                     limits = c(0, 1), option = "plasma") +
	scale_x_continuous(labels = scales::percent_format(accuracy = 1),
	                   breaks = c(0.7, 1.0, 1.3)) +
	scale_y_continuous(labels = scales::percent_format(accuracy = 1),
	                   breaks = c(0.7, 1.0, 1.3)) +
	facet_wrap(~disease, ncol = 3) +
	labs(x = "GI SD (× published value)",
	     y = "Incubation SD (× published value)",
	     title = expression("Sensitivity of " ~ hat(psi) ~ " to GI and incubation period assumptions"),
	     subtitle = "White + marks published parameter values. Contours at psi = 0.2, 0.4, 0.6, 0.8.") +
	theme_minimal(base_size = 11) +
	theme(strip.text = element_text(face = "bold"),
	      panel.spacing = unit(0.8, "lines"))

save_fig(fig_sens, "psi_sensitivity_surface",
         width = 12, height = ceiling(n_dis / 3) * 3.8)
cat("  Saved sensitivity surface figure\n")

# --------------------------------------------------------------------------
# Save numerical summary
# --------------------------------------------------------------------------

write_csv(boot_df, file.path("output", "psi_sensitivity_icc_results.csv"))
cat("  Saved results to output/psi_sensitivity_icc_results.csv\n")

cat("\n=== ICC sensitivity analysis complete ===\n")
