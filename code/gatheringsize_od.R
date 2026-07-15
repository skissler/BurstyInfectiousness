# ==============================================================================
# Gathering size restrictions: overdispersion (OD) of secondary infection counts
# ==============================================================================
#
# Directly compares OD = 1/k = Var[nu_i] / E[nu_i]^2 between the gathering-size-
# restricted (c_max = c_max_main) and unrestricted (c_max = Inf) Gamma-Poisson
# contact processes, across psi. 
#
# Per cell: generate n_cases index cases, count secondary-infection offspring,
# compute the method-of-moments OD estimator
#
#   OD_hat = (Var[N] - E[N]) / E[N]^2,
#
# (or equivalently 1 / k_hat where k_hat = E[N]^2 / (Var[N] - E[N])). Bootstrap
# standard errors are computed by resampling the per-cell offspring counts.
#
# Analytic OD curves are computed using the integral
#
#   I(psi*alpha, lambda/beta) = integral_{-Inf}^{Inf}
#                                (lambda / (pi*(lambda^2 + omega^2)))
#                                * (1 + omega^2 / beta^2)^(-psi*alpha) domega,
#
# giving OD = V_c * I / mu_c^2 with (V_c, mu_c) = (1/sigma, 1) unrestricted and
# (V_T, mu_T) restricted (V_T, mu_T derived from the truncated Gamma(sigma, sigma)).
# ==============================================================================

library(tidyverse)
source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

# --- Parameters ---
sigma           <- 1            # contact heterogeneity (Gamma shape = rate)
lambda_gp       <- 1            # Poisson switching rate
c_max_main      <- 2            # gathering cap
psi_vals_sim    <- c(0, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 1)
psi_vals_theory <- seq(0, 1, length.out = 101)
n_cases         <- 10000        # index cases per (psi, scenario) cell
B_boot          <- 200          # bootstrap reps for OD SE

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Analytic helpers
# ==============================================================================

# Spectral integral I(psi*alpha, lambda/beta).
spectral_integral <- function(psi_alpha, lambda, beta) {
	# psi_alpha == 0 -> kernel collapses to (1)^0 = 1, integral of Cauchy = 1.
	if (psi_alpha < 1e-12) return(1)
	integrand <- function(omega) {
		(lambda / (pi * (lambda^2 + omega^2))) *
		(1 + omega^2 / beta^2)^(-psi_alpha)
	}
	integrate(integrand, lower = -Inf, upper = Inf,
	          rel.tol = 1e-8, subdivisions = 1000L)$value
}

# Analytic OD = Var[c_tilde] / E[c_tilde]^2 = V_c * I(psi*alpha, lambda/beta) / mu_c^2.
od_analytic <- function(psi, alpha, beta, sigma, lambda, c_max) {
	I_psi <- spectral_integral(psi * alpha, lambda, beta)
	if (is.infinite(c_max)) {
		mu_c <- 1
		V_c  <- 1 / sigma
	} else {
		mu_c <- mu_truncated_gamma(sigma, c_max)
		V_c  <- var_truncated_gamma(sigma, c_max)
	}
	V_c * I_psi / mu_c^2
}

# extinction_prob_negbin() (Lloyd-Smith et al., Nature 2005) now lives in utils.R.
# Expected establishment probability given (R0, k) in this NegBin framework.
p_est_theory <- function(R0, k) 1 - extinction_prob_negbin(R0, k)

# ==============================================================================
# Pathogen loop
# ==============================================================================

scenario_names <- c("Restricted", "Unrestricted", "Constant")

for (idx_pathogen in seq_along(parslist)) {

pars <- parslist[[idx_pathogen]]
pathogen <- pars$pathogen
pathogen_label <- switch(pathogen,
    influenza = "Influenza",
    omicron   = "SARS-CoV-2 Omicron",
    measles   = "Measles",
    pathogen
)
Tgen     <- pars$Tgen
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

mu_T_main   <- mu_truncated_gamma(sigma, c_max_main)
R0_new_main <- R0 * mu_T_main

set.seed(42 + idx_pathogen)

cat(sprintf("\n===== Gathering size OD: %s (alpha=%.2f, beta=%.3f, R0=%.2f -> R0_eff=%.2f) =====\n",
            pathogen, alpha, beta, R0, R0_new_main))

# --------------------------------------------------------------------------
# Simulated OD per (psi, scenario) cell
# --------------------------------------------------------------------------

cells <- expand_grid(
		psi_idx      = seq_along(psi_vals_sim),
		scenario_idx = seq_along(scenario_names)) %>%
	mutate(
		psi      = psi_vals_sim[psi_idx],
		scenario = scenario_names[scenario_idx],
		seed     = 42L + 1000L * idx_pathogen + 100L * psi_idx + scenario_idx)

sim_cache <- file.path(
	"output",
	sprintf("gathering_od_sim_%s_n%d_kc%g_lam%g_cmax%g.csv",
	        pathogen, n_cases, sigma, lambda_gp, c_max_main)
)

if (file.exists(sim_cache)) {
	od_sim_df <- read_csv(sim_cache, show_col_types = FALSE)
	cat(sprintf("  %s: loaded simulated OD from %s\n",
	    pathogen, basename(sim_cache)))
} else {

	run_cell <- function(psi_val, scenario_name, seed) {
		set.seed(seed)
		# R0 is matched at R0_new_main = R0 * mu_T across all scenarios, so
		# expected offspring counts are equal (effective R0 is the same).
		# Restricted: passes raw R0, since the truncated contact process itself
		# contributes the mu_T factor (E[c_i^*(t)] = mu_T), giving expected
		# offspring R0 * mu_T = R0_new. The other two scenarios have
		# E[c_i(t)] = 1, so we pass R0_new directly.
		gfun <- switch(scenario_name,
			"Restricted"   = gen_inf_attempts_gammapoisson_contacts(
			                    Tgen, R0, alpha, psi_val, sigma, lambda_gp,
			                    c_max = c_max_main),
			"Unrestricted" = gen_inf_attempts_gammapoisson_contacts(
			                    Tgen, R0_new_main, alpha, psi_val,
			                    sigma, lambda_gp, c_max = Inf),
			"Constant"     = gen_inf_attempts_gamma(
			                    Tgen, R0_new_main, alpha, psi_val))
		noffspring <- replicate(n_cases, length(gfun(0)))

		m <- mean(noffspring)
		v <- var(noffspring)
		od_hat <- if (v > m) (v - m) / m^2 else 0

		# Bootstrap SE for OD
		od_boot <- replicate(B_boot, {
			idx <- sample.int(n_cases, n_cases, replace = TRUE)
			mb  <- mean(noffspring[idx])
			vb  <- var(noffspring[idx])
			if (vb > mb) (vb - mb) / mb^2 else 0
		})

		tibble(
			psi      = psi_val,
			scenario = scenario_name,
			od_sim   = od_hat,
			od_se    = sd(od_boot, na.rm = TRUE),
			mean_n   = m,
			var_n    = v,
			n_cases  = n_cases)
	}

	n_cores <- max(1L, parallel::detectCores() - 1L)
	cat(sprintf("  Running %d cells across %d core(s)...\n", nrow(cells), n_cores))

	od_sim_df <- parallel::mclapply(seq_len(nrow(cells)), function(i) {
			run_cell(cells$psi[i], cells$scenario[i], cells$seed[i])
		}, mc.cores = n_cores) %>%
		bind_rows()

	write_csv(od_sim_df, sim_cache)
}

cat(sprintf("\n  %s: simulated OD by (psi, scenario)\n", pathogen))
od_sim_df %>%
	select(psi, scenario, od_sim, od_se, mean_n) %>%
	print(n = Inf)

# --------------------------------------------------------------------------
# Analytic OD curves
# --------------------------------------------------------------------------

od_theory_df <- expand_grid(
		psi      = psi_vals_theory,
		scenario = scenario_names) %>%
	rowwise() %>%
	mutate(
		od_theory = case_when(
			scenario == "Restricted"   ~ od_analytic(psi, alpha, beta, sigma, lambda_gp, c_max = c_max_main),
			scenario == "Unrestricted" ~ od_analytic(psi, alpha, beta, sigma, lambda_gp, c_max = Inf),
			scenario == "Constant"     ~ 0)) %>%
	ungroup()

# --------------------------------------------------------------------------
# Figure 1: OD vs psi (three scenarios, theory + sim)
# --------------------------------------------------------------------------

pal_scenario <- c("Restricted"   = "#FF9800",
                  "Unrestricted" = "#4CAF50",
                  "Constant"     = "#2196F3")

fig_od <- ggplot() +
	geom_line(data = od_theory_df,
	          aes(x = psi, y = od_theory, color = scenario),
	          linewidth = 0.8) +
	geom_errorbar(data = od_sim_df,
	              aes(x = psi,
	                  ymin = pmax(od_sim - 1.96 * od_se, 0),
	                  ymax = od_sim + 1.96 * od_se,
	                  color = scenario),
	              width = 0.015) +
	geom_point(data = od_sim_df,
	           aes(x = psi, y = od_sim, color = scenario),
	           size = 2) +
	scale_color_manual(values = pal_scenario, name = NULL) +
	theme_classic(base_size = 16) +
	labs(x = expression(psi),
	     y = "OD = 1/k",
	     title = pathogen_label)

save_fig(fig_od, sprintf("fig_gathering_od_%s", pathogen), width = 6, height = 5)
cat(sprintf("  Saved fig_gathering_od_%s\n", pathogen))

# --------------------------------------------------------------------------
# Figure 2: absolute OD reduction (OD_Unrestricted - OD_Restricted) vs psi
# --------------------------------------------------------------------------

od_diff_theory <- od_theory_df %>%
	filter(scenario %in% c("Restricted", "Unrestricted")) %>%
	pivot_wider(names_from = scenario, values_from = od_theory) %>%
	mutate(dOD = Unrestricted - Restricted)

od_diff_sim <- od_sim_df %>%
	filter(scenario %in% c("Restricted", "Unrestricted")) %>%
	select(psi, scenario, od_sim, od_se) %>%
	pivot_wider(names_from = scenario, values_from = c(od_sim, od_se)) %>%
	mutate(
		dOD    = od_sim_Unrestricted - od_sim_Restricted,
		dOD_se = sqrt(od_se_Unrestricted^2 + od_se_Restricted^2))

fig_dod <- ggplot() +
	geom_line(data = od_diff_theory,
	          aes(x = psi, y = dOD),
	          color = "black", linewidth = 0.8) +
	geom_errorbar(data = od_diff_sim,
	              aes(x = psi,
	                  ymin = dOD - 1.96 * dOD_se,
	                  ymax = dOD + 1.96 * dOD_se),
	              width = 0.015) +
	geom_point(data = od_diff_sim,
	           aes(x = psi, y = dOD),
	           size = 2) +
	theme_classic(base_size = 16) +
	labs(x = expression(psi),
	     y = expression(Delta*"OD = OD"[Unrestricted] - "OD"[Restricted]),
	     title = pathogen_label)

save_fig(fig_dod, sprintf("fig_gathering_dod_%s", pathogen), width = 6, height = 5)
cat(sprintf("  Saved fig_gathering_dod_%s\n", pathogen))

# --------------------------------------------------------------------------
# Figure 3: expected P(establishment) vs psi (Lloyd-Smith NegBin extinction)
# All three scenarios at matched R0_eff = R0 * mu_T; only k varies.
#   Restricted:   k = mu_T^2 / (V_T * I(psi*alpha, lambda/beta))
#   Unrestricted: k = sigma  / I(psi*alpha, lambda/beta)
#   Constant:     k = Inf (Poisson offspring)
# Monte Carlo points (if cached from gatheringsize_main.R) overlaid.
# --------------------------------------------------------------------------

mu_T_main   <- mu_truncated_gamma(sigma, c_max_main)
V_T_main    <- var_truncated_gamma(sigma, c_max_main)
R0_eff_main <- R0 * mu_T_main

pest_theory_df <- expand_grid(
		psi      = psi_vals_theory,
		scenario = c("Restricted", "Unrestricted", "Constant")) %>%
	rowwise() %>%
	mutate(
		k_eff = case_when(
			scenario == "Restricted"   ~ mu_T_main^2 / (V_T_main * spectral_integral(psi * alpha, lambda_gp, beta)),
			scenario == "Unrestricted" ~ sigma            / spectral_integral(psi * alpha, lambda_gp, beta),
			scenario == "Constant"     ~ Inf),
		p_est = p_est_theory(R0_eff_main, k_eff)) %>%
	ungroup()

# Look for Monte Carlo P_est from gatheringsize_main.R cache for an overlay.
# Scenario name mapping between scripts: gatheringsize_main writes "Poisson"
# for the constant-contact scenario and "R0-only" for the unrestricted-contact
# scenario; here we relabel those as "Constant" and "Unrestricted" to match.
pest_main_csv <- file.path("output",
	sprintf("gathering_pest_table_%s_cmax%g.csv", pathogen, c_max_main))

if (file.exists(pest_main_csv)) {
	pest_sim_df <- read_csv(pest_main_csv, show_col_types = FALSE) %>%
		transmute(
			psi,
			scenario = recode(as.character(scenario),
			                  Poisson    = "Constant",
			                  `R0-only`  = "Unrestricted"),
			p_est_sim = p_est,
			p_est_se  = se)
	overlay_sim <- TRUE
} else {
	pest_sim_df <- tibble(psi = double(), scenario = character(),
	                      p_est_sim = double(), p_est_se = double())
	overlay_sim <- FALSE
	cat(sprintf("  Note: no MC overlay (%s not found)\n", basename(pest_main_csv)))
}

pal_pest <- c("Restricted"   = "#FF9800",
              "Unrestricted" = "#4CAF50",
              "Constant"     = "#2196F3")

fig_pest <- ggplot() +
	geom_line(data = pest_theory_df,
	          aes(x = psi, y = p_est, color = scenario),
	          linewidth = 0.8) +
	{if (overlay_sim) geom_errorbar(data = pest_sim_df,
	          aes(x = psi,
	              ymin = pmax(p_est_sim - 1.96 * p_est_se, 0),
	              ymax = pmin(p_est_sim + 1.96 * p_est_se, 1),
	              color = scenario),
	          width = 0.015)} +
	{if (overlay_sim) geom_point(data = pest_sim_df,
	          aes(x = psi, y = p_est_sim, color = scenario),
	          size = 2)} +
	scale_color_manual(values = pal_pest, name = NULL) +
	scale_y_continuous(limits = c(0, 1)) +
	theme_classic(base_size = 16) +
	labs(x = expression(psi),
	     y = expression(P(establishment)),
	     title = pathogen_label)

save_fig(fig_pest, sprintf("fig_gathering_pest_theory_%s", pathogen),
         width = 6, height = 5)
cat(sprintf("  Saved fig_gathering_pest_theory_%s\n", pathogen))

# Tabulate expected P_est at the simulation psi grid for the supplement.
pest_theory_tbl <- pest_theory_df %>%
	filter(psi %in% psi_vals_sim) %>%
	mutate(p_est = round(p_est, 3)) %>%
	pivot_wider(names_from = scenario, values_from = c(p_est, k_eff))

cat(sprintf("\n  %s: expected P(establishment) (theory) at simulation psi grid\n",
            pathogen))
pest_theory_tbl %>% print(n = Inf, width = Inf)

} # end pathogen loop
