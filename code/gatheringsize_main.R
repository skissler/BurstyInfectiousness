# ==============================================================================
# Gathering size restrictions: three-scenario comparison at a focal c_max
# ==============================================================================
#
# Compares establishment probability under, at c_max = c_max_main:
#   1. "Restricted" — truncated contact process (true R0_new, true k_new)
#   2. "R0-only"   — original contact process, R0 scaled down to R0_new
#   3. "Poisson"   — no contact process, Poisson(R0_new) offspring
#
# Expected ordering: P_est_R0only < P_est_restricted < P_est_Poisson,
# with gaps larger for small psi (punctuated profiles).
# ==============================================================================

library(tidyverse)
source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

# --- Parameters ---
k_c        <- 1                          # moderate contact heterogeneity
lambda_gp  <- 1                          # Poisson switching rate
psi_vals   <- c(0, 0.25, 0.5, 0.75, 1)
c_max_main <- 2                          # focal restriction for main figure

establishment_threshold <- as.integer(round(establishment_threshold))

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Start pathogen loop
# ==============================================================================

for (idx_pathogen in seq_along(parslist)) {

pars <- parslist[[idx_pathogen]]
pathogen <- pars$pathogen
Tgen     <- pars$Tgen
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

# Per-pathogen seed so results don't depend on the order of earlier pathogens.
set.seed(42 + idx_pathogen)

mu_T_main   <- mu_truncated_gamma(k_c, c_max_main)
R0_new_main <- R0 * mu_T_main

cat(sprintf("\n===== Gathering size (main): %s (T=%.2f, alpha=%.2f, R0=%.1f, c_max=%g) =====\n",
            pathogen, Tgen, alpha, R0, c_max_main))
cat(sprintf("  mu_T = %.4f, R0_new = %.3f\n", mu_T_main, R0_new_main))

# ==============================================================================
# Three-scenario epidemic simulations at the focal c_max
# ==============================================================================

scenario_names <- c("Restricted", "R0-only", "Poisson")
need_combos    <- expand_grid(psi = psi_vals, scenario = scenario_names)

sims_cache <- file.path(
	"output",
	sprintf("gathering_sims_main_%s_n%d_p%d_kc%g_lam%g_cmax%g.csv",
	        pathogen, nsim, popsize, k_c, lambda_gp, c_max_main)
)

# For backward compatibility with the original combined cache (which held all
# c_max values), fall back to it if the dedicated cache is missing.
sims_combined_cache <- file.path(
	"output",
	sprintf("gathering_sims_%s_n%d_p%d_kc%g_lam%g.csv",
	        pathogen, nsim, popsize, k_c, lambda_gp)
)

epi_df <- NULL

if (file.exists(sims_cache)) {
	cached  <- read_csv(sims_cache, show_col_types = FALSE)
	missing <- anti_join(need_combos,
	                     cached %>% distinct(psi, scenario),
	                     by = c("psi", "scenario"))
	if (nrow(missing) == 0) {
		cat(sprintf("  %s: loaded from %s\n", pathogen, basename(sims_cache)))
		epi_df <- cached
	}
}

if (is.null(epi_df) && file.exists(sims_combined_cache)) {
	cached <- read_csv(sims_combined_cache, show_col_types = FALSE) %>%
		filter(c_max == c_max_main) %>%
		select(-c_max)
	missing <- anti_join(need_combos,
	                     cached %>% distinct(psi, scenario),
	                     by = c("psi", "scenario"))
	if (nrow(missing) == 0) {
		cat(sprintf("  %s: loaded c_max = %g subset from combined cache %s\n",
		    pathogen, c_max_main, basename(sims_combined_cache)))
		epi_df <- cached
		write_csv(epi_df, sims_cache)
	}
}

if (is.null(epi_df)) {

	cells <- expand_grid(
			psi_idx      = seq_along(psi_vals),
			scenario_idx = seq_along(scenario_names)) %>%
		mutate(
			psi      = psi_vals[psi_idx],
			scenario = scenario_names[scenario_idx],
			seed     = 42L + 1000L * idx_pathogen + 100L * psi_idx + scenario_idx)

	run_cell <- function(psi_val, scenario_name, seed) {
		set.seed(seed)
		gfun <- switch(scenario_name,
			"Restricted" = gen_inf_attempts_gammapoisson_contacts(
			                  Tgen, R0, alpha, psi_val, k_c, lambda_gp, c_max = c_max_main),
			"R0-only"    = gen_inf_attempts_gammapoisson_contacts(
			                  Tgen, R0_new_main, alpha, psi_val, k_c, lambda_gp),
			"Poisson"    = gen_inf_attempts_gamma(Tgen, R0_new_main, alpha, psi_val))

		map_dfr(seq_len(nsim), function(sim) {
			tinf <- sim_stochastic_fast(n = popsize,
			                            gen_inf_attempts = gfun,
			                            maxinf = establishment_threshold)
			n_infected <- sum(is.finite(tinf))
			tibble(
				psi         = psi_val,
				scenario    = scenario_name,
				sim         = sim,
				n_infected  = n_infected,
				established = as.integer(n_infected >= establishment_threshold))
		})
	}

	n_cores <- max(1L, parallel::detectCores() - 1L)
	cat(sprintf("  Running %d cells across %d core(s)...\n", nrow(cells), n_cores))

	epi_df <- parallel::mclapply(seq_len(nrow(cells)), function(i) {
			run_cell(cells$psi[i], cells$scenario[i], cells$seed[i])
		}, mc.cores = n_cores) %>%
		bind_rows()
	write_csv(epi_df, sims_cache)
}

# ==============================================================================
# Compute P(establishment) summary
# ==============================================================================

pest_summary <- epi_df %>%
	group_by(psi, scenario) %>%
	summarise(
		p_est  = mean(established),
		n_sims = n(),
		se     = sqrt(p_est * (1 - p_est) / n_sims),
		.groups = "drop"
	) %>%
	mutate(scenario = factor(scenario, levels = c("Poisson", "Restricted", "R0-only")))

cat(sprintf("\n  %s: P(establishment) summary\n", pathogen))
pest_summary %>% select(psi, scenario, p_est, se) %>% print(n = Inf)

# ==============================================================================
# Verification: simulated mean offspring vs analytical R0_new
# ==============================================================================

cat(sprintf("\n  %s: R0 cross-check (simulated mean offspring vs analytical R0_new)\n",
    pathogen))
gfun_check  <- gen_inf_attempts_gammapoisson_contacts(
	Tgen, R0, alpha, 1, k_c, lambda_gp, c_max = c_max_main)
n_offspring <- replicate(5000, length(gfun_check(0)))
R0_sim      <- mean(n_offspring)
cat(sprintf("    c_max = %g: R0_new(theory) = %.3f, R0_new(sim) = %.3f\n",
    c_max_main, R0_new_main, R0_sim))

# ==============================================================================
# Main figure: three scenarios at fixed c_max, across psi
# ==============================================================================

fig_gathering_pest <- ggplot(pest_summary,
		aes(x = factor(psi), y = p_est, fill = scenario)) +
	geom_col(position = position_dodge(width = 0.7), width = 0.6) +
	geom_errorbar(aes(ymin = pmax(p_est - 1.96 * se, 0),
	                  ymax = pmin(p_est + 1.96 * se, 1)),
	              position = position_dodge(width = 0.7), width = 0.2) +
	scale_fill_manual(
		values = c("Poisson" = "#2196F3", "Restricted" = "#FF9800", "R0-only" = "#4CAF50"),
		name = "Scenario") +
	theme_classic() +
	labs(x = expression(psi),
	     y = expression(P(establishment)),
	     title = sprintf("%s: gathering cap c_max = %g (R0: %.1f -> %.2f)",
	                     pathogen, c_max_main, R0, R0_new_main))

save_fig(fig_gathering_pest, sprintf("fig_gathering_pest_%s", pathogen), width = 8, height = 5)
cat(sprintf("  Saved fig_gathering_pest_%s\n", pathogen))

# ==============================================================================
# Main table: establishment probabilities with standard errors (focal c_max)
# Mirrors the main figure data, written to output/ in two forms:
#   - CSV (tidy, machine-readable): one row per psi x scenario
#   - Markdown (human-readable): psi as rows, scenarios as columns, "p_est ± se"
# ==============================================================================

pest_tbl <- pest_summary %>%
	transmute(
		psi,
		scenario,
		p_est,
		se,
		ci_lower = pmax(p_est - 1.96 * se, 0),
		ci_upper = pmin(p_est + 1.96 * se, 1),
		n_sims
	) %>%
	arrange(scenario, psi)

pest_table_csv <- file.path("output",
	sprintf("gathering_pest_table_%s_cmax%g.csv", pathogen, c_max_main))
write_csv(pest_tbl, pest_table_csv)

# Human-readable wide table: cells show p_est ± se to 3 decimals
pest_wide <- pest_tbl %>%
	mutate(cell = sprintf("%.3f ± %.3f", p_est, se)) %>%
	select(psi, scenario, cell) %>%
	pivot_wider(names_from = scenario, values_from = cell) %>%
	arrange(psi)

n_per_cell <- max(pest_tbl$n_sims)
pest_table_md <- file.path("output",
	sprintf("gathering_pest_table_%s_cmax%g.md", pathogen, c_max_main))
md_lines <- c(
	sprintf("# %s: P(establishment) ± SE at gathering cap c_max = %g", pathogen, c_max_main),
	"",
	sprintf("R0: %.1f -> %.2f; %d simulations per cell; cells show p_est ± SE.",
	        R0, R0_new_main, n_per_cell),
	"",
	knitr::kable(pest_wide, format = "pipe", align = "c")
)
writeLines(md_lines, pest_table_md)
cat(sprintf("  Saved %s and %s\n", basename(pest_table_csv), basename(pest_table_md)))

} # end pathogen loop
