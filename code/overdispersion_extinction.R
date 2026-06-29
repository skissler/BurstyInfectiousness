library(tidyverse)
source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

set.seed(42)

# ==============================================================================
# Set parameters
# ==============================================================================

c_per      <- 7     # Period of the periodic contact model (days)
lambda_gp  <- 1     # Poisson switching rate (switches/day)

c_amp_vals <- c(0, 0.5, 1)         # Contact amplitude values
psi_vals   <- c(0, 0.5, 1)         # Psi values
k_c_vals   <- c(0.1, 1, 100, 1000) # Gamma/Poisson contact shape values

# Integer cast so maxinf is safe under any choice of popsize * establishment_prop.
establishment_threshold <- as.integer(round(establishment_threshold))

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Start pathogen loop
# ==============================================================================

for (idx in seq_along(parslist)) {

pars <- parslist[[idx]]
pathogen <- pars$pathogen
Tgen     <- pars$Tgen
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

p_extinct_bp <- 1 - pars$p_est  # Poisson(R0) branching-process baseline

cat(sprintf("\n===== %s: T=%.2f, alpha=%.2f, R0=%.1f (BP p_ext=%.3f) =====\n",
            pathogen, Tgen, alpha, R0, p_extinct_bp))

# ==============================================================================
# Periodic contacts: epidemic simulations
# ==============================================================================

periodic_cache <- file.path(
	"output",
	sprintf("pest_periodic_%s_n%d_p%d_cper%g.csv", pathogen, nsim, popsize, c_per)
)

need_combos <- expand_grid(c_amp = c_amp_vals, psi = psi_vals)
pest_table  <- NULL
if (file.exists(periodic_cache)) {
	cached <- read_csv(periodic_cache, show_col_types = FALSE)
	missing <- anti_join(need_combos, cached %>% distinct(c_amp, psi),
	                     by = c("c_amp", "psi"))
	if (nrow(missing) == 0) {
		cat(sprintf("  %s: loaded periodic results from %s\n", pathogen, periodic_cache))
		pest_table <- cached
	}
}

if (is.null(pest_table)) {
	epi_results <- vector("list", nsim * length(psi_vals) * length(c_amp_vals))
	epi_idx     <- 0L

	for (c_amp in c_amp_vals) {
		for (psi in psi_vals) {
			z     <- make_contact_fn_periodic(R0, c_amp, c_per)
			z_max <- R0 * (1 + c_amp)
			gfun  <- gen_inf_attempts_gamma_contacts(Tgen=Tgen, z=z, z_max=z_max,
			                                         alpha=alpha, psi=psi)
			for (sim in 1:nsim) {
				tinf <- sim_stochastic_fast(n=popsize, gen_inf_attempts=gfun,
				                            maxinf=establishment_threshold)
				n_infected <- sum(is.finite(tinf))

				epi_idx <- epi_idx + 1L
				epi_results[[epi_idx]] <- tibble(
					sim = sim, c_amp = c_amp, psi = psi,
					n_infected = n_infected,
					established = as.integer(n_infected >= establishment_threshold))

				if (sim %% 500 == 0) cat(sprintf("  [%s] c_amp=%.1f psi=%.1f: sim %d/%d\n", pathogen, c_amp, psi, sim, nsim))
			}
		}
	}

	epi_df <- bind_rows(epi_results)

	pest_table <- epi_df %>%
		group_by(c_amp, psi) %>%
		summarise(p_extinct = 1 - mean(established), .groups = "drop") %>%
		mutate(pathogen = pathogen, c_per = c_per, p_extinct_bp = p_extinct_bp)

	write_csv(pest_table, periodic_cache)
}

cat(sprintf("  %s: P(extinction) by c_amp and psi:\n", pathogen))
print(pest_table)

# ==============================================================================
# Gamma/Poisson contacts: epidemic simulations
# ==============================================================================

gp_cache <- file.path(
	"output",
	sprintf("pest_gammapoisson_%s_n%d_p%d_lam%g.csv",
	        pathogen, nsim, popsize, lambda_gp)
)

need_combos_gp <- expand_grid(psi = psi_vals, k_c = k_c_vals)
pest_table_gp <- NULL
if (file.exists(gp_cache)) {
	cached <- read_csv(gp_cache, show_col_types = FALSE)
	missing <- anti_join(need_combos_gp, cached %>% distinct(psi, k_c),
	                     by = c("psi", "k_c"))
	if (nrow(missing) == 0) {
		cat(sprintf("  %s: loaded GP results from %s\n", pathogen, gp_cache))
		pest_table_gp <- cached
	}
}

if (is.null(pest_table_gp)) {
	epi_results_gp <- vector("list", nsim * length(psi_vals) * length(k_c_vals))
	epi_idx_gp     <- 0L

	for (k_c in k_c_vals) {
		for (psi in psi_vals) {
			gfun <- gen_inf_attempts_gammapoisson_contacts(Tgen, R0, alpha, psi, k_c, lambda_gp)
			for (sim in 1:nsim) {
				tinf <- sim_stochastic_fast(n=popsize, gen_inf_attempts=gfun,
				                            maxinf=establishment_threshold)
				n_infected <- sum(is.finite(tinf))

				epi_idx_gp <- epi_idx_gp + 1L
				epi_results_gp[[epi_idx_gp]] <- tibble(
					sim = sim, psi = psi, k_c = k_c, n_infected = n_infected,
					established = as.integer(n_infected >= establishment_threshold))

				if (sim %% 500 == 0) cat(sprintf("  [%s] GP k_c=%g psi=%.1f: sim %d/%d\n",
				                                  pathogen, k_c, psi, sim, nsim))
			}
		}
	}

	epi_df_gp <- bind_rows(epi_results_gp)

	pest_table_gp <- epi_df_gp %>%
		group_by(psi, k_c) %>%
		summarise(p_extinct = 1 - mean(established), .groups = "drop") %>%
		mutate(pathogen = pathogen, lambda = lambda_gp,
		       p_extinct_bp = p_extinct_bp)

	write_csv(pest_table_gp, gp_cache)
}

cat(sprintf("  %s: P(extinction) by psi and k_c (Gamma/Poisson, lambda = %g):\n",
            pathogen, lambda_gp))
print(pest_table_gp)

} # end pathogen loop
