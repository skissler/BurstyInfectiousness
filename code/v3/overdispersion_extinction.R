# ==============================================================================
# Overdispersion: extinction probability simulations (periodic and Gamma/Poisson)
# ==============================================================================

library(tidyverse)
source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

set.seed(42)

# c_per             <- 1
c_per             <- 7
lambda_gp         <- 1     # Poisson switching rate (switches/day)

c_amp_epi_vals    <- c(0, 0.5, 1)
# psi_epi         <- seq(from=0, to=1, by=0.1)
psi_epi           <- c(0, 0.5, 1)

k_c_epi           <- 1     # moderate heterogeneity in contact levels

make_contact_fn_periodic <- function(R0, c_amp, c_per) {
	stopifnot(c_amp >= 0, c_amp <= 1)
	function(t) R0*(1 - c_amp*cos(2*pi*t/c_per))
}

# ==============================================================================
# Start pathogen loop
# ==============================================================================

for (pars in parslist) {

pathogen <- pars$pathogen
Tgen     <- pars$Tgen
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

cat(sprintf("\n===== %s: T=%.2f, alpha=%.2f, R0=%.1f =====\n",
            pathogen, Tgen, alpha, R0))

# ==============================================================================
# Periodic contacts: epidemic simulations
# ==============================================================================

epi_results <- vector("list", nsim * length(psi_epi) * length(c_amp_epi_vals))
epi_idx     <- 0L

for (c_amp_val in c_amp_epi_vals) {
	for (psi_val in psi_epi) {
		z     <- make_contact_fn_periodic(R0, c_amp_val, c_per)
		z_max <- R0 * (1 + c_amp_val)
		gfun  <- gen_inf_attempts_gamma_contacts(Tgen=Tgen, z=z, z_max=z_max,
		                                         alpha=alpha, psi=psi_val)
		for (sim in 1:nsim) {
			tinf <- sim_stochastic_fast(n=popsize, gen_inf_attempts=gfun,
			                            maxinf=establishment_threshold)
			n_infected <- sum(is.finite(tinf))

			epi_idx <- epi_idx + 1L
			epi_results[[epi_idx]] <- tibble(
				sim = sim, c_amp = c_amp_val, psi = psi_val,
				n_infected = n_infected,
				established = as.integer(n_infected >= establishment_threshold))

			if (sim %% 500 == 0) cat(sprintf("  [%s] c_amp=%.1f psi=%.1f: sim %d/%d\n",
			                                  pathogen, c_amp_val, psi_val, sim, nsim))
		}
	}
}

epi_df <- bind_rows(epi_results)

pest_table <- epi_df %>%
	group_by(c_amp, psi) %>%
	summarise(p_extinct = 1 - mean(established), .groups = "drop")

cat(sprintf("  %s: P(extinction) by c_amp and psi:\n", pathogen))
print(pest_table)

write_csv(
	pest_table %>% mutate(pathogen = pathogen),
	file.path("output", sprintf("pest_periodic_%s.csv", pathogen))
)

# ==============================================================================
# Gamma/Poisson contacts: epidemic simulations
# ==============================================================================

epi_results_gp <- vector("list", nsim * length(psi_epi))
epi_idx_gp     <- 0L

for (psi_val in psi_epi) {
	gfun <- gen_inf_attempts_gammapoisson_contacts(Tgen, R0, alpha, psi_val, k_c_epi, lambda_gp)
	for (sim in 1:nsim) {
		tinf <- sim_stochastic_fast(n=popsize, gen_inf_attempts=gfun,
		                            maxinf=establishment_threshold)
		n_infected <- sum(is.finite(tinf))

		epi_idx_gp <- epi_idx_gp + 1L
		epi_results_gp[[epi_idx_gp]] <- tibble(
			sim = sim, psi = psi_val, n_infected = n_infected,
			established = as.integer(n_infected >= establishment_threshold))

		if (sim %% 500 == 0) cat(sprintf("  [%s] GP psi=%.1f: sim %d/%d\n",
		                                  pathogen, psi_val, sim, nsim))
	}
}

epi_df_gp <- bind_rows(epi_results_gp)

pest_table_gp <- epi_df_gp %>%
	group_by(psi) %>%
	summarise(p_extinct = 1 - mean(established), .groups = "drop")

cat(sprintf("  %s: P(extinction) by psi (Gamma/Poisson, k_c = %g, lambda = %g):\n",
            pathogen, k_c_epi, lambda_gp))
print(pest_table_gp)

write_csv(
	pest_table_gp %>% mutate(pathogen = pathogen, k_c = k_c_epi, lambda = lambda_gp),
	file.path("output", sprintf("pest_gammapoisson_%s.csv", pathogen))
)

} # end pathogen loop
