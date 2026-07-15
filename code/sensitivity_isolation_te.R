library(tidyverse)
library(patchwork)

source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

# ==============================================================================
# SENSITIVITY ANALYSIS: TE curves under alternative burst models (Omicron only)
#
# Replicates the four detection scenarios from isolation_te.R under:
#   - Type-I Gamma  (original): epsilon ~ Gamma(alpha*psi, beta)
#   - Log-normal:               epsilon ~ LogNormal(mu_tau + 0.5*log(psi), sigma_tau)
#   - Type-II Gamma:            epsilon ~ Gamma(alpha, beta/sqrt(psi))
#
# Detection is anchored to the mode of each model's individual infectiousness
# profile (the distribution of epsilon_j). Produces 12 figures:
#   4 detection mechanisms x 3 burst models, all for SARS-CoV-2 omicron.
#
# Detection mechanisms (identical parameters to isolation_te.R):
#   1. Fixed-time isolation
#   2. Symptom-triggered isolation
#   3. Regular testing (no action delay)
#   4. Regular testing with exponential action delay
# ==============================================================================

# --- Omicron parameters -------------------------------------------------------
omicron  <- parslist[[2]]          # omicron is entry 2 in parslist
pathogen <- omicron$pathogen       # "omicron"
Tgen     <- omicron$Tgen
alpha    <- omicron$alpha
beta     <- omicron$beta
R0       <- omicron$R0

cat(sprintf("\n=== sensitivity_isolation_te.R — %s (R0=%g, alpha=%.2f, beta=%.3f) ===\n",
			pathogen, R0, alpha, beta))

# --- Shared detection parameters (matching isolation_te.R) --------------------
psi_vals_te <- c(0, 0.2, 0.5, 0.8, 1)
sigma_sym   <- 0.5       # SD of symptom onset relative to peak (days)
d_pre       <- 3         # days before peak at which test detects
d_post      <- 7         # days after peak at which test stops detecting
w           <- d_pre + d_post
Delta_vals  <- seq(0.5, 14, by = 0.5)
lambda_act  <- 1         # rate of exponential action delay (mean 1 day)

# x-axis range based on the Gamma(alpha,beta) tau distribution (2.5 SDs)
sd_tau           <- sqrt(alpha) / beta
tau_offset_range <- 2.5 * sd_tau
tau_offset       <- seq(-tau_offset_range, tau_offset_range, length.out = 201)
mu_sym_vals      <- seq(-tau_offset_range, tau_offset_range, length.out = 201)

# --- Log-normal parameters for tau (moment-matched to Gamma(alpha, beta)) -----
sigma2_tau_ln <- log(1 + 1 / alpha)
sigma_tau_ln  <- sqrt(sigma2_tau_ln)
mu_tau_ln     <- log(Tgen) - sigma2_tau_ln / 2

# ==============================================================================
# Burst model specifications
#
# eps_mode(psi)   — mode of epsilon_j distribution (scalar -> scalar)
# eps_cdf(t, psi) — P(epsilon_j <= t), vectorised over t
#
# For the testing-with-delay closed form (Gamma-family only):
#   E_D[1 - F_eps(t0 + D)] = (1 - pgamma(t0, a, b))
#                             - (b/(b+lam))^a * exp(lam*t0) * (1 - pgamma(t0, a, b+lam))
# where a = delay_a(psi), b = delay_b(psi).
# ==============================================================================

burst_specs <- list(

	gamma = list(
		label        = "Type-I Gamma",
		eps_mode     = function(psi) pmax(0, (alpha * psi - 1) / beta),
		eps_cdf      = function(t, psi) pgamma(t, alpha * psi, beta),
		delay_closed = TRUE,
		delay_a      = function(psi) alpha * psi,
		delay_b      = function(psi) beta
	),

	lognormal = list(
		label    = "Log-normal",
		eps_mode = function(psi) {
			if (psi < 1e-8) return(0)
			exp(mu_tau_ln + 0.5 * log(psi) - sigma2_tau_ln)
		},
		eps_cdf = function(t, psi) {
			if (psi < 1e-8) return(as.numeric(t >= 0))
			plnorm(t, mu_tau_ln + 0.5 * log(psi), sigma_tau_ln)
		},
		delay_closed = FALSE,
		delay_a      = NULL,
		delay_b      = NULL
	),

	gamma2 = list(
		label    = "Type-II Gamma",
		eps_mode = function(psi) {
			if (psi < 1e-8) return(0)
			pmax(0, (alpha - 1) * sqrt(psi) / beta)
		},
		eps_cdf = function(t, psi) {
			if (psi < 1e-8) return(as.numeric(t >= 0))
			pgamma(t, alpha, beta / sqrt(psi))
		},
		delay_closed = TRUE,
		delay_a      = function(psi) alpha,
		delay_b      = function(psi) beta / sqrt(psi)   # -> Inf when psi=0 (handled by pgamma)
	)
)

# ==============================================================================
# Closed-form testing-with-delay integrand for Gamma-family epsilon
# ==============================================================================

delay_integrand_closed <- function(u, mode_e, a, b) {
	t0 <- u - d_pre + mode_e
	if (is.infinite(b)) {
		# eps degenerate at 0 (psi=0 in gamma2): E_D[P(eps > t0+D)] analytically
		return(ifelse(t0 < 0, 1 - exp(lambda_act * t0), 0))
	}
	cr <- (b / (b + lambda_act))^a
	(1 - pgamma(t0, a, b)) -
		exp(lambda_act * t0) * cr * (1 - pgamma(t0, a, b + lambda_act))
}

# Storage for composite figures
te_fixed_list   <- list()
te_symp_list    <- list()
te_testing_list <- list()
te_delay_list   <- list()

# ==============================================================================
# Main loop: 3 burst models x 4 detection mechanisms = 12 figures
# ==============================================================================

for (mod in names(burst_specs)) {

	spec  <- burst_specs[[mod]]
	label <- spec$label

	cat(sprintf("\n  --- %s ---\n", label))

	# ---- 1. Fixed-time isolation TE ------------------------------------------
	# TE(tau_offset) = 1 - F_eps(tau_offset + mode_eps)
	# x-axis: isolation time relative to mode of individual profile

	te_fixed_df <- expand_grid(psi = psi_vals_te, tau_offset = tau_offset) %>%
		rowwise() %>%
		mutate(
			mode_e = spec$eps_mode(psi),
			TE     = 1 - spec$eps_cdf(tau_offset + mode_e, psi)
		) %>%
		ungroup()

	fig_fixed <- ggplot(te_fixed_df,
						aes(x = tau_offset, y = TE, col = factor(psi))) +
		geom_line(linewidth = 0.8, alpha = 0.8) +
		theme_classic(base_size = 14) +
		labs(x     = "Isolation time (days relative to peak infectiousness)",
			 y     = "Transmission effectiveness",
			 col   = expression(psi),
			 title = label)

	save_fig(fig_fixed, sprintf("sensitivity_te_fixed_%s_%s", mod, pathogen))
	te_fixed_list[[mod]] <- fig_fixed
	cat("    Fixed isolation done.\n")

	# ---- 2. Symptom-triggered isolation TE -----------------------------------
	# TE(mu_sym) = E_{t ~ N(mu_sym, sigma_sym)}[1 - F_eps(t + mode_eps)]
	# x-axis: mean symptom onset time relative to mode of individual profile

	te_symp_df <- expand_grid(psi = psi_vals_te, mu_sym = mu_sym_vals) %>%
		rowwise() %>%
		mutate(
			mode_e = spec$eps_mode(psi),
			TE     = integrate(
						 function(t) (1 - spec$eps_cdf(t + mode_e, psi)) *
									 dnorm(t, mu_sym, sigma_sym),
						 lower = mu_sym - 5 * sigma_sym,
						 upper = mu_sym + 5 * sigma_sym)$value
		) %>%
		ungroup()

	fig_symp <- ggplot(te_symp_df,
					   aes(x = mu_sym, y = TE, col = factor(psi))) +
		geom_line(linewidth = 0.8, alpha = 0.8) +
		theme_classic(base_size = 14) +
		labs(x     = "Mean symptom onset (days relative to peak infectiousness)",
			 y     = "Transmission effectiveness",
			 col   = expression(psi),
			 title = label)

	save_fig(fig_symp, sprintf("sensitivity_te_symp_%s_%s", mod, pathogen))
	te_symp_list[[mod]] <- fig_symp
	cat("    Symptom-triggered isolation done.\n")

	# ---- 3. Regular testing TE (no action delay) -----------------------------
	# TE(Delta) = (1/Delta) * integral_0^min(Delta,w) [1 - F_eps(u - d_pre + mode_eps)] du
	# Detection at uniform time in [0, w]; x-axis: gap between tests

	te_testing_df <- expand_grid(psi = psi_vals_te, Delta = Delta_vals) %>%
		rowwise() %>%
		mutate(
			mode_e = spec$eps_mode(psi),
			TE     = integrate(
						 function(u) 1 - spec$eps_cdf(u - d_pre + mode_e, psi),
						 lower = 0, upper = min(Delta, w))$value / Delta
		) %>%
		ungroup()

	fig_testing <- ggplot(te_testing_df,
						  aes(x = Delta, y = TE, col = factor(psi))) +
		geom_line(linewidth = 0.8, alpha = 0.8) +
		theme_classic(base_size = 14) +
		labs(x     = "Gap between tests (days)",
			 y     = "Transmission effectiveness",
			 col   = expression(psi),
			 title = label)

	save_fig(fig_testing, sprintf("sensitivity_te_testing_%s_%s", mod, pathogen))
	te_testing_list[[mod]] <- fig_testing
	cat("    Regular testing done.\n")

	# ---- 4. Regular testing TE with exponential action delay -----------------
	# For Gamma-family epsilon: closed-form integrand via shifted Gamma identity.
	# For log-normal epsilon: numerical double integration.

	if (spec$delay_closed) {
		te_delay_df <- expand_grid(psi = psi_vals_te, Delta = Delta_vals) %>%
			rowwise() %>%
			mutate(
				mode_e = spec$eps_mode(psi),
				TE     = integrate(
							 function(u) delay_integrand_closed(
											 u, mode_e,
											 spec$delay_a(psi),
											 spec$delay_b(psi)),
							 lower = 0, upper = min(Delta, w))$value / Delta
			) %>%
			ungroup()
	} else {
		# integrate() passes a vector of u values to its function; the inner
		# integrate() requires a scalar t0, so we loop over u with sapply.
		te_delay_df <- expand_grid(psi = psi_vals_te, Delta = Delta_vals) %>%
			rowwise() %>%
			mutate(
				mode_e = spec$eps_mode(psi),
				TE     = integrate(function(u) {
							 sapply(u, function(u_i) {
								 t0 <- u_i - d_pre + mode_e
								 if (psi < 1e-8) {
									 # eps = 0 (degenerate): E_D[P(eps > t0+D)]
									 # = P(D < -t0) = 1 - exp(lam*t0) if t0<0, else 0
									 if (t0 >= 0) 0 else 1 - exp(lambda_act * t0)
								 } else {
									 integrate(function(d)
										 (1 - spec$eps_cdf(t0 + d, psi)) *
										 lambda_act * exp(-lambda_act * d),
										 lower = 0, upper = 50)$value
								 }
							 })
						 }, lower = 0, upper = min(Delta, w))$value / Delta
			) %>%
			ungroup()
	}

	fig_delay <- ggplot(te_delay_df,
						aes(x = Delta, y = TE, col = factor(psi))) +
		geom_line(linewidth = 0.8, alpha = 0.8) +
		theme_classic(base_size = 14) +
		labs(x     = "Gap between tests (days)",
			 y     = "Transmission effectiveness",
			 col   = expression(psi),
			 title = label)

	save_fig(fig_delay, sprintf("sensitivity_te_delay_%s_%s", mod, pathogen))
	te_delay_list[[mod]] <- fig_delay
	cat("    Testing with action delay done.\n")

} # end burst model loop

# ==============================================================================
# Composite figures: one per detection mechanism, panels = burst models
# ==============================================================================

fig_te_fixed_composite <- wrap_plots(te_fixed_list, nrow = 1)
save_fig(fig_te_fixed_composite,
		 sprintf("sensitivity_te_fixed_%s", pathogen), width = 14, height = 5)

fig_te_symp_composite <- wrap_plots(te_symp_list, nrow = 1)
save_fig(fig_te_symp_composite,
		 sprintf("sensitivity_te_symp_%s", pathogen), width = 14, height = 5)

fig_te_testing_composite <- wrap_plots(te_testing_list, nrow = 1)
save_fig(fig_te_testing_composite,
		 sprintf("sensitivity_te_testing_%s", pathogen), width = 14, height = 5)

fig_te_delay_composite <- wrap_plots(te_delay_list, nrow = 1)
save_fig(fig_te_delay_composite,
		 sprintf("sensitivity_te_delay_%s", pathogen), width = 14, height = 5)

cat("\nSensitivity TE analysis complete. 12 individual + 4 composite figures saved.\n")
