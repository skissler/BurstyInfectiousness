library(tidyverse)
library(patchwork)

source("code/utils.R")
source("code/parameters.R")

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Global parameters (shared across pathogens)
# ==============================================================================

psi_vals_te <- c(0, 0.2, 0.5, 0.8, 1)  # for TE curves

# Symptom-triggered isolation
sigma_sym <- 0.5       # SD of symptom onset (days) relative to profile peak

# Regular screening
d_pre  <- 3            # days before peak at which test detects
d_post <- 7            # days after peak at which test stops detecting
w      <- d_pre + d_post  # detectability window
Delta_vals <- seq(from = 0.5, to = 14, by = 0.5)  # time between tests
lambda_act <- 1        # rate of exponential delay to action (mean 1 day)

# Storage for composite figures and combined CSV
te_fixed_list         <- list()
te_symp_list          <- list()
te_symp_interp_list   <- list()
te_testing_list       <- list()
te_testing_delay_list <- list()

# ==============================================================================
# Loop over pathogens
# ==============================================================================

for (pars in parslist) {

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

cat(sprintf("\n========== %s (R0=%g, alpha=%.2f, beta=%.3f) ==========\n",
            pathogen, R0, alpha, beta))

# Scale the tau_offset / mu_sym axis to cover ~2.5 SDs of the profile
sd_g             <- sqrt(alpha) / beta 
tau_offset_range <- 2.5 * sd_g
tau_offset       <- seq(-tau_offset_range, tau_offset_range, length.out = 201)
mu_sym           <- seq(-tau_offset_range, tau_offset_range, length.out = 201)

# ==============================================================================
# Fixed isolation TE curves
# ==============================================================================

te_fixed_df <- expand_grid(psi = psi_vals_te, tau_offset = tau_offset) %>%
	mutate(mode = pmax(0, ((alpha*psi - 1)/beta))) %>%
	mutate(TE = 1 - pgamma(tau_offset + mode, alpha*psi, beta))

fig_te_fixed <- te_fixed_df %>%
	ggplot(aes(x = tau_offset, y = TE, col = factor(psi))) +
		geom_line(linewidth = 0.8, alpha = 0.8) +
		theme_classic(base_size = 14) +
		labs(x = "Isolation time\n(days relative to peak infectiousness)",
		     y = "Test effectiveness",
		     col = expression(psi),
		     title = pathogen_label)

save_fig(fig_te_fixed, sprintf("fig_te_fixed_%s", pathogen))
te_fixed_list[[pathogen]] <- fig_te_fixed

# ==============================================================================
# Symptom-triggered TE curves
# ==============================================================================

te_symp_df <- expand_grid(psi = psi_vals_te, mu_sym = mu_sym) %>%
	mutate(mode = pmax(0, ((alpha*psi - 1)/beta))) %>%
	rowwise() %>%
	mutate(TE = integrate(function(t)
		(1 - pgamma(t + mode, alpha*psi, beta)) *
		dnorm(t, mu_sym, sigma_sym),
		lower = mu_sym - 5*sigma_sym, upper = mu_sym + 5*sigma_sym)$value) %>%
	ungroup()

fig_te_symp <- te_symp_df %>%
	ggplot(aes(x = mu_sym, y = TE, col = factor(psi))) +
		geom_line(linewidth = 0.8, alpha = 0.8) +
		theme_classic(base_size = 14) +
		labs(x = "Mean symptom onset time\n(days relative to peak infectiousness)",
		     y = "Test effectiveness",
		     col = expression(psi),
		     title = pathogen_label)

save_fig(fig_te_symp, sprintf("fig_te_symp_%s", pathogen))
te_symp_list[[pathogen]] <- fig_te_symp

# Interpolated TE at each integer day for each psi
day_range <- floor(min(mu_sym)):ceiling(max(mu_sym))
te_symp_interp <- map_dfr(psi_vals_te, function(p) {
  sub <- te_symp_df %>% filter(psi == p)
  tibble(psi = p, day = day_range,
         TE  = approx(sub$mu_sym, sub$TE, xout = day_range, rule = 2)$y)
}) %>% mutate(pathogen = pathogen)

te_symp_interp_list[[pathogen]] <- te_symp_interp
cat(sprintf("\n  [%s] TE (symptom-triggered) at integer days:\n", pathogen))
print(te_symp_interp %>%
        pivot_wider(names_from = psi, names_prefix = "psi=",
                    values_from = TE) %>%
        mutate(across(starts_with("psi="), ~ round(.x, 4))), n=Inf)

# ==============================================================================
# Regular testing TE (perfect sensitivity, no delay)
# ==============================================================================

te_testing_df <- expand_grid(psi = psi_vals_te, Delta = Delta_vals) %>%
	mutate(mode = pmax(0, ((alpha*psi - 1)/beta))) %>%
	rowwise() %>%
	mutate(TE = integrate(function(u)
		(1 - pgamma(u - d_pre + mode, alpha*psi, beta)),
		lower = 0, upper = min(Delta, w))$value / Delta) %>%
	ungroup()

fig_te_testing <- te_testing_df %>%
	ggplot(aes(x = Delta, y = TE, col = factor(psi))) +
		geom_line(linewidth = 0.8, alpha = 0.8) +
		theme_classic(base_size = 14) +
		labs(x = "Gap between tests (days)",
		     y = "Test effectiveness",
		     col = expression(psi),
		     title = pathogen_label)

save_fig(fig_te_testing, sprintf("fig_te_testing_%s", pathogen))
te_testing_list[[pathogen]] <- fig_te_testing

# ==============================================================================
# Regular testing TE with exponential action delay
# ==============================================================================

te_testing_delay_df <- expand_grid(psi = psi_vals_te, Delta = Delta_vals) %>%
	mutate(mode = pmax(0, ((alpha*psi - 1)/beta))) %>%
	rowwise() %>%
	mutate(TE = {
		a  <- alpha*psi
		cr <- (beta / (beta + lambda_act))^a
		integrate(function(u) {
			t0 <- u - d_pre + mode
			(1 - pgamma(t0, a, beta)) -
			exp(lambda_act*t0) * cr * (1 - pgamma(t0, a, beta + lambda_act))
		}, lower = 0, upper = min(Delta, w))$value / Delta
	}) %>%
	ungroup()

fig_te_testing_delay <- te_testing_delay_df %>%
	ggplot(aes(x = Delta, y = TE, col = factor(psi))) +
		geom_line(linewidth = 0.8, alpha = 0.8) +
		theme_classic(base_size = 14) +
		labs(x = "Gap between tests (days)",
		     y = "Test effectiveness",
		     col = expression(psi),
		     title = pathogen_label)

save_fig(fig_te_testing_delay, sprintf("fig_te_testing_delay_%s", pathogen))
te_testing_delay_list[[pathogen]] <- fig_te_testing_delay

} # end pathogen loop

# ==============================================================================
# Composite figures across pathogens
# ==============================================================================

fig_te_fixed_composite <- wrap_plots(te_fixed_list, nrow = 1)
save_fig(fig_te_fixed_composite, "fig_te_fixed", width = 14, height = 5)

fig_te_symp_composite <- wrap_plots(te_symp_list, nrow = 1)
save_fig(fig_te_symp_composite, "fig_te_symp", width = 14, height = 5)

fig_te_testing_composite <- wrap_plots(te_testing_list, nrow = 1)
save_fig(fig_te_testing_composite, "fig_te_testing", width = 14, height = 5)

fig_te_testing_delay_composite <- wrap_plots(te_testing_delay_list, nrow = 1)
save_fig(fig_te_testing_delay_composite, "fig_te_testing_delay", width = 14, height = 5)


te_symp_interp_all <- bind_rows(te_symp_interp_list) %>%
  select(pathogen, psi, day, TE) %>%
  arrange(pathogen, psi, day)

write_csv(te_symp_interp_all,
          file.path("output", "te_symp_by_day.csv"))
cat("\nTE (symptom-triggered) by pathogen, psi, and day saved to output/te_symp_by_day.csv\n")

cat("\nAll figures saved to figures/.\n")
