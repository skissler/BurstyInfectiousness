library(tidyverse)

source("code/utils.R")
source("code/parameters.R")

# ==============================================================================
# Parameters
# ==============================================================================

omicron  <- Filter(function(p) p$pathogen == "omicron", parslist)[[1]]
alpha    <- omicron$alpha
beta     <- omicron$beta
gi_mean  <- omicron$Tgen
gi_sd    <- sqrt(omicron$Tvar)

Delta    <- 2           # Deterministic mode-anchored detection offset (days)
rho_eta  <- 1           # Combined participation x effectiveness

psi_vals <- c(0, 0.5, 1)
tau_max  <- 14
tau_grid <- seq(1e-3, tau_max, length.out = 1001)

# ==============================================================================
# Analytic g*(tau) under deterministic Delta
# ==============================================================================
#
# g*(tau) = [ (1 - rho_eta) * g(tau) + rho_eta * I(tau) ] / (1 - TE)
#
# where
#   I(tau) = integral_0^{min(tau, c)} f_l(tau - eps) f_eps(eps) deps
#   c      = m_eps + Delta
#   f_l    ~ Gamma((1-psi)*alpha, beta)
#   f_eps  ~ Gamma(psi*alpha,     beta)
#   TE     = rho_eta * (1 - F_eps(c))
#
# The integral I(tau) is computed numerically.

gstar_density <- function(tau_grid, psi, alpha, beta, Delta, rho_eta) {
	g_tau <- dgamma(tau_grid, shape = alpha, rate = beta)

	# Spike limit (psi -> 0): epsilon is a point mass at 0, so all attempts from
	# index i happen at tau = l_i. Suppression depends only on the sign of Delta
	# and is independent of l_i, so the surviving GI density equals g(tau).
	# Returned directly to avoid dgamma(., shape = 0) = NaN issues.
	if (psi < 1e-9) return(g_tau)

	# Smooth limit (psi -> 1): l is a point mass at 0, so tau = eps and the
	# convolution collapses to g(tau) * 1(tau <= c_iso). Returned directly to
	# avoid dgamma(., shape = 0) = NaN issues.
	if (psi > 1 - 1e-9) {
		m_eps <- (alpha - 1) / beta
		c_iso <- m_eps + Delta
		TE    <- rho_eta * (1 - pgamma(c_iso, shape = alpha, rate = beta))
		return(((1 - rho_eta) * g_tau + rho_eta * g_tau * (tau_grid <= c_iso)) / (1 - TE))
	}

	alpha_l <- (1 - psi) * alpha
	alpha_e <- psi       * alpha
	m_eps   <- if (alpha_e >= 1) (alpha_e - 1) / beta else 0
	c_iso   <- m_eps + Delta

	I_tau <- sapply(tau_grid, function(t) {
		upper <- min(t, c_iso)
		if (upper <= 0) return(0)
		integrand <- function(eps) {
			dgamma(t - eps, shape = alpha_l, rate = beta) *
			dgamma(eps,     shape = alpha_e, rate = beta)
		}
		tryCatch(
			integrate(integrand, lower = 0, upper = upper,
			          rel.tol = 1e-6, subdivisions = 1000L,
			          stop.on.error = FALSE)$value,
			error = function(e) NA_real_
		)
	})

	TE <- rho_eta * (1 - pgamma(c_iso, shape = alpha_e, rate = beta))

	((1 - rho_eta) * g_tau + rho_eta * I_tau) / (1 - TE)
}

# ==============================================================================
# Evaluate g*(tau) for each psi
# ==============================================================================

density_df <- map_dfr(psi_vals, function(psi) {
	tibble(
		tau     = tau_grid,
		density = gstar_density(tau_grid, psi, alpha, beta, Delta, rho_eta),
		label   = sprintf("psi = %s", psi)
	)
}) %>%
	mutate(label = factor(label, levels = sprintf("psi = %s", psi_vals)))

# ==============================================================================
# Figure
# ==============================================================================

all_labels <- c("No intervention", sprintf("psi = %s", psi_vals))
pal <- setNames(
	c("grey40", "#9ecae1", "#3182bd", "#08306b"),
	all_labels
)
lty <- setNames(
	c("dashed", "solid", "solid", "solid"),
	all_labels
)

fig_gi_trunc_det <- ggplot() +
	stat_function(fun = dgamma, args = list(shape = alpha, rate = beta),
	              aes(color = "No intervention", linetype = "No intervention"),
	              linewidth = 1.5, n = 1001, xlim = c(0, tau_max)) +
	geom_line(data = density_df,
	          aes(x = tau, y = density, color = label, linetype = label),
	          linewidth = 1.5, alpha=0.6) +
	scale_color_manual(values = pal, breaks = all_labels, name = NULL) +
	scale_linetype_manual(values = lty, breaks = all_labels, name = NULL) +
	scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
	coord_cartesian(xlim = c(0, tau_max)) +
	theme_classic(base_size = 16) +
	theme(
		legend.position         = "inside",
		legend.position.inside  = c(0.78, 0.72),
		legend.background       = element_rect(fill = "white", color = NA),
		legend.key.width        = unit(1.2, "cm")
	) +
	labs(
		x        = "Generation interval (days)",
		y        = "Density"
		# title    = "Generation interval distortion under deterministic isolation",
		# subtitle = sprintf(
		# 	"SARS-CoV-2 omicron (mean=%.1f d, SD=%.1f d)  |  Delta = %.1f d from peak  |  rho*eta = %.1f",
		# 	gi_mean, gi_sd, Delta, rho_eta
		# )
	)

print(fig_gi_trunc_det)

# ==============================================================================
# Save
# ==============================================================================

fig_name <- "isolation_gi_truncation_deterministic"

save_fig(fig_gi_trunc_det, fig_name, width = 7, height = 5)
cat(sprintf("Saved %s\n", fig_name))
