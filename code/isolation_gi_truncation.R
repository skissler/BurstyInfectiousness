library(tidyverse)

source("code/utils.R")
source("code/parameters.R")

set.seed(42)

# ==============================================================================
# Define parameters
# ==============================================================================

# Pull omicron GI parameters from parslist
omicron  <- Filter(function(p) p$pathogen == "omicron", parslist)[[1]]
alpha    <- omicron$alpha
beta     <- omicron$beta
gi_mean  <- omicron$Tgen
gi_sd    <- sqrt(omicron$Tvar)

sigma_sym  <- 0.5          # SD of symptom onset around infectiousness peak (days)

psi_vals <- c(0, 0.5, 1)
# psi_vals <- c(0.01, 0.5, 0.99)
# psi_vals <- c(0.1, 0.5, 0.9)

N        <- 2e6          # Monte Carlo draws per psi
tau_max  <- 14           # x-axis limit (days)

# ==============================================================================
# Define simulation function
# ==============================================================================

# tau = l + eps, where
#   l   ~ Gamma((1-psi)*alpha, rate=beta)   [latent period, shared]
#   eps ~ Gamma(psi*alpha,     rate=beta)   [transmission timing]
#
# Isolation cutoff (in eps units): c = m_psi + D0, D0 ~ N(0, sigma_sym^2)
# Transmission survives iff: eps < c  (and c > 0)

sim_gi_post_isolation <- function(psi, alpha, beta, sigma_sym, N) {
	alpha_l <- (1 - psi) * alpha
	alpha_e <- psi       * alpha
	m_psi   <- if (alpha_e >= 1) (alpha_e - 1) / beta else 0

	l        <- if (alpha_l > 1e-9) rgamma(N, shape = alpha_l, rate = beta) else rep(0, N)
	eps      <- if (alpha_e > 1e-9) rgamma(N, shape = alpha_e, rate = beta) else rep(0, N)
	c_iso    <- m_psi + rnorm(N, mean = 0, sd = sigma_sym)
	survives <- eps < c_iso & c_iso > 0

	(l + eps)[survives]
}

# ==============================================================================
# Sample surviving generation intervals for each psi
# ==============================================================================

sim_df <- map_dfr(psi_vals, function(psi) {
	tibble(
		tau   = sim_gi_post_isolation(psi, alpha, beta, sigma_sym, N),
		label = paste0("psi = ", psi)
	)
}) %>%
	mutate(label = factor(label, levels = paste0("psi = ", psi_vals)))

# ==============================================================================
# Figure
# ==============================================================================

all_labels <- c("No intervention", paste0("psi = ", psi_vals))
pal <- setNames(
	c("grey40", "#9ecae1", "#3182bd", "#08306b"),
	all_labels
)
lty <- setNames(
	c("dashed", "solid", "solid", "solid"),
	all_labels
)

fig_gi_trunc <- ggplot() +
	# Analytic "no intervention" GI density (Gamma)
	stat_function(fun = dgamma, args = list(shape = alpha, rate = beta),
								aes(color = "No intervention", linetype = "No intervention"),
								linewidth = 0.9, n = 1001, xlim = c(0, tau_max)) +
	# Post-isolation GI densities (KDE of surviving samples)
	geom_density(data = sim_df,
							 aes(x = tau, color = label, linetype = label),
							 bw = "sj", adjust = 2, n = 1024,
							 linewidth = 0.9) +
	scale_color_manual(values = pal, breaks = all_labels, name = NULL) +
	scale_linetype_manual(values = lty, breaks = all_labels, name = NULL) +
	scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
	coord_cartesian(xlim = c(0, tau_max)) +
	theme_classic(base_size = 13) +
	theme(
		legend.position        = "inside",
		legend.position.inside  = c(0.78, 0.72),
		legend.background = element_rect(fill = "white", color = NA),
		legend.key.width  = unit(1.2, "cm")
	) +
	labs(
		x        = "Generation interval (days)",
		y        = "Density",
		title    = "Effect of detect-and-isolate on generation interval distribution",
		subtitle = sprintf(
			"SARS-CoV-2 omicron (mean=%.1f d, SD=%.1f d)  |  symptom onset SD = %.0f d from peak",
			gi_mean, gi_sd, sigma_sym
		)
	)

print(fig_gi_trunc)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

fig_name <- "isolation_gi_truncation"

save_fig(fig_gi_trunc, fig_name, width = 7, height = 5)
cat(sprintf("Saved %s\n", fig_name))
