# ==============================================================================
# Detect-and-isolate: post-intervention growth rate vs psi
# ==============================================================================
#
# Analytic post-D&I growth rate r* as a function of psi, for the deterministic
# isolation case (mode-anchored offset iota, complete participation/effectiveness
# rho = eta = 1). Derived from the Euler-Lotka equation applied to the truncated
# convolution form of the post-intervention generation interval g^*(tau):
#
#   1 = R0 * (beta / (beta + r))^alpha * F_Gamma(psi*alpha, beta+r)(m_eps + iota)
#
# (The (1 - TE) normalization in g^* cancels with the renewal-kernel rescaling,
# leaving a single one-line numerical equation for r*.)
#
# Three lines are plotted per pathogen:
#   - r0       (dashed grey):   unperturbed Lotka-Euler growth rate
#   - r_naive  (orange dotted): r if intervention only reduces R0 by (1 - TE)
#                               (no GI shape change)
#   - r*       (blue solid):    actual post-D&I r (R0 reduction + GI truncation)
#
# The gap r* - r_naive quantifies the GI-truncation acceleration effect:
# survivors of the cutoff have systematically shorter generation intervals, so
# the surviving renewal kernel runs faster than the naive prediction suggests.
# ==============================================================================

library(tidyverse)
source("code/utils.R")
source("code/parameters.R")

# --- Parameters ---
iota_main       <- 2                    # deterministic mode-anchored isolation offset (days)
psi_vals_theory <- seq(0, 1, length.out = 101)

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Helpers
# ==============================================================================

# Mode of f_epsilon = Gamma(psi*alpha, beta).
m_eps_fn <- function(psi, alpha, beta) {
	shape <- psi * alpha
	if (shape >= 1) (shape - 1) / beta else 0
}

# Testing effectiveness under deterministic iota, rho = eta = 1:
#   TE = P(epsilon > m_eps + iota), epsilon ~ Gamma(psi*alpha, beta).
te_fn <- function(psi, alpha, beta, iota) {
	shape <- psi * alpha
	if (shape < 1e-9) return(0)  # spike limit: epsilon = 0 a.s.
	c_iso <- m_eps_fn(psi, alpha, beta) + iota
	pgamma(c_iso, shape = shape, rate = beta, lower.tail = FALSE)
}

# Post-D&I growth rate r* from the Euler-Lotka equation
#   1 = R0 (beta/(beta+r))^alpha F_Gamma(psi*alpha, beta+r)(m_eps + iota).
# The integrand is monotone decreasing in r (Laplace transform).
r_post_fn <- function(R0, alpha, beta, psi, iota) {
	shape <- psi * alpha
	# Spike limit: f_epsilon -> point mass at 0, no GI truncation, r* = r0.
	if (shape < 1e-6) return(beta * (R0^(1/alpha) - 1))

	c_iso <- m_eps_fn(psi, alpha, beta) + iota
	el_eq <- function(r) {
		R0 * (beta / (beta + r))^alpha *
		     pgamma(c_iso, shape = shape, rate = beta + r) - 1
	}
	# Search r > -beta (Laplace transform is finite there).
	uniroot(el_eq, interval = c(-beta * 0.9, 5), extendInt = "downX")$root
}

# Naive post-D&I growth rate: R0 reduction only, no GI shape change.
#   1 = R0(1 - TE) (beta/(beta+r))^alpha => r = beta((R0(1-TE))^(1/alpha) - 1).
# Returns 0 if subcritical (R_eff <= 1).
r_naive_fn <- function(R0, TE, alpha, beta) {
	R_eff <- R0 * (1 - TE)
	if (R_eff <= 1) return(0)
	beta * (R_eff^(1/alpha) - 1)
}

# ==============================================================================
# Compute curves
# ==============================================================================

growth_df <- map_dfr(seq_along(parslist), function(idx) {
	pars     <- parslist[[idx]]
	pathogen <- pars$pathogen
	alpha    <- pars$alpha
	beta     <- pars$beta
	R0       <- pars$R0
	r0       <- pars$r

	TE <- sapply(psi_vals_theory, function(psi)
		te_fn(psi, alpha, beta, iota_main))
	r_post  <- sapply(psi_vals_theory, function(psi)
		r_post_fn(R0, alpha, beta, psi, iota_main))
	r_naive <- sapply(seq_along(TE), function(i)
		r_naive_fn(R0, TE[i], alpha, beta))

	tibble(
		pathogen = pathogen,
		psi      = psi_vals_theory,
		r0       = r0,
		r_post   = r_post,
		r_naive  = r_naive,
		TE       = TE)
})

# Sanity print at canonical psi values
cat("\nGrowth rates at psi = 0, 0.5, 1:\n")
growth_df %>%
	filter(psi %in% c(0, 0.5, 1)) %>%
	mutate(across(c(r0, r_post, r_naive, TE), ~ round(.x, 4))) %>%
	select(pathogen, psi, r0, r_naive, r_post, TE) %>%
	print(n = Inf)

# ==============================================================================
# Per-pathogen figures
# ==============================================================================

pal <- c(
	"No intervention"               = "grey40",
	"R0 reduction only (no GI change)"   = "#FF9800",
	"Post-D&I (with GI truncation)"  = "#08306b"
)

lty <- c(
	"No intervention"               = "dashed",
	"R0 reduction only (no GI change)"   = "dotted",
	"Post-D&I (with GI truncation)"  = "solid"
)

scenario_levels <- names(pal)

for (pathogen_name in unique(growth_df$pathogen)) {

	pathogen_label <- switch(pathogen_name,
		influenza = "Influenza",
		omicron   = "SARS-CoV-2 Omicron",
		measles   = "Measles",
		pathogen_name
	)

	this_long <- growth_df %>%
		filter(pathogen == pathogen_name) %>%
		select(psi, r0, r_post, r_naive) %>%
		pivot_longer(cols = -psi, names_to = "scenario", values_to = "r") %>%
		mutate(scenario = recode(scenario,
			"r0"      = "No intervention",
			"r_naive" = "R0 reduction only (no GI change)",
			"r_post"  = "Post-D&I (with GI truncation)")) %>%
		mutate(scenario = factor(scenario, levels = scenario_levels))

	fig <- ggplot(this_long,
				  aes(x = psi, y = r, color = scenario, linetype = scenario)) +
		geom_line(linewidth = 1.5) +
		scale_color_manual(values = pal, name = NULL) +
		scale_linetype_manual(values = lty, name = NULL) +
		theme_classic(base_size = 16) +
		theme(legend.position = "bottom",
		      legend.direction = "vertical",
		      legend.key.width = unit(1.0, "cm")) +
		labs(
			x     = expression(psi),
			y     = expression("Growth rate " * italic(r) * " (1/day)"),
			title = pathogen_label)

	save_fig(fig, sprintf("fig_di_growth_rate_%s", pathogen_name),
	         width = 5.5, height = 5.5)
	cat(sprintf("  Saved fig_di_growth_rate_%s\n", pathogen_name))
}

# Save the underlying data for reuse / verification.
write_csv(growth_df,
          file.path("output",
                    sprintf("di_growth_rate_iota%g.csv", iota_main)))
cat(sprintf("\nSaved di_growth_rate_iota%g.csv\n", iota_main))
