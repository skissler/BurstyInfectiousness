library(tidyverse)
source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

# ==============================================================================
# Parameters
# ==============================================================================

n_draws  <- 100000  # for the theoretical distributions
psivals  <- c(0, 0.5, 1)

psi_colors <- setNames(c("red", "blue", "black")[seq_along(psivals)],
                       as.character(psivals))

establishment_threshold <- as.integer(round(establishment_threshold))

dir.create("output", showWarnings = FALSE, recursive = TRUE)

set.seed(42)

# ==============================================================================
# Helper functions
# ==============================================================================

#' Compute Var(W) for Gamma generation interval with punctuation psi
#'
#' Uses the rho_1/rho_2 parameterization from the writeup:
#'   z = R0^(1/alpha) - 1
#'   rho_k = (1+z)^k / (1+2z)
#'   Var(W) = (rho_1^alpha + rho_2^((1-psi)*alpha) - 1) / (1 - rho_1^alpha)
#'
#' @param R0    Basic reproduction number
#' @param alpha Shape parameter of Gamma generation interval
#' @param psi   Punctuation parameter (0 = spike, 1 = smooth)
compute_var_W <- function(R0, alpha, psi) {
	z <- R0^(1 / alpha) - 1
	rho_1 <- (1 + z) / (1 + 2*z)
	rho_2 <- (1 + z)^2 / (1 + 2*z)

	rho_1_a <- rho_1^alpha
	(rho_1_a + rho_2^((1 - psi) * alpha) - 1) / (1 - rho_1_a)
}

#' Fit 2-moment Gamma to W | survival and return shape/rate
#'
#' Unconditional: E[W] = 1, E[W^2] = 1 + Var(W)
#' Conditional on survival (W > 0): E_cond = 1/(1-q), Var_cond from above
fit_W_gamma <- function(R0, alpha, psi) {
	q <- extinction_prob(R0)
	var_W <- compute_var_W(R0, alpha, psi)
	EW2 <- 1 + var_W  # E[W^2] = Var(W) + (E[W])^2 = Var(W) + 1

	E_cond   <- 1 / (1 - q)
	E2_cond  <- EW2 / (1 - q)
	var_cond <- E2_cond - E_cond^2

	list(
		shape = E_cond^2 / var_cond,
		rate  = E_cond   / var_cond,
		E_cond = E_cond,
		var_cond = var_cond,
		var_W = var_W,
		q = q
	)
}

# ==============================================================================
# Loop over pathogens
# ==============================================================================

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
r        <- pars$r

cat(sprintf("\n===== %s: Tgen=%.2f, alpha=%.2f, R0=%.1f =====\n",
            pathogen, Tgen, alpha, R0))

# --------------------------------------------------------------------------
# Deterministic time to threshold from renewal equation
# --------------------------------------------------------------------------

# Extend tmax so slow pathogens still reach threshold within the integration
# horizon. 15*Tgen comfortably covers most epidemics at popsize/threshold scales
# used here.
ren_out <- renewal_epidemic(R0, alpha, Tgen, popsize, tmax = max(100, 15 * Tgen))

t_det <- approx(
	x = ren_out$cuminf * popsize,
	y = ren_out$t,
	xout = establishment_threshold,
	ties = mean
)$y

if (is.na(t_det)) {
	stop(sprintf(
		"%s: deterministic curve never reaches threshold %d by tmax (max cuminf*popsize = %.1f).",
		pathogen, establishment_threshold, max(ren_out$cuminf) * popsize
	))
}

# --------------------------------------------------------------------------
# Theoretical W distributions and survival curves for each psi
# --------------------------------------------------------------------------

theory_surv_list <- list()

for (psi in psivals) {

	fit <- fit_W_gamma(R0, alpha, psi)

	# Sample W | survival from fitted 2-moment Gamma
	W_draws <- rgamma(n_draws, shape = fit$shape, rate = fit$rate)
	t_stoch <- t_det - log(W_draws) / r

	theory_surv_list[[as.character(psi)]] <- tibble(
		psi = factor(psi), t_stoch = t_stoch)
}

theory_draws_df <- bind_rows(theory_surv_list)

# --------------------------------------------------------------------------
# Load empirical simulation data and build empirical survival curves
# --------------------------------------------------------------------------

cache <- load_cache(pathogen, nsim, popsize, psivals)

if (is.null(cache)) {
	cat(sprintf("  WARNING: cache not found for %s (n=%d, s=%d), skipping\n",
	    pathogen, popsize, nsim))
	next
}

# Extract time to establishment from summary (pre-computed in episims.R)
time_to_establishment <- cache$summary %>%
	filter(established == 1, !is.na(establishment_time)) %>%
	transmute(sim, psi, tinf = establishment_time)

# Build a t-grid that spans both empirical observations and the right tail of
# the theoretical W distribution, so theory curves aren't visually clipped.
t_max <- max(c(
	time_to_establishment$tinf,
	quantile(theory_draws_df$t_stoch, 0.999, na.rm = TRUE)
))
t_grid <- seq(0, t_max, length.out = 500)

empirical_surv <- time_to_establishment %>%
	group_by(psi) %>%
	summarise(
		t = list(t_grid),
		surv = list(sapply(t_grid, function(tt) mean(tinf > tt))),
		.groups = "drop") %>%
	unnest(cols = c(t, surv))

# Build theoretical survival curves on the same grid (2-moment gamma)
theory_surv <- theory_draws_df %>%
	group_by(psi) %>%
	summarise(
		t = list(t_grid),
		surv = list(sapply(t_grid, function(tt) mean(t_stoch > tt))),
		.groups = "drop") %>%
	unnest(cols = c(t, surv))

# --------------------------------------------------------------------------
# Overlay figures
# --------------------------------------------------------------------------

breakvals <- if (t_max > 120) seq(0, 365, by = 15) else seq(0, 365, by = 7)
surv_ylab <- paste0("Proportion not yet reaching ", establishment_threshold, " cases")

fig_survival <- ggplot() +
	geom_line(data = empirical_surv,
	          aes(x = t, y = surv, col = psi),
	          alpha = 0.6, linewidth = 0.8) +
	geom_line(data = theory_surv,
	          aes(x = t, y = surv, col = psi),
	          linetype = "dashed", linewidth = 0.8) +
	scale_x_continuous(breaks = breakvals) +
	scale_color_manual(values = psi_colors) +
	labs(x = "Days", y = surv_ylab, col = expression(psi),
	     title = pathogen_label) +
	theme_classic(base_size=14)

save_fig(fig_survival, paste0("fig_survival_", pathogen))

# Output the survival at weekly intervals:
empirical_surv_summary <- empirical_surv %>%
	mutate(week = floor(t / 7)) %>%
	group_by(psi, week) %>%
	filter(t == max(t))
cat(sprintf("  %s: Weekly survival values:\n", pathogen))
print(empirical_surv_summary, n = Inf)

# Output the time each curve takes to reach 5% survival:
empirical_surv_summary_pct <- empirical_surv %>%
	filter(surv < 0.05) %>%
	group_by(psi) %>%
	slice(1)
cat(sprintf("  %s: Time to reach 5pct survival:\n", pathogen))
print(empirical_surv_summary_pct, n = Inf)

# Find the date of the max deviation between psi=0 and psi=1 survival curves:
empirical_surv_summary_max <- empirical_surv %>%
	filter(psi %in% c("0", "1")) %>%
	mutate(day = floor(t)) %>%
	group_by(psi, day) %>%
	summarise(surv = max(surv), .groups = "drop") %>%
	pivot_wider(names_from = psi, values_from = surv) %>%
	mutate(diff = abs(`1` - `0`)) %>%
	filter(diff == max(diff))

cat(sprintf("  %s: Day of max difference:\n", pathogen))
print(empirical_surv_summary_max, n = Inf)

cat(sprintf("  %s: survival overlays saved.\n", pathogen))

# --------------------------------------------------------------------------
# W distribution: empirical histogram with moment-matched Gamma overlay
# --------------------------------------------------------------------------

# Empirical W | survival: W = exp(r * (t_det - t_emp))
W_empirical <- time_to_establishment %>%
	mutate(W = exp(r * (t_det - tinf)))

# Theoretical density functions for each psi
W_density_list <- list()
w_range <- sapply(psivals, function(p) {
	fit <- fit_W_gamma(R0, alpha, p)
	qgamma(c(0.001, 0.999), shape = fit$shape, rate = fit$rate)
})
w_grid <- seq(max(1e-4, min(w_range)), max(w_range), length.out = 500)

for (psi in psivals) {
	fit <- fit_W_gamma(R0, alpha, psi)
	dens <- dgamma(w_grid, shape = fit$shape, rate = fit$rate)
	W_density_list[[as.character(psi)]] <- tibble(
		psi = factor(psi), w = w_grid, density = dens)
}

W_density_df <- bind_rows(W_density_list)

fig_W_hist <- ggplot(W_empirical, aes(x = W)) +
	geom_histogram(aes(y = after_stat(density)), bins = 60,
	               fill = "white", col = "darkgrey") +
	geom_line(data = W_density_df,
	          aes(x = w, y = density),
	          col = "blue", linetype = "dashed", linewidth = 0.8) +
	facet_wrap(~psi, nrow = 1, scales = "free_y",
	           labeller = label_bquote(psi == .(as.numeric(as.character(psi))))) +
	labs(x = "W | survival", y = "Density", title = pathogen) +
	theme_classic()

save_fig(fig_W_hist, paste0("fig_W_hist_", pathogen))

cat(sprintf("  %s: W histogram saved.\n", pathogen))

# --------------------------------------------------------------------------
# Overlaid W distributions by psi
# --------------------------------------------------------------------------

fig_W_overlay <- ggplot() +
	geom_histogram(data = W_empirical,
	               aes(x = W, y = after_stat(density), fill = psi),
	               bins = 60, alpha = 0.25, position = "identity") +
	geom_line(data = W_density_df,
	          aes(x = w, y = density, col = psi),
	          linewidth = 0.8) +
	scale_fill_manual(values = psi_colors) +
	scale_color_manual(values = psi_colors) +
	coord_cartesian(xlim = c(0, quantile(W_empirical$W, 0.995))) +
	labs(x = "W | survival", y = "Density",
	     col = expression(psi), fill = expression(psi),
	     title = pathogen) +
	theme_classic()

save_fig(fig_W_overlay, paste0("fig_W_overlay_", pathogen))

cat(sprintf("  %s: W overlay saved.\n", pathogen))

# --------------------------------------------------------------------------
# Overlaid varsigma distributions by psi
# --------------------------------------------------------------------------

# Empirical varsigma = log(W) / r
varsigma_empirical <- W_empirical %>%
	mutate(varsigma = log(W) / r)

# Theoretical density of varsigma = log(W)/r where W ~ Gamma(k, lambda):
#   f_varsigma(s) = f_W(e^{rs}) * r * e^{rs}
varsigma_density_list <- list()
# Derive grid from theoretical distributions so tails aren't clipped by
# the empirical range of any particular psi value.
s_range <- sapply(psivals, function(p) {
	fit    <- fit_W_gamma(R0, alpha, p)
	w_edge <- qgamma(c(0.001, 0.999), shape = fit$shape, rate = fit$rate)
	log(w_edge) / r
})
s_grid <- seq(min(s_range), max(s_range), length.out = 500)

for (psi in psivals) {
	fit    <- fit_W_gamma(R0, alpha, psi)
	w_vals <- exp(r * s_grid)
	dens   <- dgamma(w_vals, shape = fit$shape, rate = fit$rate) * r * w_vals

	varsigma_density_list[[as.character(psi)]] <- tibble(
		psi = factor(psi), s = s_grid, density = dens)
}

varsigma_density_df <- bind_rows(varsigma_density_list)

# Compute empirical and theoretical means of varsigma by psi
varsigma_emp_means <- varsigma_empirical %>%
	group_by(psi) %>%
	summarise(mean_varsigma = mean(varsigma), .groups = "drop")

varsigma_theo_means <- tibble(
	psi = factor(psivals),
	mean_varsigma = sapply(psivals, function(p) {
		fit <- fit_W_gamma(R0, alpha, p)
		(digamma(fit$shape) - log(fit$rate)) / r
	}))

fig_varsigma_overlay <- ggplot() +
	geom_histogram(data = varsigma_empirical,
	               aes(x = varsigma, y = after_stat(density), fill = psi),
	               bins = 60, alpha = 0.25, position = "identity") +
	geom_line(data = varsigma_density_df,
	          aes(x = s, y = density, col = psi),
	          linewidth = 0.8) +
	geom_vline(data = varsigma_emp_means,
	           aes(xintercept = mean_varsigma, col = psi),
	           linetype = "dashed", linewidth = 0.6) +
	geom_vline(data = varsigma_theo_means,
	           aes(xintercept = mean_varsigma, col = psi),
	           linetype = "solid", linewidth = 0.6) +
	scale_fill_manual(values = psi_colors) +
	scale_color_manual(values = psi_colors) +
	labs(x = expression(varsigma ~~ "(days)"), y = "Density",
	     col = expression(psi), fill = expression(psi),
	     title = pathogen_label) +
	theme_classic(base_size=16)

save_fig(fig_varsigma_overlay, paste0("fig_varsigma_overlay_", pathogen))

cat(sprintf("  %s: varsigma overlay saved.\n", pathogen))

# --------------------------------------------------------------------------
# Combined cuminf + varsigma overlay (Omicron only)
# --------------------------------------------------------------------------

if (pathogen == "omicron") {

	plot_df <- cache$plot
	lastday <- ceiling(max(plot_df$tinf))

	# Transform varsigma density to time: t = t_det - s (mirror flip so
	# negative s -> times after t_det, matching the trajectory spray).
	# Scale so the tallest peak reaches 50% of popsize on the y-axis.
	density_height <- popsize * 0.50
	max_density    <- max(varsigma_density_df$density)

	time_density_df <- varsigma_density_df %>%
		mutate(
			t              = t_det - s,
			scaled_density = density / max_density * density_height
		)

	# Replicate all three distributions into every panel by crossing with
	# psi_panel (the faceting variable), keeping psi for fill/color.
	psi_panels <- levels(plot_df$psi)
	time_density_all <- time_density_df %>%
		tidyr::crossing(psi_panel = factor(psi_panels, levels = psi_panels))

	# Replicate the renewal curve across panels the same way.
	ren_all <- filter(ren_out, t <= lastday) %>%
		mutate(cuminf_n = cuminf * popsize) %>%
		tidyr::crossing(psi_panel = factor(psi_panels, levels = psi_panels))

	fig_cuminf_varsigma <- ggplot() +
		geom_line(data = plot_df %>% mutate(psi_panel = psi),
		          aes(x = tinf, y = cuminf, group = sim),
		          alpha = 0.2, col = "grey") +
		geom_line(data = ren_all,
		          aes(x = t, y = cuminf_n),
		          alpha = 0.8, linewidth = 1, col = "black") +
		geom_hline(yintercept = establishment_threshold,
		           linetype = "dashed", linewidth = 0.5, col = "black") +
		geom_ribbon(data = time_density_all,
		            aes(x = t, ymin = 0, ymax = scaled_density,
		                fill = psi, col = psi),
		            alpha = 0.3, linewidth = 0.5) +
		scale_fill_manual(values = psi_colors) +
		scale_color_manual(values = psi_colors) +
		facet_wrap(~psi_panel, nrow = 1,
		           labeller = label_bquote(psi == .(as.numeric(as.character(psi_panel))))) +
		labs(x = "Days", y = "Cumulative infections",
		     fill = expression(psi), col = expression(psi),
		     title = "Omicron") +
		theme_classic() +
		theme(legend.position = "bottom")

	save_fig(fig_cuminf_varsigma, "fig_cuminf_varsigma_omicron")
	cat("  omicron: cuminf + varsigma overlay saved.\n")
}

} # end pathogen loop
