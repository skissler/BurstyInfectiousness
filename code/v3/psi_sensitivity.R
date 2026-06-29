# ==============================================================================
# psi_sensitivity.R — incubation period SD sensitivity analysis
#
# Requires the psi_empirical_v2.R environment to already be loaded:
#   pathogen_params, all_disease_clusters, psi_grid, min_multi_clusters,
#   precompute_engines_litgi, compute_psi_posterior_litgi, save_fig
#
# Usage:
#   source("code/psi_empirical_v2.R")
#   source("code/psi_sensitivity.R")
# ==============================================================================

library(tidyverse)

stopifnot(
  exists("pathogen_params"),
  exists("all_disease_clusters"),
  exists("psi_grid"),
  exists("min_multi_clusters"),
  exists("precompute_engines_litgi"),
  exists("compute_psi_posterior_litgi")
)

# ------------------------------------------------------------------------------
# Helper: run full psi estimation for a given pathogen_params list
# ------------------------------------------------------------------------------

run_psi_estimation <- function(pp_list, disease_clusters, psi_grid, min_multi) {
	results <- list()
	for (disease in names(disease_clusters)) {
		dc      <- disease_clusters[[disease]]
		multi   <- dc$clusters[sapply(dc$clusters, length) >= 2]
		if (length(multi) < min_multi) next
		pp      <- pp_list[[disease]]
		engines <- tryCatch(
			precompute_engines_litgi(psi_grid, pp$alpha_gi, pp$beta_gi, pp$a_obs, pp$b_obs),
			error = function(e) NULL
		)
		if (is.null(engines)) next
		post     <- compute_psi_posterior_litgi(multi, engines, psi_grid)
		post_cdf <- cumsum(post)
		results[[disease]] <- tibble(
			disease   = disease,
			post_mode = psi_grid[which.max(post)],
			post_mean = sum(psi_grid * post),
			ci_lo     = psi_grid[which.min(abs(post_cdf - 0.025))],
			ci_hi     = psi_grid[which.min(abs(post_cdf - 0.975))]
		)
	}
	bind_rows(results)
}

# ------------------------------------------------------------------------------
# Run across inc_sd scaling factors
# ------------------------------------------------------------------------------

inc_sd_scales <- c(0.5, 0.75, 1.0, 1.25, 1.5)

cat("\n--- Incubation SD sensitivity analysis ---\n")

sensitivity_results <- map_dfr(inc_sd_scales, function(scale) {
	cat(sprintf("  inc_sd scale = %.2f\n", scale))
	pp_scaled <- lapply(pathogen_params, function(pp) {
		inc_sd_s  <- pp$inc_sd * scale
		inc_var_s <- inc_sd_s^2
		pp$inc_sd <- inc_sd_s
		pp$a_obs  <- pp$inc_mean^2 / inc_var_s
		pp$b_obs  <- pp$inc_mean   / inc_var_s
		pp
	})
	run_psi_estimation(pp_scaled, all_disease_clusters, psi_grid, min_multi_clusters) %>%
		mutate(inc_sd_scale = scale)
})

# Order diseases by psi mode at baseline (scale = 1.0)
disease_order_sens <- sensitivity_results %>%
	filter(inc_sd_scale == 1.0) %>%
	arrange(post_mode) %>%
	pull(disease)
sensitivity_results <- sensitivity_results %>%
	mutate(disease = factor(disease, levels = disease_order_sens))

# ------------------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------------------

fig_sensitivity <- ggplot(
	sensitivity_results,
	aes(x = inc_sd_scale, y = post_mode, color = disease, group = disease)
) +
	geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = disease),
	            alpha = 0.10, color = NA) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	scale_x_continuous(breaks = inc_sd_scales,
	                   labels = sprintf("%.2fx", inc_sd_scales)) +
	scale_y_continuous(limits = c(0, 1)) +
	geom_vline(xintercept = 1.0, linetype = "dashed",
	           color = "grey50", linewidth = 0.4) +
	facet_wrap(~ disease, nrow = 2) +
	theme_classic(base_size = 11) +
	theme(legend.position = "none",
	      strip.text     = element_text(size = 9),
	      axis.text.x    = element_text(angle = 30, hjust = 1, size = 8)) +
	labs(
		x        = "Incubation period SD scaling factor",
		y        = "Posterior ψ mode",
		title    = "Sensitivity of ψ estimates to incubation period SD",
		subtitle = "Shaded band = 95% CI; vertical line at 1.0 = verified published values"
	)

print(fig_sensitivity)

if (exists("save_fig")) {
	save_fig(fig_sensitivity, "psi_sensitivity_incsd", width = 12, height = 6)
	cat("Saved psi_sensitivity_incsd figure.\n")
}
