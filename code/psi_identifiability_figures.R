# ==============================================================================
# Figures for the psi identifiability analysis
# ==============================================================================
#
# Three figures:
#   A. Per-pathogen identifiability ribbon: simulated 95% CI envelope as a
#      function of psi_true, with the empirical posterior CI overlaid.
#   B. Ascertainment sensitivity forest plot: psi_hat distribution across
#      p_asc values, faceted by pathogen.
#   C. Incubation misspecification forest plot: psi_hat distribution across
#      0.5x, 1x, 2x literature incubation variance.
#
# Reads:
#   output/psi_empirical_results.csv
#   output/psi_identifiability_calibration.csv
#   output/psi_identifiability_ascertainment.csv
#   output/psi_identifiability_incubation.csv
#
# Writes:
#   figures/psi_identifiability_calibration.pdf (Figure A)
#   figures/psi_identifiability_ascertainment.pdf (Figure B)
#   figures/psi_identifiability_incubation.pdf (Figure C)
# ==============================================================================

cat("\n=== Psi identifiability figures ===\n")

# Theme used across all three figures, with explicit white background
theme_id <- function() {
	theme_classic(base_size = 18) +
		theme(
			plot.background  = element_rect(fill = "white", color = NA),
			panel.background = element_rect(fill = "white", color = NA),
			strip.background = element_rect(fill = "white", color = NA),
			strip.text       = element_text(face = "bold", size = 14, color = "black"),
			plot.title       = element_text(face = "bold", size = 18, color = "black"),
			plot.subtitle    = element_text(size = 8, color = "grey40"),
			axis.text        = element_text(color = "grey30"),
			axis.title       = element_text(color = "black"),
			legend.text      = element_text(color = "black"),
			legend.title     = element_text(color = "black"),
			panel.grid.minor = element_blank()
		)
}

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

emp <- read_csv("output/psi_empirical_results.csv", show_col_types = FALSE) %>%
	distinct(disease, post_mean, ci_lo, ci_hi, n_multi, n_si)

cal <- read_csv("output/psi_identifiability_calibration.csv", show_col_types = FALSE)
asc <- if (file.exists("output/psi_identifiability_ascertainment.csv")) {
	read_csv("output/psi_identifiability_ascertainment.csv", show_col_types = FALSE)
} else NULL
inc <- if (file.exists("output/psi_identifiability_incubation.csv")) {
	read_csv("output/psi_identifiability_incubation.csv", show_col_types = FALSE)
} else NULL

# Order pathogens by empirical posterior mean (matches the empirical figure)
pathogen_order <- emp %>% arrange(post_mean) %>% pull(disease)

# Helper: factor with consistent ordering
order_disease <- function(df) {
	df %>% mutate(disease = factor(disease, levels = pathogen_order))
}

emp <- order_disease(emp)
cal <- order_disease(cal)
if (!is.null(asc)) asc <- order_disease(asc)
if (!is.null(inc)) inc <- order_disease(inc)

# ------------------------------------------------------------------------------
# Figure A: identifiability ribbon
# ------------------------------------------------------------------------------

cat("  Figure A: per-pathogen identifiability ribbon\n")

# For each (pathogen, psi_true), compute the empirical distribution of
# (post_mean, ci_lo, ci_hi) across replicates.
ribbon_df <- cal %>%
	group_by(disease, psi_true) %>%
	summarise(
		mean_post_mean = mean(post_mean, na.rm = TRUE),
		q05_post_mean  = quantile(post_mean, 0.05, na.rm = TRUE),
		q95_post_mean  = quantile(post_mean, 0.95, na.rm = TRUE),
		mean_ci_lo     = mean(ci_lo, na.rm = TRUE),
		mean_ci_hi     = mean(ci_hi, na.rm = TRUE),
		.groups = "drop"
	)

fig_a <- ggplot(ribbon_df, aes(x = psi_true)) +
	# Identity reference (perfect recovery)
	geom_abline(slope = 1, intercept = 0, linetype = "dotted",
	            color = "grey60", linewidth = 0.4) +
	# Replicate-to-replicate spread of post_mean
	geom_ribbon(aes(ymin = q05_post_mean, ymax = q95_post_mean),
	            fill = "steelblue", alpha = 0.25) +
	# Mean posterior mean across replicates
	geom_line(aes(y = mean_post_mean), color = "steelblue", linewidth = 0.7) +
	# Typical 95% CI envelope (averaged across replicates)
	geom_line(aes(y = mean_ci_lo), color = "steelblue",
	          linewidth = 0.4, linetype = "dashed") +
	geom_line(aes(y = mean_ci_hi), color = "steelblue",
	          linewidth = 0.4, linetype = "dashed") +
	# Empirical CI overlay
	geom_rect(data = emp,
	          aes(xmin = -Inf, xmax = Inf, ymin = ci_lo, ymax = ci_hi),
	          fill = "firebrick", alpha = 0.15, inherit.aes = FALSE) +
	geom_hline(data = emp, aes(yintercept = post_mean),
	           color = "firebrick", linewidth = 0.5) +
	facet_wrap(~ disease, ncol = 5) +
	coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
	scale_x_continuous(breaks = c(0, 0.5, 1)) +
	scale_y_continuous(breaks = c(0, 0.5, 1)) +
	labs(
		x = expression("True" ~ psi),
		y = expression("Recovered" ~ psi ~ "(posterior mean and 95% CI)")
		# title = expression("Identifiability of " ~ psi ~ "given the empirical sample size"),
		# subtitle = "Blue: simulated mean (line) and 5-95% range (band) of recovered posterior mean. Dashed: mean 95% CI bounds. Red: empirical posterior mean (line) and 95% CI (band)."
	) +
	theme_id()

save_fig(fig_a, "psi_identifiability_calibration", width = 14, height = 6)

# ------------------------------------------------------------------------------
# Figure B: ascertainment sensitivity
# ------------------------------------------------------------------------------

if (!is.null(asc)) {
	cat("  Figure B: ascertainment sensitivity\n")

	# For each (pathogen, psi_true, p_asc), get mean recovered post_mean
	asc_summary <- asc %>%
		group_by(disease, psi_true, p_asc) %>%
		summarise(
			mean_post_mean = mean(post_mean, na.rm = TRUE),
			q05_post_mean  = quantile(post_mean, 0.05, na.rm = TRUE),
			q95_post_mean  = quantile(post_mean, 0.95, na.rm = TRUE),
			.groups = "drop"
		) %>%
		mutate(p_asc_label = sprintf("p_asc = %.1f", p_asc))

	fig_b <- ggplot(asc_summary, aes(x = psi_true, y = mean_post_mean,
	                                  color = factor(p_asc), group = p_asc)) +
		geom_abline(slope = 1, intercept = 0, linetype = "dotted",
		            color = "grey60", linewidth = 0.4) +
		geom_ribbon(aes(ymin = q05_post_mean, ymax = q95_post_mean,
		                fill = factor(p_asc)), alpha = 0.12, color = NA) +
		geom_line(linewidth = 0.7) +
		facet_wrap(~ disease, ncol = 5) +
		scale_color_viridis_d(name = "Ascertainment", end = 0.85,
		                      labels = c("0.3", "0.5", "0.7", "1.0")) +
		scale_fill_viridis_d(end = 0.85, guide = "none") +
		coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
		labs(
			x = expression("True" ~ psi),
			y = expression("Recovered" ~ psi ~ "(posterior mean)")
			# title = expression("Sensitivity of" ~ hat(psi) ~ "to ascertainment bias"),
			# subtitle = "For each p_asc value, total offspring R0 chosen so expected observed cluster size matches the empirical multi-cluster mean."
		) +
		theme_id() +
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

	save_fig(fig_b, "psi_identifiability_ascertainment", width = 14, height = 6)
}

# ------------------------------------------------------------------------------
# Figure C: incubation misspecification
# ------------------------------------------------------------------------------

if (!is.null(inc)) {
	cat("  Figure C: incubation misspecification\n")

	# For each (pathogen, incub_factor), distribution of post_mean.
	# Reference line is psi_hat under correctly specified incubation (factor=1),
	# so we isolate the misspecification effect from baseline posterior-mean
	# shrinkage toward the prior interior.
	inc_summary <- inc %>%
		group_by(disease, incub_factor) %>%
		summarise(
			mean_post_mean = mean(post_mean, na.rm = TRUE),
			q05_post_mean  = quantile(post_mean, 0.05, na.rm = TRUE),
			q95_post_mean  = quantile(post_mean, 0.95, na.rm = TRUE),
			.groups        = "drop"
		) %>%
		group_by(disease) %>%
		mutate(reference_psi = mean_post_mean[incub_factor == 1.0]) %>%
		ungroup()

	fig_c <- ggplot(inc_summary, aes(x = factor(incub_factor), y = mean_post_mean)) +
		geom_hline(aes(yintercept = reference_psi),
		           linetype = "dotted", color = "grey50", linewidth = 0.4) +
		geom_errorbar(aes(ymin = q05_post_mean, ymax = q95_post_mean),
		              width = 0.2, color = "steelblue", linewidth = 0.5) +
		geom_point(size = 2.5, color = "steelblue") +
		facet_wrap(~ disease, ncol = 5) +
		coord_cartesian(ylim = c(0, 1)) +
		labs(
			x = expression("Incubation variance multiplier (truth = 1" * x * ")"),
			y = expression(hat(psi) ~ "(posterior mean across reps)")
			# title = expression("Bias in" ~ hat(psi) ~ "under incubation-variance misspecification"),
			# subtitle = "Dotted line: posterior mean under correctly specified incubation (factor=1), the within-pathogen reference. Bars: 5-95% range across replicates."
		) +
		theme_id()

	save_fig(fig_c, "psi_identifiability_incubation", width = 14, height = 6)
}

cat("\n=== Psi identifiability figures complete ===\n")
