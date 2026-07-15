#!/usr/bin/env Rscript
# =============================================================================
# psi_cluster_figures.R
#
# Builds Figure 2 panels A-C: for each of measles, MERS, and norovirus, the six
# largest transmission clusters from OutbreakTrees, drawn as points (serial
# intervals) along one horizontal line per infector, with:
#   (a) a red bar spanning each cluster's within-cluster range,
#   (b) the incubation-period distribution as light-blue shading, shifted so its
#       lower tail sits at the cluster's earliest case. If all of an infector's
#       secondary infections occurred at a single instant (psi = 0), onset times
#       would spread over exactly one incubation-period distribution, so this is
#       the serial-interval spread consistent with psi = 0,
#   (c) the generation-interval distribution g(tau) as a dashed reference curve.
#
# Reads output/psi_cluster_data.RDS (written by psi_inference.R), so it runs
# immediately after psi_inference.R without recomputing anything.
# =============================================================================

suppressMessages({library(dplyr); library(readr); library(purrr); library(tibble); library(ggplot2); library(patchwork)})

CLUST <- readRDS("output/psi_cluster_data.RDS")
adc   <- CLUST$all_disease_clusters
pp    <- CLUST$pathogen_params

# psi posteriors from psi_inference.R: full grid (panel D densities) + per-disease
# summaries (panel A-C subtitles and panel D points/CIs/labels)
psi_full <- read_csv("output/psi_empirical_results.csv", show_col_types = FALSE)
psi_res  <- psi_full %>% distinct(disease, post_mean, ci_lo, ci_hi, n_si, n_multi)

# viral-kinetics panel (E): plotting helpers + stored posterior (no MCMC re-run)
source("code/sc2_kinetics_utils.R")
kfeats <- readRDS("output/posterior_features_single_category.rds")

# --- config ------------------------------------------------------------------
DISEASES  <- c("Measles", "MERS", "Norovirus")
N_SHOW    <- 6          # clusters (infectors) per panel
ANCHOR_Q  <- 0.025      # incubation shading anchored so this quantile sits at earliest case
FIGDIR    <- "figures/"
BASE_SIZE <- 7          # base font (pt) -- tuned for a 6in-wide, 3-panel figure
PT_SIZE   <- 1.0        # infection marker size
RANGE_LW  <- 0.5        # cluster-range bar width
FIG_W     <- 6; FIG_H <- 2.7

# generation-interval params for the dashed curve. "data" uses the data-derived
# GI actually used for psi inference (pp[[disease]]$alpha_gi/beta_gi in the RDS),
# so the figure is consistent with the estimation; set to "literature" to instead
# use the published values in gi_lit below.
GI_SOURCE <- "data"
gi_lit <- list(  # mean, sd (days)
	Measles   = list(mean = 12.2, sd = 3.6),
	MERS      = list(mean = 6.8,  sd = 4.1),
	Norovirus = list(mean = 3.6,  sd = 2.0)
)
gi_params <- function(disease) {
	if (GI_SOURCE == "data") {
		list(alpha = pp[[disease]]$alpha_gi, beta = pp[[disease]]$beta_gi)
	} else {
		m <- gi_lit[[disease]]$mean; s <- gi_lit[[disease]]$sd
		list(alpha = m^2 / s^2, beta = m / s^2)
	}
}

# GI dashed-reference style:
#   "envelope" (default): per row, anchored at the cluster's earliest case, draw
#     the generation interval CONVOLVED with the incubation distribution -- the
#     within-infector serial-interval spread expected under psi = 1 -- at the same
#     scale as the incubation (psi = 0) shading, for a direct width comparison.
#   "single": one g(tau) curve spanning the panel height (original look).
#   "none":   omit the GI reference entirely (rely on incubation shading + points).
GI_STYLE <- "none"

# numerical convolution of the GI and incubation gammas (unit-peak-scaled) =
# the psi = 1 within-infector serial-interval spread
conv_env <- function(a_gi, b_gi, a_inc, b_inc, dx = 0.02) {
	xs <- seq(0, qgamma(0.9995, a_gi, b_gi) + qgamma(0.9995, a_inc, b_inc), dx)
	fc <- convolve(dgamma(xs, a_gi, b_gi), rev(dgamma(xs, a_inc, b_inc)), type = "open") * dx
	list(x = seq(0, by = dx, length.out = length(fc)), d = fc / max(fc))
}

COL_PT <- "grey15"; COL_RANGE <- "#d1362f"; COL_INC <- "#9ecae1"; COL_GI <- "grey30"

# --- one panel ---------------------------------------------------------------
build_panel <- function(disease, show_y = FALSE) {
	cl  <- adc[[disease]]$clusters
	ord <- order(sapply(cl, length), decreasing = TRUE)[seq_len(N_SHOW)]  # six largest clusters
	ord <- ord[order(sapply(ord, function(i) min(cl[[i]])))]              # then earliest-first (top)
	a_inc <- pp[[disease]]$a_obs; b_inc <- pp[[disease]]$b_obs
	gp <- gi_params(disease)

	# rows: earliest first infection at top -> latest at bottom
	rows <- map_dfr(seq_along(ord), function(r) {
		tibble(rank = r, y = N_SHOW - r + 1, si = cl[[ord[r]]])
	})
	ranges <- rows %>% group_by(rank, y) %>%
		summarise(lo = min(si), hi = max(si), .groups = "drop")

	xmax <- max(rows$si) * 1.05
	xpad <- max(0.15, 0.025 * xmax)   # left margin so day-0 points aren't clipped
	ybot <- 0.45; ytop <- N_SHOW + 0.7

	# stack overlapping (integer) points upward; step adapts so the tallest
	# stack fits within the row band
	pts <- rows %>% mutate(sib = round(si)) %>%
		group_by(rank, sib) %>% mutate(k = row_number() - 1L) %>% ungroup()
	step <- min(0.11, 0.60 / max(pts$k + 1L))
	pts  <- pts %>% mutate(yp = y + k * step)

	# (b) incubation shading, per row, shifted so ANCHOR_Q quantile aligns w/ earliest case
	tinc <- seq(0, qgamma(0.999, a_inc, b_inc), length.out = 200)
	dinc <- dgamma(tinc, a_inc, b_inc); dinc <- dinc / max(dinc)
	t_lo <- qgamma(ANCHOR_Q, a_inc, b_inc)
	inc_df <- map_dfr(seq_len(N_SHOW), function(r) {
		yb <- ranges$y[ranges$rank == r]; x0 <- ranges$lo[ranges$rank == r]
		tibble(rank = r, x = x0 + (tinc - t_lo), ymin = yb, ymax = yb + dinc * 0.62)
	}) %>% filter(x >= 0)

	# (c) generation-interval dashed reference (style-dependent)
	if (GI_STYLE == "none") {
		gi_geom <- NULL
	} else if (GI_STYLE == "envelope") {
		ce  <- conv_env(gp$alpha, gp$beta, a_inc, b_inc)
		clo <- ce$x[which(cumsum(ce$d) / sum(ce$d) >= ANCHOR_Q)[1]]   # anchor lower tail at earliest case
		gi_df <- map_dfr(seq_len(N_SHOW), function(r) {
			yb <- ranges$y[ranges$rank == r]; x0 <- ranges$lo[ranges$rank == r]
			tibble(rank = r, x = x0 + (ce$x - clo), y = yb + ce$d * 0.62)
		}) %>% filter(x >= 0)
		gi_geom <- geom_line(data = gi_df, aes(x, y, group = rank,
								 linetype = "Generation interval distribution"), colour = COL_GI, linewidth = 0.35)
	} else {
		tau <- seq(0, xmax, length.out = 300); dgi <- dgamma(tau, gp$alpha, gp$beta); dgi <- dgi / max(dgi)
		gi_df <- tibble(x = tau, y = ybot + dgi * (ytop - ybot) * 0.95)
		gi_geom <- geom_line(data = gi_df, aes(x, y,
								 linetype = "Generation interval distribution"), colour = COL_GI, linewidth = 0.4)
	}

	p <- ggplot() +
		geom_hline(data = ranges, aes(yintercept = y), colour = "grey85", linewidth = 0.3) +
		geom_ribbon(data = inc_df, aes(x = x, ymin = ymin, ymax = ymax, group = rank,
																	 fill = "Incubation period distribution"), alpha = 0.55) +
		gi_geom +
		geom_segment(data = ranges, aes(x = lo, xend = hi, y = y, yend = y, colour = "Cluster range"),
								 linewidth = RANGE_LW) +
		geom_point(data = pts, aes(x = si, y = yp, colour = "Infections"), shape = 18, size = PT_SIZE) +
		scale_colour_manual(NULL, values = c("Infections" = COL_PT, "Cluster range" = COL_RANGE),
												breaks = c("Infections", "Cluster range")) +
		scale_fill_manual(NULL, values = c("Incubation period distribution" = COL_INC)) +
		(if (GI_STYLE != "none")
			 scale_linetype_manual(NULL, values = c("Generation interval distribution" = "dashed"))) +
		guides(colour = guide_legend(order = 1, override.aes = list(
											shape = c(18, NA), linetype = c(NA, 1), linewidth = c(NA, RANGE_LW))),
					 linetype = if (GI_STYLE != "none") guide_legend(order = 2, override.aes = list(colour = COL_GI)),
					 fill = guide_legend(order = 3)) +
		scale_y_continuous(breaks = NULL) +
		coord_cartesian(xlim = c(-xpad, xmax), ylim = c(ybot, ytop), expand = FALSE) +
		labs(x = "Serial interval (days)", y = if (show_y) "Infector" else NULL,
				 title = disease, subtitle = psi_label(disease)) +
		theme_classic(base_size = BASE_SIZE) +
		theme(plot.title = element_text(face = "bold", size = BASE_SIZE + 1),
					plot.subtitle = element_text(size = BASE_SIZE, colour = "grey30"),
					axis.line.y = element_blank(), legend.position = "bottom",
					legend.box = "horizontal", legend.margin = margin(0, 0, 0, 0),
					legend.key.size = unit(9, "pt"),
					plot.margin = margin(t = 2, r = 4, b = 2, l = 16))  # left room for the flush tag
	p
}

# "psi = 0.03 (0.00-0.12)" subtitle from the stored posterior summary
psi_label <- function(disease) {
	r <- psi_res[psi_res$disease == disease, ]
	sprintf("ψ = %.2f (%.2f–%.2f)", r$post_mean, r$ci_lo, r$ci_hi)
}

# --- panel D: psi posterior across pathogens ---------------------------------
# ridgeline of the posterior density per pathogen (light grey), with black
# posterior-mean point + 95% CI whisker, ordered least- to most-bursty (lowest
# psi at top). n-labels sit in the right margin, clear of the data.
build_panel_D <- function() {
	lev <- psi_res %>% arrange(desc(post_mean)) %>% pull(disease)   # highest psi first -> bottom
	pf  <- psi_full %>% mutate(disease = factor(disease, levels = lev)) %>%
		group_by(disease) %>%
		mutate(ypos = as.integer(disease), dn = posterior / max(posterior) * 0.42) %>% ungroup()
	sm  <- psi_res %>% mutate(disease = factor(disease, levels = lev),
														ypos = as.integer(disease))
	txt_sz <- BASE_SIZE / 3.5   # n-labels a bit smaller than the pathogen names
	ggplot() +
		geom_ribbon(data = pf, aes(x = psi_grid, ymin = ypos - dn, ymax = ypos + dn, group = disease),
								fill = "grey80") +
		geom_errorbarh(data = sm, aes(y = ypos, xmin = ci_lo, xmax = ci_hi),
									 height = 0.22, linewidth = 0.4, colour = "black") +
		geom_point(data = sm, aes(x = post_mean, y = ypos), size = 1.3, colour = "black") +
		# n-labels: "n=x," right-aligned (commas line up) + "x clusters" left-aligned
		geom_text(data = sm, aes(x = 1.13, y = ypos, label = paste0("n=", n_si, ",")),
							hjust = 1, size = txt_sz, colour = "grey20") +
		geom_text(data = sm, aes(x = 1.16, y = ypos, label = paste0(n_multi, " clusters")),
							hjust = 0, size = txt_sz, colour = "grey20") +
		scale_y_continuous(breaks = seq_along(lev), labels = lev, expand = expansion(add = 0.7)) +
		scale_x_continuous(breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
		coord_cartesian(xlim = c(0, 1), clip = "off") +
		labs(x = "ψ (posterior mean with 95% CI)", y = NULL) +
		theme_classic(base_size = BASE_SIZE) +
		theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
					plot.margin = margin(5, 74, 5, 5))
}

# --- panel E: viral-kinetics infectious-window figure ------------------------
build_panel_E <- function() {
	plot_window_vs_ct_pub(kfeats, base_size = BASE_SIZE, gi_short = TRUE)
}

panels <- list(build_panel("Measles", show_y = TRUE),
							 build_panel("MERS"), build_panel("Norovirus"))
# NB: panel letters (A-E) are intentionally omitted -- added in post (Inkscape)
fig <- wrap_plots(panels, nrow = 1, guides = "collect") &
	theme(legend.position = "bottom")
ggsave(file.path(FIGDIR, "psi_cluster_panels.pdf"), fig, width = FIG_W, height = FIG_H, device = cairo_pdf)
cat(sprintf("Saved: %spsi_cluster_panels.pdf  (GI source: %s)\n", FIGDIR, GI_SOURCE))

# panel D standalone (for review; assembled with A-C + E below)
ggsave(file.path(FIGDIR, "psi_panelD.pdf"), build_panel_D(), width = 4, height = 3.4, device = cairo_pdf)
cat(sprintf("Saved: %spsi_panelD.pdf\n", FIGDIR))

# --- full Figure 2: A B C across the top, D + E across the bottom ------------
# legend collected within the top row (so it tucks under A-C, not the whole fig).
# free() detaches the top row from D's wide y-axis labels so A-C span the full width.
top    <- (panels[[1]] | panels[[2]] | panels[[3]]) + plot_layout(guides = "collect")
bottom <- (build_panel_D() | build_panel_E()) + plot_layout(widths = c(1.1, 1))
full   <- (free(top, side = "l") / bottom) + plot_layout(heights = c(0.95, 1.2)) &
	theme(legend.position = "bottom", legend.box.spacing = unit(3, "pt"))
ggsave(file.path(FIGDIR, "fig2_full.pdf"), full, width = 7, height = 5.5, device = cairo_pdf)
cat(sprintf("Saved: %sfig2_full.pdf\n", FIGDIR))
