# =============================================================================
# make_extended_data.R
# -----------------------------------------------------------------------------
# Build Nature-formatted Extended Data figures from the plot objects stashed by
# save_fig() during run_analysis.R (in figures/_objects/*.rds).
#
# WHY re-render instead of resize the exported PDFs:
#   Text size *relative to the plot* is fixed at ggsave() time. Resizing an
#   exported PDF/PNG scales the text along with everything else, so it can't fix
#   text that's too small. Re-rendering the plot object at the target size (with
#   an appropriate base font size) is the only thing that works.
#
# WORKFLOW:
#   1. Run the full pipeline once:      source("code/run_analysis.R")
#      (this now also writes figures/_objects/<name>.rds via save_fig)
#   2. Build the Extended Data figures:  Rscript code/make_extended_data.R
#      (or: source("code/make_extended_data.R") from the repo root)
#
# Originals in figures/ are never touched. Outputs go to figures/ExtendedData/,
# and a copy of each is placed in the manuscript figure directories (COPY_DIRS).
# Run from the repository root (same working directory as run_analysis.R).
# =============================================================================

library(tidyverse)
library(patchwork)

OBJ_DIR   <- file.path("figures", "_objects")
OUT_DIR   <- file.path("figures", "ExtendedData")
COPY_DIRS <- c(file.path("writeup", "v7_nature", "figures"))  # manuscript copies
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Restore pipeline globals (e.g. popsize) that some saved-plot aes() reference at
# render time -- these are lost when a ggplot object is read back from .rds.
source(file.path("code", "global_parameters.R"))

# --- Nature figure widths (mm) -----------------------------------------------
# Single column 89 mm; 1.5 column ~120 mm; double column 183 mm.
NATURE_WIDTH <- c(single = 89, onehalf = 120, double = 183)

# --- Shared Extended Data theme ----------------------------------------------
# Applied on TOP of each plot's existing theme (via `&`), so it only resets text
# sizes / margins -- it does not clobber backgrounds, colors, or facet strips.
# Nature wants ~5-9 pt text at final size and a sans-serif face (Arial/Helvetica).
theme_ed <- function(base_size = 7, axis_text = NULL, strip_size = NULL,
										 strip_clip = "inherit", margin_pt = c(2, 4, 2, 4),
										 legend_tight = TRUE, legend_position = "bottom",
										 base_family = "") {
	if (is.null(axis_text))  axis_text  <- base_size - 1
	if (is.null(strip_size)) strip_size <- base_size
	th <- theme(
		text            = element_text(size = base_size,  family = base_family),
		axis.title      = element_text(size = base_size),
		axis.text       = element_text(size = axis_text),
		plot.title      = element_text(size = base_size, hjust = 0.5),
		legend.title    = element_text(size = base_size),
		legend.text     = element_text(size = axis_text),
		legend.position = legend_position,   # Nature-friendly: legends below by default
		legend.box      = "vertical",        # stack multiple legends so none overflow the width
		legend.spacing.y = unit(3, "pt"),    # tighten the gap between stacked legends
		strip.text      = element_text(size = strip_size),
		plot.tag        = element_text(size = base_size + 2, face = "bold"),
		plot.margin     = margin(margin_pt[1], margin_pt[2], margin_pt[3], margin_pt[4], "pt")
	)
	# strip.clip = "off" lets long facet labels (e.g. "Pneumonic plague") overflow
	# rather than being clipped to the panel width.
	if (!identical(strip_clip, "inherit")) th <- th + theme(strip.clip = strip_clip)
	# Pull a moved (e.g. bottom) legend snug against its panel.
	if (legend_tight) th <- th + theme(legend.box.spacing = unit(2, "pt"),
																		 legend.margin = margin(1, 1, 1, 1, "pt"))
	th
}

# --- Helpers -----------------------------------------------------------------
load_fig <- function(name) {
	path <- file.path(OBJ_DIR, paste0(name, ".rds"))
	if (!file.exists(path)) {
		stop("No saved plot object for '", name, "' at ", path,
				 ".\n  Has run_analysis.R been run since save_fig() was updated? ",
				 "Mathematica-generated figures (e.g. 'contacts') won't have one.",
				 call. = FALSE)
	}
	readRDS(path)
}

# Shrink (or otherwise resize) the points in a saved ggplot's GeomPoint layer(s).
set_point_size <- function(p, size) {
	for (j in which(vapply(p$layers, function(L) inherits(L$geom, "GeomPoint"), logical(1))))
		p$layers[[j]]$aes_params$size <- size
	p
}

# Thicken lines in the given geom classes (default: GeomLine).
set_linewidth <- function(p, lw, geoms = "GeomLine") {
	for (j in which(vapply(p$layers, function(L) class(L$geom)[1] %in% geoms, logical(1))))
		p$layers[[j]]$aes_params$linewidth <- lw
	p
}

# Drop specified layers (by index) from the legend, so they don't spawn guides.
hide_legend_layers <- function(p, idx) {
	for (k in idx) p$layers[[k]]$show.legend <- FALSE
	p
}

# Remove all layers of a given geom class from a saved ggplot (e.g. "GeomVline").
drop_geom <- function(p, geom) {
	p$layers <- p$layers[!vapply(p$layers, function(L) inherits(L$geom, geom), logical(1))]
	p
}

# --- Builder -----------------------------------------------------------------
# Spec fields (all but id/figs optional):
#   id         : output basename (e.g. "ED_psi_identifiability")
#   figs       : character vector of source figure names (1 = single; >1 = grid)
#   prep       : optional function(list_of_plots) -> list, for per-panel edits
#                (relabel axes, move legends, resize points, ...)
#   layout     : list of wrap_plots() args, e.g. list(ncol = 1, heights = c(1,1.1,1))
#   tags       : tag_levels for panel letters (e.g. "A"); NULL for none
#   width      : "single"|"onehalf"|"double" (mm), or a number (units below)
#   height     : height (units below)
#   units      : "mm" (default) or "in"
#   base_size  : base font size in pt at final size (default 7)
#   axis_text  : tick-label size (default base_size - 1)
#   strip_size : facet-strip size (default base_size)
#   strip_clip : "off" to stop clipping long strip labels
#   margin_pt  : c(t,r,b,l) plot margin in pt (default c(2,4,2,4))
#   rasterize  : TRUE to flatten to 300-dpi raster PDF (needs ImageMagick)
build_ed_figure <- function(spec) {
	plots <- lapply(spec$figs, load_fig)
	if (!is.null(spec$prep)) plots <- spec$prep(plots)

	fig <- do.call(wrap_plots, c(plots, spec$layout %||% list()))
	if (!is.null(spec$tags)) fig <- fig + plot_annotation(tag_levels = spec$tags)

	fig <- fig & theme_ed(
		base_size       = spec$base_size %||% 7,
		axis_text       = spec$axis_text,
		strip_size      = spec$strip_size,
		strip_clip      = spec$strip_clip %||% "inherit",
		margin_pt       = spec$margin_pt %||% c(2, 4, 2, 4),
		legend_position = spec$legend_position %||% "bottom"
	)

	units    <- spec$units %||% "mm"
	width    <- if (is.numeric(spec$width)) spec$width else NATURE_WIDTH[[spec$width]]
	out_pdf  <- file.path(OUT_DIR, paste0(spec$id, ".pdf"))
	ggsave(out_pdf, fig, width = width, height = spec$height, units = units)

	if (isTRUE(spec$rasterize) && nzchar(Sys.which("magick"))) {
		system2("magick", c("-density", "300", "-background", "white",
												shQuote(out_pdf), "-flatten", "-compress", "jpeg",
												"-quality", "90", shQuote(out_pdf)))
	}

	for (d in COPY_DIRS)
		if (dir.exists(d)) file.copy(out_pdf, file.path(d, basename(out_pdf)), overwrite = TRUE)

	message(sprintf("  %-28s %4.1f x %4.1f %s  (%d panel%s)  -> %s + copies",
									spec$id, width, spec$height, units,
									length(plots), if (length(plots) > 1) "s" else "", OUT_DIR))
	invisible(out_pdf)
}

# =============================================================================
# EXTENDED DATA SPEC  --  one completed, tuned entry so far; add the rest of the
# triage as they're dialed in. (ED5/ED9 etc. templates are in the block below.)
# =============================================================================
ed_spec <- list(

	# Identifiability & robustness of psi estimation:
	#   (A) calibration/recovery, (B) ascertainment, (C) incubation misspecification.
	list(
		id   = "ED_psi_identifiability",
		figs = c("psi_identifiability_calibration",
						 "psi_identifiability_ascertainment",
						 "psi_identifiability_incubation"),
		prep = function(p) {
			# (A) break the long y-axis label so the parenthetical drops to a new line
			p[[1]] <- p[[1]] + labs(y = expression(atop("Recovered" ~ psi,
																									"(posterior mean and 95% CI)")))
			# (B) legend below the panel (also aligns B's columns with A and C)
			p[[2]] <- p[[2]] + theme(legend.position = "bottom")
			# (C) shrink points so the error-bar whiskers stay visible
			p[[3]] <- set_point_size(p[[3]], 1.0)
			p
		},
		layout     = list(ncol = 1, heights = c(1, 1.1, 1)),  # give B room for its legend
		tags       = "A",
		width      = 6, height = 7.6, units = "in",
		base_size  = 9, axis_text = 8, strip_size = 7, strip_clip = "off",
		margin_pt  = c(2, 4, 2, 4)
	),

	# Overdispersion heatmaps (coincidence superspreading):
	#   Row 1 (A-C): pathogens under the stochastic contact model (Type-I Gamma burst).
	#   Row 2 (D-F): omicron, burst models, PERIODIC contact model.
	#   Row 3 (G-I): omicron, burst models, STOCHASTIC contact model.
	# Shared k colorbar + shared extinction-probability size legend (guides collected).
	list(
		id   = "ED_overdispersion_heatmaps",
		figs = c("fig_heatmap_gammapoisson_influenza",
						 "fig_heatmap_gammapoisson_omicron",
						 "fig_heatmap_gammapoisson_measles",
						 "sensitivity_od_periodic_gamma_omicron",
						 "sensitivity_od_periodic_lognormal_omicron",
						 "sensitivity_od_periodic_gamma2_omicron",
						 "sensitivity_od_gp_gamma_omicron",
						 "sensitivity_od_gp_lognormal_omicron",
						 "sensitivity_od_gp_gamma2_omicron"),
		prep = function(p) {
			# Row 1: common [0,1] limits so the 3 extinction-prob legends collapse to one;
			# max_size 7 (was 10) keeps the largest discs from dominating panel A.
			for (i in 1:3)
				p[[i]] <- p[[i]] +
					scale_size_area(max_size = 7, name = "Extinction prob.", limits = c(0, 1))
			# Rows 2-3: shorten the long titles (the "(periodic/stochastic contacts)" part
			# overflowed); the caption carries the row meaning.
			bt <- c("Gamma (Type I)", "Log-normal", "Gamma (Type II)")
			for (i in 4:6) p[[i]] <- p[[i]] + ggtitle(bt[i - 3])
			for (i in 7:9) p[[i]] <- p[[i]] + ggtitle(bt[i - 6])
			p
		},
		layout    = list(ncol = 3, guides = "collect"),
		tags      = "A",
		width     = 6, height = 7.5, units = "in",   # taller -> panels closer to square
		base_size = 8, axis_text = 7,
		margin_pt = c(2, 3, 2, 3)
	),

	# Epidemic latent period:
	#   Row 1 (A-C): time-shift varsigma histograms; Row 2 (D-F): establishment
	#   survival curves; Row 3 (G, full width): omicron cumulative-infection
	#   trajectories across burst models x psi. Rasterized (trajectory-heavy).
	#   Single psi legend: only the survival panels emit it (identical -> merge);
	#   the histograms' extra fill scale would otherwise spawn a second legend.
	list(
		id   = "ED_latent_period",
		figs = c("fig_varsigma_overlay_influenza","fig_varsigma_overlay_omicron","fig_varsigma_overlay_measles",
						 "fig_survival_influenza","fig_survival_omicron","fig_survival_measles",
						 "sensitivity_episims_cuminf_omicron"),
		prep = function(p) {
			psi_scale <- scale_colour_manual(
				name = expression(psi), values = c(`0` = "red", `0.5` = "blue", `1` = "black"),
				guide = guide_legend(override.aes = list(alpha = 1, linewidth = 1.0, linetype = 1)))
			for (i in 1:6) p[[i]] <- set_linewidth(p[[i]], 0.4)          # thin curves (half width)
			for (i in 1:3) p[[i]] <- drop_geom(p[[i]], "GeomVline")      # remove cluttered mean bars
			for (i in 1:3) p[[i]] <- hide_legend_layers(p[[i]], seq_along(p[[i]]$layers))  # varsigma: no legend
			for (i in 4:6) p[[i]] <- hide_legend_layers(p[[i]], 2)       # survival: drop dashed key
			for (i in 1:6) p[[i]] <- p[[i]] + psi_scale                  # unified psi colour scale
			# Survival panels: zoom x to the informative part, and wrap the long y-label.
			for (i in 4:6) p[[i]] <- p[[i]] + labs(y = "Proportion not yet\nreaching 500 cases")
			p[[4]] <- p[[4]] + coord_cartesian(xlim = c(7, 49))          # influenza
			p[[5]] <- p[[5]] + coord_cartesian(xlim = c(7, 49))          # omicron (drop flat tail)
			p[[6]] <- p[[6]] + coord_cartesian(xlim = c(14, NA))         # measles (drop long lead-in)
			# Episims (G): drop the ugly boxes around the facet-strip labels.
			p[[7]] <- p[[7]] + theme(strip.background = element_blank())
			p
		},
		layout     = list(design = "ABC\nDEF\nGGG", heights = c(0.7, 0.7, 1.4), guides = "collect"),
		tags       = "A",
		width      = 6, height = 6.9, units = "in",   # flatter top rows
		base_size  = 8, axis_text = 7, strip_size = 6,
		strip_clip = "off",   # episims' rotated right-side strip labels overflow otherwise
		margin_pt  = c(2, 3, 2, 3),
		rasterize  = TRUE   # episims trajectories -> flatten to 300 dpi to keep size sane
	),

	# Detect-and-isolate interventions (15 panels, A-O):
	#   Row 1 (A-C): symptom-based TE, pathogens.   Row 2 (D-F): screening TE, pathogens.
	#   Row 3 (G-I): post-D&I growth rate, pathogens (own legend).
	#   Row 4 (J-L): symptom-based TE across burst models.  Row 5 (M-O): screening TE, burst models.
	#   Two collected legends: shared psi (rows 1/2/4/5) and intervention scenario (row 3).
	list(
		id   = "ED_di_interventions",
		figs = c("fig_te_symp_influenza","fig_te_symp_omicron","fig_te_symp_measles",
						 "fig_te_testing_delay_influenza","fig_te_testing_delay_omicron","fig_te_testing_delay_measles",
						 "fig_di_growth_rate_influenza","fig_di_growth_rate_omicron","fig_di_growth_rate_measles",
						 "sensitivity_te_symp_gamma_omicron","sensitivity_te_symp_lognormal_omicron","sensitivity_te_symp_gamma2_omicron",
						 "sensitivity_te_delay_gamma_omicron","sensitivity_te_delay_lognormal_omicron","sensitivity_te_delay_gamma2_omicron"),
		prep = function(p) {
			symp_x <- "Mean symptom onset\n(days rel. to peak)"
			for (i in c(1, 2, 3, 10, 11, 12)) p[[i]] <- p[[i]] + labs(x = symp_x)  # shorten long x-label
			for (i in c(1:6, 10:15))          p[[i]] <- p[[i]] + labs(y = "Test effectiveness")  # unify y-label
			# Row 3: thin lines, tighter dash/dot spacing, R_0 subscript in the labels,
			# and the scenario legend on a single row to save vertical space.
			scen <- c("No intervention", "R0 reduction only (no GI change)", "Post-D&I (with GI truncation)")
			pal  <- setNames(c("grey40", "#FF9800", "#08306b"), scen)
			lty  <- setNames(c("22", "12", "solid"), scen)
			labs <- c(expression("No intervention"),
								expression(R[0] ~ "reduction only (no GI change)"),
								expression("Post-D&I (with GI truncation)"))
			for (i in 7:9)
				p[[i]] <- set_linewidth(p[[i]], 0.8) +
					scale_colour_manual(values = pal, breaks = scen, labels = labs, name = NULL,
															guide = guide_legend(nrow = 1)) +
					scale_linetype_manual(values = lty, breaks = scen, labels = labs, name = NULL,
																guide = guide_legend(nrow = 1))
			p
		},
		layout    = list(ncol = 3, guides = "collect", heights = c(0.85, 0.85, 1, 0.85, 0.85)),
		tags      = "A",
		width     = 6, height = 9, units = "in",
		base_size = 8, axis_text = 7,
		margin_pt = c(2, 3, 2, 3)
	),

	# Gathering-size restrictions (9 panels, A-I; pathogens in columns):
	#   Row 1 (A-C): overdispersion OD = 1/k under Constant/Restricted/Unrestricted contacts.
	#   Row 2 (D-F): OD reduction from restrictions (single series, no legend).
	#   Row 3 (G-I): establishment probability. Rows 1 & 3 share one contact-scenario legend.
	list(
		id   = "ED_gathering_restrictions",
		figs = c("fig_gathering_od_influenza","fig_gathering_od_omicron","fig_gathering_od_measles",
						 "fig_gathering_dod_influenza","fig_gathering_dod_omicron","fig_gathering_dod_measles",
						 "fig_gathering_pest_theory_influenza","fig_gathering_pest_theory_omicron","fig_gathering_pest_theory_measles"),
		prep = function(p) { for (i in seq_along(p)) p[[i]] <- set_point_size(p[[i]], 1.0); p },  # shrink dots off the whiskers
		layout    = list(ncol = 3, guides = "collect"),
		tags      = "A",
		width     = 6, height = 6.5, units = "in",
		base_size = 8, axis_text = 7,
		margin_pt = c(2, 3, 2, 3)
	),

	# Growth-rate estimation under burstiness (9 panels, A-I):
	#   Rows 1-3 (A, B, C): full-width daily-incidence trajectories (faceted by psi)
	#     for influenza / omicron / measles.
	#   Row 4 (D-F): empirical growth-rate distributions by pathogen.
	#   Row 5 (G-I): the same across burst models (omicron). Rows 4-5 share the psi legend.
	list(
		id   = "ED_growthrate_estimation",
		figs = c("fig_growthrate_infpop_lines_influenza",
						 "fig_growthrate_infpop_lines_omicron",
						 "fig_growthrate_infpop_lines_measles",
						 "fig_growthrate_infpop_hists_overlay_influenza","fig_growthrate_infpop_hists_overlay_omicron","fig_growthrate_infpop_hists_overlay_measles",
						 "sensitivity_gr_overlay_gamma_omicron","sensitivity_gr_overlay_lognormal_omicron","sensitivity_gr_overlay_gamma2_omicron"),
		prep = function(p) {
			for (i in 1:3) p[[i]] <- set_linewidth(p[[i]], 0.4, geoms = "GeomAbline")                    # halve dashed fit lines
			for (i in 4:9) {
				p[[i]] <- set_linewidth(p[[i]], 0.15, geoms = "GeomBar")                                   # thin histogram bar borders
				p[[i]] <- set_linewidth(p[[i]], 0.6,  geoms = c("GeomDensity", "GeomVline"))               # thicker density outlines / mean markers
			}
			p
		},
		layout    = list(design = "AAA\nBBB\nCCC\nDEF\nGHI", heights = c(0.7, 0.7, 0.7, 1, 1), guides = "collect"),
		tags      = "A",
		width     = 6, height = 8.5, units = "in",
		base_size = 8, axis_text = 7,
		margin_pt = c(2, 3, 2, 3)
	),

	# Generation-interval inference (3 stacked full-width panels, A-C):
	#   A: marginal posteriors for alpha (GI shape); B: for beta (GI rate)
	#   -- both 3 pathogens x 3 psi, squashed. C: empirical 95% CI coverage.
	#   Single Parameter (alpha/beta) legend from C; strip boxes removed.
	list(
		id   = "ED_gi_inference",
		figs = c("g_identifiability_alpha_posteriors",
						 "g_identifiability_beta_posteriors",
						 "g_identifiability_coverage"),
		prep = function(p) {
			for (i in seq_along(p)) p[[i]] <- p[[i]] + theme(strip.background = element_blank())
			p[[3]] <- set_point_size(p[[3]], 1.0)   # shrink coverage points off the whiskers
			p
		},
		layout    = list(ncol = 1, heights = c(1, 1, 0.75), guides = "collect"),
		tags      = "A",
		width     = 6, height = 6.5, units = "in",
		base_size = 8, axis_text = 7,
		margin_pt = c(2, 3, 2, 3)
	)

	# ---- TEMPLATES for the remaining triage (activate as tuned) ----------------
	# , list(id = "ED_testing_effectiveness",
	#        figs = c("fig_te_symp", "fig_te_testing"),
	#        layout = list(nrow = 2), tags = "A",
	#        width = "double", height = 140, base_size = 7)
	# , list(id = "ED_gi_inference",
	#        figs = c("g_identifiability_coverage",
	#                 "g_identifiability_alpha_posteriors",
	#                 "g_identifiability_beta_posteriors"),
	#        layout = list(nrow = 1), tags = "A",
	#        width = "double", height = 70, base_size = 7)
	# ED1  burst-model validation  : marginal_gi_hist_* + marginal_gi_qq_* (merge)
	# ED3  overdispersion robustness: heatmap_gammapoisson + alternative_burst_od
	# ED4  epidemic latent period   : varsigma_overlay + survival_full
	# ED6  D&I secondary effects    : gi_truncation_deterministic + di_growth_rate
	# ED7  gathering-size           : gathering_od + gathering_dod + gathering_pest
	# ED8  growth-rate estimation   : growthrate_hists + growthrate_daily_counts
	# ED10 alt-burst dynamics       : alt-burst cuminf/daily (set rasterize = TRUE)
)

# --- Build all ---------------------------------------------------------------
message("Building Extended Data figures -> ", OUT_DIR)
invisible(lapply(ed_spec, build_ed_figure))
message("Done.")
