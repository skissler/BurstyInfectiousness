#!/usr/bin/env Rscript
# =============================================================================
# fig5_figures.R
#
# Assembles Figure 5 (detect-and-isolate interventions) panels B-E from the
# saved ggplot objects (figures/_objects/*.rds, written by save_fig in the
# isolation_* scripts) -- re-rendered at a common size/style, NOT re-computed.
# Panel A is a Mathematica schematic added in post; a blank slot is reserved.
#
# Layout: [A slot | B | C] on top, [D | E] on the bottom.
#   B: symptom-based test effectiveness   C: screening-based test effectiveness
#   D: realized generation interval        E: individual reproduction number nu_i
#
# Strengthening vs. the old Inkscape composite:
#   - psi is encoded by ONE light->dark purple gradient in every panel, with a
#     single shared psi legend (bursty psi~0 light, sustained psi~1 dark)
#   - B/C titles tinted to their detection mechanism (symptom = red, screen = blue)
#   - E spikes at nu=0 / nu=R0 annotated (all-or-nothing detection under burst)
# =============================================================================

suppressMessages({library(tidyverse); library(patchwork)})

OBJ <- "figures/_objects"; FIGDIR <- "figures/"
load_obj <- function(nm) readRDS(file.path(OBJ, paste0(nm, ".rds")))
BASE <- 7

# Lightness encodes psi in every panel: dark = bursty (psi~0) -> light = sustained
# (psi~1). Hue encodes mechanism: symptom = reds (B), screening = blues (C),
# downstream consequences = purples (D, E). The shared psi legend uses neutral
# greys so intensity reads as psi regardless of hue. (ColorBrewer Reds/Blues/Purples.)
psi_lvls   <- c("0", "0.2", "0.5", "0.8", "1")
psi_red    <- setNames(c("#67000d", "#a50f15", "#ef3b2c", "#fc9272", "#fcbba1"), psi_lvls)
psi_blue   <- setNames(c("#08306b", "#08519c", "#4292c6", "#9ecae1", "#c6dbef"), psi_lvls)
psi_purple <- setNames(c("#3f007d", "#54278f", "#756bb1", "#9e9ac8", "#cbc9e2"), psi_lvls)
psi_grey   <- setNames(c("#2b2b2b", "#5a5a5a", "#8c8c8c", "#b5b5b5", "#d9d9d9"), psi_lvls)
d_cols <- c("psi = 0" = psi_purple[["0"]], "psi = 0.5" = psi_purple[["0.5"]],
						"psi = 1" = psi_purple[["1"]], "No intervention" = "grey30")
e_cols <- c("psi = 0" = psi_purple[["0"]], "psi = 1" = psi_purple[["1"]])
COL_SYMP <- psi_red[["0.5"]]; COL_SCREEN <- psi_blue[["0.5"]]   # title accents (match the gradients)

theme_fig5 <- theme_classic(base_size = BASE) +
	theme(plot.title = element_text(size = BASE, face = "bold"),
				plot.subtitle = element_blank(),
				legend.key.size = unit(9, "pt"),
				plot.margin = margin(3, 4, 3, 4))

pct_y <- scale_y_continuous(labels = function(x) round(100 * x))  # TE stored as 0-1

# --- Panel B: symptom-based test effectiveness (reds) ------------------------
pB <- load_obj("fig_te_symp_omicron") +
	scale_colour_manual(values = psi_red, name = "ψ",
											guide = guide_legend(override.aes = list(colour = unname(psi_grey)))) + pct_y +
	labs(title = "Symptom-based detection", subtitle = NULL,
			 x = "Mean symptom onset\n(days from peak)", y = "Test effectiveness (%)") +
	theme_fig5 + theme(plot.title = element_text(colour = COL_SYMP, size = BASE, face = "bold"))

# --- Panel C: screening-based test effectiveness (blues) ---------------------
pC <- load_obj("fig_te_testing_omicron") +
	scale_colour_manual(values = psi_blue, guide = "none") + pct_y +
	labs(title = "Screening-based detection", subtitle = NULL,
			 x = "Gap between tests (days)", y = "Test effectiveness (%)") +
	theme_fig5 + theme(plot.title = element_text(colour = COL_SCREEN, size = BASE, face = "bold"))

# --- Panel D: realized generation interval -----------------------------------
# draw the dashed "uncontrolled" reference on top (it coincides with bursty psi=0)
dobj <- load_obj("isolation_gi_truncation")
dobj$layers <- dobj$layers[c(2, 1)]
pD <- dobj +
	scale_colour_manual(values = d_cols, guide = "none") +
	scale_linetype_manual(values = c("psi = 0" = "solid", "psi = 0.5" = "solid",
																	 "psi = 1" = "solid", "No intervention" = "dashed"),
												name = NULL, breaks = "No intervention", labels = "uncontrolled") +
	labs(title = "Realized generation interval", subtitle = NULL,
			 x = "Generation interval (days)", y = "Density") +
	theme_fig5

# --- Panel E: individual reproduction number ---------------------------------
# saved object is faceted (psi=0 spikes vs psi=1 bell, free-y); overlay them in a
# single panel with the tall bursty spikes clipped so both are visible.
eobj <- load_obj("isolation_overdispersion")
eobj$facet <- ggplot2::facet_null()
eobj$layers[[2]] <- NULL                       # drop the per-group mean vlines (clutter)
eobj$layers[[1]]$position <- position_identity()
eobj$layers[[1]]$aes_params$alpha <- 0.85
pE <- eobj +
	scale_fill_manual(values = e_cols, guide = "none") +
	coord_cartesian(ylim = c(0, 1.4)) +
	annotate("text", x = 0.15, y = 1.36, label = "fully\nsuppressed", hjust = 0, vjust = 1,
					 size = BASE / 3.2, colour = "grey25", lineheight = 0.9) +
	annotate("text", x = 5.85, y = 1.36, label = "fully\nmissed", hjust = 1, vjust = 1,
					 size = BASE / 3.2, colour = "grey25", lineheight = 0.9) +
	labs(title = "Individual reproduction number", subtitle = NULL, y = "Density") +
	theme_fig5

# --- assemble: [A slot | B | C] / [D | E] ------------------------------------
top    <- plot_spacer() | pB | pC
bottom <- pD | pE
fig5 <- (top / bottom) + plot_layout(guides = "collect", heights = c(1, 1)) &
	theme(legend.position = "bottom")

ggsave(file.path(FIGDIR, "fig5_full.pdf"), fig5, width = 7, height = 5, device = cairo_pdf)
cat(sprintf("Saved: %sfig5_full.pdf  (panel A slot reserved for the Mathematica schematic)\n", FIGDIR))
