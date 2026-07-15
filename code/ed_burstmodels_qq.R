#!/usr/bin/env Rscript
# =============================================================================
# ed_burstmodels_qq.R
#
# Extended Data (standalone, 7in): QQ plots showing that each alternative burst
# model reproduces the target marginal generation interval across psi, for the
# three benchmark pathogens. For each pathogen, a 3x3 grid of QQ plots:
#   rows    = burst model (Type-I Gamma, Type-II Gamma, Log-normal)
#   columns = psi in {0, 0.5, 1}
# with the three pathogens side by side. Purple line = simulated-vs-theoretical
# quantiles; grey diagonal = perfect agreement.
#
# Recomposed from the saved QQ plot objects (figures/_objects/, written by
# sensitivity_marginal_gi.R) -- no recomputation. Meant to sit below the
# hand-made burst-shape panels in post (Inkscape).
# =============================================================================

suppressMessages({library(tidyverse); library(patchwork)})

OBJ <- "figures/_objects"; FIGDIR <- "figures/"
models   <- c("gamma", "gamma2", "lognormal")
mod_lab  <- c(gamma = "Type-I Gamma", gamma2 = "Type-II Gamma", lognormal = "Log-normal")
paths    <- c("measles", "omicron", "influenza")
path_lab <- c(measles = "Measles", omicron = "SARS-CoV-2 omicron", influenza = "Influenza")
PSI      <- c(0, 0.5, 1)

qq <- map_dfr(models, function(m) map_dfr(paths, function(pth) {
	d <- readRDS(file.path(OBJ, sprintf("sensitivity_marginal_gi_qq_%s_%s.rds", m, pth)))$data
	tibble(model = m, pathogen = pth,
				 psi = as.numeric(sub("psi == ", "", as.character(d$psi_label))),
				 theoretical = d$theoretical, empirical = d$empirical)
})) %>% filter(psi %in% PSI) %>%
	mutate(model = factor(mod_lab[model], levels = mod_lab),
				 psi_lab = factor(paste0("psi==", psi)))

# subsample the 300 quantile pairs to 50 per cell for a conventional QQ scatter
qq_pts <- qq %>% group_by(pathogen, model, psi_lab) %>%
	slice(round(seq(1, n(), length.out = 50))) %>% ungroup()

panel <- function(pth, leftmost = FALSE) {
	p <- ggplot(filter(qq_pts, pathogen == pth), aes(theoretical, empirical)) +
		geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.3) +
		geom_point(colour = "black", size = 0.25, alpha = 0.9) +
		facet_grid(model ~ psi_lab, labeller = labeller(psi_lab = label_parsed), switch = "y") +
		coord_equal() +
		labs(title = path_lab[pth], x = "Theoretical quantile",
				 y = if (leftmost) "Empirical quantile" else NULL) +
		theme_bw(base_size = 6) +
		theme(plot.title = element_text(face = "bold", size = 7, hjust = 0.5),
					panel.grid = element_blank(),
					axis.text = element_blank(), axis.ticks = element_blank(),
					strip.text.x = element_text(size = 5.5), strip.background = element_blank(),
					panel.spacing = unit(1.5, "pt"))
	# model (row) strips only on the leftmost pathogen; drop the redundant repeats
	if (leftmost) p + theme(strip.text.y.left = element_text(size = 5.5, angle = 90))
	else          p + theme(strip.text.y = element_blank())
}

fig <- wrap_plots(panel("measles", leftmost = TRUE),
									panel("omicron"), panel("influenza"), nrow = 1)
ggsave(file.path(FIGDIR, "ED_burstmodels_qq.pdf"), fig, width = 7, height = 3, device = cairo_pdf)
cat(sprintf("Saved: %sED_burstmodels_qq.pdf\n", FIGDIR))
