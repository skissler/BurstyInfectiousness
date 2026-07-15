#!/usr/bin/env Rscript
# =============================================================================
# sc2_kinetics_figures.R
#
# Regenerates the SARS-CoV-2 viral-kinetics figures from the stored posterior
# features (output/posterior_features_single_category.rds) WITHOUT re-running the
# MCMC. All plotting logic lives in code/sc2_kinetics_utils.R.
#
# Config via env vars (all optional):
#   CT_FEATS   path to the posterior-features rds
#              (default: output/posterior_features_single_category.rds)
#   CT_FIGDIR  output dir for figures            (default: figures/)
# =============================================================================

suppressMessages(library(tidyverse))
source("code/sc2_kinetics_utils.R")

FEATS  <- Sys.getenv("CT_FEATS",  "output/posterior_features_single_category.rds")
FIGDIR <- Sys.getenv("CT_FIGDIR", "figures/")
OUTDIR <- Sys.getenv("CT_OUT",    "output/")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

feats <- readRDS(FEATS)
cat(sprintf("Loaded %d infections x %d draws from %s\n",
						ncol(feats$dp), nrow(feats$dp), FEATS))

## window-width summary table (Ct 40..20) -- rewritten from the stored posterior
## so column/threshold changes don't require a re-fit
write_csv(window_summary_by_ct(feats, cts = 40:20),
					file.path(OUTDIR, "window_summary_by_Ct.csv"))
cat(sprintf("Saved: %swindow_summary_by_Ct.csv\n", OUTDIR))

## per-infection window features at the reference thresholds (Ct 32 / 26)
write_csv(infection_window_features(feats, thresholds = c(32, 26)),
					file.path(OUTDIR, "infection_windows_single_category.csv"))
cat(sprintf("Saved: %sinfection_windows_single_category.csv\n", OUTDIR))

## window width vs. Ct threshold ----------------------------------------------
p_win <- plot_window_vs_ct(feats)
ggsave(file.path(FIGDIR, "window_vs_ct.pdf"), p_win, width = 6, height = 4)
cat(sprintf("Saved: %swindow_vs_ct.pdf\n", FIGDIR))

## publication-quality version (main Fig. 2 subplot), wide aspect
p_win_pub <- plot_window_vs_ct_pub(feats)
ggsave(file.path(FIGDIR, "window_vs_ct_pub.pdf"), p_win_pub, width = 4.5, height = 2.8)  # ~1.6:1
cat(sprintf("Saved: %swindow_vs_ct_pub.pdf\n", FIGDIR))

## individual windows composing the generation interval -----------------------
p_gi <- plot_gi_composition(feats)
ggsave(file.path(FIGDIR, "gi_composition.pdf"), p_gi, width = 6, height = 3.4)
cat(sprintf("Saved: %sgi_composition.pdf\n", FIGDIR))

## fitted trajectories vs raw Ct (25 infections) ------------------------------
p_traj <- plot_trajectories(feats, n = 25, seed = 1)
ggsave(file.path(FIGDIR, "trajectories_ct26.png"), p_traj, width = 10, height = 8, dpi = 150)
cat(sprintf("Saved: %strajectories_ct26.png\n", FIGDIR))
