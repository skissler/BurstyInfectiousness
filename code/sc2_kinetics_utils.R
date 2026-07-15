# =============================================================================
# sc2_kinetics_utils.R
# Shared helpers for the SARS-CoV-2 viral-kinetics ("tent") analysis. Every
# function operates on a `feats` list (as saved by sc2_kinetics_fit.R):
#   feats$dp, feats$wp, feats$wr, feats$tp  : [ndraws x n_infection] posterior draws
#   feats$lod                               : Ct limit of detection
#   feats$raw                               : long data (InfectionEvent, id_clean, t, y)
#   feats$key                               : per-infection covariates
# so the window analysis and every figure can be regenerated from the saved
# .rds alone -- no MCMC re-run. Assumes the tidyverse is attached.
# =============================================================================

# --- Window (time below a Ct threshold) computations -------------------------

# Per-draw window widths (days a trajectory spends below `ct`): [ndraws x n_id].
ct_window_draws <- function(feats, ct) {
	D <- feats$lod - ct
	(feats$wp + feats$wr) * pmax(feats$dp - D, 0) / feats$dp
}

# Summaries at one Ct threshold.
ct_window <- function(feats, ct) {
	W <- ct_window_draws(feats, ct)
	list(ct     = ct,
			 perinf = colMeans(W),                          # per-infection posterior-mean window
			 reach  = mean(colMeans(feats$dp) > (feats$lod - ct)),  # fraction whose mean peak reaches
			 popmed = apply(W, 1, median))                  # per-draw population median
}

# Window-width table across a range of Ct thresholds.
window_summary_by_ct <- function(feats, cts = 40:20) {
	n <- ncol(feats$dp)
	map_dfr(cts, function(ctv) {
		r <- ct_window(feats, ctv); q <- quantile(r$perinf, c(.25, .5, .75))
		tibble(Ct = ctv,
					 n_exceed        = sum(colMeans(feats$dp) > (feats$lod - ctv)),
					 n_total         = n,
					 pct_exceed      = round(100 * r$reach, 1),
					 window_median_d = round(q[2], 3),
					 window_iqr_lo   = round(q[1], 3),
					 window_iqr_hi   = round(q[3], 3),
					 window_q05      = round(quantile(r$perinf, .05), 3),
					 window_q95      = round(quantile(r$perinf, .95), 3),
					 popmed_mean_d   = round(mean(r$popmed), 3),
					 popmed_cri_lo   = round(quantile(r$popmed, .025), 3),
					 popmed_cri_hi   = round(quantile(r$popmed, .975), 3))
	})
}

# Per-infection feature table (window at each threshold + peak params + covariates).
# Reference Ct thresholds: 32 = "bare" infectiousness (virus rarely culturable
# above this Ct), 26 = ~50% culture positivity (see Methods/Supplement).
infection_window_features <- function(feats, thresholds = c(32, 26)) {
	out <- tibble(InfectionEvent = as.integer(colnames(feats$dp)),
								dp = colMeans(feats$dp), wp = colMeans(feats$wp), wr = colMeans(feats$wr))
	for (ctv in thresholds) out[[paste0("window_ct", ctv)]] <- ct_window(feats, ctv)$perinf
	if (!is.null(feats$key)) out <- left_join(out, feats$key, by = "InfectionEvent")
	out
}

# Console report of windows at a few thresholds + the crude-psi bound.
report_ct_windows <- function(feats, thresholds = c(35, 32, 26), gi_sd = 2.0) {
	cat(sprintf("\n=== Duration below Ct threshold, n = %d ===\n", ncol(feats$dp)))
	for (ctv in thresholds) {
		r <- ct_window(feats, ctv); q <- quantile(r$perinf, c(.25, .5, .75))
		cat(sprintf("Ct<%d : %.0f%% reach it | window (d): median %.1f (IQR %.1f-%.1f) | pop-median %.1f (95%% CrI %.1f-%.1f)\n",
								ctv, 100 * r$reach, q[2], q[1], q[3],
								mean(r$popmed), quantile(r$popmed, .025), quantile(r$popmed, .975)))
	}
	# crude psi bound from the "bare infectiousness" window (Ct<32)
	r_bare <- ct_window(feats, 32); sd_bare <- median(r_bare$perinf) / sqrt(12)
	cat(sprintf("\nCt<32 within-infector timing SD ~%.1f d (uniform) vs GI SD ~%.1f d => crude psi (upper bound) ~%.2f\n",
							sd_bare, gi_sd, min((sd_bare / gi_sd)^2, 1)))
	cat("NOTE: time-below-Ct32 over-states the transmission window (viability/emission narrow it) -> conservative psi.\n")
	invisible(NULL)
}

# --- Figures ------------------------------------------------------------------

# Fitted "tent" trajectories vs. raw Ct for a random subset of infections.
plot_trajectories <- function(feats, n = 25, seed = 1, ncol = 5, ct_line = 26) {
	lod  <- feats$lod
	tp_m <- colMeans(feats$tp); dp_m <- colMeans(feats$dp)
	wp_m <- colMeans(feats$wp); wr_m <- colMeans(feats$wr)  # all named by InfectionEvent
	set.seed(seed)
	elig <- feats$raw %>% count(InfectionEvent, name = "npt") %>% filter(npt >= 4) %>% pull(InfectionEvent)
	sel  <- sort(sample(elig, min(n, length(elig))))
	tent <- function(tg, tp, wp, wr, dp) ifelse(tg <= tp, dp/wp*(tg-(tp-wp)), dp - dp/wr*(tg-tp))
	curve_df <- map_dfr(as.character(sel), function(ie) {
		tg <- seq(tp_m[ie] - wp_m[ie] - 1, tp_m[ie] + wr_m[ie] + 1, length.out = 120)
		tibble(InfectionEvent = as.integer(ie), t = tg,
					 ct = lod - pmax(tent(tg, tp_m[ie], wp_m[ie], wr_m[ie], dp_m[ie]), 0))
	})
	raw_df <- feats$raw %>% filter(InfectionEvent %in% sel)
	ggplot(curve_df, aes(t, ct)) +
		geom_hline(yintercept = ct_line, linetype = 2, colour = "#d1682e", linewidth = 0.4) +
		geom_line(colour = "#2f6db0", linewidth = 0.6) +
		geom_point(data = raw_df, aes(t, y), size = 0.9, alpha = 0.7) +
		facet_wrap(~ InfectionEvent, scales = "free_x", ncol = ncol) +
		scale_y_reverse(limits = c(lod, 10)) +
		labs(x = "day", y = "Ct") +
		theme_bw(base_size = 8)
}

# Per-infection window cloud (`long`) + per-Ct median/5-95% summary (`summ`) used
# by the window-vs-Ct figures.
window_vs_ct_data <- function(feats, cts = 20:40) {
	long <- map_dfr(cts, ~ tibble(ct = .x, window = ct_window(feats, .x)$perinf)) %>%
		filter(window > 0)
	summ <- long %>% group_by(ct) %>%
		summarise(med = median(window), lo = quantile(window, .05), hi = quantile(window, .95),
							n = n(), .groups = "drop")
	list(long = long, summ = summ)
}

# Window width vs. Ct threshold: median + 5/95% of the per-infection posterior
# means, over a horizontally-jittered cloud of the individual means. Ct axis is
# reversed (40 at origin -> 20), so leftward = higher viral load. Optionally marks
# the 5-95% span of the *population* generation interval g(tau) as a horizontal
# reference (how broad the population GI is, vs. how narrow an individual window
# is). gi defaults = SARS-CoV-2 omicron Gamma(2.39, 0.339).
plot_window_vs_ct <- function(feats, cts = 20:40, point_alpha = 0.05, jitter_w = 0.35,
															mark_gi = TRUE, gi_shape = 2.39, gi_rate = 0.339,
															gi_span = c(.05, .95)) {
	d <- window_vs_ct_data(feats, cts); long <- d$long; summ <- d$summ
	p <- ggplot()
	if (mark_gi) {
		q <- qgamma(gi_span, gi_shape, gi_rate); gi_w <- diff(q)
		p <- p +
			geom_hline(yintercept = gi_w, linetype = 2, colour = "grey40", linewidth = 0.5) +
			annotate("text", x = max(cts), y = gi_w, vjust = -0.6, hjust = 0, size = 2.8, colour = "grey30",
							 label = sprintf("population generation interval, %g-%g%% span (%.1f d)",
															 100 * gi_span[1], 100 * gi_span[2], gi_w))
	}
	p +
		geom_jitter(data = long, aes(ct, window), width = jitter_w, height = 0,
								size = 0.12, alpha = point_alpha, colour = "#2f6db0") +
		geom_linerange(data = summ, aes(ct, ymin = lo, ymax = hi), linewidth = 0.5) +
		geom_line(data = summ, aes(ct, med), linewidth = 0.5) +
		geom_point(data = summ, aes(ct, med), size = 1.4) +
		scale_x_reverse(breaks = seq(min(cts), max(cts), 2)) +
		labs(x = "Ct threshold  (leftward = higher viral load)", y = "Time below Ct threshold (days)") +
		theme_bw(base_size = 9)
}

# Publication-quality version of the window-vs-Ct panel (intended as a subplot of
# main Fig. 2). Same data as plot_window_vs_ct, but theme_classic, y clipped to
# [0, ymax], and terser labels ("Ct threshold" / "Time above threshold (days)" --
# "above" = viral load above that threshold; explained in the caption).
# reach_overlay adds a floating, monotone-declining stepped area at the top of the
# panel: the fraction of infections whose peak viral load reaches (crosses) each
# Ct threshold (= pct_exceed, 100% at Ct 40 down to ~6.5% at Ct 20), with a small
# floating 0/50/100% axis tucked into the open right-side space.
plot_window_vs_ct_pub <- function(feats, cts = 20:40, point_alpha = 0.05, jitter_w = 0.35,
																	ymax = 25, base_size = 9, x_break = 4,
																	mark_gi = TRUE, gi_shape = 2.39, gi_rate = 0.339, gi_span = c(.05, .95),
																	gi_short = FALSE,
																	reach_overlay = TRUE, overlay_baseline = 22, overlay_band = 3) {
	d <- window_vs_ct_data(feats, cts); long <- d$long; summ <- d$summ
	p <- ggplot()
	if (mark_gi) {
		gi_w <- diff(qgamma(gi_span, gi_shape, gi_rate))
		gi_lbl <- if (gi_short) sprintf("GI %g–%g%% span (%.1f d)", 100 * gi_span[1], 100 * gi_span[2], gi_w)
							else sprintf("generation interval, %g-%g%% span (%.1f d)", 100 * gi_span[1], 100 * gi_span[2], gi_w)
		p <- p +
			geom_hline(yintercept = gi_w, linetype = 2, colour = "grey40", linewidth = 0.5) +
			# label right-justified at the right (low-Ct) edge, clear of the plotted data
			annotate("text", x = min(cts), y = gi_w, vjust = -0.6, hjust = 1, size = 2.6, colour = "grey30",
							 label = gi_lbl)
	}
	p <- p +
		geom_jitter(data = long, aes(ct, window), width = jitter_w, height = 0,
								size = 0.12, alpha = point_alpha, colour = "#2f6db0") +
		geom_linerange(data = summ, aes(ct, ymin = lo, ymax = hi), linewidth = 0.5) +
		geom_line(data = summ, aes(ct, med), linewidth = 0.5) +
		geom_point(data = summ, aes(ct, med), size = 1.4)
	xlim_r <- NULL
	if (reach_overlay) {
		y0 <- overlay_baseline; band <- overlay_band
		reach <- map_dfr(cts, ~ tibble(ct = .x, pct = 100 * ct_window(feats, .x)$reach)) %>%
			mutate(top = y0 + band * pct / 100)
		xr <- min(cts) - 0.4                                        # floating % axis, in the right margin
		ticks <- tibble(p = c(0, 50, 100), y = y0 + band * c(0, 50, 100) / 100)
		p <- p +
			geom_rect(data = reach, inherit.aes = FALSE,                # contiguous rects -> stepped area
								aes(xmin = ct - 0.5, xmax = ct + 0.5, ymin = y0, ymax = top),
								fill = "grey55", alpha = 0.55) +
			geom_step(data = reach, inherit.aes = FALSE, aes(ct, top), direction = "mid",
								colour = "grey35", linewidth = 0.3) +
			annotate("segment", x = xr, xend = xr, y = y0, yend = y0 + band, colour = "grey35", linewidth = 0.3) +
			annotate("segment", x = xr, xend = xr - 0.22, y = ticks$y, yend = ticks$y, colour = "grey35", linewidth = 0.3) +
			annotate("text", x = xr - 0.35, y = ticks$y, label = paste0(ticks$p, "%"),
							 hjust = 0, size = 1.9, colour = "grey35") +
			annotate("text", x = min(cts) + 4.5, y = y0 + band - 0.4, label = "reaching Ct",
							 hjust = 0.5, vjust = 1, size = 2, colour = "grey20")
		xlim_r <- c(max(cts) + 0.8, min(cts) - 2)                    # widen right for the floating axis
	}
	p +
		scale_x_reverse(breaks = seq(min(cts), max(cts), x_break)) +
		coord_cartesian(ylim = c(0, ymax), xlim = xlim_r) +
		labs(x = "Ct threshold", y = "Time above threshold (days)") +
		theme_classic(base_size = base_size)
}

# Narrow individual infectious windows composing the broad population generation
# interval g(tau). Widths are the real fitted medians; window placement is
# illustrative (centers drawn from g(tau)). gi defaults = SARS-CoV-2 omicron.
plot_gi_composition <- function(feats, gi_shape = 2.39, gi_rate = 0.339,
																ct_ind = 32, ct_viable = 26, n_windows = 8, seed = 3) {
	med_win <- function(ct) { w <- ct_window(feats, ct)$perinf; median(w[w > 0]) }
	w_ind <- med_win(ct_ind); w_via <- med_win(ct_viable)
	x  <- seq(0, 16, length.out = 600); g <- dgamma(x, gi_shape, gi_rate); pk <- max(g)
	bump <- function(center, width) { y <- dnorm(x, center, width / 3.92); y / max(y) }
	set.seed(seed); centers <- rgamma(n_windows, gi_shape, gi_rate)
	gi  <- tibble(x = x, d = g)
	ind <- imap_dfr(centers, ~ tibble(id = .y, x = x, d = 0.55 * pk * bump(.x, w_ind)))
	via <- tibble(x = x, d = 0.90 * pk * bump(gi_shape / gi_rate, w_via))
	ggplot() +
		geom_area(data = gi,  aes(x, d), fill = "grey80", colour = "grey40", linewidth = 0.6) +
		geom_line(data = ind, aes(x, d, group = id), colour = "#2f6db0", linewidth = 0.5, alpha = 0.7) +
		geom_line(data = via, aes(x, d), colour = "#d1682e", linewidth = 0.8) +
		annotate("text", x = 12.5, y = 0.92 * pk, label = "generation interval",       colour = "grey30",  hjust = 0, size = 3) +
		annotate("text", x = 12.5, y = 0.58 * pk, label = sprintf("individual windows\n(Ct<%d)", ct_ind), colour = "#2f6db0", hjust = 0, size = 3) +
		annotate("text", x = gi_shape/gi_rate - 0.4, y = 0.98 * pk, label = sprintf("viable (Ct<%d)", ct_viable), colour = "#d1682e", hjust = 0, size = 3) +
		labs(x = expression("Days since infection (" * tau * ")"), y = "Infectiousness (arb.)") +
		theme_classic(base_size = 9) +
		theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}
