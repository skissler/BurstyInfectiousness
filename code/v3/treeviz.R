# ==============================================================================
# treeviz.R — within-cluster vs. random serial interval pair scatterplot
#
# Requires all_disease_clusters and pathogen_params in the environment
# (produced by psi_empirical_litgi.R). Run that script first, or source it.
# ==============================================================================

library(ggplot2)
library(dplyr)
library(purrr)
library(hexbin)

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------

target_disease  <- "MERS"
n_random_pairs  <- 500   # number of random (S1, S2) draws for the background
set.seed(42)

# ------------------------------------------------------------------------------
# 1. Pull MERS clusters (require >= 2 offspring so we have real pairs)
# ------------------------------------------------------------------------------

all_clusters <- all_disease_clusters[[target_disease]]$clusters
multi_clusters <- all_clusters[sapply(all_clusters, length) >= 2]

si_pool <- unlist(multi_clusters)

# ------------------------------------------------------------------------------
# 2. Random pairs — shuffle the pool and pair up
# ------------------------------------------------------------------------------

# Draw with replacement so we always get exactly n_random_pairs pairs
idx1 <- sample(seq_along(si_pool), n_random_pairs, replace = TRUE)
idx2 <- sample(seq_along(si_pool), n_random_pairs, replace = TRUE)

random_df <- tibble(
  s1   = si_pool[idx1],
  s2   = si_pool[idx2],
  type = "random"
)

# ------------------------------------------------------------------------------
# 3. Cluster-aware pairs — all ordered pairs within each cluster
# Both (a, b) and (b, a) included so the scatter is symmetric
# ------------------------------------------------------------------------------

cluster_df <- map_dfr(multi_clusters, function(cl) {
  if (length(cl) < 2) return(tibble())
  pairs <- combn(cl, 2, simplify = FALSE)
  map_dfr(pairs, function(p) {
    bind_rows(
      tibble(s1 = p[1], s2 = p[2]),
      tibble(s1 = p[2], s2 = p[1])
    )
  })
}) %>% mutate(type = "cluster")

# ------------------------------------------------------------------------------
# 4. Plot
# ------------------------------------------------------------------------------

all_vals <- c(random_df$s1, random_df$s2, cluster_df$s1, cluster_df$s2)
pad      <- diff(range(all_vals)) * 0.05
ax_lim   <- range(all_vals) + c(-pad, pad)

fig_treeviz <- ggplot() +
  geom_jitter(
    data  = random_df,
    aes(x = s1, y = s2),
    color = "grey80", size = 1.5, alpha = 0.3
  ) +
  # geom_density_2d(
  #   data  = random_df,
  #   aes(x = s1, y = s2),
  #   color = "grey50", linewidth = 0.5, alpha = 0.8, breaks=c(0.0001, 0.0002)
  # ) +
  geom_jitter(
    data  = cluster_df,
    aes(x = s1, y = s2),
    color = "darkorange", size = 2, alpha = 0.8
  ) +
  geom_density_2d(
    data  = cluster_df,
    aes(x = s1, y = s2),
    color = "black", linewidth = 0.5, alpha = 0.8, breaks=c(0.00001, 0.0001, 0.001, 0.01)
  ) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", linewidth = 0.5, color = "grey40") +
  coord_fixed(xlim = ax_lim, ylim = ax_lim) +
  theme_classic(base_size = 13) +
  labs(
    x     = expression(S[1] ~ "(days)"),
    y     = expression(S[2] ~ "(days)"),
    title = sprintf("%s: within-cluster vs. random SI pairs", target_disease),
    subtitle = sprintf(
      "Grey: %d random pairs drawn from pool  |  Orange: %d within-cluster pairs (%d clusters)",
      n_random_pairs, nrow(cluster_df) / 2, length(multi_clusters)
    )
  )

print(fig_treeviz)

if (exists("save_fig")) {
  save_fig(fig_treeviz, "treeviz_mers", width = 6, height = 6)
  cat("Saved treeviz_mers figure.\n")
}

# ==============================================================================
# COMBINED MANUSCRIPT FIGURE
# Panel A: normalised (S1/ḡ, S2/ḡ) scatter for 3 selected pathogens
# Panel B: ECDF of |S1−S2|/ḡ for all pathogens, faceted and ordered by ψ
# ==============================================================================

library(patchwork)

set.seed(42)

# Pathogens for Panel A — pick one low-ψ, one mid-ψ, one high-ψ
scatter_diseases <- c("Measles", "Norovirus", "MERS")

# Order all pathogens by posterior ψ mode (uses posterior_results if present)
if (exists("posterior_results")) {
  psi_modes         <- sapply(posterior_results, function(r) r$post_mode[1])
  all_diseases_ord  <- names(sort(psi_modes))
} else {
  all_diseases_ord  <- c("Measles", "Pneumonic plague", "COVID-19",
                          "Hepatitis A", "Norovirus", "MERS", "Ebola", "Smallpox")
  psi_modes         <- NULL
}

# ------------------------------------------------------------------------------
# Helper: normalised pairs for one disease
# ------------------------------------------------------------------------------

make_norm_pairs <- function(disease, n_rand = 500) {
  mu_d     <- pathogen_params[[disease]]$gi_mean
  clusters <- all_disease_clusters[[disease]]$clusters
  multi    <- clusters[sapply(clusters, length) >= 2]

  pool <- unlist(multi) / mu_d
  if (length(pool) < 2) return(NULL)

  rand_df <- tibble(
    s1      = sample(pool, n_rand, replace = TRUE),
    s2      = sample(pool, n_rand, replace = TRUE),
    group   = "Random",
    disease = disease
  )

  cl_df <- map_dfr(multi, function(cl) {
    cl <- cl / mu_d
    if (length(cl) < 2) return(tibble())
    map_dfr(combn(cl, 2, simplify = FALSE), function(p) {
      bind_rows(tibble(s1 = p[1], s2 = p[2]),
                tibble(s1 = p[2], s2 = p[1]))
    })
  }) %>% mutate(group = "Within-cluster", disease = disease)

  list(random = rand_df, cluster = cl_df)
}

# Helper: theoretical ellipse for (S1/ḡ, S2/ḡ) under the Gamma burst model.
#
# Under the model, S_j = l + ε_j + d_j − d_inf where l ~ Gamma((1-ψ)α, β) is
# shared across siblings and ε_j ~ Gamma(ψα, β) is individual jitter. This
# gives Var[S_j/ḡ] = (σ²_GI + 2σ²_inc)/ḡ² and
# Corr[S1/ḡ, S2/ḡ] = ((1−ψ)σ²_GI + σ²_inc) / (σ²_GI + 2σ²_inc).
# The bivariate normal ellipse axes always align with and against y = x.
make_si_ellipse <- function(disease, psi, level = 0.68, n_pts = 300) {
  pp    <- pathogen_params[[disease]]
  mu_gi <- pp$gi_mean
  alpha <- pp$alpha_gi
  beta  <- pp$beta_gi
  a_obs <- pp$a_obs
  b_obs <- pp$b_obs

  sigma2_gi  <- alpha / beta^2
  sigma2_inc <- a_obs / b_obs^2

  sigma2 <- (sigma2_gi + 2 * sigma2_inc) / mu_gi^2
  rho    <- ((1 - psi) * sigma2_gi + sigma2_inc) / (sigma2_gi + 2 * sigma2_inc)

  s_diag <- sqrt(sigma2 * (1 + rho))   # std dev along y = x
  s_perp <- sqrt(sigma2 * (1 - rho))   # std dev perpendicular to y = x
  c_val  <- sqrt(qchisq(level, df = 2))
  theta  <- seq(0, 2 * pi, length.out = n_pts + 1)

  tibble(
    x   = 1 + c_val * (s_diag * cos(theta) + s_perp * sin(theta)) / sqrt(2),
    y   = 1 + c_val * (s_diag * cos(theta) - s_perp * sin(theta)) / sqrt(2),
    psi = factor(psi)
  )
}

psi_ref_vals   <- c(0.25, 0.5, 0.75, 1.0)
psi_ref_colors <- c("0.25" = "#c6dbef", "0.5" = "#6baed6",
                    "0.75" = "#2171b5", "1"   = "#084594")

# ------------------------------------------------------------------------------
# Panel A: scatter for 3 pathogens on a shared normalised axis
# ------------------------------------------------------------------------------

scatter_pairs <- map(scatter_diseases, make_norm_pairs)
names(scatter_pairs) <- scatter_diseases

scatter_vals <- unlist(map(scatter_pairs, function(x) {
  c(x$random$s1, x$random$s2, x$cluster$s1, x$cluster$s2)
}))
pad_s    <- diff(range(scatter_vals)) * 0.05
ax_lim_s <- range(scatter_vals) + c(-pad_s, pad_s)

scatter_panels <- map(scatter_diseases, function(disease) {
  rd <- scatter_pairs[[disease]]$random
  cd <- scatter_pairs[[disease]]$cluster

  psi_str <- if (!is.null(psi_modes) && disease %in% names(psi_modes)) {
    sprintf("ψ ≈ %.2f", psi_modes[[disease]])
  } else { "" }

  ggplot() +
    geom_point(data = rd, aes(x = s1, y = s2),
               color = "grey80", size = 1, alpha = 0.4) +
    geom_point(data = cd, aes(x = s1, y = s2),
               color = "darkorange", size = 1.8, alpha = 0.75) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", linewidth = 0.5, color = "grey40") +
    coord_fixed(xlim = ax_lim_s, ylim = ax_lim_s, expand = FALSE) +
    theme_classic(base_size = 11) +
    labs(
      x        = expression(S[1] / bar(g)),
      y        = expression(S[2] / bar(g)),
      title    = disease,
      subtitle = psi_str
    )
})

panel_a <- wrap_plots(scatter_panels, nrow = 1) +
  plot_annotation(title = "A")

# ------------------------------------------------------------------------------
# Panel B: ECDF of |S1−S2|/ḡ for all pathogens, faceted by ψ
# ------------------------------------------------------------------------------

ecdf_data <- map_dfr(all_diseases_ord, function(disease) {
  pairs <- make_norm_pairs(disease, n_rand = 1000)
  if (is.null(pairs)) return(tibble())
  bind_rows(pairs$random, pairs$cluster) %>%
    mutate(abs_diff = abs(s1 - s2)) %>%
    select(disease, group, abs_diff)
})

if (!is.null(psi_modes)) {
  ecdf_data <- ecdf_data %>%
    mutate(
      fac_label = sprintf("%s (ψ=%.2f)", disease, psi_modes[disease]),
      fac_label = factor(fac_label,
                         levels = sprintf("%s (ψ=%.2f)",
                                          all_diseases_ord,
                                          psi_modes[all_diseases_ord]))
    )
} else {
  ecdf_data <- ecdf_data %>%
    mutate(fac_label = factor(disease, levels = all_diseases_ord))
}

panel_b <- ggplot(ecdf_data, aes(x = abs_diff, color = group, linetype = group)) +
  stat_ecdf(linewidth = 0.8) +
  scale_color_manual(
    values = c("Random" = "grey50", "Within-cluster" = "darkorange"),
    name   = NULL
  ) +
  scale_linetype_manual(
    values = c("Random" = "dashed", "Within-cluster" = "solid"),
    name   = NULL
  ) +
  facet_wrap(~ fac_label, nrow = 2) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", strip.text = element_text(size = 9)) +
  labs(
    x     = expression("|" * S[1] - S[2] * "| / " * bar(g)),
    y     = "Cumulative proportion",
    title = "B"
  )

# ------------------------------------------------------------------------------
# Combine and save
# ------------------------------------------------------------------------------

fig_combined <- (panel_a / panel_b) +
  plot_annotation(
    title = "Within-cluster serial interval correlations across pathogens",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  ) +
  plot_layout(heights = c(1, 1.4))

print(fig_combined)

if (exists("save_fig")) {
  save_fig(fig_combined, "treeviz_combined", width = 12, height = 9)
  cat("Saved treeviz_combined figure.\n")
}

# ==============================================================================
# Per-pathogen scatter plots: S2/ḡ vs S1/ḡ
# ==============================================================================

for (disease in names(all_disease_clusters)) {
  pairs <- make_norm_pairs(disease, n_rand = 500)
  if (is.null(pairs)) {
    cat(sprintf("  Skipping %s (insufficient data)\n", disease))
    next
  }

  rd <- pairs$random
  cd <- pairs$cluster

  all_vals_d <- c(rd$s1, rd$s2, cd$s1, cd$s2)
  pad_d      <- diff(range(all_vals_d)) * 0.05
  ax_lim_d   <- range(all_vals_d) + c(-pad_d, pad_d)

  psi_str <- if (!is.null(psi_modes) && disease %in% names(psi_modes)) {
    sprintf("ψ ≈ %.2f  |  %d within-cluster pairs from %d clusters",
            psi_modes[[disease]], nrow(cd) / 2,
            sum(sapply(all_disease_clusters[[disease]]$clusters,
                       function(cl) length(cl) >= 2)))
  } else {
    sprintf("%d within-cluster pairs", nrow(cd) / 2)
  }

  sd_rand  <- sd(rd$s1 - rd$s2)
  sd_clust <- sd(cd$s1 - cd$s2)
  x_seq    <- seq(ax_lim_d[1], ax_lim_d[2], length.out = 200)
  band_rand  <- tibble(x = x_seq, ymin = x - sd_rand,  ymax = x + sd_rand)
  band_clust <- tibble(x = x_seq, ymin = x - sd_clust, ymax = x + sd_clust)

  ellipse_df_d <- map_dfr(psi_ref_vals, function(pv) {
    tryCatch(make_si_ellipse(disease, pv), error = function(e) tibble())
  })

  fig_d <- ggplot() +
    # geom_ribbon(data = band_rand,  aes(x = x, ymin = ymin, ymax = ymax),
                # fill = "grey70", alpha = 0.25, inherit.aes = FALSE) +
    # geom_ribbon(data = band_clust, aes(x = x, ymin = ymin, ymax = ymax),
                # fill = "darkorange", alpha = 0.25, inherit.aes = FALSE) +
    # geom_path(data = ellipse_df_d,
              # aes(x = x, y = y, group = psi, color = psi),
              # linewidth = 0.7) +
    scale_color_manual(values = psi_ref_colors, name = "ψ") +
    geom_jitter(data = rd, aes(x = s1, y = s2),
               color = "grey80", size = 1.2, alpha = 0.4) +
    geom_jitter(data = cd, aes(x = s1, y = s2),
               color = "darkorange", size = 2, alpha = 0.75) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", linewidth = 0.5, color = "grey40") +
    coord_fixed(xlim = ax_lim_d, ylim = ax_lim_d, expand = FALSE) +
    theme_classic(base_size = 13) +
    labs(
      x        = expression(S[1] / bar(g)),
      y        = expression(S[2] / bar(g)),
      title    = disease,
      subtitle = psi_str
    )

  print(fig_d)

  if (exists("save_fig")) {
    slug <- gsub("[^A-Za-z0-9]", "_", tolower(disease))
    save_fig(fig_d, sprintf("treeviz_scatter_%s", slug), width = 6, height = 6)
    cat(sprintf("  Saved treeviz_scatter_%s\n", slug))
  }
}

# ==============================================================================
# Per-pathogen hex-grid heat maps: S2/ḡ vs S1/ḡ  (own axis limits)
#
# Within-cluster pairs shown as hexbins coloured by count; random pairs as
# light grey points.  Saves as treeviz_hex_{slug}.
# ==============================================================================

for (disease in names(all_disease_clusters)) {
  pairs <- make_norm_pairs(disease, n_rand = 500)
  if (is.null(pairs)) next

  rd <- pairs$random
  cd <- pairs$cluster

  all_vals_d <- c(rd$s1, rd$s2, cd$s1, cd$s2)
  pad_d      <- diff(range(all_vals_d)) * 0.05
  ax_lim_d   <- range(all_vals_d) + c(-pad_d, pad_d)

  psi_str <- if (!is.null(psi_modes) && disease %in% names(psi_modes)) {
    sprintf("ψ ≈ %.2f  |  %d within-cluster pairs from %d clusters",
            psi_modes[[disease]], nrow(cd) / 2,
            sum(sapply(all_disease_clusters[[disease]]$clusters,
                       function(cl) length(cl) >= 2)))
  } else {
    sprintf("%d within-cluster pairs", nrow(cd) / 2)
  }

  bins_d <- max(5L, min(25L, round(nrow(cd)^(1/3) * 2)))

  fig_hex_d <- ggplot() +
    geom_point(data = rd, aes(x = s1, y = s2),
               color = "grey85", size = 0.8, alpha = 0.3) +
    geom_hex(data = cd, aes(x = s1, y = s2), bins = bins_d, alpha = 0.9) +
    scale_fill_gradient(low = "#fff5e6", high = "#d94801", name = "Count") +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", linewidth = 0.5, color = "grey40") +
    coord_fixed(xlim = ax_lim_d, ylim = ax_lim_d, expand = FALSE) +
    theme_classic(base_size = 13) +
    labs(
      x        = expression(S[1] / bar(g)),
      y        = expression(S[2] / bar(g)),
      title    = disease,
      subtitle = psi_str
    )

  print(fig_hex_d)

  if (exists("save_fig")) {
    slug <- gsub("[^A-Za-z0-9]", "_", tolower(disease))
    save_fig(fig_hex_d, sprintf("treeviz_hex_%s", slug), width = 6, height = 6)
    cat(sprintf("  Saved treeviz_hex_%s\n", slug))
  }
}

# ==============================================================================
# Per-pathogen scatter plots: capped shared axis
#
# Axis limits are shared across all pathogens and capped at the 99th percentile
# of all normalised SI values, so the panels are directly comparable.
# ==============================================================================

# Collect all normalised values across every pathogen to set shared limits
all_norm_vals <- unlist(map(names(all_disease_clusters), function(disease) {
  pairs <- make_norm_pairs(disease, n_rand = 500)
  if (is.null(pairs)) return(NULL)
  c(pairs$random$s1, pairs$random$s2, pairs$cluster$s1, pairs$cluster$s2)
}))

cap      <- quantile(all_norm_vals, 0.99)
pad_cap  <- cap * 0.05
ax_lim_cap <- c(0, cap + pad_cap)   # start at 0 since SIs are positive

for (disease in names(all_disease_clusters)) {
  pairs <- make_norm_pairs(disease, n_rand = 500)
  if (is.null(pairs)) {
    cat(sprintf("  Skipping %s (insufficient data)\n", disease))
    next
  }

  rd <- pairs$random
  cd <- pairs$cluster

  psi_str <- if (!is.null(psi_modes) && disease %in% names(psi_modes)) {
    sprintf("ψ ≈ %.2f  |  %d within-cluster pairs from %d clusters",
            psi_modes[[disease]], nrow(cd) / 2,
            sum(sapply(all_disease_clusters[[disease]]$clusters,
                       function(cl) length(cl) >= 2)))
  } else {
    sprintf("%d within-cluster pairs", nrow(cd) / 2)
  }

  sd_rand  <- sd(rd$s1 - rd$s2)
  sd_clust <- sd(cd$s1 - cd$s2)
  x_seq    <- seq(ax_lim_cap[1], ax_lim_cap[2], length.out = 200)
  band_rand  <- tibble(x = x_seq, ymin = x - sd_rand,  ymax = x + sd_rand)
  band_clust <- tibble(x = x_seq, ymin = x - sd_clust, ymax = x + sd_clust)

  ellipse_df_cap <- map_dfr(psi_ref_vals, function(pv) {
    tryCatch(make_si_ellipse(disease, pv), error = function(e) tibble())
  })

  fig_cap <- ggplot() +
    # geom_ribbon(data = band_rand,  aes(x = x, ymin = ymin, ymax = ymax),
                # fill = "grey70", alpha = 0.25, inherit.aes = FALSE) +
    # geom_ribbon(data = band_clust, aes(x = x, ymin = ymin, ymax = ymax),
                # fill = "darkorange", alpha = 0.25, inherit.aes = FALSE) +
    # geom_path(data = ellipse_df_cap,
              # aes(x = x, y = y, group = psi, color = psi),
              # linewidth = 0.7) +
    scale_color_manual(values = psi_ref_colors, name = "ψ") +
    geom_jitter(data = rd, aes(x = s1, y = s2),
               color = "grey80", size = 2, alpha = 0.4) +
    geom_jitter(data = cd, aes(x = s1, y = s2),
               color = "darkorange", size = 4, alpha = 0.75) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", linewidth = 0.5, color = "grey40") +
    coord_fixed(xlim = ax_lim_cap, ylim = ax_lim_cap, expand = FALSE) +
    theme_classic(base_size = 13) +
    labs(
      x        = expression(S[1] / bar(g)),
      y        = expression(S[2] / bar(g)),
      title    = disease,
      subtitle = psi_str
    )

  print(fig_cap)

  if (exists("save_fig")) {
    slug <- gsub("[^A-Za-z0-9]", "_", tolower(disease))
    save_fig(fig_cap, sprintf("treeviz_scatter_%s_capped", slug), width = 6, height = 6)
    cat(sprintf("  Saved treeviz_scatter_%s_capped\n", slug))
  }
}

# ==============================================================================
# Per-pathogen hex-grid heat maps: S2/ḡ vs S1/ḡ  (capped shared axis)
#
# Same layout as the capped scatter but using hexbins for within-cluster pairs.
# Saves as treeviz_hex_{slug}_capped.
# ==============================================================================

for (disease in names(all_disease_clusters)) {
  pairs <- make_norm_pairs(disease, n_rand = 500)
  if (is.null(pairs)) {
    cat(sprintf("  Skipping %s (insufficient data)\n", disease))
    next
  }

  rd <- pairs$random
  cd <- pairs$cluster

  psi_str <- if (!is.null(psi_modes) && disease %in% names(psi_modes)) {
    sprintf("ψ ≈ %.2f  |  %d within-cluster pairs from %d clusters",
            psi_modes[[disease]], nrow(cd) / 2,
            sum(sapply(all_disease_clusters[[disease]]$clusters,
                       function(cl) length(cl) >= 2)))
  } else {
    sprintf("%d within-cluster pairs", nrow(cd) / 2)
  }

  bins_cap <- max(5L, min(25L, round(nrow(cd)^(1/3) * 2)))

  fig_hex_cap <- ggplot() +
    geom_point(data = rd, aes(x = s1, y = s2),
               color = "grey85", size = 1, alpha = 0.3) +
    geom_hex(data = cd, aes(x = s1, y = s2), bins = bins_cap, alpha = 0.9) +
    scale_fill_gradient(low = "#fff5e6", high = "#d94801", name = "Count") +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", linewidth = 0.5, color = "grey40") +
    coord_fixed(xlim = ax_lim_cap, ylim = ax_lim_cap, expand = FALSE) +
    theme_classic(base_size = 13) +
    labs(
      x        = expression(S[1] / bar(g)),
      y        = expression(S[2] / bar(g)),
      title    = disease,
      subtitle = psi_str
    )

  print(fig_hex_cap)

  if (exists("save_fig")) {
    slug <- gsub("[^A-Za-z0-9]", "_", tolower(disease))
    save_fig(fig_hex_cap, sprintf("treeviz_hex_%s_capped", slug), width = 6, height = 6)
    cat(sprintf("  Saved treeviz_hex_%s_capped\n", slug))
  }
}

# ==============================================================================
# Standalone ECDF of |S1 - S2| / gbar
# ==============================================================================

fig_ecdf <- ggplot(ecdf_data, aes(x = abs_diff, color = group, linetype = group)) +
  stat_ecdf(linewidth = 0.8) +
  scale_color_manual(
    values = c("Random" = "grey50", "Within-cluster" = "darkorange"),
    name   = NULL
  ) +
  scale_linetype_manual(
    values = c("Random" = "dashed", "Within-cluster" = "solid"),
    name   = NULL
  ) +
  facet_wrap(~ fac_label, nrow = 2) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", strip.text = element_text(size = 9)) +
  labs(
    x     = expression("|" * S[1] - S[2] * "| / " * bar(g)),
    y     = "Cumulative proportion"
  )

print(fig_ecdf)

if (exists("save_fig")) {
  save_fig(fig_ecdf, "treeviz_ecdf", width = 12, height = 5)
  cat("Saved treeviz_ecdf figure.\n")
}

# ==============================================================================
# ICC analysis: observed vs. model-predicted intraclass correlation
#
# ICC = (MS_between - MS_within) / (MS_between + (n0 - 1) * MS_within)
# where n0 is the harmonic-ish effective cluster size for unequal group sizes.
# Negative ICC = within-cluster variance exceeds between-cluster (no clustering).
# Bootstrap resamples clusters (the independent unit), not individual SIs.
# ==============================================================================

compute_icc <- function(clusters) {
  clusters  <- clusters[sapply(clusters, length) >= 2]
  if (length(clusters) < 2) return(NA_real_)
  all_si    <- unlist(clusters)
  grand_mean <- mean(all_si)
  total_var  <- var(all_si)
  if (total_var == 0) return(NA_real_)
  # Mean of (Si - mu)(Sj - mu) over all within-cluster ordered pairs i != j
  cross_prods <- unlist(lapply(clusters, function(cl) {
    centered <- cl - grand_mean
    outer(centered, centered)[lower.tri(matrix(0, length(cl), length(cl)))]
  }))
  mean(cross_prods) / total_var
}

icc_model_pred <- function(psi, disease) {
  pp         <- pathogen_params[[disease]]
  sigma2_gi  <- pp$alpha_gi / pp$beta_gi^2
  sigma2_inc <- pp$a_obs    / pp$b_obs^2
  ((1 - psi) * sigma2_gi + sigma2_inc) / (sigma2_gi + 2 * sigma2_inc)
}

set.seed(123)
B <- 1000

icc_results <- map_dfr(names(all_disease_clusters), function(disease) {
  multi <- all_disease_clusters[[disease]]$clusters
  multi <- multi[sapply(multi, length) >= 2]
  if (length(multi) < 2) return(tibble())

  obs_icc   <- compute_icc(multi)
  boot_iccs <- replicate(B, compute_icc(sample(multi, length(multi), replace = TRUE)))
  psi_est   <- if (!is.null(psi_modes) && disease %in% names(psi_modes))
                 psi_modes[[disease]] else NA_real_

  tibble(
    disease  = disease,
    obs_icc  = obs_icc,
    ci_lo    = quantile(boot_iccs, 0.025, na.rm = TRUE),
    ci_hi    = quantile(boot_iccs, 0.975, na.rm = TRUE),
    pred_icc = if (!is.na(psi_est)) icc_model_pred(psi_est, disease) else NA_real_,
    psi_est  = psi_est
  )
})

icc_results <- icc_results %>%
  arrange(psi_est) %>%
  mutate(label = factor(sprintf("%s\n(ψ=%.2f)", disease, psi_est),
                        levels = sprintf("%s\n(ψ=%.2f)", disease, psi_est)))

fig_icc <- ggplot(icc_results, aes(x = label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.25, color = "darkorange", linewidth = 0.8) +
  geom_point(aes(y = obs_icc), color = "darkorange", size = 3) +
  geom_point(aes(y = pred_icc), color = "steelblue", size = 3, shape = 18) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9)) +
  labs(
    x        = NULL,
    y        = "Intraclass correlation (ICC)",
    title    = "Within-cluster SI correlation: observed vs. model-predicted",
    subtitle = "Orange circle: observed ICC with 95% bootstrap CI  |  Blue diamond: model-predicted ICC at estimated ψ"
  )

print(fig_icc)

if (exists("save_fig")) {
  save_fig(fig_icc, "treeviz_icc", width = 10, height = 5)
  cat("Saved treeviz_icc figure.\n")
}
