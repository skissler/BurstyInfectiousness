# ==============================================================================
# isolation_overdispersion.R
#
# Distribution of the individual reproduction number
#   nu_i = R0 * pgamma(max(0, c); psi*alpha, beta)
# under symptom-based detect-and-isolate, as a function of psi.
#
# Detection model (same as isolation_gi_truncation.R)
# ---------------
# Symptom onset: t_sx = l_i + m_psi + D0,  D0 ~ N(0, sigma_d^2)
# Key cancellation: l_i drops out, so survival fraction depends only on D0.
#   c = m_psi + D0;  nu_i = R0 * pgamma(max(0, c); psi*alpha, beta)
#
# Interpretation by psi
# ---------------------
# psi -> 0 (spike profile): m_psi -> 0, so c = D0 ~ N(0,1).
#   ~50% of individuals have c <= 0 -> nu_i = 0 (isolated before transmitting).
#   ~50% have c > 0 -> nu_i ~ R0 (nearly all mass survives the tiny spike).
#   Result: bimodal distribution near {0, R0} -> high overdispersion.
#
# psi -> 1 (smooth profile): m_psi large, c rarely <= 0.
#   nu_i = R0 * pgamma(c; alpha, beta), c ~ N(m_psi, sigma_d).
#   Result: bell-shaped continuous distribution below R0 -> low overdispersion.
#
# GI / R0 parameters: SARS-CoV-2 omicron from code/parameters.R
# ==============================================================================

library(tidyverse)

source("code/utils.R")
source("code/parameters.R")

set.seed(42)

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------

omicron  <- Filter(function(p) p$pathogen == "omicron", parslist)[[1]]
alpha_gi <- omicron$alpha
beta_gi  <- omicron$beta
R0       <- omicron$R0

sigma_d  <- 1.0          # SD of symptom onset around infectiousness peak (days)
psi_vals <- c(0.01, 0.5, 0.99)
N        <- 2e6
nu_grid  <- seq(0, R0, length.out = 500)

# ------------------------------------------------------------------------------
# Simulation
# ------------------------------------------------------------------------------

simulate_nu <- function(psi, alpha_gi, beta_gi, R0, sigma_d, N) {
  alpha_e <- psi * alpha_gi
  m_psi   <- if (alpha_e >= 1) (alpha_e - 1) / beta_gi else 0
  D0      <- rnorm(N, 0, sigma_d)
  R0 * pgamma(pmax(0, m_psi + D0), shape = alpha_e, rate = beta_gi)
}

# ------------------------------------------------------------------------------
# Build density data frame
# ------------------------------------------------------------------------------

# Use nrd0 (Silverman) bandwidth: more stable than SJ for the bimodal
# psi=0.1 distribution which has a genuine 50% point mass at nu=0.

nu_df <- map_dfr(psi_vals, function(psi) {
  nu_i        <- simulate_nu(psi, alpha_gi, beta_gi, R0, sigma_d, N)
  dens        <- density(nu_i, from = 0, to = R0, n = 1024,
                         bw = "nrd0", adjust = 1)
  dens_interp <- pmax(0, approx(dens$x, dens$y, xout = nu_grid, rule = 2)$y)
  tibble(nu = nu_grid, density = dens_interp, label = paste0("psi = ", psi))
}) %>%
  mutate(label = factor(label, levels = paste0("psi = ", psi_vals)))

# ------------------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------------------

pal <- setNames(
  colorRampPalette(c("#9ecae1", "#08306b"))(length(psi_vals)),
  paste0("psi = ", psi_vals)
)

# Arrow height for the "no intervention" delta function at R0:
# set to 80% of the tallest KDE peak so it reads as large but doesn't
# force the y-axis to expand further.
arrow_height <- max(nu_df$density) * 0.8

fig_nu <- ggplot(nu_df, aes(x = nu, y = density, color = label)) +
  geom_line(linewidth = 0.9) +
  # "No intervention": delta function at R0, represented as an upward arrow
  annotate("segment",
           x = R0, xend = R0, y = 0, yend = arrow_height,
           color = "grey40", linewidth = 0.9,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  annotate("text",
           x = R0 - 0.15, y = arrow_height * 0.55,
           label = "No\nintervention",
           hjust = 1, color = "grey40", size = 3.5) +
  scale_color_manual(values = pal, name = NULL) +
  scale_x_continuous(
    breaks = 0:R0,
    labels = parse(text = c(as.character(0:(R0 - 1)),
                             sprintf("R[0] == %g", R0))),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, NA)
  ) +
  coord_cartesian(xlim = c(0, R0 + 0.2)) +
  theme_classic(base_size = 13) +
  theme(
    legend.position   = c(0.55, 0.80),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.width  = unit(1.2, "cm")
  ) +
  labs(
    x        = expression(nu[i] ~ "(individual reproduction number)"),
    y        = "Density",
    title    = "Detect-and-isolate introduces overdispersion in individual reproduction numbers",
    subtitle = sprintf(
      "SARS-CoV-2 omicron (R₀ = %g)  |  symptom onset SD = %.0f d from peak infectiousness",
      R0, sigma_d
    )
  )

print(fig_nu)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

fig_name <- "isolation_overdispersion"

if (exists("save_fig")) {
  save_fig(fig_nu, fig_name, width = 7, height = 5)
  cat(sprintf("Saved %s\n", fig_name))
} else {
  ggsave(file.path("figures", paste0(fig_name, ".pdf")),
         fig_nu, width = 7, height = 5)
  ggsave(file.path("figures", paste0(fig_name, ".png")),
         fig_nu, width = 7, height = 5, dpi = 150)
  cat(sprintf("Saved figures/%s\n", fig_name))
}

# ------------------------------------------------------------------------------
# Histogram version — shows point masses at 0 and R0 honestly without KDE
# smoothing.  Uses a fresh sample; faceted so each distribution is legible.
# ------------------------------------------------------------------------------

set.seed(42)
N_hist <- 2e5

raw_df <- map_dfr(psi_vals, function(psi) {
  tibble(nu    = simulate_nu(psi, alpha_gi, beta_gi, R0, sigma_d, N_hist),
         label = paste0("psi = ", psi))
}) %>%
  mutate(label = factor(label, levels = paste0("psi = ", psi_vals)))

mean_df <- raw_df %>%
  group_by(label) %>%
  summarise(mean_nu = mean(nu), .groups = "drop")

fig_nu_hist <- ggplot(raw_df, aes(x = nu, fill = label)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.05, boundary = 0, color = NA) +
  geom_vline(data = mean_df, aes(xintercept = mean_nu),
             color = "black", linewidth = 0.8, linetype = "dashed") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_x_continuous(
    breaks = 0:R0,
    labels = parse(text = c(as.character(0:(R0 - 1)),
                             sprintf("R[0] == %g", R0))),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(xlim = c(0, R0 + 0.05)) +
  facet_wrap(~ label, ncol = 1, scales = "free_y",
             strip.position = "right") +
  theme_classic(base_size = 13) +
  theme(
    legend.position  = "none",
    strip.background = element_blank(),
    strip.text       = element_text(size = 11)
  ) +
  labs(
    x        = expression(nu[i] ~ "(individual reproduction number)"),
    y        = "Density",
    title    = "Detect-and-isolate introduces overdispersion in individual reproduction numbers",
    subtitle = sprintf(
      "SARS-CoV-2 omicron (R₀ = %g)  |  symptom onset SD = %.0f d from peak infectiousness",
      R0, sigma_d
    )
  )

print(fig_nu_hist)

fig_name_hist <- "isolation_overdispersion_hist"

if (exists("save_fig")) {
  save_fig(fig_nu_hist, fig_name_hist, width = 7, height = 7)
  cat(sprintf("Saved %s\n", fig_name_hist))
} else {
  ggsave(file.path("figures", paste0(fig_name_hist, ".pdf")),
         fig_nu_hist, width = 7, height = 7)
  ggsave(file.path("figures", paste0(fig_name_hist, ".png")),
         fig_nu_hist, width = 7, height = 7, dpi = 150)
  cat(sprintf("Saved figures/%s\n", fig_name_hist))
}
