# ==============================================================================
# isolation_gi_truncation.R
#
# Generation interval distribution under detect-and-isolate, as a function
# of psi (individual infectiousness profile burstiness).
#
# Detection model
# ---------------
# Symptom onset occurs at individual infectiousness peak + Gaussian noise:
#   t_sx = l_i + m_psi + D0,   D0 ~ N(0, sigma_d^2)
# where m_psi = mode of Gamma(psi*alpha, beta) = max(0, (psi*alpha-1)/beta).
# Perfect isolation at symptom onset.
#
# Key cancellation (from math_output.html derivation): l_i drops out of the
# survival condition, so a transmission at eps survives iff eps < m_psi + D0.
#
# GI parameters
# -------------
# SARS-CoV-2 omicron: sourced from code/parameters.R (parslist entry "omicron").
#   Lancet Europe (2022). doi:10.1016/j.lanepe.2022.100446
#   Gamma shape = 2.39, scale = 2.95  →  mean ≈ 7.05 d, SD ≈ 4.56 d
# ==============================================================================

library(tidyverse)

source("code/utils.R")
source("code/parameters.R")

set.seed(42)

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------

# Pull omicron GI parameters from parslist
omicron  <- Filter(function(p) p$pathogen == "omicron", parslist)[[1]]
alpha_gi <- omicron$alpha
beta_gi  <- omicron$beta
gi_mean  <- omicron$Tgen
gi_sd    <- sqrt(omicron$Tvar)

sigma_d  <- 1.0          # SD of symptom onset around infectiousness peak (days)

psi_vals <- c(0.1, 0.5, 0.9)

N        <- 2e6          # Monte Carlo draws per psi
tau_max  <- 14           # x-axis limit (days)
tau_grid <- seq(0, tau_max, length.out = 500)

# ------------------------------------------------------------------------------
# Simulation function
# ------------------------------------------------------------------------------

# tau = l + eps, where
#   l   ~ Gamma((1-psi)*alpha_gi, rate=beta_gi)   [latent period, shared]
#   eps ~ Gamma(psi*alpha_gi,     rate=beta_gi)   [transmission timing]
#
# Isolation cutoff (in eps units): c = m_psi + D0, D0 ~ N(0, sigma_d^2)
# Transmission survives iff: eps < c  (and c > 0)

sim_gi_post_isolation <- function(psi, alpha_gi, beta_gi, sigma_d, N) {
  alpha_l <- (1 - psi) * alpha_gi
  alpha_e <- psi       * alpha_gi
  m_psi   <- if (alpha_e >= 1) (alpha_e - 1) / beta_gi else 0

  l        <- if (alpha_l > 1e-9) rgamma(N, shape = alpha_l, rate = beta_gi) else rep(0, N)
  eps      <- rgamma(N, shape = alpha_e, rate = beta_gi)
  c_iso    <- m_psi + rnorm(N, mean = 0, sd = sigma_d)
  survives <- eps < c_iso & c_iso > 0

  (l + eps)[survives]
}

# ------------------------------------------------------------------------------
# Build density data frame
# ------------------------------------------------------------------------------

# Natural GI — analytic Gamma density
gi_df <- tibble(
  tau     = tau_grid,
  density = dgamma(tau_grid, shape = alpha_gi, rate = beta_gi),
  label   = "No intervention"
)

# Post-isolation GI — kernel density estimate from Monte Carlo samples
for (psi in psi_vals) {
  tau_star    <- sim_gi_post_isolation(psi, alpha_gi, beta_gi, sigma_d, N)
  dens        <- density(tau_star, from = 0, to = tau_max + 2, n = 1024, bw = "SJ", adjust = 2)
  dens_interp <- pmax(0, approx(dens$x, dens$y, xout = tau_grid, rule = 2)$y)

  gi_df <- bind_rows(gi_df, tibble(
    tau     = tau_grid,
    density = dens_interp,
    label   = paste0("psi = ", psi)
  ))
}

gi_df <- gi_df %>%
  mutate(label = factor(label, levels = c(
    "No intervention", paste0("psi = ", psi_vals)
  )))

# ------------------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------------------

pal <- c(
  "No intervention" = "grey40",
  "psi = 0.1"       = "#9ecae1",
  "psi = 0.5"       = "#3182bd",
  "psi = 0.9"       = "#08306b"
)
lty <- c(
  "No intervention" = "dashed",
  "psi = 0.1"       = "solid",
  "psi = 0.5"       = "solid",
  "psi = 0.9"       = "solid"
)

fig_gi_trunc <- ggplot(gi_df, aes(x = tau, y = density,
                                   color = label, linetype = label)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = pal, name = NULL) +
  scale_linetype_manual(values = lty, name = NULL) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
  coord_cartesian(xlim = c(0, tau_max)) +
  theme_classic(base_size = 13) +
  theme(
    legend.position   = c(0.78, 0.72),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.width  = unit(1.2, "cm")
  ) +
  labs(
    x        = "Generation interval (days)",
    y        = "Density",
    title    = "Effect of detect-and-isolate on generation interval distribution",
    subtitle = sprintf(
      "SARS-CoV-2 omicron (mean=%.1f d, SD=%.1f d)  |  symptom onset SD = %.0f d from peak",
      gi_mean, gi_sd, sigma_d
    )
  )

print(fig_gi_trunc)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

fig_name <- "isolation_gi_truncation"

if (exists("save_fig")) {
  save_fig(fig_gi_trunc, fig_name, width = 7, height = 5)
  cat(sprintf("Saved %s\n", fig_name))
} else {
  ggsave(file.path("figures", paste0(fig_name, ".pdf")),
         fig_gi_trunc, width = 7, height = 5)
  ggsave(file.path("figures", paste0(fig_name, ".png")),
         fig_gi_trunc, width = 7, height = 5, dpi = 150)
  cat(sprintf("Saved figures/%s\n", fig_name))
}
