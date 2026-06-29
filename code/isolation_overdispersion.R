library(tidyverse)

source("code/utils.R")
source("code/parameters.R")

set.seed(42)

# ==============================================================================
# Define parameters
# ==============================================================================

omicron  <- Filter(function(p) p$pathogen == "omicron", parslist)[[1]]
alpha    <- omicron$alpha
beta     <- omicron$beta
R0       <- omicron$R0

sigma_sym  <- 0.5    # SD of symptom onset around infectiousness peak (days)
# psi_vals   <- c(0.01, 0.5, 0.99)
psi_vals   <- c(0, 1)
N          <- 2e5

# ==============================================================================
# Simulate individual reproduction numbers 
# ==============================================================================

simulate_nu <- function(psi, alpha, beta, R0, sigma_sym, N) {
  alpha_e <- psi * alpha
  m_psi   <- if (alpha_e >= 1) (alpha_e - 1) / beta else 0
  D0      <- rnorm(N, 0, sigma_sym)
  cutoff  <- pmax(0, m_psi + D0)
  surv    <- if (alpha_e > 1e-9) pgamma(cutoff, shape = alpha_e, rate = beta) else as.numeric(cutoff > 0)
  R0 * surv
}

# ------------------------------------------------------------------------------
# Plot histogram
# ------------------------------------------------------------------------------

pal <- setNames(
  colorRampPalette(c("#9ecae1", "#08306b"))(length(psi_vals)),
  paste0("psi = ", psi_vals)
)

raw_df <- map_dfr(psi_vals, function(psi) {
  tibble(nu    = simulate_nu(psi, alpha, beta, R0, sigma_sym, N),
         label = paste0("psi = ", psi))
}) %>%
  mutate(label = factor(label, levels = paste0("psi = ", psi_vals)))

mean_df <- raw_df %>%
  group_by(label) %>%
  summarise(mean_nu = mean(nu), .groups = "drop")

int_breaks <- seq(0, floor(R0))
int_breaks <- int_breaks[int_breaks != R0]
x_breaks   <- c(int_breaks, R0)
x_labels   <- parse(text = c(as.character(int_breaks),
                             sprintf("R[0] == %g", R0)))

fig_nu <- ggplot(raw_df, aes(x = nu, fill = label)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.05, boundary = 0, color = NA) +
  geom_vline(data = mean_df, aes(xintercept = mean_nu),
             color = "black", linewidth = 0.8, linetype = "dashed") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
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
      R0, sigma_sym
    )
  )

print(fig_nu)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

fig_name <- "isolation_overdispersion"

save_fig(fig_nu, fig_name, width = 7, height = 7)
cat(sprintf("Saved %s\n", fig_name))
