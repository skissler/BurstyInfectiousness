library(tidyverse)

source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

# ==============================================================================
# SENSITIVITY ANALYSIS: marginal GI convergence checks (all pathogens)
#
# Verifies that the Type-II Gamma and Log-normal burst models produce the
# correct marginal generation interval distribution across a range of psi
# values.  For each model and psi, n_samples draws of tau = l + eps are
# compared to the target distribution via:
#   - K-S test (D statistic + p-value)
#   - Density histogram with theoretical curve overlaid
#   - QQ plot against theoretical quantiles
#
# Target distributions:
#   - Type-I Gamma:  Gamma(alpha, beta)   [exact by additive-Gamma property]
#   - Type-II Gamma: Gamma(alpha, beta)   [exact by Lévy decomposition]
#   - Log-normal:    LogNormal(mu_tau, sigma_tau)  [approximate via FW]
#
# Note on K-S tests: with n_samples = 50000 the test has high power. A large D
# statistic indicates a meaningful discrepancy; a small D with a nominally
# significant p-value is an artefact of sample size.
#
# Produces 6 figures per pathogen (3 burst models x 2 figure types).
# ==============================================================================

psi_check_vals    <- c(0, 0.25, 0.5, 0.75, 1)
n_samples         <- 50000L
psi_factor_levels <- paste0("psi == ", psi_check_vals)

models <- c("gamma", "lognormal", "gamma2")

model_labels <- c(gamma     = "Type-I Gamma",
                  lognormal = "Log-normal",
                  gamma2    = "Type-II Gamma")

target_labels <- c(gamma     = "Gamma(alpha, beta)",
                   lognormal = "LogNormal(mu_tau, sigma_tau)",
                   gamma2    = "Gamma(alpha, beta)")

target_xlabs <- list(
    gamma     = expression("Theoretical quantiles —" ~ "Gamma(" * alpha * ", " * beta * ")"),
    lognormal = expression("Theoretical quantiles —" ~ "LogNormal(" * mu[tau] * ", " * sigma[tau] * ")"),
    gamma2    = expression("Theoretical quantiles —" ~ "Gamma(" * alpha * ", " * beta * ")")
)

hist_title_suffixes <- list(
    gamma     = quote("Simulated histogram vs." ~ "Gamma(" * alpha * ", " * beta * ")"),
    lognormal = quote("Simulated histogram vs." ~ "LogNormal(" * mu[tau] * ", " * sigma[tau] * ")"),
    gamma2    = quote("Simulated histogram vs." ~ "Gamma(" * alpha * ", " * beta * ")")
)

qq_title_suffixes <- list(
    gamma     = quote("QQ plot vs." ~ "Gamma(" * alpha * ", " * beta * ")"),
    lognormal = quote("QQ plot vs." ~ "LogNormal(" * mu[tau] * ", " * sigma[tau] * ")"),
    gamma2    = quote("QQ plot vs." ~ "Gamma(" * alpha * ", " * beta * ")")
)

set.seed(42L)

ks_all <- list()

# ==============================================================================
# Outer loop: pathogens
# ==============================================================================

for (idx_pathogen in seq_along(parslist)) {

    pars     <- parslist[[idx_pathogen]]
    pathogen <- pars$pathogen
    pathogen_label <- switch(pathogen,
        influenza = "Influenza",
        omicron   = "SARS-CoV-2 Omicron",
        measles   = "Measles",
        pathogen
    )
    Tgen     <- pars$Tgen
    alpha    <- pars$alpha
    beta     <- pars$beta

    cat(sprintf("\n=== %s (alpha=%.2f, beta=%.3f) ===\n", pathogen, alpha, beta))

    # --- Log-normal tau parameters (pathogen-specific) ------------------------
    sigma2_tau_ln <- log(1 + 1 / alpha)
    sigma_tau_ln  <- sqrt(sigma2_tau_ln)
    mu_tau_ln     <- log(Tgen) - sigma2_tau_ln / 2

    # Precompute fl_objs for all intermediate psi values
    fl_objs_ln <- list()
    for (psi in psi_check_vals[psi_check_vals > 0 & psi_check_vals < 1]) {
        fl_objs_ln[[as.character(psi)]] <- precompute_fl_lognormal(
            psi, mu_tau = mu_tau_ln, sigma_tau = sigma_tau_ln)
    }

    # ---- Sampler -------------------------------------------------------------
    # Captures alpha, beta, mu_tau_ln, sigma_tau_ln from the current pathogen.

    sample_tau <- function(n, model, psi) {
        if (model == "gamma") {
            if (psi < 1e-8 || psi > 1 - 1e-8) return(rgamma(n, alpha, beta))
            rgamma(n, alpha * (1 - psi), beta) + rgamma(n, alpha * psi, beta)

        } else if (model == "gamma2") {
            if (psi < 1e-8 || psi > 1 - 1e-8) return(rgamma(n, alpha, beta))
            rL_gamma2(n, alpha, beta, psi) + rgamma(n, alpha, beta / sqrt(psi))

        } else if (model == "lognormal") {
            if (psi < 1e-8 || psi > 1 - 1e-8) return(rlnorm(n, mu_tau_ln, sigma_tau_ln))
            fl_obj <- fl_objs_ln[[as.character(psi)]]
            sample_l_lognormal(n, fl_obj) +
                rlnorm(n, mu_tau_ln + 0.5 * log(psi), sigma_tau_ln)
        }
    }

    # ---- Theoretical CDF / PDF / quantile function ---------------------------

    theo_cdf <- function(x, model) {
        if (model %in% c("gamma", "gamma2")) pgamma(x, alpha, beta)
        else                                  plnorm(x, mu_tau_ln, sigma_tau_ln)
    }

    theo_pdf <- function(x, model) {
        if (model %in% c("gamma", "gamma2")) dgamma(x, alpha, beta)
        else                                  dlnorm(x, mu_tau_ln, sigma_tau_ln)
    }

    theo_quant <- function(p, model) {
        if (model %in% c("gamma", "gamma2")) qgamma(p, alpha, beta)
        else                                  qlnorm(p, mu_tau_ln, sigma_tau_ln)
    }

    # ==========================================================================
    # Inner loop: burst models
    # ==========================================================================

    for (mod in models) {

        label        <- model_labels[mod]
        target_label <- target_labels[mod]

        cat(sprintf("\n  --- %s ---\n", label))

        # Draw samples for all psi values
        all_tau <- map_dfr(psi_check_vals, function(psi) {
            tibble(psi       = psi,
                   psi_label = factor(paste0("psi == ", psi), levels = psi_factor_levels),
                   tau       = sample_tau(n_samples, mod, psi))
        })

        # --- K-S tests --------------------------------------------------------
        ks_rows <- map_dfr(psi_check_vals, function(psi) {
            tau  <- all_tau$tau[all_tau$psi == psi]
            test <- ks.test(tau, function(x) theo_cdf(x, mod))
            tibble(pathogen = pathogen,
                   model    = label,
                   psi      = psi,
                   D        = round(test$statistic, 5),
                   p_value  = signif(test$p.value, 3))
        })
        ks_all[[paste(pathogen, mod, sep = "_")]] <- ks_rows

        cat(sprintf("    K-S results for %s:\n", label))
        print(as.data.frame(ks_rows), row.names = FALSE)

        # --- Figure 1: density histograms -------------------------------------
        tau_upper  <- quantile(all_tau$tau, 0.999)
        tau_grid   <- seq(1e-3, tau_upper, length.out = 500)
        theo_curve <- expand_grid(
            tibble(tau = tau_grid, density = theo_pdf(tau_grid, mod)),
            psi_label = factor(psi_factor_levels, levels = psi_factor_levels)
        )

        fig_hist <- ggplot(all_tau, aes(x = tau)) +
            geom_histogram(aes(y = after_stat(density)),
                           bins = 60, fill = "grey70", col = "white", linewidth = 0.2) +
            geom_line(data        = theo_curve,
                      aes(x = tau, y = density),
                      col         = "firebrick",
                      linewidth   = 0.8,
                      inherit.aes = FALSE) +
            facet_wrap(~ psi_label, nrow = 1,
                       labeller = label_parsed) +
            scale_x_continuous(limits = c(0, tau_upper)) +
            theme_classic(base_size = 14) +
            theme(strip.background = element_blank()) +
            labs(x     = "Generation interval (days)",
                 y     = "Density",
                 title = bquote(atop(.(paste0(label, " - ", pathogen_label)),
                                     .(hist_title_suffixes[[mod]]))))

        save_fig(fig_hist,
                 sprintf("sensitivity_marginal_gi_hist_%s_%s", mod, pathogen),
                 width = 14, height = 4)

        # --- Figure 2: QQ plots -----------------------------------------------
        qq_probs <- seq(0.005, 0.995, length.out = 300)

        qq_df <- map_dfr(psi_check_vals, function(psi) {
            tau <- all_tau$tau[all_tau$psi == psi]
            tibble(
                psi_label   = factor(paste0("psi == ", psi), levels = psi_factor_levels),
                empirical   = quantile(tau, qq_probs),
                theoretical = theo_quant(qq_probs, mod)
            )
        })

        ref_df <- qq_df %>%
            group_by(psi_label) %>%
            summarise(lo = min(theoretical), hi = max(theoretical), .groups = "drop")

        fig_qq <- ggplot(qq_df, aes(x = theoretical, y = empirical)) +
            geom_segment(data        = ref_df,
                         aes(x = lo, xend = hi, y = lo, yend = hi),
                         col         = "firebrick",
                         linewidth   = 0.7,
                         linetype    = "dashed",
                         inherit.aes = FALSE) +
            geom_point(size = 0.8, alpha = 0.5) +
            facet_wrap(~ psi_label, nrow = 1, scales = "free",
                       labeller = label_parsed) +
            theme_classic(base_size = 14) +
            theme(strip.background = element_blank()) +
            labs(x     = target_xlabs[[mod]],
                 y     = "Empirical quantiles (simulated)",
                 title = bquote(atop(.(paste0(label, " - ", pathogen_label)),
                                     .(qq_title_suffixes[[mod]]))))

        save_fig(fig_qq,
                 sprintf("sensitivity_marginal_gi_qq_%s_%s", mod, pathogen),
                 width = 14, height = 4)

        cat(sprintf("    2 figures saved.\n"))

    } # end model loop

} # end pathogen loop

# ==============================================================================
# Combined K-S summary table
# ==============================================================================

cat("\n=== K-S test summary (all pathogens and models) ===\n")
ks_summary <- bind_rows(ks_all)
print(as.data.frame(ks_summary), row.names = FALSE)

cat("\nMarginal GI convergence checks complete.\n")

# ==============================================================================
# Write K-S results as a LaTeX table
# ==============================================================================

ks_wide <- ks_summary %>%
    mutate(model = factor(model,
                          levels = c("Type-I Gamma", "Log-normal", "Type-II Gamma"))) %>%
    select(pathogen, psi, model, D) %>%
    pivot_wider(names_from = model, values_from = D)

D_star <- 1.36 / sqrt(n_samples)   # KS critical value at alpha = 0.05

fmt_D <- function(d) {
    s <- sprintf("%.4f", d)
    if (!is.na(d) && d > D_star) paste0(s, "$^{*}$") else s
}

pathogen_order  <- c("influenza", "omicron", "measles")
pathogen_labels <- c(influenza = "Influenza",
                     omicron   = "SARS-CoV-2 omicron",
                     measles   = "Measles")

tex_lines <- c(
    "% Auto-generated by sensitivity_marginal_gi.R -- do not edit by hand",
    "\\begin{table}[h]",
    "\\centering",
    paste0(
        "\\caption{{\\bf Kolmogorov--Smirnov $D$ statistics for the marginal generation interval ",
        "distribution under alternative burst models.} ",
        "For each pathogen and burst model, $n = 50{,}000$ generation intervals were ",
        "simulated across five values of $\\psi$; $D$ is the maximum absolute deviation ",
        "between the empirical and target cumulative distribution functions ",
        "(target: $\\text{Gamma}(\\alpha, \\beta)$ for the Gamma burst models; ",
        "$\\text{LogNormal}(\\mu_\\tau, \\sigma_\\tau)$ for the log-normal model). ",
        sprintf(
            "Values marked $^*$ are nominally significant at $\\alpha = 0.05$ (critical value $D^* \\approx %.4f$ for $n = %d$). ",
            D_star, n_samples
        ),
        "}"
    ),
    "\\label{tab:ks_marginal_gi}",
    "\\begin{tabular}{llrrr}",
    "{\\bf Pathogen} & {\\bf $\\psi$} & {\\bf Type-I Gamma} & {\\bf Log-normal} & {\\bf Type-II Gamma} \\\\",
    "\\hline"
)

for (p in pathogen_order) {
    sub    <- ks_wide %>% filter(pathogen == p) %>% arrange(psi)
    plabel <- pathogen_labels[p]
    for (j in seq_len(nrow(sub))) {
        row     <- sub[j, ]
        pname   <- if (j == 1) plabel else ""
        psi_str <- as.character(row$psi)
        dg      <- fmt_D(row[["Type-I Gamma"]])
        dln     <- fmt_D(row[["Log-normal"]])
        dg2     <- fmt_D(row[["Type-II Gamma"]])
        tex_lines <- c(tex_lines,
                       sprintf("%s & %s & %s & %s & %s \\\\",
                               pname, psi_str, dg, dln, dg2))
    }
    tex_lines <- c(tex_lines, "\\hline")
}

tex_lines <- c(tex_lines,
               "\\end{tabular}",
               "\\end{table}")

out_path <- file.path("writeup", "v5_nature", "ks_table_marginal_gi.tex")
writeLines(tex_lines, out_path)
cat(sprintf("  K-S table written to %s\n", out_path))
