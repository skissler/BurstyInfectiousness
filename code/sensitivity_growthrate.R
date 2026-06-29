library(tidyverse)
library(parallel)

source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

# ==============================================================================
# SENSITIVITY ANALYSIS: growth rate inference under alternative burst models
# (Omicron only)
#
# Replicates the three figures from growthrate.R under:
#   - Type-I Gamma  (original): gen_inf_attempts_gamma
#   - Log-normal:               gen_inf_attempts_lognormal
#   - Type-II Gamma:            gen_inf_attempts_gamma2
#
# psi in {0, 0.5, 1}. Produces 9 figures: 3 burst models x 3 figure types.
#
# Theoretical r:
#   Gamma and Type-II Gamma share the same Gamma(alpha, beta) marginal GI, so
#   r is the same (pars$r). Log-normal has a different marginal so r is solved
#   numerically from the Euler-Lotka equation.
# ==============================================================================

psivals              <- c(0, 0.5, 1)
min_growth_threshold <- 100
growth_window_days   <- 7

n_cores <- getOption("mc.cores",
                     as.integer(Sys.getenv("MC_CORES",
                                           max(1L, detectCores() - 1L))))

# --- Omicron parameters -------------------------------------------------------
omicron  <- parslist[[2]]
pathogen <- omicron$pathogen
Tgen     <- omicron$Tgen
Tvar     <- omicron$Tvar
alpha    <- omicron$alpha
beta     <- omicron$beta
R0       <- omicron$R0

cat(sprintf("\n=== sensitivity_growthrate.R — %s (R0=%g, alpha=%.2f, beta=%.3f) ===\n",
            pathogen, R0, alpha, beta))

# --- Shared plotting aesthetics -----------------------------------------------
psi_colors <- setNames(c("red", "blue", "black"), as.character(psivals))

# --- Log-normal tau parameters (moment-matched to Gamma(alpha, beta)) ---------
sigma2_tau_ln <- log(1 + 1 / alpha)
sigma_tau_ln  <- sqrt(sigma2_tau_ln)
mu_tau_ln     <- log(Tgen) - sigma2_tau_ln / 2

# --- Theoretical growth rates -------------------------------------------------
# Both Gamma models have marginal GI ~ Gamma(alpha, beta), so r = pars$r.
r_gamma <- omicron$r

# Log-normal: solve Euler-Lotka numerically (lognormal tail differs from Gamma).
r_lognormal <- uniroot(function(r) {
    R0 * integrate(function(tau)
                       exp(-r * tau) * dlnorm(tau, mu_tau_ln, sigma_tau_ln),
                   lower = 0, upper = Inf)$value - 1
}, lower = 0, upper = 10)$root

cat(sprintf("  Theoretical r: Gamma/Gamma2 = %.4f, Log-normal = %.4f\n",
            r_gamma, r_lognormal))

# --- max_cases (using Gamma r as reference for the buffer calculation) --------
max_cases <- ceiling(min_growth_threshold *
                     exp(r_gamma * (growth_window_days + 1)) * 3)

# --- Precompute lognormal fl_obj for psi = 0.5 --------------------------------
fl_obj_ln_05 <- precompute_fl_lognormal(0.5,
                                        mu_tau    = mu_tau_ln,
                                        sigma_tau = sigma_tau_ln)

# --- Burst model specs --------------------------------------------------------
models <- c("gamma", "lognormal", "gamma2")

model_labels <- c(gamma     = "Type-I Gamma",
                  lognormal = "Log-normal",
                  gamma2    = "Type-II Gamma")

r_by_model <- c(gamma     = r_gamma,
                lognormal = r_lognormal,
                gamma2    = r_gamma)

# --- Cache helpers ------------------------------------------------------------
sensitivity_gr_cache_summary <- function(model, nsim, max_cases) {
    file.path("output",
              sprintf("sensitivity_gr_summary_%s_%s_n%d_m%d.csv",
                      model, pathogen, nsim, max_cases))
}

sensitivity_gr_cache_plot <- function(model, nsim, max_cases) {
    file.path("output",
              sprintf("sensitivity_gr_plot_%s_%s_n%d_m%d.csv",
                      model, pathogen, nsim, max_cases))
}

# ==============================================================================
# Main loop: 3 burst models
# ==============================================================================

for (mod in models) {

    label        <- model_labels[mod]
    r_malthusian <- r_by_model[mod]

    cat(sprintf("\n  --- %s (r_theory = %.4f) ---\n", label, r_malthusian))

    f_sum  <- sensitivity_gr_cache_summary(mod, nsim_small, max_cases)
    f_plot <- sensitivity_gr_cache_plot(mod, nsim_small, max_cases)

    # ------------------------------------------------------------------
    # Simulate or load from cache
    # ------------------------------------------------------------------
    if (file.exists(f_sum) && file.exists(f_plot)) {

        infpop_summary_df <- read_csv(f_sum, show_col_types = FALSE)
        infpop_plot_df    <- read_csv(f_plot, show_col_types = FALSE)
        cat(sprintf("    Loading from cache: %s\n", basename(f_sum)))

    } else {

        task_grid <- expand_grid(sim = seq_len(nsim_small), psi = psivals)

        RNGkind("L'Ecuyer-CMRG")
        set.seed(sum(utf8ToInt(mod)) + sum(utf8ToInt(pathogen)) + 42L)

        run_one <- function(i) {
            sim_i <- task_grid$sim[i]
            psi_i <- task_grid$psi[i]

            gen <- switch(mod,
                gamma     = gen_inf_attempts_gamma(Tgen, R0, alpha, psi_i),
                lognormal = gen_inf_attempts_lognormal(
                                Tgen, Tvar, R0, psi_i,
                                fl_obj = if (psi_i > 0 && psi_i < 1) fl_obj_ln_05
                                         else NULL),
                gamma2    = gen_inf_attempts_gamma2(Tgen, R0, alpha, beta, psi_i)
            )

            infection_times <- sim_infinite_pop(max_cases        = max_cases,
                                                gen_inf_attempts = gen)
            n_infected <- length(infection_times)

            growthrate <- suppressWarnings(
                compute_growth_rate(infection_times,
                                    min_growth_threshold,
                                    growth_window_days)
            )

            summary_row <- tibble(sim        = sim_i,
                                  psi        = psi_i,
                                  n_infected = n_infected,
                                  growthrate = growthrate)

            plot_row <- if (sim_i <= max_plot_sims && n_infected > 0) {
                tibble(tinf   = infection_times,
                       cuminf = seq_along(infection_times),
                       sim    = sim_i,
                       psi    = psi_i)
            } else NULL

            list(summary = summary_row, plot = plot_row)
        }

        results <- mclapply(seq_len(nrow(task_grid)), run_one,
                            mc.cores = n_cores, mc.set.seed = TRUE)

        errs <- vapply(results, function(x) inherits(x, "try-error"), logical(1))
        if (any(errs))
            stop(sprintf("%s/%s: %d task(s) failed; first error:\n%s",
                         mod, pathogen, sum(errs),
                         as.character(results[[which(errs)[1]]])))

        infpop_summary_df <- bind_rows(lapply(results, `[[`, "summary"))
        infpop_plot_df    <- bind_rows(
            Filter(Negate(is.null), lapply(results, `[[`, "plot")))

        write_csv(infpop_summary_df, f_sum)
        write_csv(infpop_plot_df,    f_plot)
        cat(sprintf("    Simulations saved to %s\n", basename(f_sum)))
    }

    infpop_summary_df <- infpop_summary_df %>% mutate(psi = factor(psi))
    infpop_plot_df    <- infpop_plot_df    %>% mutate(psi = factor(psi))

    # ------------------------------------------------------------------
    # Figure 1: faceted growth-rate histograms (one panel per psi)
    # ------------------------------------------------------------------
    infpop_growthrate_df <- infpop_summary_df %>% filter(!is.na(growthrate))

    infpop_growthrate_table <- infpop_growthrate_df %>%
        group_by(psi) %>%
        summarise(mean   = mean(growthrate),
                  sd     = sd(growthrate),
                  lwr95  = quantile(growthrate, 0.025),
                  upr95  = quantile(growthrate, 0.975),
                  .groups = "drop")

    cat(sprintf("    %s growth rate summary:\n", label))
    print(infpop_growthrate_table)

    fig_hists <- ggplot(infpop_growthrate_df, aes(x = growthrate)) +
        geom_histogram(aes(y = after_stat(density)), bins = 40,
                       fill = "white", col = "darkgrey") +
        geom_density(adjust = 2) +
        geom_vline(xintercept = r_malthusian,
                   col = "blue", lty = "dashed", linewidth = 0.8) +
        theme_classic() +
        facet_wrap(~psi, nrow = 1) +
        labs(x     = "Empirical growth rate (1/day)",
             y     = "Density",
             title = label)

    save_fig(fig_hists,
             sprintf("sensitivity_gr_hists_%s_%s", mod, pathogen))

    # ------------------------------------------------------------------
    # Figure 2: overlaid growth-rate histograms (psi as colour)
    # ------------------------------------------------------------------
    fig_overlay <- ggplot(infpop_growthrate_df,
                          aes(x = growthrate, fill = psi, col = psi)) +
        geom_histogram(aes(y = after_stat(density)), bins = 40,
                       alpha = 0.3, position = "identity") +
        geom_density(adjust = 2, linewidth = 0.8, fill = NA) +
        geom_vline(xintercept = r_malthusian,
                   col = "black", lty = "dashed", linewidth = 0.8) +
        scale_color_manual(values = psi_colors) +
        scale_fill_manual(values  = psi_colors) +
        theme_classic(base_size = 16) +
        labs(x     = "Empirical growth rate (1/day)",
             y     = "Density",
             fill  = expression(psi),
             col   = expression(psi),
             title = label)

    save_fig(fig_overlay,
             sprintf("sensitivity_gr_overlay_%s_%s", mod, pathogen))

    # ------------------------------------------------------------------
    # Figure 3: daily incidence trajectories on log scale
    # ------------------------------------------------------------------
    start_days <- infpop_plot_df %>%
        filter(cuminf == min_growth_threshold) %>%
        transmute(sim, psi, start_day = floor(tinf) + 1L)

    daily_counts <- infpop_plot_df %>%
        mutate(day = floor(tinf)) %>%
        inner_join(start_days, by = c("sim", "psi")) %>%
        filter(day >= start_day, day < start_day + growth_window_days) %>%
        group_by(sim, psi, day) %>%
        summarise(count = n(), .groups = "drop")

    infpop_growth_incidence <- start_days %>%
        group_by(sim, psi) %>%
        reframe(day = seq(start_day, start_day + growth_window_days - 1)) %>%
        left_join(daily_counts, by = c("sim", "psi", "day")) %>%
        replace_na(list(count = 0)) %>%
        group_by(sim, psi) %>%
        mutate(day0 = day - min(day)) %>%
        ungroup()

    refline_df <- infpop_growth_incidence %>%
        group_by(psi) %>%
        summarise(
            intercept = coef(glm(count ~ day0, family = poisson))[1] / log(10),
            slope     = r_malthusian / log(10),
            .groups   = "drop")

    fig_lines <- infpop_growth_incidence %>%
        filter(count > 0) %>%
        ggplot(aes(x = day0, y = count, group = factor(sim))) +
            geom_line(alpha = 0.1, linewidth = 0.3, col = "grey") +
            geom_point(alpha = 0.2, size = 0.3, col = "grey") +
            geom_abline(data        = refline_df,
                        aes(intercept = intercept, slope = slope),
                        col         = "blue",
                        linewidth   = 0.8,
                        lty         = "dashed",
                        inherit.aes = FALSE) +
            scale_y_log10() +
            theme_classic(base_size = 14) +
            theme(strip.background = element_blank()) +
            facet_wrap(~psi, nrow = 1,
                       labeller = as_labeller(function(x) paste0("psi == ", x), label_parsed)) +
            labs(x     = sprintf("Days since case %d", min_growth_threshold),
                 y     = "Daily incidence",
                 title = label)

    save_fig(fig_lines,
             sprintf("sensitivity_gr_lines_%s_%s", mod, pathogen))

    cat(sprintf("    3 figures saved.\n"))

} # end burst model loop

cat("\nSensitivity growth rate analysis complete. 9 figures saved.\n")
