library(tidyverse)
source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

# ==============================================================================
# SENSITIVITY ANALYSIS: overdispersion heatmaps under alternative burst models
#
# Replicates the periodic-contact and stochastic-contact (Gamma/Poisson)
# overdispersion heatmaps from overdispersion_heatmaps.R for SARS-CoV-2 omicron
# under three burst models:
#   - gamma:      original shared-rate Gamma burst model
#   - lognormal:  Log-normal burst model (Fenton-Wilkinson f_l)
#   - gamma2:     equal-shape Gamma burst model (exact compound-Poisson f_l)
#
# Output: six heatmaps (2 contact models × 3 burst models).
# No theory contours, no extinction overlays.
# Cache files prefixed "sensitivity_od_" to avoid overwriting main analysis.
# ==============================================================================

c_per          <- 7
c_amp_vals_sim <- seq(0, 1, length.out = 25)
psi_vals_sim   <- seq(0, 1, length.out = 25)

lambda_gp      <- 1
k_c_vals_sim   <- exp(seq(log(0.1), log(1000), length.out = 25))

n_index <- 10000
k_cap   <- 50
od_cap  <- 1

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# Omicron only
pars     <- parslist[[which(sapply(parslist, `[[`, "pathogen") == "omicron")]]
pathogen <- pars$pathogen
pathogen_label <- switch(pathogen,
    influenza = "Influenza",
    omicron   = "SARS-CoV-2 Omicron",
    measles   = "Measles",
    pathogen
)
Tgen     <- pars$Tgen
Tvar     <- pars$Tvar
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

cat(sprintf("\n===== SENSITIVITY OD: %s T=%.2f alpha=%.2f R0=%.1f =====\n",
            pathogen, Tgen, alpha, R0))

set.seed(99L)

# ------------------------------------------------------------------------------
# Lognormal burst parameters (derived once)
# ------------------------------------------------------------------------------
sigma_tau_ln <- sqrt(log(1 + 1 / alpha))
mu_tau_ln    <- log(Tgen) - sigma_tau_ln^2 / 2

# Precompute fl_obj for each interior psi value (moment-matched lognormal).
# Outer psi boundary cases (0, 1) are handled analytically in the generators.
fl_objs_ln <- setNames(
    lapply(psi_vals_sim, function(psi) {
        if (psi > 1e-8 && psi < 1 - 1e-8)
            precompute_fl_lognormal(psi, mu_tau_ln, sigma_tau_ln)
        else
            NULL
    }),
    as.character(psi_vals_sim)
)

# ------------------------------------------------------------------------------
# Contact-aware generator factories for alternative burst models
#
# Both use the thinning principle: propose N ~ Poisson(z_max) candidate times
# from the burst structure (l_i + eps_j), then accept each with prob z(t)/z_max.
# This is burst-model-agnostic; only the timing distribution changes.
# ------------------------------------------------------------------------------

# --- Periodic contacts --------------------------------------------------------

make_gfun_lognormal_periodic <- function(z, z_max, psi, fl_obj) {
    mu_eps <- mu_tau_ln + 0.5 * log(max(psi, .Machine$double.xmin))
    if (psi < 1e-8) {
        function(tinf) {
            n <- rpois(1L, z_max)
            if (n == 0L) return(numeric(0))
            times <- tinf + rep(rlnorm(1L, mu_tau_ln, sigma_tau_ln), n)
            times[runif(n) < z(times) / z_max]
        }
    } else if (psi > 1 - 1e-8) {
        function(tinf) {
            n <- rpois(1L, z_max)
            if (n == 0L) return(numeric(0))
            times <- tinf + rlnorm(n, mu_tau_ln, sigma_tau_ln)
            times[runif(n) < z(times) / z_max]
        }
    } else {
        force(fl_obj); force(mu_eps)
        function(tinf) {
            n <- rpois(1L, z_max)
            if (n == 0L) return(numeric(0))
            l_i   <- sample_l_lognormal(1L, fl_obj)
            times <- tinf + l_i + rlnorm(n, mu_eps, sigma_tau_ln)
            times[runif(n) < z(times) / z_max]
        }
    }
}

make_gfun_gamma2_periodic <- function(z, z_max, psi) {
    beta_eps <- beta / sqrt(max(psi, .Machine$double.xmin))
    if (psi < 1e-8) {
        function(tinf) {
            n <- rpois(1L, z_max)
            if (n == 0L) return(numeric(0))
            times <- tinf + rep(rgamma(1L, shape = alpha, rate = beta), n)
            times[runif(n) < z(times) / z_max]
        }
    } else if (psi > 1 - 1e-8) {
        function(tinf) {
            n <- rpois(1L, z_max)
            if (n == 0L) return(numeric(0))
            times <- tinf + rgamma(n, shape = alpha, rate = beta)
            times[runif(n) < z(times) / z_max]
        }
    } else {
        force(beta_eps)
        function(tinf) {
            n <- rpois(1L, z_max)
            if (n == 0L) return(numeric(0))
            l_i   <- rL_gamma2(1L, alpha, beta, psi)
            times <- tinf + l_i + rgamma(n, shape = alpha, rate = beta_eps)
            times[runif(n) < z(times) / z_max]
        }
    }
}

# --- Stochastic (Gamma/Poisson) contacts --------------------------------------

make_gfun_lognormal_gp <- function(k_c, psi, fl_obj) {
    mu_eps        <- mu_tau_ln + 0.5 * log(max(psi, .Machine$double.xmin))
    traj_duration <- 5 * Tgen
    if (psi < 1e-8) {
        function(tinf) {
            n_sw    <- rpois(1L, lambda_gp * traj_duration)
            offsets <- if (n_sw > 0L) sort(runif(n_sw, 0, traj_duration)) else numeric(0)
            levels  <- rgamma(n_sw + 1L, k_c, k_c)
            breaks  <- tinf + c(0, offsets)
            zm      <- R0 * max(levels)
            n       <- rpois(1L, zm)
            if (n == 0L) return(numeric(0))
            times   <- tinf + rep(rlnorm(1L, mu_tau_ln, sigma_tau_ln), n)
            seg     <- pmax(1L, pmin(findInterval(times, breaks), length(levels)))
            times[runif(n) < R0 * levels[seg] / zm]
        }
    } else if (psi > 1 - 1e-8) {
        function(tinf) {
            n_sw    <- rpois(1L, lambda_gp * traj_duration)
            offsets <- if (n_sw > 0L) sort(runif(n_sw, 0, traj_duration)) else numeric(0)
            levels  <- rgamma(n_sw + 1L, k_c, k_c)
            breaks  <- tinf + c(0, offsets)
            zm      <- R0 * max(levels)
            n       <- rpois(1L, zm)
            if (n == 0L) return(numeric(0))
            times   <- tinf + rlnorm(n, mu_tau_ln, sigma_tau_ln)
            seg     <- pmax(1L, pmin(findInterval(times, breaks), length(levels)))
            times[runif(n) < R0 * levels[seg] / zm]
        }
    } else {
        force(fl_obj); force(mu_eps)
        function(tinf) {
            n_sw    <- rpois(1L, lambda_gp * traj_duration)
            offsets <- if (n_sw > 0L) sort(runif(n_sw, 0, traj_duration)) else numeric(0)
            levels  <- rgamma(n_sw + 1L, k_c, k_c)
            breaks  <- tinf + c(0, offsets)
            zm      <- R0 * max(levels)
            n       <- rpois(1L, zm)
            if (n == 0L) return(numeric(0))
            l_i     <- sample_l_lognormal(1L, fl_obj)
            times   <- tinf + l_i + rlnorm(n, mu_eps, sigma_tau_ln)
            seg     <- pmax(1L, pmin(findInterval(times, breaks), length(levels)))
            times[runif(n) < R0 * levels[seg] / zm]
        }
    }
}

make_gfun_gamma2_gp <- function(k_c, psi) {
    beta_eps      <- beta / sqrt(max(psi, .Machine$double.xmin))
    traj_duration <- 5 * Tgen
    if (psi < 1e-8) {
        function(tinf) {
            n_sw    <- rpois(1L, lambda_gp * traj_duration)
            offsets <- if (n_sw > 0L) sort(runif(n_sw, 0, traj_duration)) else numeric(0)
            levels  <- rgamma(n_sw + 1L, k_c, k_c)
            breaks  <- tinf + c(0, offsets)
            zm      <- R0 * max(levels)
            n       <- rpois(1L, zm)
            if (n == 0L) return(numeric(0))
            times   <- tinf + rep(rgamma(1L, shape = alpha, rate = beta), n)
            seg     <- pmax(1L, pmin(findInterval(times, breaks), length(levels)))
            times[runif(n) < R0 * levels[seg] / zm]
        }
    } else if (psi > 1 - 1e-8) {
        function(tinf) {
            n_sw    <- rpois(1L, lambda_gp * traj_duration)
            offsets <- if (n_sw > 0L) sort(runif(n_sw, 0, traj_duration)) else numeric(0)
            levels  <- rgamma(n_sw + 1L, k_c, k_c)
            breaks  <- tinf + c(0, offsets)
            zm      <- R0 * max(levels)
            n       <- rpois(1L, zm)
            if (n == 0L) return(numeric(0))
            times   <- tinf + rgamma(n, shape = alpha, rate = beta)
            seg     <- pmax(1L, pmin(findInterval(times, breaks), length(levels)))
            times[runif(n) < R0 * levels[seg] / zm]
        }
    } else {
        force(beta_eps)
        function(tinf) {
            n_sw    <- rpois(1L, lambda_gp * traj_duration)
            offsets <- if (n_sw > 0L) sort(runif(n_sw, 0, traj_duration)) else numeric(0)
            levels  <- rgamma(n_sw + 1L, k_c, k_c)
            breaks  <- tinf + c(0, offsets)
            zm      <- R0 * max(levels)
            n       <- rpois(1L, zm)
            if (n == 0L) return(numeric(0))
            l_i     <- rL_gamma2(1L, alpha, beta, psi)
            times   <- tinf + l_i + rgamma(n, shape = alpha, rate = beta_eps)
            seg     <- pmax(1L, pmin(findInterval(times, breaks), length(levels)))
            times[runif(n) < R0 * levels[seg] / zm]
        }
    }
}

# ------------------------------------------------------------------------------
# Helper: compute k (negative binomial dispersion) from offspring counts
# ------------------------------------------------------------------------------
compute_k <- function(noffspring) {
    m <- mean(noffspring)
    v <- var(noffspring)
    if (v > m) m^2 / (v - m) else Inf
}

# ==============================================================================
# Periodic contact heatmaps
# ==============================================================================

need_periodic <- expand_grid(psi = psi_vals_sim, c_amp = c_amp_vals_sim)

for (model in c("gamma", "lognormal", "gamma2")) {

    cache_file <- file.path(
        "output",
        sprintf("sensitivity_od_periodic_%s_%s_n%d_cper%g.csv",
                model, pathogen, n_index, c_per)
    )

    sim_grid <- NULL
    if (file.exists(cache_file)) {
        cached  <- read_csv(cache_file, show_col_types = FALSE)
        missing <- anti_join(need_periodic, cached %>% distinct(psi, c_amp),
                             by = c("psi", "c_amp"))
        if (nrow(missing) == 0) {
            cat(sprintf("  [%s periodic] loaded from cache\n", model))
            sim_grid <- cached
        }
    }

    if (is.null(sim_grid)) {
        sim_grid <- need_periodic %>% mutate(k_sim = NA_real_)

        for (idx in seq_len(nrow(sim_grid))) {
            psi   <- sim_grid$psi[idx]
            c_amp <- sim_grid$c_amp[idx]
            z     <- make_contact_fn_periodic(R0, c_amp, c_per)
            z_max <- R0 * (1 + c_amp)
            fl    <- fl_objs_ln[[as.character(psi)]]

            gfun <- switch(model,
                gamma     = gen_inf_attempts_gamma_contacts(
                                Tgen, z, z_max, alpha, psi),
                lognormal = make_gfun_lognormal_periodic(z, z_max, psi, fl),
                gamma2    = make_gfun_gamma2_periodic(z, z_max, psi)
            )

            tinfs       <- c_per * runif(n_index)
            noffspring  <- lengths(lapply(tinfs, gfun))
            sim_grid$k_sim[idx] <- compute_k(noffspring)

            if (idx %% 100 == 0)
                cat(sprintf("  [%s periodic] %d / %d\n", model, idx, nrow(sim_grid)))
        }

        write_csv(sim_grid, cache_file)
        cat(sprintf("  [%s periodic] saved to %s\n", model, cache_file))
    }

    sim_grid <- sim_grid %>%
        mutate(k_capped  = pmin(k_sim, k_cap),
               od_capped = pmin(1 / k_sim, od_cap))

    model_label <- switch(model,
        gamma     = "Gamma (Type I)",
        lognormal = "Log-normal",
        gamma2    = "Gamma (Type II)"
    )

    fig <- sim_grid %>%
        ggplot(aes(x = psi, y = c_amp, fill = k_capped)) +
            geom_tile() +
            scale_fill_viridis_c(option = "inferno",
                                 name   = sprintf("k\n(capped\nat %d)", k_cap),
                                 limits = c(0, k_cap)) +
            theme_classic(base_size = 14) +
            labs(x     = expression(psi),
                 y     = expression("Contact amplitude (" * zeta * ")"),
                 title = sprintf("%s (periodic contacts)", model_label))

    save_fig(fig,
             sprintf("sensitivity_od_periodic_%s_%s", model, pathogen),
             width = 5, height = 4)
    cat(sprintf("  [%s periodic] figure saved\n", model))
}

# ==============================================================================
# Gamma/Poisson (stochastic) contact heatmaps
# ==============================================================================

need_gp <- expand_grid(psi = psi_vals_sim, k_c = k_c_vals_sim)

for (model in c("gamma", "lognormal", "gamma2")) {

    cache_file <- file.path(
        "output",
        sprintf("sensitivity_od_gp_%s_%s_n%d_lam%g.csv",
                model, pathogen, n_index, lambda_gp)
    )

    sim_grid <- NULL
    if (file.exists(cache_file)) {
        cached     <- read_csv(cache_file, show_col_types = FALSE)
        have       <- cached %>% distinct(psi, k_c)
        missing_psi <- !(psi_vals_sim %in% have$psi)
        missing_kc  <- sapply(k_c_vals_sim,
                               function(v) !any(abs(have$k_c - v) < 1e-8 * max(1, v)))
        if (!any(missing_psi) && !any(missing_kc)) {
            cat(sprintf("  [%s GP] loaded from cache\n", model))
            sim_grid <- cached
        }
    }

    if (is.null(sim_grid)) {
        sim_grid <- need_gp %>% mutate(k_sim = NA_real_)

        for (idx in seq_len(nrow(sim_grid))) {
            psi <- sim_grid$psi[idx]
            k_c <- sim_grid$k_c[idx]
            fl  <- fl_objs_ln[[as.character(psi)]]

            gfun <- switch(model,
                gamma     = gen_inf_attempts_gammapoisson_contacts(
                                Tgen, R0, alpha, psi, k_c, lambda_gp),
                lognormal = make_gfun_lognormal_gp(k_c, psi, fl),
                gamma2    = make_gfun_gamma2_gp(k_c, psi)
            )

            noffspring  <- replicate(n_index, length(gfun(0)))
            sim_grid$k_sim[idx] <- compute_k(noffspring)

            if (idx %% 100 == 0)
                cat(sprintf("  [%s GP] %d / %d\n", model, idx, nrow(sim_grid)))
        }

        write_csv(sim_grid, cache_file)
        cat(sprintf("  [%s GP] saved to %s\n", model, cache_file))
    }

    sim_grid <- sim_grid %>%
        mutate(k_capped  = pmin(k_sim, k_cap),
               od_capped = pmin(1 / k_sim, od_cap))

    model_label <- switch(model,
        gamma     = "Gamma (Type I)",
        lognormal = "Log-normal",
        gamma2    = "Gamma (Type II)"
    )

    fig <- sim_grid %>%
        ggplot(aes(x = psi, y = k_c, fill = k_capped)) +
            geom_tile() +
            scale_y_log10() +
            scale_fill_viridis_c(option = "inferno",
                                 name   = sprintf("k\n(capped\nat %d)", k_cap),
                                 limits = c(0, k_cap)) +
            theme_classic(base_size = 14) +
            labs(x     = expression(psi),
                 y     = expression("Contact shape (" * sigma * ")"),
                 title = sprintf("%s (stochastic contacts)", model_label))

    save_fig(fig,
             sprintf("sensitivity_od_gp_%s_%s", model, pathogen),
             width = 5, height = 4)
    cat(sprintf("  [%s GP] figure saved\n", model))
}

cat("\nSensitivity overdispersion heatmaps complete.\n")
