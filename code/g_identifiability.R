# ==============================================================================
# GI parameter identifiability under punctuated infectiousness (all pathogens)
# ==============================================================================
#
# Demonstrates that the standard iid serial-interval likelihood produces
# overconfident posteriors for GI shape (alpha) and rate (beta) when
# infectiousness is punctuated (psi < 1), because it ignores within-cluster
# correlation among sibling serial intervals.
#
# Approach:
#   1. Compute marginal serial interval density f_S(s; alpha, beta, a_obs, b_obs)
#      via numerical convolution (same for all psi)
#   2. Evaluate the iid posterior on a 2D grid over (alpha, beta)
#   3. Repeat across simulated datasets at different psi values
#   4. Show coverage of 95% CIs drops as psi -> 0
#
# Depends on: parslist (from parameters.R), save_fig, simulate_clusters (from utils.R)
# ==============================================================================

cat("=== GI parameter identifiability (all pathogens) ===\n")

# Defensive check: this script assumes run_analysis.R has already sourced
# utils.R and the tidyverse. Bail out early with a clear message if a key
# dependency is missing.
if (!"dplyr" %in% loadedNamespaces()) {
	stop("g_identifiability.R expects tidyverse and code/utils.R to be loaded ",
	     "first (typically via run_analysis.R).")
}

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Helper functions: SI density, grid posterior, posterior summaries, ICC.
#
# Gamma sum/difference convolution helpers (dgamma_sum, dgamma_diff) live in
# utils.R, sourced via run_analysis.R.
# ==============================================================================

#' Precompute the serial-interval density machinery for fast evaluation
#'
#' The serial interval is S = tau + (d_j - d_0), where:
#'   tau ~ Gamma(alpha, beta)         (generation interval)
#'   d_j, d_0 ~ Gamma(a_obs, b_obs)  (detection delays)
#'
#' The density is: f_S(s) = integral f_tau(tau) * f_{d-d0}(s - tau) dtau
#'
#' Strategy: pre-tabulate f_{d-d0} on a fine grid (fixed, known params),
#' then for each (alpha, beta), form g(tau_k) = dgamma(tau_k; alpha, beta)
#' and compute f_S(s_j) = dtau * sum_k D[j,k] * g[k], where
#' D[j,k] = f_{d-d0}(s_j - tau_k) is pre-computed once. This mirrors the
#' precompute-then-evaluate pattern used in psi_inference.R.
#'
#' @param a_obs Detection delay Gamma shape
#' @param b_obs Detection delay Gamma rate
#' @param alpha_range Range of alpha values to support
#' @param beta_range Range of beta values to support
#' @return List with D_matrix, tau_grid, s_grid, dtau, and a loglik() closure
precompute_si_density <- function(a_obs, b_obs, alpha_range, beta_range) {
	sd_d <- sqrt(a_obs) / b_obs

	# Determine grid ranges to cover all (alpha, beta) combinations
	# Max GI mean and SD occur at max alpha, min beta
	max_gi_mean <- max(alpha_range) / min(beta_range)
	max_gi_sd   <- sqrt(max(alpha_range)) / min(beta_range)

	# tau grid: covers Gamma(alpha, beta) for all (alpha, beta)
	# Start at small positive value to avoid Inf from dgamma(0, shape<1, rate)
	tau_grid <- seq(1e-4, max_gi_mean + 8 * max_gi_sd, length.out = 501)
	dtau <- diff(tau_grid)[1]

	# s grid: covers f_S support for all (alpha, beta)
	s_lo <- -8 * sd_d
	s_hi <- max_gi_mean + 8 * (max_gi_sd + sd_d)
	s_grid <- seq(s_lo, s_hi, length.out = 501)

	# Pre-tabulate f_{d_j - d_0} on a fine grid
	ddiff_extent <- max(abs(s_lo), s_hi) + max(tau_grid)
	ddiff_grid <- seq(-ddiff_extent, ddiff_extent, length.out = 2001)
	ddiff_vals <- dgamma_diff(ddiff_grid, a_obs, b_obs, a_obs, b_obs)
	ddiff_interp <- approxfun(ddiff_grid, ddiff_vals, rule = 2, yleft = 0, yright = 0)

	# Pre-compute D matrix: D[j, k] = f_{d-d0}(s_j - tau_k)
	D_matrix <- outer(s_grid, tau_grid, function(s, t) ddiff_interp(s - t))

	list(
		D_matrix = D_matrix,
		tau_grid = tau_grid,
		s_grid   = s_grid,
		dtau     = dtau,
		loglik   = function(alpha, beta, obs) {
			# Compute serial interval density for given (alpha, beta)
			g_vec <- dgamma(tau_grid, shape = alpha, rate = beta)
			# Guard against Inf/NaN from dgamma (e.g. shape < 1 near tau=0)
			g_vec[!is.finite(g_vec)] <- 0
			f_s_vals <- as.numeric(D_matrix %*% g_vec) * dtau
			f_s <- approxfun(s_grid, pmax(f_s_vals, 0), rule = 2, yleft = 0, yright = 0)

			f_vals <- f_s(obs)
			f_vals[f_vals < .Machine$double.xmin] <- .Machine$double.xmin
			sum(log(f_vals))
		}
	)
}

#' Compute the iid posterior on a 2D grid over (alpha, beta)
compute_ab_posterior <- function(serial_intervals, alpha_grid, beta_grid, si_density) {
	grid <- expand.grid(alpha = alpha_grid, beta = beta_grid)
	n_cores <- max(1L, parallel::detectCores(logical = FALSE))
	grid$loglik <- unlist(parallel::mclapply(seq_len(nrow(grid)), function(idx) {
		si_density$loglik(grid$alpha[idx], grid$beta[idx], serial_intervals)
	}, mc.cores = n_cores))

	# Normalize to posterior (flat prior, log-sum-exp)
	grid$loglik[is.nan(grid$loglik)] <- -Inf
	finite_ll <- grid$loglik[is.finite(grid$loglik)]
	if (length(finite_ll) == 0) {
		grid$posterior <- 0
		return(grid[, c("alpha", "beta", "loglik", "posterior")])
	}
	max_ll <- max(finite_ll)
	grid$posterior <- exp(grid$loglik - max_ll)
	grid$posterior[!is.finite(grid$posterior)] <- 0
	grid$posterior <- grid$posterior / sum(grid$posterior)

	grid[, c("alpha", "beta", "loglik", "posterior")]
}

#' Extract summary statistics from 2D posterior
posterior_summary_ab <- function(post_df) {
	# Marginal for alpha
	alpha_marg <- aggregate(posterior ~ alpha, data = post_df, FUN = sum)
	alpha_marg$posterior <- alpha_marg$posterior / sum(alpha_marg$posterior)
	alpha_cdf <- cumsum(alpha_marg$posterior)
	alpha_mean <- sum(alpha_marg$alpha * alpha_marg$posterior)
	alpha_ci_lo <- alpha_marg$alpha[which.min(abs(alpha_cdf - 0.025))]
	alpha_ci_hi <- alpha_marg$alpha[which.min(abs(alpha_cdf - 0.975))]

	# Marginal for beta
	beta_marg <- aggregate(posterior ~ beta, data = post_df, FUN = sum)
	beta_marg$posterior <- beta_marg$posterior / sum(beta_marg$posterior)
	beta_cdf <- cumsum(beta_marg$posterior)
	beta_mean <- sum(beta_marg$beta * beta_marg$posterior)
	beta_ci_lo <- beta_marg$beta[which.min(abs(beta_cdf - 0.025))]
	beta_ci_hi <- beta_marg$beta[which.min(abs(beta_cdf - 0.975))]

	list(
		alpha_mean  = alpha_mean,
		alpha_ci_lo = alpha_ci_lo,
		alpha_ci_hi = alpha_ci_hi,
		beta_mean   = beta_mean,
		beta_ci_lo  = beta_ci_lo,
		beta_ci_hi  = beta_ci_hi
	)
}

#' Intraclass correlation coefficient for sibling serial intervals
compute_si_icc <- function(psi, alpha, beta, a_obs, b_obs) {
	var_gi <- alpha / beta^2
	var_d  <- a_obs / b_obs^2
	((1 - psi) * var_gi + var_d) / (var_gi + 2 * var_d)
}

# ==============================================================================
# 1. Coverage simulation study across all pathogens
# ==============================================================================

cat("--- Coverage simulation study ---\n")

psi_vals     <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.0)
K            <- 100
n_replicates <- 500
a_obs        <- 4
b_obs        <- 1

cache_file <- file.path(
	"output",
	sprintf("g_identifiability_K%d_reps%d_aobs%g_bobs%g.csv",
	        K, n_replicates, a_obs, b_obs)
)

if (file.exists(cache_file)) {
	cat("  Loading cached results from", cache_file, "\n")
	results_df <- read_csv(cache_file, show_col_types = FALSE)
} else {

	results_list <- list()
	result_idx   <- 0

	for (pathogen_idx in seq_along(parslist)) {
		pars       <- parslist[[pathogen_idx]]
		alpha_true <- pars$alpha
		beta_true  <- pars$beta
		R0         <- pars$R0

		cat(sprintf("  Pathogen: %s (alpha=%.3f, beta=%.4f, R0=%g)\n",
		    pars$pathogen, alpha_true, beta_true, R0))

		# Grid for (alpha, beta) — centered on true values
		# Alpha lower bound >= 1.0: avoids Gamma density pole at tau=0 (shape < 1)
		# which causes spurious numerical integration artifacts.
		# Range (0.5, 1.8) × true value with 100 points gives adequate resolution
		# even for measles (alpha~11, narrow posterior with ~1200 serial intervals).
		alpha_grid <- seq(max(1.0, alpha_true * 0.5), alpha_true * 1.8, length.out = 100)
		beta_grid  <- seq(max(0.05, beta_true * 0.5), beta_true * 1.8, length.out = 100)

		# Precompute serial-interval density machinery for this pathogen
		cat("    Precomputing SI density...\n")
		si_density <- precompute_si_density(a_obs, b_obs, range(alpha_grid), range(beta_grid))
		cat("    SI density ready.\n")

		for (psi_true in psi_vals) {
			icc <- compute_si_icc(psi_true, alpha_true, beta_true, a_obs, b_obs)
			cat(sprintf("    psi = %.1f (ICC = %.3f)\n", psi_true, icc))

			for (rep in seq_len(n_replicates)) {
				set.seed(30000 * pathogen_idx + 10000 * which(psi_vals == psi_true) + rep)

				clusters <- simulate_clusters(
					K, psi_true, R0, alpha_true, beta_true, a_obs, b_obs, p_asc = 1
				)

				all_si <- unlist(clusters)
				n_obs <- length(all_si)
				n_clusters <- length(clusters)

				if (n_obs < 5) {
					result_idx <- result_idx + 1
					results_list[[result_idx]] <- data.frame(
						pathogen = pars$pathogen,
						psi_true = psi_true, replicate = rep,
						n_obs = n_obs, n_clusters = n_clusters,
						icc = icc,
						alpha_mean = NA, alpha_ci_lo = NA, alpha_ci_hi = NA,
						alpha_covers = NA,
						beta_mean = NA, beta_ci_lo = NA, beta_ci_hi = NA,
						beta_covers = NA
					)
					next
				}

				post <- compute_ab_posterior(all_si, alpha_grid, beta_grid, si_density)
				summ <- posterior_summary_ab(post)

				result_idx <- result_idx + 1
				results_list[[result_idx]] <- data.frame(
					pathogen    = pars$pathogen,
					psi_true    = psi_true,
					replicate   = rep,
					n_obs       = n_obs,
					n_clusters  = n_clusters,
					icc         = icc,
					alpha_mean  = summ$alpha_mean,
					alpha_ci_lo = summ$alpha_ci_lo,
					alpha_ci_hi = summ$alpha_ci_hi,
					alpha_covers = (alpha_true >= summ$alpha_ci_lo & alpha_true <= summ$alpha_ci_hi),
					beta_mean   = summ$beta_mean,
					beta_ci_lo  = summ$beta_ci_lo,
					beta_ci_hi  = summ$beta_ci_hi,
					beta_covers = (beta_true >= summ$beta_ci_lo & beta_true <= summ$beta_ci_hi)
				)

				if (rep %% 50 == 0) {
					cat(sprintf("      psi=%.1f, rep %d/%d done\n", psi_true, rep, n_replicates))
				}
			}
		}
	}

	results_df <- bind_rows(results_list)
	write_csv(results_df, cache_file)
	cat(sprintf("  Saved %d results to %s\n", nrow(results_df), cache_file))
}

# ==============================================================================
# 2. Figures
# ==============================================================================

cat("--- Generating figures ---\n")

pathogen_order <- sapply(parslist, function(p) p$pathogen)

# --------------------------------------------------------------------------
# Figure 1: Coverage vs psi, faceted by pathogen
# --------------------------------------------------------------------------

cat("  Figure 1: Coverage vs psi\n")

coverage_df <- results_df %>%
	filter(!is.na(alpha_covers)) %>%
	group_by(pathogen, psi_true) %>%
	summarise(
		alpha_coverage = mean(alpha_covers),
		beta_coverage  = mean(beta_covers),
		n_reps    = n(),
		mean_nobs = mean(n_obs),
		mean_icc  = mean(icc),
		.groups = "drop"
	)

coverage_long <- coverage_df %>%
	pivot_longer(
		cols = c(alpha_coverage, beta_coverage),
		names_to = "parameter",
		values_to = "coverage"
	) %>%
	mutate(parameter = case_when(
		parameter == "alpha_coverage" ~ "alpha",
		parameter == "beta_coverage"  ~ "beta"
	)) %>%
	mutate(se = sqrt(coverage * (1 - coverage) / n_reps))

coverage_long$pathogen <- factor(coverage_long$pathogen, levels = pathogen_order)

fig1 <- ggplot(coverage_long, aes(x = psi_true, y = coverage, color = parameter)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2.5) +
	geom_errorbar(aes(ymin = coverage - 1.96 * se,
	                  ymax = pmin(1, coverage + 1.96 * se)),
	              width = 0.03, linewidth = 0.5) +
	geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50") +
	facet_wrap(~ pathogen, nrow = 1,
	           labeller = as_labeller(c(
	               influenza = "Influenza",
	               omicron   = "SARS-CoV-2 Omicron",
	               measles   = "Measles"
	           ))) +
	scale_x_continuous(breaks = psi_vals) +
	scale_color_manual(
		values = c(alpha = "#F8766D", beta = "#00BFC4"),
		labels = c(alpha = expression(alpha ~ "(GI shape)"),
		           beta  = expression(beta  ~ "(GI rate)"))) +
	ylim(0.4, 1.0) +
	labs(
		x = expression(psi),
		y = "Empirical coverage of 95% CI",
		color = "Parameter"
		# title = sprintf("IID posterior coverage vs psi (K = %d, %d replicates)", K, n_replicates),
		# subtitle = "Dashed line = nominal 95%. Coverage drops at low psi due to ignored within-cluster correlation"
	) +
	theme_classic(base_size = 20) +
	theme(legend.position = "bottom",
	      strip.background = element_blank())
save_fig(fig1, "g_identifiability_coverage", width = 14, height = 6)

# --------------------------------------------------------------------------
# Figures 2 & 3: Marginal posterior distributions for alpha and beta
# --------------------------------------------------------------------------

cat("  Figures 2-3: Marginal posterior distributions\n")

poster_psi_vals <- c(0.1, 0.5, 1.0)
poster_reps     <- 1:25

posterior_curves <- list()
curve_idx <- 0

for (pathogen_idx in seq_along(parslist)) {
	pars       <- parslist[[pathogen_idx]]
	alpha_true <- pars$alpha
	beta_true  <- pars$beta
	R0         <- pars$R0

	cat(sprintf("    Posteriors for %s\n", pars$pathogen))

	alpha_grid <- seq(max(1.0, alpha_true * 0.5), alpha_true * 1.8, length.out = 100)
	beta_grid  <- seq(max(0.05, beta_true * 0.5), beta_true * 1.8, length.out = 100)

	si_density <- precompute_si_density(a_obs, b_obs, range(alpha_grid), range(beta_grid))

	for (psi_true in poster_psi_vals) {
		for (rep in poster_reps) {
			set.seed(30000 * pathogen_idx + 10000 * which(psi_vals == psi_true) + rep)

			clusters <- simulate_clusters(
				K, psi_true, R0, alpha_true, beta_true, a_obs, b_obs, p_asc = 1
			)
			all_si <- unlist(clusters)
			if (length(all_si) < 5) next

			post <- compute_ab_posterior(all_si, alpha_grid, beta_grid, si_density)

			# Marginal for alpha
			alpha_marg <- aggregate(posterior ~ alpha, data = post, FUN = sum)
			alpha_marg$posterior <- alpha_marg$posterior / sum(alpha_marg$posterior)

			# Marginal for beta
			beta_marg <- aggregate(posterior ~ beta, data = post, FUN = sum)
			beta_marg$posterior <- beta_marg$posterior / sum(beta_marg$posterior)

			curve_idx <- curve_idx + 1
			posterior_curves[[curve_idx]] <- rbind(
				data.frame(
					pathogen  = pars$pathogen,
					psi       = psi_true,
					replicate = rep,
					parameter = "alpha",
					value     = alpha_marg$alpha,
					density   = alpha_marg$posterior
				),
				data.frame(
					pathogen  = pars$pathogen,
					psi       = psi_true,
					replicate = rep,
					parameter = "beta",
					value     = beta_marg$beta,
					density   = beta_marg$posterior
				)
			)

			cat(sprintf("      psi=%.1f, rep=%d done\n", psi_true, rep))
		}
	}
}

posterior_df <- bind_rows(posterior_curves)
posterior_df$pathogen <- factor(posterior_df$pathogen, levels = pathogen_order)
posterior_df$psi_label <- factor(
	paste0("psi == ", posterior_df$psi),
	levels = paste0("psi == ", poster_psi_vals)
)

# True parameter values for vertical lines
true_vals <- do.call(rbind, lapply(parslist, function(p) {
	data.frame(
		pathogen = p$pathogen,
		alpha_true = p$alpha,
		beta_true  = p$beta
	)
}))
true_vals$pathogen <- factor(true_vals$pathogen, levels = pathogen_order)

# Figure 2: Alpha posteriors
alpha_df <- posterior_df %>% filter(parameter == "alpha")

fig2 <- ggplot(alpha_df, aes(x = value, y = density,
                              group = replicate)) +
	geom_line(alpha = 0.5, linewidth = 0.6, color = "#2c7bb6") +
	geom_vline(data = true_vals, aes(xintercept = alpha_true),
	           linetype = "dashed", color = "grey30") +
	facet_grid(psi_label ~ pathogen, scales = "free",
	           labeller = labeller(
	               psi_label = label_parsed,
	               pathogen  = as_labeller(c(
	                   influenza = "Influenza",
	                   omicron   = "SARS-CoV-2 Omicron",
	                   measles   = "Measles"
	               )))) +
	labs(
		x = expression(alpha ~ "(GI shape)"),
		y = "Marginal posterior density"
		# title = expression("Marginal posterior for" ~ alpha ~ "(GI shape) under IID likelihood"),
		# subtitle = "Dashed line = true value. Low psi yields overconfident, shifted posteriors."
	) +
	theme_classic(base_size = 18)
save_fig(fig2, "g_identifiability_alpha_posteriors", width = 14, height = 10)

# Figure 3: Beta posteriors
beta_df <- posterior_df %>% filter(parameter == "beta")

fig3 <- ggplot(beta_df, aes(x = value, y = density,
                              group = replicate)) +
	geom_line(alpha = 0.5, linewidth = 0.6, color = "#2c7bb6") +
	geom_vline(data = true_vals, aes(xintercept = beta_true),
	           linetype = "dashed", color = "grey30") +
	facet_grid(psi_label ~ pathogen, scales = "free",
	           labeller = labeller(
	               psi_label = label_parsed,
	               pathogen  = as_labeller(c(
	                   influenza = "Influenza",
	                   omicron   = "SARS-CoV-2 Omicron",
	                   measles   = "Measles"
	               )))) +
	labs(
		x = expression(beta ~ "(GI rate)"),
		y = "Marginal posterior density"
		# title = expression("Marginal posterior for" ~ beta ~ "(GI rate) under IID likelihood"),
		# subtitle = "Dashed line = true value. Low psi yields overconfident, shifted posteriors."
	) +
	theme_classic(base_size = 18)
save_fig(fig3, "g_identifiability_beta_posteriors", width = 14, height = 10)

# --------------------------------------------------------------------------
# Print summary table
# --------------------------------------------------------------------------

cat("\n--- Coverage summary ---\n")
for (path in pathogen_order) {
	cat(sprintf("\n  %s:\n", path))
	cat(sprintf("  %-6s  %-8s  %-10s  %-10s  %-10s\n",
	    "psi", "ICC", "alpha cov", "beta cov", "n_obs"))
	sub <- coverage_df %>% filter(pathogen == path)
	for (i in seq_len(nrow(sub))) {
		r <- sub[i, ]
		cat(sprintf("  %-6.1f  %-8.3f  %-10.3f  %-10.3f  %-10.1f\n",
		    r$psi_true, r$mean_icc,
		    r$alpha_coverage, r$beta_coverage,
		    r$mean_nobs))
	}
}

cat("\n=== GI parameter identifiability (all pathogens) complete ===\n")
