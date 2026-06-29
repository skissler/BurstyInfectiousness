# ==============================================================================
# Shared likelihood machinery for psi inference and identifiability analysis
# ==============================================================================
#
# Functions exported:
#   make_density_interp(dfun, grid, ...)
#       Tabulate a density on a grid; return an approxfun interpolator.
#
#   loglik_cluster(si_vec, f_eps_tilde_interp, f_l_tilde_interp, u_grid)
#       Marginal log-likelihood for one cluster of serial intervals, integrating
#       over the shared infector-level offset l_tilde = l_i - d_i. See
#       supplement Eq. for L(psi).
#
#   precompute_densities(psi_grid, alpha, beta, a_obs, b_obs)
#       For each psi in psi_grid, build interpolators for the two marginal
#       densities (f_eps_tilde and f_l_tilde) that loglik_cluster needs. Returns
#       a list with parallel slots f_eps_tilde_interps, f_l_tilde_interps, and
#       the shared integration grid u_grid.
#
#   compute_psi_posterior(clusters, densities, psi_grid)
#       Discrete posterior over psi_grid given a list of serial-interval
#       clusters and the precomputed density bank.
#
#   compute_tail_prob(clusters, densities, psi_grid, threshold, side)
#   compute_interval_prob(clusters, densities, psi_grid, lo, hi)
#       Posterior probability summaries used by the identifiability analysis.
#
#   compute_icc_psi(psi, alpha, beta, a_obs, b_obs)
#       Closed-form intraclass correlation coefficient for sibling serial
#       intervals under the Gamma burst model (SI eq. ICC).
#
# Depends on: dgamma_sum, dgamma_diff (utils.R); parallel package.
# Used by:    psi_inference.R, psi_identifiability.R
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Low-level helpers
# ------------------------------------------------------------------------------

#' Pre-tabulate a density on a grid and return an approxfun interpolator
make_density_interp <- function(dfun, grid, ...) {
	vals <- dfun(grid, ...)
	approxfun(grid, vals, rule = 2, yleft = 0, yright = 0)
}

#' Log-likelihood for one cluster via Riemann sum over the latent offset grid
#'
#' Implements
#'   L_i(psi) = integral f_{l_tilde}(u; psi) * prod_j f_{eps_tilde}(s_j - u; psi) du
#' using log-sum-exp for numerical stability.
loglik_cluster <- function(si_vec, f_eps_tilde_interp, f_l_tilde_interp, u_grid) {
	du <- diff(u_grid)[1]
	chi_i <- length(si_vec)

	f_l_tilde_vals <- f_l_tilde_interp(u_grid)

	log_integrand <- rep(0, length(u_grid))
	for (j in seq_len(chi_i)) {
		f_eps_tilde_vals <- f_eps_tilde_interp(si_vec[j] - u_grid)
		f_eps_tilde_vals[f_eps_tilde_vals < .Machine$double.xmin] <- .Machine$double.xmin
		log_integrand <- log_integrand + log(f_eps_tilde_vals)
	}
	log_integrand <- log_integrand + log(pmax(f_l_tilde_vals, .Machine$double.xmin))

	max_li <- max(log_integrand)
	if (is.infinite(max_li) && max_li < 0) return(-Inf)
	log(sum(exp(log_integrand - max_li))) + max_li + log(du)
}


# ------------------------------------------------------------------------------
# 2. Precompute density bank for a full psi grid
# ------------------------------------------------------------------------------

#' Build f_eps_tilde and f_l_tilde interpolators for every psi in psi_grid
#'
#' For a given (alpha, beta, a_obs, b_obs), the marginal densities depend
#' only on psi. We compute them once and reuse across all replicates /
#' clusters. The two boundary cases (psi = 0, psi = 1) are handled
#' analytically rather than by numerical convolution.
#'
#' @return list(f_eps_tilde_interps, f_l_tilde_interps, u_grid)
precompute_densities <- function(psi_grid, alpha, beta, a_obs, b_obs) {
	mean_gi <- alpha / beta
	sd_gi   <- sqrt(alpha) / beta
	mean_d  <- a_obs / b_obs
	sd_d    <- sqrt(a_obs) / b_obs

	# u is the integration variable; its grid spans the support of l_tilde.
	u_lo <- -(mean_d + 5 * sd_d)
	u_hi <- mean_gi + 5 * sd_gi
	u_grid <- seq(u_lo, u_hi, length.out = 501)

	eps_tilde_hi <- mean_gi + mean_d + 5 * (sd_gi + sd_d)
	eps_tilde_grid <- seq(1e-6, eps_tilde_hi, length.out = 501)

	n_cores <- max(1L, parallel::detectCores(logical = FALSE))

	per_psi <- parallel::mclapply(psi_grid, function(psi) {
		if (psi < 1e-10) {
			# psi = 0: eps_j is a point mass at 0, so eps_tilde = d_j ~ Gamma(a_obs, b_obs).
			f_eps_tilde <- approxfun(
				eps_tilde_grid,
				dgamma(eps_tilde_grid, shape = a_obs, rate = b_obs),
				rule = 2, yleft = 0, yright = 0
			)
			f_l_tilde <- make_density_interp(
				dgamma_diff, u_grid,
				shape1 = alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
		} else if (psi > 1 - 1e-10) {
			# psi = 1: l_i is a point mass at 0, so l_tilde = -d_i and
			# f_{l_tilde}(u) = f_{d}(- u) for u < 0, zero otherwise.
			f_eps_tilde <- make_density_interp(
				dgamma_sum, eps_tilde_grid,
				shape1 = alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
			f_l_tilde <- approxfun(
				u_grid,
				dgamma(-u_grid, shape = a_obs, rate = b_obs),
				rule = 2, yleft = 0, yright = 0
			)
		} else {
			f_eps_tilde <- make_density_interp(
				dgamma_sum, eps_tilde_grid,
				shape1 = psi * alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
			f_l_tilde <- make_density_interp(
				dgamma_diff, u_grid,
				shape1 = (1 - psi) * alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
		}
		list(f_eps_tilde = f_eps_tilde, f_l_tilde = f_l_tilde)
	}, mc.cores = n_cores)

	list(
		f_eps_tilde_interps = lapply(per_psi, `[[`, "f_eps_tilde"),
		f_l_tilde_interps   = lapply(per_psi, `[[`, "f_l_tilde"),
		u_grid              = u_grid
	)
}


# ------------------------------------------------------------------------------
# 3. Posterior given precomputed densities
# ------------------------------------------------------------------------------

#' Discrete posterior weights over psi_grid (flat prior)
compute_psi_posterior <- function(clusters, densities, psi_grid) {
	f_eps_tilde_interps <- densities$f_eps_tilde_interps
	f_l_tilde_interps   <- densities$f_l_tilde_interps
	u_grid              <- densities$u_grid

	logliks <- sapply(seq_along(psi_grid), function(idx) {
		total <- 0
		for (k in seq_along(clusters)) {
			ll <- loglik_cluster(
				clusters[[k]],
				f_eps_tilde_interps[[idx]],
				f_l_tilde_interps[[idx]],
				u_grid
			)
			total <- total + ll
		}
		total
	})

	max_ll <- max(logliks)
	log_post <- logliks - max_ll
	post <- exp(log_post)
	post / sum(post)
}


# ------------------------------------------------------------------------------
# 4. Posterior-summary helpers (identifiability wrappers)
# ------------------------------------------------------------------------------

#' P(psi > threshold | data) or P(psi < threshold | data)
#'
#' @param side One of "upper" (default; returns P(psi > threshold)) or
#'   "lower" (returns P(psi < threshold)).
compute_tail_prob <- function(clusters, densities, psi_grid, threshold,
                              side = c("upper", "lower")) {
	side <- match.arg(side)
	post <- compute_psi_posterior(clusters, densities, psi_grid)
	if (side == "upper") sum(post[psi_grid > threshold])
	else                 sum(post[psi_grid < threshold])
}

#' P(lo <= psi <= hi | data)
compute_interval_prob <- function(clusters, densities, psi_grid, lo, hi) {
	post <- compute_psi_posterior(clusters, densities, psi_grid)
	sum(post[psi_grid >= lo & psi_grid <= hi])
}


# ------------------------------------------------------------------------------
# 5. Closed-form ICC under the Gamma burst model
# ------------------------------------------------------------------------------

#' Intraclass correlation coefficient for sibling serial intervals
#'
#' Closed form from SI:
#'   ICC = [(1 - psi) * alpha / beta^2 + sigma_d^2]
#'         / [alpha / beta^2 + 2 * sigma_d^2]
#'   sigma_d^2 = a_obs / b_obs^2 = Var[d].
compute_icc_psi <- function(psi, alpha, beta, a_obs, b_obs) {
	var_gi <- alpha / beta^2
	var_d  <- a_obs / b_obs^2
	((1 - psi) * var_gi + var_d) / (var_gi + 2 * var_d)
}
