# ==============================================================================
# Identifiability calibration for empirical psi estimates
# ==============================================================================
#
# For each OutbreakTrees pathogen, given the actual data we have
# (cluster-size distribution, generation interval parameters, incubation
# parameters), which values of psi could we have identified?
#
# Analyses run by this script:
#   Part 1 -- per-pathogen calibration: sweep psi_true on a fine
#               grid; report the distribution of recovered psi_hat
#               and CI given the empirical sample size.
#   Part 2 -- ascertainment sweep: re-run calibration across
#               p_asc in {0.3, 0.5, 0.7, 1.0}. 
#   Part 3 -- incubation misspecification: fit with 0.5x and 2x
#               the literature incubation variance. 
#
# Depends on: psi_likelihood.R, utils.R (dgamma_sum / dgamma_diff /
#             simulate_clusters_from_sizes / save_fig), output from
#             psi_inference.R (output/psi_empirical_results.csv,
#             output/psi_cluster_data.RDS).
#
# Caches results to output/psi_identifiability_*.csv. Delete those files to
# force recomputation.
# ==============================================================================

source("code/psi_likelihood.R")

cat("\n=== Psi identifiability calibration ===\n")

# ------------------------------------------------------------------------------
# 1. Parameters
# ------------------------------------------------------------------------------

psi_grid       <- seq(0, 1, length.out = 101)   # for posterior estimation
psi_true_grid  <- seq(0, 1, length.out = 21)    # for the calibration sweep
n_reps_part1   <- 500
n_reps_part2   <- 200    # ascertainment sweep -- fewer reps to keep compute in check
n_reps_part3   <- 100    # incubation misspecification -- one psi_true value, fewer reps
p_asc_grid     <- c(0.3, 0.5, 0.7, 1.0)
incub_factors  <- c(0.5, 1.0, 2.0)     # multipliers applied to incubation variance

run_ascertainment    <- TRUE
run_misspecification <- TRUE

cache_part1 <- file.path("output", "psi_identifiability_calibration.csv")
cache_part2 <- file.path("output", "psi_identifiability_ascertainment.csv")
cache_part3 <- file.path("output", "psi_identifiability_incubation.csv")

# ------------------------------------------------------------------------------
# 2. Load empirical inputs
# ------------------------------------------------------------------------------

emp_path  <- "output/psi_empirical_results.csv"
clust_path <- "output/psi_cluster_data.RDS"

if (!file.exists(emp_path) || !file.exists(clust_path)) {
	stop("Missing required inputs from psi_inference.R: ",
	     emp_path, " or ", clust_path,
	     ". Run psi_inference.R first.")
}

emp_raw <- read_csv(emp_path, show_col_types = FALSE)
emp_summary <- emp_raw %>%
	distinct(disease, n_clusters, n_multi, n_si,
	         alpha_gi, beta_gi, a_obs, b_obs,
	         post_mode, post_mean, ci_lo, ci_hi)

clust_data <- readRDS(clust_path)
all_disease_clusters <- clust_data$all_disease_clusters

# Build a per-pathogen object with everything we need
get_pathogen_record <- function(disease) {
	row <- emp_summary[emp_summary$disease == disease, ]
	if (nrow(row) == 0) stop("No empirical record for ", disease)
	clusters <- all_disease_clusters[[disease]]$clusters
	multi_sizes <- sapply(clusters, length)
	multi_sizes <- multi_sizes[multi_sizes >= 2]
	list(
		disease     = disease,
		alpha_gi    = row$alpha_gi,
		beta_gi     = row$beta_gi,
		a_obs       = row$a_obs,
		b_obs       = row$b_obs,
		multi_sizes = multi_sizes,
		n_multi     = length(multi_sizes),
		emp_mode    = row$post_mode,
		emp_mean    = row$post_mean,
		emp_ci_lo   = row$ci_lo,
		emp_ci_hi   = row$ci_hi
	)
}

pathogen_records <- lapply(emp_summary$disease, get_pathogen_record)
names(pathogen_records) <- emp_summary$disease

cat(sprintf("  Loaded empirical records for %d pathogens.\n", length(pathogen_records)))

# ------------------------------------------------------------------------------
# 3. Helper: posterior summary from a single simulated dataset
# ------------------------------------------------------------------------------

#' Given clusters + precomputed densities, return mode/mean/95%CI of the
#' posterior over psi_grid.
posterior_summary <- function(clusters, densities, psi_grid) {
	if (length(clusters) < 2) {
		return(list(post_mode = NA_real_, post_mean = NA_real_,
		            ci_lo = NA_real_, ci_hi = NA_real_))
	}
	post <- compute_psi_posterior(clusters, densities, psi_grid)
	post_cdf <- cumsum(post)
	list(
		post_mode = psi_grid[which.max(post)],
		post_mean = sum(psi_grid * post),
		ci_lo     = psi_grid[which.min(abs(post_cdf - 0.025))],
		ci_hi     = psi_grid[which.min(abs(post_cdf - 0.975))]
	)
}

#' Build a deterministic seed for a (pathogen-index, phase, psi_idx, p_asc_idx,
#' incub_idx, rep) tuple. Decimal-place packing that fits within R's 32-bit
#' signed integer range (max ~2.15e9) for all sane grid sizes used here:
#' path_idx <= 10, phase <= 9, psi_idx <= 99, p_asc_idx <= 9, incub_idx <= 9,
#' rep <= 9999. Max packed value ~1.1e9 << 2^31 - 1.
make_seed <- function(path_idx, phase, psi_idx = 0, p_asc_idx = 0,
                      incub_idx = 0, rep = 0) {
	as.integer(path_idx)   * 100000000L +
	as.integer(phase)      * 10000000L  +
	as.integer(psi_idx)    * 100000L    +
	as.integer(p_asc_idx)  * 10000L     +
	as.integer(incub_idx)  * 1000L      +
	as.integer(rep)
}

# ------------------------------------------------------------------------------
# 4. Part 1: per-pathogen identifiability calibration
# ------------------------------------------------------------------------------

cat("\n--- Part 1: per-pathogen calibration ---\n")

if (file.exists(cache_part1)) {
	cat("  Loading cached results from", cache_part1, "\n")
	results_part1 <- read_csv(cache_part1, show_col_types = FALSE)
} else {
	n_cores <- max(1L, parallel::detectCores(logical = FALSE))
	cat(sprintf("  Running %d reps per (pathogen, psi_true) cell on %d cores.\n",
	    n_reps_part1, n_cores))

	results_part1_list <- list()

	for (path_idx in seq_along(pathogen_records)) {
		pp <- pathogen_records[[path_idx]]
		cat(sprintf("  [%2d/%d] %-18s  K_obs=%d  ",
		    path_idx, length(pathogen_records), pp$disease, pp$n_multi))

		# Precompute densities ONCE per pathogen
		densities <- precompute_densities(psi_grid, pp$alpha_gi, pp$beta_gi,
		                                  pp$a_obs, pp$b_obs)

		for (psi_idx in seq_along(psi_true_grid)) {
			psi_true <- psi_true_grid[psi_idx]
			rep_rows <- parallel::mclapply(seq_len(n_reps_part1), function(rep) {
				set.seed(make_seed(path_idx, phase = 3L, psi_idx = psi_idx, rep = rep))
				# Bootstrap empirical cluster-size distribution
				sampled_sizes <- sample(pp$multi_sizes, pp$n_multi, replace = TRUE)
				clusters <- simulate_clusters_from_sizes(
					sampled_sizes, psi_true,
					pp$alpha_gi, pp$beta_gi, pp$a_obs, pp$b_obs
				)
				s <- posterior_summary(clusters, densities, psi_grid)
				data.frame(
					disease   = pp$disease,
					psi_true  = psi_true,
					replicate = rep,
					post_mode = s$post_mode,
					post_mean = s$post_mean,
					ci_lo     = s$ci_lo,
					ci_hi     = s$ci_hi
				)
			}, mc.cores = n_cores)
			results_part1_list[[length(results_part1_list) + 1]] <- bind_rows(rep_rows)
		}
		cat("done.\n")
	}

	results_part1 <- bind_rows(results_part1_list)
	write_csv(results_part1, cache_part1)
	cat(sprintf("  Saved %d rows to %s\n", nrow(results_part1), cache_part1))
}

# Diagnostic table: mean recovered post_mean and 95% CI of post_mean across reps
part1_summary <- results_part1 %>%
	group_by(disease, psi_true) %>%
	summarise(
		mean_post_mean = mean(post_mean, na.rm = TRUE),
		q05_post_mean  = quantile(post_mean, 0.05, na.rm = TRUE),
		q95_post_mean  = quantile(post_mean, 0.95, na.rm = TRUE),
		mean_ci_lo     = mean(ci_lo, na.rm = TRUE),
		mean_ci_hi     = mean(ci_hi, na.rm = TRUE),
		mean_ci_width  = mean(ci_hi - ci_lo, na.rm = TRUE),
		.groups = "drop"
	)

cat("\n--- Part 1 diagnostic (mean CI width by pathogen, averaged over psi_true) ---\n")
part1_summary %>%
	group_by(disease) %>%
	summarise(mean_ci_width = mean(mean_ci_width), .groups = "drop") %>%
	arrange(mean_ci_width) %>%
	print(n = Inf)

# ------------------------------------------------------------------------------
# 5. Part 2: ascertainment sweep
# ------------------------------------------------------------------------------

if (run_ascertainment) {
	cat("\n--- Part 2: ascertainment sweep ---\n")

	if (file.exists(cache_part2)) {
		cat("  Loading cached results from", cache_part2, "\n")
		results_part2 <- read_csv(cache_part2, show_col_types = FALSE)
	} else {
		# For the ascertainment sweep we use simulate_clusters (Poisson offspring
		# + binomial thinning) so that p_asc has its natural interpretation.
		# Effective R0 is set so the expected post-thinning mean cluster size
		# matches the empirical mean of multi-cluster sizes -- this keeps the
		# data structure comparable as p_asc varies.

		n_cores <- max(1L, parallel::detectCores(logical = FALSE))
		cat(sprintf("  Running %d reps per (pathogen, psi_true, p_asc) cell on %d cores.\n",
		    n_reps_part2, n_cores))

		results_part2_list <- list()

		for (path_idx in seq_along(pathogen_records)) {
			pp <- pathogen_records[[path_idx]]
			# Mean observed multi-cluster size for this pathogen
			emp_mean_size <- mean(pp$multi_sizes)
			cat(sprintf("  [%2d/%d] %-18s  emp_mean_size=%.1f\n",
			    path_idx, length(pathogen_records), pp$disease, emp_mean_size))

			densities <- precompute_densities(psi_grid, pp$alpha_gi, pp$beta_gi,
			                                  pp$a_obs, pp$b_obs)

			for (p_asc_idx in seq_along(p_asc_grid)) {
				p_asc <- p_asc_grid[p_asc_idx]
				# Choose R0 so that R0 * p_asc = emp_mean_size, i.e., the
				# expected observed cluster size matches the empirical mean.
				R0_eff <- emp_mean_size / p_asc

				for (psi_idx in seq_along(psi_true_grid)) {
					psi_true <- psi_true_grid[psi_idx]
					rep_rows <- parallel::mclapply(seq_len(n_reps_part2), function(rep) {
						set.seed(make_seed(path_idx, phase = 4L, psi_idx = psi_idx,
						                   p_asc_idx = p_asc_idx, rep = rep))
						# Oversample K so that on average we get ~n_multi multi-clusters
						K_attempt <- ceiling(pp$n_multi * 3)
						raw <- simulate_clusters(
							K_attempt, psi_true, R0_eff,
							pp$alpha_gi, pp$beta_gi,
							pp$a_obs, pp$b_obs, p_asc = p_asc
						)
						multi <- raw[sapply(raw, length) >= 2]
						# Cap at n_multi to match empirical sample size
						if (length(multi) > pp$n_multi) {
							multi <- multi[seq_len(pp$n_multi)]
						}
						s <- posterior_summary(multi, densities, psi_grid)
						data.frame(
							disease   = pp$disease,
							psi_true  = psi_true,
							p_asc     = p_asc,
							replicate = rep,
							n_multi_obs = length(multi),
							post_mode = s$post_mode,
							post_mean = s$post_mean,
							ci_lo     = s$ci_lo,
							ci_hi     = s$ci_hi
						)
					}, mc.cores = n_cores)
					results_part2_list[[length(results_part2_list) + 1]] <- bind_rows(rep_rows)
				}
			}
		}

		results_part2 <- bind_rows(results_part2_list)
		write_csv(results_part2, cache_part2)
		cat(sprintf("  Saved %d rows to %s\n", nrow(results_part2), cache_part2))
	}
}

# ------------------------------------------------------------------------------
# 6. Part 3: incubation misspecification bias
# ------------------------------------------------------------------------------

if (run_misspecification) {
	cat("\n--- Part 3: incubation misspecification ---\n")

	if (file.exists(cache_part3)) {
		cat("  Loading cached results from", cache_part3, "\n")
		results_part3 <- read_csv(cache_part3, show_col_types = FALSE)
	} else {
		# For each pathogen, fix psi_true at the empirical posterior mean.
		# Generate data with the literature (a_obs, b_obs); fit with scaled
		# variance. Report shift in post_mean.

		n_cores <- max(1L, parallel::detectCores(logical = FALSE))
		cat(sprintf("  Running %d reps per (pathogen, incub_factor) cell on %d cores.\n",
		    n_reps_part3, n_cores))

		results_part3_list <- list()

		for (path_idx in seq_along(pathogen_records)) {
			pp <- pathogen_records[[path_idx]]
			psi_true <- pp$emp_mean
			mean_d   <- pp$a_obs / pp$b_obs
			var_d    <- pp$a_obs / pp$b_obs^2
			cat(sprintf("  [%2d/%d] %-18s  psi_true=%.2f (empirical mean)\n",
			    path_idx, length(pathogen_records), pp$disease, psi_true))

			for (incub_idx in seq_along(incub_factors)) {
				factor <- incub_factors[incub_idx]
				# Mis-specified Gamma: same mean, variance scaled by `factor`.
				var_d_mis <- var_d * factor
				b_obs_mis <- mean_d / var_d_mis
				a_obs_mis <- mean_d * b_obs_mis

				# Precompute densities at the *misspecified* incubation params
				densities_mis <- precompute_densities(
					psi_grid, pp$alpha_gi, pp$beta_gi, a_obs_mis, b_obs_mis
				)

				rep_rows <- parallel::mclapply(seq_len(n_reps_part3), function(rep) {
					set.seed(make_seed(path_idx, phase = 5L,
					                   incub_idx = incub_idx, rep = rep))
					# Generate data with TRUE incubation params, fit with MISSPECIFIED
					sampled_sizes <- sample(pp$multi_sizes, pp$n_multi, replace = TRUE)
					clusters <- simulate_clusters_from_sizes(
						sampled_sizes, psi_true,
						pp$alpha_gi, pp$beta_gi, pp$a_obs, pp$b_obs   # true
					)
					s <- posterior_summary(clusters, densities_mis, psi_grid)
					data.frame(
						disease       = pp$disease,
						psi_true      = psi_true,
						incub_factor  = factor,
						replicate     = rep,
						post_mode     = s$post_mode,
						post_mean     = s$post_mean,
						ci_lo         = s$ci_lo,
						ci_hi         = s$ci_hi
					)
				}, mc.cores = n_cores)
				results_part3_list[[length(results_part3_list) + 1]] <- bind_rows(rep_rows)
			}
		}

		results_part3 <- bind_rows(results_part3_list)
		write_csv(results_part3, cache_part3)
		cat(sprintf("  Saved %d rows to %s\n", nrow(results_part3), cache_part3))
	}

	# Quick diagnostic table
	cat("\n--- Part 3 diagnostic (mean post_mean by incub_factor) ---\n")
	results_part3 %>%
		group_by(disease, incub_factor) %>%
		summarise(
			mean_post_mean = mean(post_mean, na.rm = TRUE),
			.groups = "drop"
		) %>%
		pivot_wider(names_from = incub_factor, values_from = mean_post_mean,
		            names_prefix = "incub_x") %>%
		print(n = Inf)
}

cat("\n=== Psi identifiability calibration complete ===\n")
