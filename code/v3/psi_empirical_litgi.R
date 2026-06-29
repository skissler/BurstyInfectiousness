# ==============================================================================
# Empirical estimation of psi using literature GI values
# ==============================================================================
#
# Same as psi_empirical.R but uses published generation interval parameters
# rather than estimating them from the data. This avoids circularity since
# we're estimating psi from the same serial intervals used to characterize
# the GI.
#
# Restricted to pathogens with well-characterized GI distributions.
#
# Data source: OutbreakTrees database (Taube et al., PLOS Biology 2022)
# Depends on: save_fig (from utils.R)
# ==============================================================================

cat("=== Empirical psi estimation (literature GI) from OutbreakTrees ===\n")

library(igraph)

# ==============================================================================
# 1. Published GI and incubation period parameters
# ==============================================================================
#
# Each entry contains:
#   gi_mean, gi_sd: generation interval mean and SD (days), Gamma-distributed
#   inc_mean, inc_sd: incubation period mean and SD (days), Gamma-distributed
#
# GI sources:
#   COVID-19: Ganyani et al. (2020) Eurosurveillance, mean 5.2d, sd 1.7d
#             (Singapore clusters; Gamma fit to infector-infectee pairs)
#   MERS:     Cauchemez et al. (2016) PNAS, SI mean 6.8d, sd 4.1d
#             (used as GI proxy; largest published MERS transmission dataset)
#   Measles:  parameters.R values: mean 12.2d, sd 3.6d
#             (Klinkenberg & Nishiura 2011, J Theor Biol)
#   Ebola:    Faye et al. (2015) Lancet ID, mean 14.2d, sd 7.1d
#             (Guinea 2014 outbreak, 152 transmission pairs)
#   Smallpox: Eichner & Dietz (2003), SI mean 16.0d, sd 4.0d
#             (5 outbreaks, 223 intervals; used as GI proxy)
#   Pneumonic plague: Gani & Leach (2004) EID, SI mean 5.7d, sd 3.7d
#             (Manchuria 1910-11 outbreak)
#   Norovirus: Gastanaduy et al. / Lee et al., SI mean 3.6d, sd 2.0d
#             (childcare outbreaks, Gamma fit)
#   Hepatitis A: Zhang et al. (2018) PLOS ONE, SI mean 23.9d, sd 20.9d
#             (Chinese elementary school outbreak, Gamma fit)
#
# Incubation period sources: same as psi_empirical.R

pathogen_params <- list(
	"COVID-19" = list(
		gi_mean = 5.2, gi_sd = 1.7,
		inc_mean = 5.2, inc_sd = 2.8,
		gi_source = "Ganyani et al. (2020) Eurosurveillance"
	),
	"MERS" = list(
		gi_mean = 6.8, gi_sd = 4.1,
		inc_mean = 5.2, inc_sd = 2.5,
		gi_source = "Cauchemez et al. (2016) PNAS"
	),
	"Measles" = list(
		gi_mean = 12.2, gi_sd = 3.6,
		inc_mean = 11.5, inc_sd = 2.2,
		gi_source = "Klinkenberg & Nishiura (2011) J Theor Biol"
	),
	"Ebola" = list(
		gi_mean = 14.2, gi_sd = 7.1,
		inc_mean = 11.4, inc_sd = 6.5,
		gi_source = "Faye et al. (2015) Lancet ID"
	),
	"Smallpox" = list(
		gi_mean = 16.0, gi_sd = 4.0,
		inc_mean = 12.5, inc_sd = 2.2,
		gi_source = "Eichner & Dietz (2003), 223 intervals"
	),
	"Pneumonic plague" = list(
		gi_mean = 5.7, gi_sd = 3.7,
		inc_mean = 3.5, inc_sd = 1.5,
		gi_source = "Gani & Leach (2004) EID"
	),
	"Norovirus" = list(
		gi_mean = 3.6, gi_sd = 2.0,
		inc_mean = 1.2, inc_sd = 0.5,
		gi_source = "Lee et al. (2013), childcare outbreaks"
	),
	"Hepatitis A" = list(
		gi_mean = 23.9, gi_sd = 20.9,
		inc_mean = 28.0, inc_sd = 7.0,
		gi_source = "Zhang et al. (2018) PLOS ONE"
	)
)

# Moment-match to Gamma parameters
for (d in names(pathogen_params)) {
	pp <- pathogen_params[[d]]
	pathogen_params[[d]]$alpha_gi <- pp$gi_mean^2 / pp$gi_sd^2
	pathogen_params[[d]]$beta_gi  <- pp$gi_mean / pp$gi_sd^2
	pathogen_params[[d]]$a_obs    <- pp$inc_mean^2 / pp$inc_sd^2
	pathogen_params[[d]]$b_obs    <- pp$inc_mean / pp$inc_sd^2
}

cat("\n--- Literature GI parameters ---\n")
cat(sprintf("  %-20s  %8s  %8s  %8s  %8s  %s\n",
    "Disease", "GI mean", "GI sd", "alpha", "beta", "Source"))
for (d in names(pathogen_params)) {
	pp <- pathogen_params[[d]]
	cat(sprintf("  %-20s  %8.1f  %8.1f  %8.2f  %8.3f  %s\n",
	    d, pp$gi_mean, pp$gi_sd, pp$alpha_gi, pp$beta_gi, pp$gi_source))
}

# ==============================================================================
# 2. Extract serial interval clusters from OutbreakTrees
# ==============================================================================

cat("\n--- Loading and parsing OutbreakTrees data ---\n")

dat <- readRDS("data/outbreaktrees_data.RDS")
has_onset <- grepl("symptom_onset", dat$Attributes)
onset_dat <- dat[has_onset, ]

cat(sprintf("  %d trees with symptom onset data (out of %d total)\n",
    sum(has_onset), nrow(dat)))

#' Parse onset values to numeric days (relative to earliest onset in tree)
parse_onset <- function(onset_raw) {
	n <- length(onset_raw)
	result <- rep(NA_real_, n)

	skip <- is.na(onset_raw) |
		grepl("asymptomatic|unclear|before|after|largely|diagnosed",
		      onset_raw, ignore.case = TRUE)
	working <- onset_raw
	working[skip] <- NA

	# Try 1: already numeric (relative days)
	nums <- suppressWarnings(as.numeric(working))
	if (sum(!is.na(nums)) > 0) {
		result[!is.na(nums)] <- nums[!is.na(nums)]
		return(result)
	}

	# Try 2: "Day X" format
	day_pattern <- "^Day ([0-9]+)$"
	day_matches <- grepl(day_pattern, working)
	if (any(day_matches, na.rm = TRUE)) {
		for (i in which(!is.na(working) & day_matches)) {
			result[i] <- as.numeric(sub(day_pattern, "\\1", working[i]))
		}
		if (sum(!is.na(result)) > 0) return(result)
	}

	# Try 3: M/D date format (handle Dec->Jan year boundary)
	dates <- as.Date(working, format = "%m/%d")
	if (sum(!is.na(dates)) > sum(!is.na(working)) * 0.3) {
		month_nums <- as.numeric(format(dates, "%m"))
		has_late <- any(month_nums >= 10, na.rm = TRUE)
		has_early <- any(month_nums <= 3, na.rm = TRUE)
		if (has_late && has_early) {
			early_idx <- which(!is.na(dates) & month_nums <= 3)
			dates[early_idx] <- dates[early_idx] + 365
		}
		ref <- min(dates, na.rm = TRUE)
		result[!is.na(dates)] <- as.numeric(dates[!is.na(dates)] - ref)
		return(result)
	}

	result
}

#' Extract clusters of serial intervals from a transmission tree
extract_clusters <- function(tree) {
	tree <- upgrade_graph(tree)
	onset_raw <- vertex_attr(tree, "onset")
	onset_num <- parse_onset(onset_raw)
	names(onset_num) <- V(tree)$name

	el <- as_edgelist(tree)
	infectors <- unique(el[, 1])

	clusters <- list()
	for (inf in infectors) {
		infectees <- el[el[, 1] == inf, 2]
		onset_inf <- onset_num[inf]
		onset_infectees <- onset_num[infectees]
		if (is.na(onset_inf)) next
		valid <- !is.na(onset_infectees)
		if (sum(valid) < 1) next
		si <- onset_infectees[valid] - onset_inf
		clusters[[length(clusters) + 1]] <- si
	}
	clusters
}

# Process all trees, keep only diseases we have literature GI for
all_disease_clusters <- list()
for (i in seq_len(nrow(onset_dat))) {
	disease <- as.character(onset_dat$Disease[i])
	if (!(disease %in% names(pathogen_params))) next

	tree_id <- as.character(onset_dat$id[i])

	clusters <- tryCatch(
		extract_clusters(onset_dat$tree[[i]]),
		error = function(e) list()
	)
	if (length(clusters) == 0) next

	si_all <- unlist(clusters)
	if (abs(median(si_all)) > 100) next

	if (!(disease %in% names(all_disease_clusters))) {
		all_disease_clusters[[disease]] <- list(clusters = list(), tree_ids = character(0))
	}
	all_disease_clusters[[disease]]$clusters <- c(
		all_disease_clusters[[disease]]$clusters, clusters
	)
	all_disease_clusters[[disease]]$tree_ids <- c(
		all_disease_clusters[[disease]]$tree_ids, tree_id
	)
}

# ----------------------------------------------------------------------
# 2b. Supplement with additional published datasets
# ----------------------------------------------------------------------

# --- Tianjin COVID-19 clusters (Ganyani et al. 2020 Eurosurveillance) ---
# cat("  Supplementing with Tianjin COVID-19 data (Ganyani et al. 2020)...\n")
# 
# tj <- tryCatch(read.csv("data/tianjin_covid_data.csv", stringsAsFactors = FALSE),
#                error = function(e) NULL)
# if (!is.null(tj)) {
# 	# Parse infection sources: extract TJ## case IDs
# 	tj_pairs <- list()
# 	for (i in seq_len(nrow(tj))) {
# 		if (is.na(tj$Infection_source[i])) next
# 		infector_ids <- regmatches(tj$Infection_source[i],
# 		                           gregexpr("TJ[0-9]+", tj$Infection_source[i]))[[1]]
# 		for (inf_id in infector_ids) {
# 			tj_pairs[[length(tj_pairs) + 1]] <- data.frame(
# 				infector = inf_id, infectee = tj$case_id[i],
# 				stringsAsFactors = FALSE
# 			)
# 		}
# 	}
# 	if (length(tj_pairs) > 0) {
# 		tj_pair_df <- do.call(rbind, tj_pairs)
# 		onset_map <- setNames(tj$symptom_onset, tj$case_id)
# 		tj_pair_df$onset_inf <- as.Date(onset_map[tj_pair_df$infector], format = "%d/%m/%Y")
# 		tj_pair_df$onset_ife <- as.Date(onset_map[tj_pair_df$infectee], format = "%d/%m/%Y")
# 		tj_pair_df$si <- as.numeric(tj_pair_df$onset_ife - tj_pair_df$onset_inf)
# 		tj_valid <- tj_pair_df[!is.na(tj_pair_df$si), ]
# 
# 		tj_clusters <- split(tj_valid$si, tj_valid$infector)
# 		tj_clusters <- unname(tj_clusters)
# 
# 		if (!("COVID-19" %in% names(all_disease_clusters))) {
# 			all_disease_clusters[["COVID-19"]] <- list(clusters = list(), tree_ids = character(0))
# 		}
# 		all_disease_clusters[["COVID-19"]]$clusters <- c(
# 			all_disease_clusters[["COVID-19"]]$clusters, tj_clusters
# 		)
# 		all_disease_clusters[["COVID-19"]]$tree_ids <- c(
# 			all_disease_clusters[["COVID-19"]]$tree_ids, "tianjin_ganyani_2020"
# 		)
# 		cat(sprintf("    Added %d Tianjin clusters (%d with >=2 offspring), %d SI\n",
# 		    length(tj_clusters), sum(sapply(tj_clusters, length) >= 2),
# 		    sum(sapply(tj_clusters, length))))
# 	}
# }

# --- Hart et al. (2022) eLife: Zhang superspreader cluster ---
# cat("  Supplementing with Hart et al. (2022) transmission pairs...\n")
# 
# hart <- tryCatch(readxl::read_xlsx("data/hart_transmission_pairs.xlsx"),
#                  error = function(e) NULL)
# if (!is.null(hart)) {
# 	# Zhang dataset: 35 pairs from one superspreader — one big cluster
# 	zhang <- hart[hart$Dataset == "Zhang", ]
# 	zhang_si <- zhang$t_s2 - zhang$t_s1
# 	zhang_cluster <- list(zhang_si)
# 
# 	# Other datasets: treat each unique t_s1 as a putative infector
# 	# This is imperfect (two infectors could share onset day) but
# 	# provides additional data, especially for large groups
# 	other_clusters <- list()
# 	for (ds in setdiff(unique(hart$Dataset), "Zhang")) {
# 		sub <- hart[hart$Dataset == ds, ]
# 		sub$si <- sub$t_s2 - sub$t_s1
# 		by_infector <- split(sub$si, sub$t_s1)
# 		other_clusters <- c(other_clusters, unname(by_infector))
# 	}
# 
# 	hart_clusters <- c(zhang_cluster, other_clusters)
# 
# 	if (!("COVID-19" %in% names(all_disease_clusters))) {
# 		all_disease_clusters[["COVID-19"]] <- list(clusters = list(), tree_ids = character(0))
# 	}
# 	all_disease_clusters[["COVID-19"]]$clusters <- c(
# 		all_disease_clusters[["COVID-19"]]$clusters, hart_clusters
# 	)
# 	all_disease_clusters[["COVID-19"]]$tree_ids <- c(
# 		all_disease_clusters[["COVID-19"]]$tree_ids, "hart_2022_elife"
# 	)
# 	n_hart_multi <- sum(sapply(hart_clusters, length) >= 2)
# 	cat(sprintf("    Added %d Hart clusters (%d with >=2 offspring), %d SI\n",
# 	    length(hart_clusters), n_hart_multi,
# 	    sum(sapply(hart_clusters, length))))
# }

cat("\n--- Extracted clusters by disease ---\n")
cat(sprintf("  %-20s  %5s  %5s  %5s  %8s\n",
    "Disease", "Trees", "Clust", "Multi", "TotalSI"))
for (d in names(all_disease_clusters)) {
	dc <- all_disease_clusters[[d]]
	n_trees <- length(unique(dc$tree_ids))
	n_clust <- length(dc$clusters)
	n_multi <- sum(sapply(dc$clusters, length) >= 2)
	n_si <- sum(sapply(dc$clusters, length))
	cat(sprintf("  %-20s  %5d  %5d  %5d  %8d\n",
	    d, n_trees, n_clust, n_multi, n_si))
}

# ==============================================================================
# 3. Density functions for the serial interval likelihood
# ==============================================================================

#' Density of X + Y where X ~ Gamma(shape1, rate1), Y ~ Gamma(shape2, rate2)
dgamma_sum_litgi <- function(x, shape1, rate1, shape2, rate2) {
	if (abs(rate1 - rate2) < 1e-10 * max(rate1, rate2)) {
		return(dgamma(x, shape = shape1 + shape2, rate = rate1))
	}
	sapply(x, function(xi) {
		if (xi <= 0) return(0)
		integrand <- function(t) {
			dgamma(t, shape1, rate1) * dgamma(xi - t, shape2, rate2)
		}
		tryCatch(
			integrate(integrand, lower = 0, upper = xi,
			          rel.tol = 1e-8, abs.tol = 1e-12)$value,
			error = function(e) 0
		)
	})
}

#' Density of X - Y where X ~ Gamma(shape1, rate1), Y ~ Gamma(shape2, rate2)
dgamma_diff_litgi <- function(d, shape1, rate1, shape2, rate2) {
	sapply(d, function(di) {
		integrand <- function(y) {
			dgamma(di + y, shape1, rate1) * dgamma(y, shape2, rate2)
		}
		lower <- max(0, -di)
		tryCatch(
			integrate(integrand, lower = lower, upper = Inf,
			          rel.tol = 1e-8, abs.tol = 1e-12)$value,
			error = function(e) 0
		)
	})
}

make_density_interp_litgi <- function(dfun, grid, ...) {
	vals <- dfun(grid, ...)
	approxfun(grid, vals, rule = 2, yleft = 0, yright = 0)
}

#' Log-likelihood for one cluster via delta quadrature
#' Integrate over delta, the shared infector-level offset
loglik_cluster_litgi <- function(s_vec, f_nu_interp, f_delta_interp, delta_grid) {
	dd <- diff(delta_grid)[1]
	m <- length(s_vec)

	f_delta_vals <- f_delta_interp(delta_grid)

	log_integrand <- rep(0, length(delta_grid))
	for (j in seq_len(m)) {
		f_nu_vals <- f_nu_interp(s_vec[j] - delta_grid)
		f_nu_vals[f_nu_vals < .Machine$double.xmin] <- .Machine$double.xmin
		log_integrand <- log_integrand + log(f_nu_vals)
	}
	log_integrand <- log_integrand + log(pmax(f_delta_vals, .Machine$double.xmin))

	max_li <- max(log_integrand)
	if (is.infinite(max_li) && max_li < 0) return(-Inf)
	log(sum(exp(log_integrand - max_li))) + max_li + log(dd)
}

# ==============================================================================
# 4. Precompute density engines and estimate psi for each disease
# ==============================================================================

precompute_engines_litgi <- function(psi_grid, alpha, beta, a_obs, b_obs) {
	mean_gi <- alpha / beta
	sd_gi   <- sqrt(alpha) / beta
	mean_d  <- a_obs / b_obs
	sd_d    <- sqrt(a_obs) / b_obs

	delta_lo <- -(mean_d + 5 * sd_d)
	delta_hi <- mean_gi + 5 * sd_gi
	delta_grid <- seq(delta_lo, delta_hi, length.out = 401)

	nu_hi <- mean_gi + mean_d + 5 * (sd_gi + sd_d)
	nu_grid <- seq(1e-6, nu_hi, length.out = 501)

	n_cores <- max(1L, parallel::detectCores(logical = FALSE))

	engines <- parallel::mclapply(psi_grid, function(psi) {
		if (psi < 1e-10) {
			f_nu <- approxfun(nu_grid, dgamma(nu_grid, shape = a_obs, rate = b_obs),
			                  rule = 2, yleft = 0, yright = 0)
			f_delta <- make_density_interp_litgi(
				dgamma_diff_litgi, delta_grid,
				shape1 = alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
		} else if (psi > 1 - 1e-10) {
			f_nu <- make_density_interp_litgi(
				dgamma_sum_litgi, nu_grid,
				shape1 = alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
			f_delta <- approxfun(delta_grid, dgamma(-delta_grid, shape = a_obs, rate = b_obs),
			                     rule = 2, yleft = 0, yright = 0)
		} else {
			f_nu <- make_density_interp_litgi(
				dgamma_sum_litgi, nu_grid,
				shape1 = psi * alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
			f_delta <- make_density_interp_litgi(
				dgamma_diff_litgi, delta_grid,
				shape1 = (1 - psi) * alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
		}
		list(f_nu = f_nu, f_delta = f_delta)
	}, mc.cores = n_cores)

	list(
		f_nu_interps    = lapply(engines, `[[`, "f_nu"),
		f_delta_interps = lapply(engines, `[[`, "f_delta"),
		delta_grid      = delta_grid
	)
}

compute_psi_posterior_litgi <- function(clusters, engines, psi_grid) {
	f_nu_interps    <- engines$f_nu_interps
	f_delta_interps <- engines$f_delta_interps
	delta_grid      <- engines$delta_grid

	logliks <- sapply(seq_along(psi_grid), function(idx) {
		total <- 0
		for (k in seq_along(clusters)) {
			ll <- loglik_cluster_litgi(
				clusters[[k]],
				f_nu_interps[[idx]],
				f_delta_interps[[idx]],
				delta_grid
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

# ==============================================================================
# 5. Main estimation loop
# ==============================================================================

cat("\n--- Estimating psi using literature GI values ---\n")

psi_grid <- seq(0, 1, length.out = 101)
min_multi_clusters <- 5

posterior_results <- list()

for (disease in names(all_disease_clusters)) {
	dc <- all_disease_clusters[[disease]]
	clusters <- dc$clusters
	multi_clusters <- clusters[sapply(clusters, length) >= 2]

	if (length(multi_clusters) < min_multi_clusters) {
		cat(sprintf("  Skipping %s (only %d multi-offspring clusters)\n",
		    disease, length(multi_clusters)))
		next
	}

	pp <- pathogen_params[[disease]]
	alpha_gi <- pp$alpha_gi
	beta_gi  <- pp$beta_gi
	a_obs    <- pp$a_obs
	b_obs    <- pp$b_obs

	si_all <- unlist(clusters)

	cat(sprintf("\n  %s:\n", disease))
	cat(sprintf("    Data: %d clusters (%d with >=2 offspring), %d serial intervals\n",
	    length(clusters), length(multi_clusters), length(si_all)))
	cat(sprintf("    Literature GI: mean=%.1f, sd=%.1f => alpha=%.2f, beta=%.3f (%s)\n",
	    pp$gi_mean, pp$gi_sd, alpha_gi, beta_gi, pp$gi_source))
	cat(sprintf("    Incubation: mean=%.1f, sd=%.1f\n", pp$inc_mean, pp$inc_sd))

	cat("    Building density engines...\n")
	engines <- precompute_engines_litgi(psi_grid, alpha_gi, beta_gi, a_obs, b_obs)
	cat("    Computing posterior...\n")

	post <- compute_psi_posterior_litgi(clusters, engines, psi_grid)

	# Summary statistics
	post_mean <- sum(psi_grid * post)
	post_cdf  <- cumsum(post)
	ci_lo <- psi_grid[which.min(abs(post_cdf - 0.025))]
	ci_hi <- psi_grid[which.min(abs(post_cdf - 0.975))]
	post_mode <- psi_grid[which.max(post)]

	cat(sprintf("    Posterior: mode=%.2f, mean=%.2f, 95%% CI=[%.2f, %.2f]\n",
	    post_mode, post_mean, ci_lo, ci_hi))

	posterior_results[[disease]] <- data.frame(
		disease    = disease,
		psi_grid   = psi_grid,
		posterior  = post,
		n_clusters = length(clusters),
		n_multi    = length(multi_clusters),
		n_si       = length(si_all),
		alpha_gi   = alpha_gi,
		beta_gi    = beta_gi,
		a_obs      = a_obs,
		b_obs      = b_obs,
		gi_source  = pp$gi_source,
		post_mode  = post_mode,
		post_mean  = post_mean,
		ci_lo      = ci_lo,
		ci_hi      = ci_hi
	)
}

# ==============================================================================
# 6. Figures
# ==============================================================================

cat("\n--- Generating figures ---\n")

posterior_df <- bind_rows(posterior_results)

# Order diseases by posterior mode
disease_order <- posterior_df %>%
	group_by(disease) %>%
	summarise(mode = psi_grid[which.max(posterior)], .groups = "drop") %>%
	arrange(mode) %>%
	pull(disease)
posterior_df$disease <- factor(posterior_df$disease, levels = disease_order)

# Annotation data
annot_df <- posterior_df %>%
	group_by(disease) %>%
	summarise(
		post_mode = psi_grid[which.max(posterior)],
		post_mean = first(post_mean),
		ci_lo     = first(ci_lo),
		ci_hi     = first(ci_hi),
		n_multi   = first(n_multi),
		n_si      = first(n_si),
		gi_source = first(gi_source),
		.groups   = "drop"
	)
annot_df$disease <- factor(annot_df$disease, levels = disease_order)
annot_df$label <- sprintf("n=%d SI, %d clusters", annot_df$n_si, annot_df$n_multi)

# --------------------------------------------------------------------------
# Figure 1: Posterior distributions for psi, one panel per disease
# --------------------------------------------------------------------------

cat("  Figure 1: Posterior distributions\n")

fig1 <- ggplot(posterior_df, aes(x = psi_grid, y = posterior)) +
	geom_line(linewidth = 0.9, color = "steelblue") +
	geom_ribbon(aes(ymin = 0, ymax = posterior), fill = "steelblue", alpha = 0.2) +
	geom_vline(data = annot_df, aes(xintercept = ci_lo),
	           linetype = "dashed", color = "grey50", linewidth = 0.4) +
	geom_vline(data = annot_df, aes(xintercept = ci_hi),
	           linetype = "dashed", color = "grey50", linewidth = 0.4) +
	geom_text(data = annot_df,
	          aes(x = 0.95, y = Inf, label = label),
	          hjust = 1, vjust = 1.5, size = 3, color = "grey40") +
	facet_wrap(~ disease, scales = "free_y", ncol = 3) +
	xlim(0, 1) +
	labs(
		x = expression(psi),
		y = "Posterior density",
		title = expression("Posterior distribution of" ~ psi ~ "(literature GI parameters)"),
		subtitle = "Dashed lines = 95% CI. GI parameters fixed from published estimates."
	) +
	theme_minimal(base_size = 12) +
	theme(strip.text = element_text(face = "bold"))

n_diseases <- length(disease_order)
fig_height <- ceiling(n_diseases / 3) * 3.5
save_fig(fig1, "psi_empirical_litgi_posteriors", width = 14, height = fig_height)

# --------------------------------------------------------------------------
# Figure 2: Summary dot plot with credible intervals
# --------------------------------------------------------------------------

cat("  Figure 2: Summary dot plot\n")

annot_df2 <- annot_df
annot_df2$disease <- factor(annot_df2$disease, levels = rev(disease_order))

violin_df <- posterior_df %>%
	mutate(disease = factor(disease, levels = rev(disease_order))) %>%
	group_by(disease) %>%
	mutate(
		y_pos     = as.integer(disease),
		post_norm = posterior / max(posterior) * 0.4
	) %>%
	ungroup()

fig2 <- ggplot(annot_df2, aes(y = disease, x = post_mean)) +
	geom_ribbon(data = violin_df,
	            aes(x = psi_grid, ymin = y_pos - post_norm, ymax = y_pos + post_norm,
	                group = disease),
	            fill = "steelblue", alpha = 0.2, inherit.aes = FALSE) +
	geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
	               height = 0.3, linewidth = 0.6, color = "steelblue") +
	geom_point(size = 3, color = "steelblue") +
	geom_text(aes(x = 0.99, label = label),
	          hjust = 1, size = 3, color = "grey40") +
	xlim(0, 1) +
	labs(
		x = expression(psi ~ "(posterior mean with 95% CI)"),
		y = NULL,
		title = expression("Estimated" ~ psi ~ "across diseases (literature GI)"),
		subtitle = "GI parameters fixed from published estimates. Bars = 95% credible intervals."
	) +
	theme_minimal(base_size = 13) +
	theme(panel.grid.major.y = element_blank())
save_fig(fig2, "psi_empirical_litgi_summary", width = 10, height = max(4, n_diseases * 0.6))

# --------------------------------------------------------------------------
# Figure 3: Serial interval dot plots per disease
# --------------------------------------------------------------------------

cat("  Figure 3: Serial interval dot plots\n")

for (disease in names(posterior_results)) {
	clusters      <- all_disease_clusters[[disease]]$clusters
	cluster_sizes <- sapply(clusters, length)
	order_idx     <- order(cluster_sizes, decreasing = TRUE)

	cluster_df <- map_dfr(seq_along(order_idx), function(row) {
		tibble(row = row, si = clusters[[order_idx[row]]])
	})

	n_clusters <- length(clusters)

	fig_si <- ggplot() +
		geom_hline(data = tibble(row = seq_len(n_clusters)),
		           aes(yintercept = row),
		           color = "grey80", linewidth = 0.3) +
		geom_point(data = cluster_df, aes(x = si, y = row),
		           size = 1.5, color = "steelblue", alpha = 0.8) +
		scale_y_reverse(breaks = NULL) +
		theme_classic() +
		theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
		labs(
			x = "Serial interval (days)",
			y = "Infector",
			title = sprintf("%s: serial intervals by infector", disease),
			subtitle = sprintf("%d infectors, %d total SI", n_clusters, sum(cluster_sizes))
		)

	save_fig(fig_si,
	         sprintf("psi_empirical_litgi_si_%s", gsub("[^A-Za-z0-9]", "_", disease)),
	         width = 8, height = max(3, n_clusters * 0.25))
	cat(sprintf("    Saved SI figure for %s (%d infectors)\n", disease, n_clusters))
}

# --------------------------------------------------------------------------
# Figure 4: Infectiousness kernel illustrations (Measles and COVID-19)
# --------------------------------------------------------------------------
#
# For each of 5 selected clusters per disease, plots:
#   - background: overall g(τ) = Gamma(alpha_gi, beta_gi), gray filled shape
#   - per infector: kernel = Gamma(psi*alpha, beta) hanging downward from a
#     horizontal line, positioned so its mean aligns with the cluster mean SI
#   - dots on the line at the observed serial intervals
#
# Kernel position: τ_shared ≈ mean(SI) - psi * mu_GI, so the kernel's mean
# (τ_shared + psi*mu_GI) equals mean(SI). This uses SI ≈ GI since incubation
# periods approximately cancel in expectation.
# --------------------------------------------------------------------------

cat("  Figure 4: Infectiousness kernel illustrations\n")

for (disease in c("Measles", "COVID-19")) {
	if (!(disease %in% names(posterior_results))) next

	dc      <- all_disease_clusters[[disease]]
	pp      <- pathogen_params[[disease]]
	pr_rows <- posterior_results[[disease]]

	psi_mean  <- pr_rows$post_mean[1]
	alpha_gi  <- pp$alpha_gi
	beta_gi   <- pp$beta_gi
	psi_alpha <- psi_mean * alpha_gi
	mu_gi     <- alpha_gi / beta_gi

	# Select 5 largest clusters
	clusters      <- dc$clusters
	cluster_sizes <- sapply(clusters, length)
	n_show        <- min(5, length(clusters))
	top_idx       <- order(cluster_sizes, decreasing = TRUE)[seq_len(n_show)]
	selected      <- clusters[top_idx]

	all_si <- unlist(selected)
	x_lo   <- min(c(0, all_si)) - 1
	x_hi   <- max(c(qgamma(0.995, alpha_gi, beta_gi), all_si)) + 1
	x_seq  <- seq(x_lo, x_hi, length.out = 1000)

	row_height <- 1.0
	kern_scale <- 0.85   # kernel peak height as fraction of row_height
	bg_peak    <- 0.75   # background peak height as fraction of row_height

	# Background: g(τ) = Gamma(alpha_gi, beta_gi)
	bg_x   <- x_seq[x_seq >= 0]
	bg_raw <- dgamma(bg_x, alpha_gi, beta_gi)
	bg_df  <- tibble(x = bg_x, ymin = 0,
	                 ymax = bg_raw / max(bg_raw) * bg_peak * row_height)

	# Per-infector kernels and dots
	# Row i sits at y_line = (i + 1) * row_height, leaving a gap above the background
	kernel_list <- list()
	dot_list    <- list()

	for (i in seq_len(n_show)) {
		si_i         <- selected[[i]]
		y_line       <- (i + 1) * row_height
		tau_shared_i <- mean(si_i) - psi_alpha / beta_gi

		kern_raw <- ifelse(x_seq > tau_shared_i,
		                   dgamma(x_seq - tau_shared_i, psi_alpha, beta_gi),
		                   0)
		max_kern <- max(kern_raw)
		kern_ht  <- if (max_kern > 0) kern_raw / max_kern * kern_scale * row_height else kern_raw

		kernel_list[[i]] <- tibble(x = x_seq,
		                           ymin = y_line - kern_ht, ymax = y_line, row = i)
		dot_list[[i]]    <- tibble(x = si_i, y = y_line, row = i)
	}

	kernel_df <- bind_rows(kernel_list)
	dot_df    <- bind_rows(dot_list)
	line_df   <- tibble(y = (seq_len(n_show) + 1) * row_height)

	fig_kern <- ggplot() +
		geom_ribbon(data = bg_df,
		            aes(x = x, ymin = ymin, ymax = ymax),
		            fill = "grey75", color = NA) +
		geom_ribbon(data = kernel_df,
		            aes(x = x, ymin = ymin, ymax = ymax, group = row),
		            fill = "grey20", color = NA) +
		geom_hline(data = line_df, aes(yintercept = y),
		           color = "black", linewidth = 0.4) +
		geom_point(data = dot_df, aes(x = x, y = y), size = 2) +
		scale_y_continuous(breaks = NULL) +
		theme_classic() +
		labs(x = "Serial interval (days)", y = "Infectiousness",
		     title = sprintf("%s  (ψ = %.2f)", disease, psi_mean))

	save_fig(fig_kern,
	         sprintf("psi_empirical_litgi_kernels_%s", gsub("[^A-Za-z0-9]", "_", disease)),
	         width = 7, height = n_show + 2)
	cat(sprintf("    Saved kernel figure for %s\n", disease))
}

# --------------------------------------------------------------------------
# Figure 5: Panel B — Within-cluster SI pairs at three reference ψ values
#
# Simulates (S1, S2) pairs from clusters of size 2 at ψ = 0, 0.5, 1 (gray),
# then overlays observed pairs from Measles and MERS.
# SIs are normalised by each disease's mean GI so both fall on a common scale.
# --------------------------------------------------------------------------

cat("  Figure 5: Within-cluster SI pairs\n")

set.seed(42)

n_sim_pairs  <- 600
psi_ref_vals <- c(0, 0.5, 1)
psi_labels   <- c("ψ = 0", "ψ = 0.5", "ψ = 1")

# MERS GI and incubation parameters for the simulated background
ref_pp      <- pathogen_params[["MERS"]]
alpha_sim   <- ref_pp$alpha_gi
beta_sim    <- ref_pp$beta_gi
a_obs_sim   <- ref_pp$a_obs
b_obs_sim   <- ref_pp$b_obs
mu_sim      <- ref_pp$gi_mean

# Simulate (S1, S2) pairs at each reference psi, normalised by mean GI
sim_pairs_df <- map_dfr(seq_along(psi_ref_vals), function(k) {
	psi   <- psi_ref_vals[k]
	psi_s <- pmax(pmin(psi, 1 - 1e-6), 1e-6)

	d_inf  <- rgamma(n_sim_pairs, a_obs_sim, b_obs_sim)
	g_sh   <- rgamma(n_sim_pairs, (1 - psi_s) * alpha_sim, beta_sim)
	g_ind1 <- rgamma(n_sim_pairs, psi_s * alpha_sim, beta_sim)
	g_ind2 <- rgamma(n_sim_pairs, psi_s * alpha_sim, beta_sim)
	d_ife1 <- rgamma(n_sim_pairs, a_obs_sim, b_obs_sim)
	d_ife2 <- rgamma(n_sim_pairs, a_obs_sim, b_obs_sim)

	tibble(
		s1        = (g_sh + g_ind1 + d_ife1 - d_inf) / mu_sim,
		s2        = (g_sh + g_ind2 + d_ife2 - d_inf) / mu_sim,
		psi_label = psi_labels[k]
	)
})

# Extract observed (S1, S2) pairs from all clusters with >= 2 secondaries.
# Both orderings (a,b) and (b,a) are included so the scatter is symmetric,
# matching the simulated background which is symmetric by construction.
show_diseases <- c("Measles", "MERS")
obs_pairs_df  <- map_dfr(show_diseases, function(disease) {
	if (!(disease %in% names(all_disease_clusters))) return(tibble())
	mu_d  <- pathogen_params[[disease]]$gi_mean
	multi <- Filter(function(cl) length(cl) >= 2,
	                all_disease_clusters[[disease]]$clusters)
	map_dfr(multi, function(cl) {
		pairs <- combn(cl, 2, simplify = FALSE)
		map_dfr(pairs, function(p) {
			bind_rows(
				tibble(s1 = p[1] / mu_d, s2 = p[2] / mu_d),
				tibble(s1 = p[2] / mu_d, s2 = p[1] / mu_d)
			)
		})
	}) %>% mutate(disease = disease)
})

# Replicate observed data across all three facets
obs_rep_df <- map_dfr(psi_labels, function(lbl) mutate(obs_pairs_df, psi_label = lbl))

sim_pairs_df$psi_label <- factor(sim_pairs_df$psi_label, levels = psi_labels)
obs_rep_df$psi_label   <- factor(obs_rep_df$psi_label,   levels = psi_labels)

all_vals <- c(sim_pairs_df$s1, sim_pairs_df$s2, obs_pairs_df$s1, obs_pairs_df$s2)
ax_lim   <- c(quantile(all_vals, 0.01, na.rm = TRUE),
              quantile(all_vals, 0.99, na.rm = TRUE))

fig_panel_b <- ggplot() +
	geom_point(data = sim_pairs_df, aes(x = s1, y = s2),
	           color = "grey75", size = 1, alpha = 0.6) +
	geom_abline(slope = 1, intercept = 0,
	            linetype = "dashed", linewidth = 0.5, color = "black") +
	geom_point(data = obs_rep_df, aes(x = s1, y = s2, color = disease),
	           size = 2.5, alpha = 0.85) +
	scale_color_manual(values = c("Measles" = "steelblue", "MERS" = "darkorange"),
	                   name = NULL) +
	coord_fixed(xlim = ax_lim, ylim = ax_lim) +
	facet_wrap(~ psi_label, nrow = 1) +
	theme_classic() +
	theme(legend.position = "bottom") +
	labs(x = expression(S[1] / bar(g)),
	     y = expression(S[2] / bar(g)),
	     title = "Within-cluster serial interval pairs (normalised by mean GI)")

save_fig(fig_panel_b, "psi_empirical_litgi_panel_b", width = 10, height = 4.5)
cat("    Saved Panel B\n")

# --------------------------------------------------------------------------
# Print summary table
# --------------------------------------------------------------------------

cat("\n--- Psi estimation summary (literature GI) ---\n")
cat(sprintf("  %-20s  %5s  %5s  %8s  %8s  %15s  %10s  %s\n",
    "Disease", "Multi", "TotSI", "GI_alpha", "GI_beta", "psi 95% CI", "psi mode", "GI source"))
cat(paste(rep("-", 110), collapse = ""), "\n")

for (d in disease_order) {
	r <- annot_df[annot_df$disease == d, ]
	p <- posterior_results[[d]]
	cat(sprintf("  %-20s  %5d  %5d  %8.2f  %8.3f  [%4.2f, %4.2f]      %5.2f  %s\n",
	    d, r$n_multi, r$n_si,
	    p$alpha_gi[1], p$beta_gi[1],
	    r$ci_lo, r$ci_hi, r$post_mode,
	    as.character(r$gi_source)))
}

# Save results
results_out <- posterior_df %>%
	select(disease, psi_grid, posterior, n_clusters, n_multi, n_si,
	       alpha_gi, beta_gi, a_obs, b_obs, gi_source,
	       post_mode, post_mean, ci_lo, ci_hi)
write_csv(results_out, file.path("output", "psi_empirical_litgi_results.csv"))
cat(sprintf("\n  Saved results to output/psi_empirical_litgi_results.csv\n"))

# Save cluster data for downstream scripts (e.g. psi_sensitivity_icc.R)
saveRDS(
	list(all_disease_clusters = all_disease_clusters,
	     pathogen_params      = pathogen_params),
	file.path("output", "psi_cluster_data.RDS")
)
cat("  Saved cluster data to output/psi_cluster_data.RDS\n")

cat("\n=== Empirical psi estimation (literature GI) complete ===\n")


# ⏺ Everything runs cleanly. Here's the summary with literature GI values:
#   Disease: Measles
#   psi mode: 0.00
#   95% CI: [0.00, 0.11]
#   Data: 72 SI, 8 clusters
#   GI source: Klinkenberg & Nishiura (2011)
#   ────────────────────────────────────────
#   Disease: Pneumonic plague
#   psi mode: 0.02
#   95% CI: [0.00, 0.20]
#   Data: 77 SI, 20 clusters
#   GI source: Gani & Leach (2004)
#   ────────────────────────────────────────
#   Disease: COVID-19
#   psi mode: 0.13
#   95% CI: [0.00, 0.75]
#   Data: 191 SI, 31 clusters
#   GI source: Ganyani et al. (2020)
#   ────────────────────────────────────────
#   Disease: Hepatitis A
#   psi mode: 0.15
#   95% CI: [0.04, 0.44]
#   Data: 72 SI, 18 clusters
#   GI source: Zhang et al. (2018)
#   ────────────────────────────────────────
#   Disease: Norovirus
#   psi mode: 0.66
#   95% CI: [0.35, 0.96]
#   Data: 59 SI, 14 clusters
#   GI source: Lee et al. (2013)
#   ────────────────────────────────────────
#   Disease: MERS
#   psi mode: 0.71
#   95% CI: [0.59, 0.82]
#   Data: 265 SI, 23 clusters
#   GI source: Cauchemez et al. (2016)
#   ────────────────────────────────────────
#   Disease: Ebola
#   psi mode: 0.77
#   95% CI: [0.05, 0.96]
#   Data: 43 SI, 11 clusters
#   GI source: Faye et al. (2015)
#   ────────────────────────────────────────
#   Disease: Smallpox
#   psi mode: 0.83
#   95% CI: [0.56, 0.96]
#   Data: 46 SI, 5 clusters
#   GI source: Eichner & Dietz (2003)
#   Compared to the version that estimated GI from the data, the results shift noticeably:

#   - MERS moved from 0.43 to 0.71 — using the Cauchemez GI (mean 6.8d, sd 4.1d) instead of the
#   moment-matched estimate (mean 11.5d) pushes psi higher. The literature GI is shorter than what the
#   raw SI suggests, so more of the SI variance gets attributed to individual-level variation (higher
#   psi).
#   - Measles and pneumonic plague remain firmly near 0 — these are robust to the GI specification.
#   Both suggest strongly punctuated infectiousness.
#   - COVID-19 is wide (CI [0.00, 0.75]) — despite having 191 SI and 31 clusters, most COVID trees are
#   small. The Ganyani GI (mean 5.2d, narrow sd 1.7d) is quite peaked, which limits how much the model
#   can discriminate.
#   - MERS is the tightest (CI width 0.23) — again the best-powered disease thanks to the Korean
#   superspreader outbreak.



  # COVID-19 before supplements: 69 clusters, 31 multi, 191 SI → mode=0.13, CI=[0.00, 0.75]
  # COVID-19 after supplements: 188 clusters, 84 multi, 455 SI → mode=0.97, CI=[0.46, 0.99]

  # The additional 43 Tianjin clusters and 76 Hart clusters (including the Zhang superspreader with 35
  # secondaries) more than doubled the data. The posterior shifted dramatically — the Zhang
  # superspreader cluster is particularly influential because 35 sibling serial intervals provide
  # strong evidence about within-cluster correlation.

  # The other diseases are unchanged since the supplemental data is all COVID-19.

  # Updated summary:
  # ┌──────────────────┬──────────┬──────────────┬─────────────────────┐
  # │     Disease      │ psi mode │    95% CI    │        Data         │
  # ├──────────────────┼──────────┼──────────────┼─────────────────────┤
  # │ Measles          │ 0.00     │ [0.00, 0.11] │ 72 SI, 8 clusters   │
  # ├──────────────────┼──────────┼──────────────┼─────────────────────┤
  # │ Pneumonic plague │ 0.02     │ [0.00, 0.20] │ 77 SI, 20 clusters  │
  # ├──────────────────┼──────────┼──────────────┼─────────────────────┤
  # │ Hepatitis A      │ 0.15     │ [0.04, 0.44] │ 72 SI, 18 clusters  │
  # ├──────────────────┼──────────┼──────────────┼─────────────────────┤
  # │ Norovirus        │ 0.66     │ [0.35, 0.96] │ 59 SI, 14 clusters  │
  # ├──────────────────┼──────────┼──────────────┼─────────────────────┤
  # │ MERS             │ 0.71     │ [0.59, 0.82] │ 265 SI, 23 clusters │
  # ├──────────────────┼──────────┼──────────────┼─────────────────────┤
  # │ Ebola            │ 0.77     │ [0.05, 0.96] │ 43 SI, 11 clusters  │
  # ├──────────────────┼──────────┼──────────────┼─────────────────────┤
  # │ Smallpox         │ 0.83     │ [0.56, 0.96] │ 46 SI, 5 clusters   │
  # ├──────────────────┼──────────┼──────────────┼─────────────────────┤
  # │ COVID-19         │ 0.97     │ [0.46, 0.99] │ 455 SI, 84 clusters │
  # └──────────────────┴──────────┴──────────────┴─────────────────────┘
  # The COVID-19 estimate is driven heavily by the Zhang cluster — worth noting that a single
  # superspreader event dominates the posterior. The sensitivity to the GI specification is also worth
  # flagging: the Ganyani GI (mean 5.2d, sd 1.7d) is quite narrow, and the estimate is conditional on
  # that being correct.

#   ❯ Ok -- now can we think through some of the limitations of this analysis? How much can we trust
#    these posteriors? For measles and pneumonic plague, where we estimate very low psi values:
#   could this come from biases in the data? Maybe we just captured superspreading events, but
#   didn't look at a person's entire infectiousness period? what other biases/limitations should
#   we be aware of, here?

# ⏺ This is a really important question. There are several biases that could affect these estimates,
#   and most of them push psi downward — which is concerning for exactly the diseases you mention.

#   Biases that inflate apparent punctuation (push psi toward 0):

#   1. Ascertainment of superspreading events (your point). Contact tracing preferentially captures
#   clusters where many secondaries appeared in a short window — these are conspicuous and easier to
#   trace. An infector who produces 1-2 secondaries spread over a week may never be identified as a
#   source. This systematically over-represents temporally clustered transmission, mimicking low psi.
#   2. Incomplete observation of the infectious period. Related but distinct: even when an infector is
#   identified, contact tracing may stop after the initial cluster is found. Later secondaries are
#   missed. For measles, where infectiousness spans ~8 days, this is a real concern — we may be
#   capturing the first burst and missing the tail.
#   3. Event-based transmission vs. biology. For measles (school exposure) and pneumonic plague
#   (household/close-contact), transmission often happens at discrete events. All secondaries infected
#   at the same school day or household gathering will naturally have correlated exposure times. This
#   is real clustering, but it reflects the social context of transmission rather than the biological
#   infectiousness profile. Our model doesn't distinguish these — it attributes all within-cluster SI
#   correlation to psi.
#   4. Date discretization. Symptom onset is recorded to the nearest day. For diseases with serial
#   intervals of ~10-15 days, a 1-day resolution might seem fine, but within a cluster, secondaries
#   infected at the same event will have onset dates that differ only by incubation period variation —
#   which could be sub-day. This artificially compresses within-cluster variance.
#   5. Publication/selection bias in OutbreakTrees. Published transmission trees tend to be from
#   notable outbreaks — often notable because of superspreading or unusual clustering. This is a
#   meta-level ascertainment bias that compounds point 1.
#   6. Right-censoring. In ongoing outbreaks, recent infectors haven't produced all their secondaries
#   yet. We observe the early ones, missing later ones — reducing within-cluster variance.

#   Biases that push psi upward (toward 1):

#   7. Transmission direction uncertainty. If who-infected-whom is misassigned, clusters get corrupted
#   — adding noise that looks like high within-cluster variance, inflating psi.
#   8. Heterogeneous incubation periods. If incubation periods vary across individuals within a cluster
#    (which they do), this adds within-cluster SI variance beyond what the GI alone would produce,
#   inflating apparent psi.

#   Other limitations:

#   9. Confounding with GI heterogeneity across infectors. If different infectors have systematically
#   different GI distributions (e.g., due to age, severity, viral load), between-cluster variance
#   increases for reasons unrelated to psi. This could bias psi in either direction depending on the
#   structure.
#   10. Incubation period parameter sensitivity. In the circular version, GI variance = SI variance - 2
#    × incubation variance. If published incubation parameters are wrong, GI parameters are wrong, and
#   psi estimates shift accordingly. In the literature-GI version, we depend entirely on published GI
#   estimates being correct.
#   11. Small cluster sizes. Many clusters have only 2-3 secondaries. Within-cluster variance is poorly
#    estimated from 2 observations. The posterior is being driven by a modest amount of information.

#   Net assessment: The dominant concern for measles and pneumonic plague is that biases 1-5 all push
#   psi downward, and these diseases are particularly susceptible because their transmission is often
#   event-mediated and well-traced. The low psi estimates could be partly or largely artifactual.

#   For COVID-19, where we estimate high psi, biases 7-8 could be inflating the estimate — but COVID
#   also has genuinely prolonged infectiousness with substantial pre-symptomatic transmission, which is
#    biologically consistent with high psi.

#   A key diagnostic would be: do the cluster sizes in our data look representative of the offspring
#   distribution, or are we over-sampling large clusters? If median cluster size in the data is much
#   larger than expected from the disease's reproduction number, that's evidence of ascertainment bias.