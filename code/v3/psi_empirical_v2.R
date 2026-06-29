library(igraph)

# ==============================================================================
# Define incubation periods
# ==============================================================================
# 
# Sources:
#   COVID-19: Li et al. (2020) NEJM
#	  https://www.nejm.org/doi/full/10.1056/NEJMoa2001316
#	  Mean 5.2 days, 95% percentile 12.5 days; moment matched SD = 3.9
#   MERS: Assiri et al. (2013) Lancet ID
#	  https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(13)70204-4/fulltext
#	  Median 5.2, 95% percentile 12.4; 
#   Measles: Lessler et al. (2009) Lancet ID
#	  https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(09)70069-6/fulltext
#	  Median: 12.5; 95% percentile, 17.7
#   Ebola: WHO Ebola Response Team (2014) NEJM
#	  https://www.nejm.org/doi/full/10.1056/NEJMoa1411100
#	  Mean: 9.4, sd: 7.4
#   Pneumonic plague: Gani & Leach (2004) EID
#	  https://wwwnc.cdc.gov/eid/article/10/4/03-0509_article
#	  Mean: 4.3, SD: 1.8
#   Norovirus: Lee et al. (2013) BMC ID
#	  https://link.springer.com/article/10.1186/1471-2334-13-446
#	  Median 1.2, 95th percentile 2.6; 
#   Nipah: Nikolay et al. (2019) NEJM
#	  https://www.nejm.org/doi/full/10.1056/NEJMoa1805376
#	  https://www.nejm.org/doi/suppl/10.1056/NEJMoa1805376/suppl_file/nejmoa1805376_appendix.pdf
#	  Mean: 9.7, SD: 2.2 
#   Smallpox: Nishiura & Eichner (2007) Int J Hyg
#	  https://www.jstor.org/stable/4621176?if_data=e30%3D&seq=2
#	  Mean(log(t)): 2.47; SD(log(t)):0.17
#   Hepatitis A: CDC
#	  https://www.cdc.gov/hepatitis-a/hcp/clinical-overview/index.html
#	  Mean: 28, Range: 15-50; treat range as 95% interval.
#   Influenza A: Lessler et al. (2009), mean 1.4d, sd 0.5d
#	  https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(09)70069-6/fulltext
#	  Median: 1.4; 95th percentile: 2.8

incubation_params <- list(
	"COVID-19"         = list(mean = 5.2,  sd = 3.9),  # VERIFIED 
	"MERS"             = list(mean = 6.0,  sd = 3.4),  # VERIFIED
	"Measles"          = list(mean = 12.8, sd = 2.7),  # VERIFIED 
	"Ebola"            = list(mean = 9.4,  sd = 7.4),  # VERIFIED
	"Pneumonic plague" = list(mean = 4.3,  sd = 1.8),  # VERIFIED
	"Norovirus"        = list(mean = 1.3,  sd = 0.67), # VERIFIED
	"Nipah virus"      = list(mean = 9.7,  sd = 2.2),  # VERIFIED
	"Smallpox"         = list(mean = 12.0, sd = 2.1),  # VERIFIED
	"Hepatitis A"      = list(mean = 28.0, sd = 9.3),  # VERIFIED
	"Influenza"        = list(mean = 1.5,  sd = 0.67)  # VERIFIED
)

# ==============================================================================
# Define generation interval distributions
# ==============================================================================
# 
# Sources:
# COVID-19: 
#	https://www.thelancet.com/journals/lanepe/article/PIIS2666-7762(22)00140-5/fulltext
#	Shape mean: 2.39; Scale mean: 2.95
# MERS: 
#	https://www.pnas.org/doi/10.1073/pnas.1519235113#st03
#	Mean: 6.8; SD: 4.1
# Measles: 
#	https://www.sciencedirect.com/science/article/abs/pii/S0022519311003146?via%3Dihub
#	
# Ebola: 
#	https://www.nejm.org/doi/10.1056/NEJMoa1411100#APPNEJMoa1411100SUP
#	Mean: 15.3, SD: 9.1
# Pneumonic plague: 
#	https://pmc.ncbi.nlm.nih.gov/articles/PMC2566243/
#	Mean: 5.1, SD: 2.3
# Norovirus: 
#	https://pmc.ncbi.nlm.nih.gov/articles/instance/2660689/bin/08-0299_Techapp1-s1.pdf
#	alpha: 3.35; beta: 1/1.09
# Nipah virus: 
#	https://www.nejm.org/doi/suppl/10.1056/NEJMoa1805376/suppl_file/nejmoa1805376_appendix.pdf
#	Gamma mean: 12.7; Gamma sd: 3.0
# Smallpox: 
#	https://www.cambridge.org/core/journals/epidemiology-and-infection/article/infectiousness-of-smallpox-relative-to-disease-age-estimates-based-on-transmission-network-and-incubation-period/F8449950FBEFCEED57994B50D18F96FB
#	Mean: 16.0; SD: 4.0
# Hepatitis A: 
#	https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0204201
#	Mean: 23.9, SD: 20.9
# Influenza 
#	https://pmc.ncbi.nlm.nih.gov/articles/PMC11370535/
#	Mean: 3.2; Sd: 2.1

generation_params <- list(
	"COVID-19"         = list(mean = 7.1,  sd = 4.6),   # VERIFIED 
	"MERS"             = list(mean = 6.8,  sd = 4.1),   # VERIFIED
	"Measles"          = list(mean = 12.2, sd = 3.6),   # VERIFIED 
	"Ebola"            = list(mean = 15.3, sd = 9.1),   # VERIFIED
	"Pneumonic plague" = list(mean = 5.1,  sd = 2.3),   # VERIFIED
	"Norovirus"        = list(mean = 3.6,  sd = 2.0),   # VERIFIED
	"Nipah virus"      = list(mean = 12.7, sd = 3.0),   # VERIFIED
	"Smallpox"         = list(mean = 16.0, sd = 4.0),   # VERIFIED
	"Hepatitis A"      = list(mean = 23.9, sd = 20.9),  # VERIFIED
	"Influenza"        = list(mean = 3.2,  sd = 2.1)    # VERIFIED
	)

pathogen_params <- Map(function(gi, inc) {
	list(
		gi_mean  = gi$mean,  gi_sd  = gi$sd,
		inc_mean = inc$mean, inc_sd = inc$sd,
		alpha_gi = gi$mean^2  / gi$sd^2,
		beta_gi  = gi$mean    / gi$sd^2,
		a_obs    = inc$mean^2 / inc$sd^2,
		b_obs    = inc$mean   / inc$sd^2
	)
}, generation_params, incubation_params)

# ==============================================================================
# Import and clean OutbreakTrees data
# ==============================================================================

dat <- readRDS("data/data_tibble_trees.RDS")
has_onset <- grepl("symptom_onset", dat$Attributes)
onset_dat <- dat[has_onset, ]

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
	if (disease %in% c("H1N1", "H1N1, H3N2")) disease <- "Influenza"
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

# ==============================================================================
# Optional: replace published GI parameters with data-derived estimates
#
# The marginal SI distribution satisfies E[SI] = E[GI] and
# Var(SI) = Var(GI) + 2*Var(inc), both independently of psi. So given the
# (verified) incubation parameters, we can back out consistent GI estimates
# from the observed SI pool without touching the cluster structure.
#
# Set use_data_gi <- TRUE to use data-derived GI; FALSE uses published values.
# ==============================================================================

use_data_gi <- TRUE

if (use_data_gi) {
	cat("\n--- Estimating GI parameters from observed SIs ---\n")
	for (disease in names(all_disease_clusters)) {
		si_all <- unlist(all_disease_clusters[[disease]]$clusters)
		if (length(si_all) < 10) next

		inc_var <- pathogen_params[[disease]]$inc_sd^2

		gi_mean_hat <- mean(si_all)
		gi_var_hat  <- var(si_all) - 2 * inc_var

		if (gi_var_hat <= 0) {
			cat(sprintf(
				"  %-20s  obs SI var (%.2f) <= 2*inc_var (%.2f); keeping published GI\n",
				disease, var(si_all), 2 * inc_var
			))
			next
		}

		gi_sd_hat <- sqrt(gi_var_hat)

		cat(sprintf(
			"  %-20s  published: mean=%.1f sd=%.1f  ->  data: mean=%.1f sd=%.1f\n",
			disease,
			pathogen_params[[disease]]$gi_mean,
			pathogen_params[[disease]]$gi_sd,
			gi_mean_hat, gi_sd_hat
		))

		pathogen_params[[disease]]$gi_mean <- gi_mean_hat
		pathogen_params[[disease]]$gi_sd   <- gi_sd_hat
		pathogen_params[[disease]]$alpha_gi <- gi_mean_hat^2 / gi_var_hat
		pathogen_params[[disease]]$beta_gi  <- gi_mean_hat  / gi_var_hat
	}
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
	cat(sprintf("    Literature GI: mean=%.1f, sd=%.1f => alpha=%.2f, beta=%.3f \n",
	    pp$gi_mean, pp$gi_sd, alpha_gi, beta_gi))
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
		post_mode  = post_mode,
		post_mean  = post_mean,
		ci_lo      = ci_lo,
		ci_hi      = ci_hi
	)
}


# Sensitivity analysis (incubation SD scaling) has been moved to
# code/psi_sensitivity.R — source that script after this one to run it.

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

	# Stack dots that fall in the same 1-day bin within each infector row,
	# offsetting each successive dot upward by stack_step row-units.
	row_spacing <- 5.0   # vertical distance between infector baselines
	bin_width   <- 1.0
	stack_step  <- 0.05 * row_spacing

	cluster_df <- cluster_df %>%
		mutate(si_bin = round(si / bin_width) * bin_width) %>%
		group_by(row, si_bin) %>%
		mutate(stack_idx = row_number() - 1L) %>%
		ungroup() %>%
		mutate(y_stacked = row * row_spacing + stack_idx * stack_step)

	fig_si <- ggplot() +
		geom_hline(data = tibble(y = seq_len(n_clusters) * row_spacing),
		           aes(yintercept = y),
		           color = "grey80", linewidth = 0.3) +
		geom_point(data = cluster_df, aes(x = si, y = y_stacked),
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
# Print summary table
# --------------------------------------------------------------------------

cat("\n--- Psi estimation summary (literature GI) ---\n")
cat(sprintf("  %-20s  %5s  %5s  %8s  %8s  %15s  %10s  %s\n",
    "Disease", "Multi", "TotSI", "GI_alpha", "GI_beta", "psi 95% CI", "psi mode", "GI source"))
cat(paste(rep("-", 110), collapse = ""), "\n")

for (d in disease_order) {
	r <- annot_df[annot_df$disease == d, ]
	p <- posterior_results[[d]]
	cat(sprintf("  %-20s  %5d  %5d  %8.2f  %8.3f  [%4.2f, %4.2f]      %5.2f\n",
	    d, r$n_multi, r$n_si,
	    p$alpha_gi[1], p$beta_gi[1],
	    r$ci_lo, r$ci_hi, r$post_mode))
}

# Save results
results_out <- posterior_df %>%
	select(disease, psi_grid, posterior, n_clusters, n_multi, n_si,
	       alpha_gi, beta_gi, a_obs, b_obs, 
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


 # For most pathogens the observed SI standard deviation is smaller than
 #  the model predicts — in some cases dramatically so (Influenza 0.36×, COVID 0.50×, Ebola 0.59×). For
 #   a few it's larger (MERS 1.26×, Smallpox 1.66×). The model-predicted ICC is too high for a specific
 #   reason in each direction:

 #  When observed SD < model SD (Influenza, COVID, Ebola, Nipah...): The GI and incubation parameters —
 #   which come from large population-level studies — predict more SI spread than is actually present
 #  in these specific OutbreakTrees clusters. These clusters are often from small, well-investigated
 #  outbreaks (households, wards) where transmission was tightly timed. The model's assumed variance is
 #   too large, so it predicts high ICC, but the actual within-cluster covariance doesn't materialise
 #  at that magnitude.

 #  When observed SD > model SD (MERS, Smallpox): The actual SI spread is larger than the model expects
 #   — driven here by unusual outbreak dynamics (Korean hospital chains for MERS, large
 #  historically-documented smallpox outbreaks). Yet Smallpox is the one case where blue ≈ orange,
 #  suggesting that even with excess total variance, the within-cluster covariance structure is roughly
 #   what the model expects.

 #  The deeper issue is a mismatch between the population studied and the parameter source. The GI and
 #  incubation parameters were estimated from population-representative datasets, but the OutbreakTrees
 #   clusters are a highly selected subset — they were published precisely because they had clear
 #  transmission chains and reliable symptom onset data. That selection process distorts both the total
 #   variance and the within-cluster structure relative to what the population-level parameters
 #  predict.

 #  So to directly answer your question: the blue diamonds being systematically high is telling you
 #  that the model's assumed variance components ($\sigma^2_{GI}$ and $\sigma^2_{inc}$) don't closely
 #  match the specific clusters in the data, rather than that ψ is wrong per se. The ECDF and ICC are
 #  both showing you this mismatch. It argues for caution in interpreting the absolute ψ values, while
 #  the relative ordering across pathogens (low vs. high ψ) may still be informative.

#  ❯ Is there any value in calculating the generation interval distribution parameters directly from
#   the OutbreakTrees data, rather than using published values? I imagine that's a bit circular ---
#   but since this is already so hacked together, would that be a worthwhile thing to try? If so,
#   let's maybe add it as an option in psi_empirical_v2.R, rather than getting rid of the generation
#   interval parameter values that are already there (which I painstakingly verified)

# ⏺ There's real value in it, and the circularity is less severe than it first appears. Here's why:

#   The marginal SI distribution satisfies:
#   $$E[SI] = E[GI], \quad \text{Var}(SI) = \sigma^2_{GI} + 2\sigma^2_{inc}$$

#   Both of these hold regardless of ψ — ψ affects the within-cluster covariance structure but not the
#   marginal mean or variance. So if you know the incubation parameters (which you've verified), you
#   can back out:

#   $$\hat{\mu}{GI} = \overline{SI}, \quad \hat{\sigma}{GI} = \sqrt{\max(0,; s^2_{SI} -
#   2\sigma^2_{inc})}$$

#   from the observed SI pool without touching the cluster structure at all. That's not circular —
#   you're using the marginal distribution to pin down the GI, then using the cluster structure to
#   estimate ψ. The only dependency is that the incubation parameters are taken as given.

#   The main thing it buys you is consistency: the ICC comparison you just made showed that the
#   model-predicted total SI variance often doesn't match the data. Data-derived GI parameters would
#   fix that by construction, making the ICC comparison a fairer test of the ψ estimate rather than
#   also testing the GI parameters.

#   For implementation, I'd add it as a preprocessing option that slots into pathogen_params before the
#    likelihood estimation, keeping the verified published values untouched:



  # Tier 1 — Robust estimates (trustworthy regardless of incubation SD assumption)

  # - MERS (ψ ≈ 0.4–0.5, gradual slight decline): The most stable mid-range estimate. The signal in the
  #  MERS clusters is strong enough that the incubation SD assumption barely moves the posterior. This
  # is the estimate I'd trust most.
  # - Smallpox (ψ ≈ 0.4–0.5, similarly stable): Same story as MERS. Historical outbreak data with
  # well-defined transmission chains gives enough within-cluster signal to dominate the incubation
  # assumption.
  # - Norovirus (ψ ≈ 0.75–0.8, very stable, only slight drift at 1.5x): The most robust high-ψ
  # estimate. Norovirus genuinely has a diffuse infectiousness profile, and that signal is clearly
  # present in the data. Worth highlighting.
  # - Nipah virus (ψ ≈ 1.0 through 1.25x, drops slightly at 1.5x): Mostly robust. The drop at 1.5x is
  # notable but even then ψ ≈ 0.6 — still high. The high-ψ estimate is well-supported.

  # Tier 2 — Consistent but boundary-adjacent estimates

  # - Measles (ψ ≈ 0.1 at 0.5x, → 0 at 1.0x+): Consistently low. The mode reaches the boundary at
  # published parameters but is never large even at 0.5x inc_sd. Given the large Measles clusters in
  # OutbreakTrees, this is probably a real signal, not just boundary pileup. Measles has a tight,
  # bursty infectiousness profile.
  # - Influenza (ψ ≈ 0.25 at 0.5x, gradual decline to ~0 at 1.5x): Also consistently low, just with
  # more sensitivity. The estimate never crosses into the moderate range. Probably genuinely low ψ, but
  #  the published incubation SD is right at the edge of identifiability.
  # - COVID-19 (ψ ≈ 0.25 at 0.5x, → 0 at 1.0x): Similar pattern to Influenza, slightly more sensitive.
  # Even at 0.5x, the mode is only 0.25 — the data genuinely doesn't support high ψ. But the fact that
  # 1.0x gives ψ ≈ 0 with a pileup suggests this is the unidentifiability scenario we discussed: the
  # within-cluster SI variation is entirely explained by the incubation period at published parameters,
  #  leaving no room to estimate ψ. At smaller inc_sd there's slightly more room, and the data supports
  #  ψ ≈ 0.25. The honest interpretation: COVID-19 ψ is low but exact value is hard to pin down.

  # Tier 3 — Unreliable (dominated by parameterization, not data)

  # - Hepatitis A (ψ ≈ 0.35 at 0.5x, → 0 at 1.5x, wide CI throughout): Few clusters. The CIs span
  # almost the entire [0,1] range regardless of scale. The point estimate drifts with the assumption.
  # Can't be trusted.
  # - Ebola (ψ ≈ 0.8 at 0.5x, → 0 at 1.5x, huge CI at every scale): The estimate sweeps from high to
  # low across plausible incubation SDs. This is the clearest case of unidentifiability in your dataset
  #  — there's just not enough cluster structure to separate the ψ signal from the incubation period.
  # The answer depends almost entirely on what you assume for the incubation SD. I would not report an
  # Ebola ψ estimate.
  # - Pneumonic plague (ψ ≈ 1.0 at 0.5x, → 0 at 1.0x+): Hits both extremes of the [0,1] range depending
  #  on inc_sd. With very few and small clusters, the likelihood is nearly flat and the posterior is
  # driven entirely by whichever boundary the parameterization happens to favor. Essentially
  # unidentifiable — discard this one.

  # The general mechanism driving the downward trend:

  # Across most pathogens, increasing inc_sd → lower ψ mode. This is because a larger inc_sd attributes
  #  more of the observed within-cluster SI similarity to shared incubation period variation, leaving
  # less signal for ψ. Formally, the identifiable range of ICC narrows as σ²_inc grows (range = σ²_GI /
  #  (σ²_GI + 2σ²_inc)), so at large inc_sd the posterior is pushed toward the lower ψ boundary
  # regardless of what the data actually contain. The pathogens that resist this trend (Norovirus,
  # Nipah) are those where the within-cluster ICC signal is genuinely strong enough to overcome the
  # parameter-induced compression.

  # Bottom line for the paper: The three most defensible results are MERS (ψ ≈ 0.4–0.5), Smallpox (ψ ≈
  # 0.4–0.5), and Norovirus (ψ ≈ 0.75–0.8). Nipah is also defensible as high-ψ. Measles, Influenza, and
  #  COVID-19 are defensibly low-ψ but with the caveat that they sit near the identifiability boundary
  # at published parameters. Ebola, Pneumonic plague, and Hepatitis A should probably be reported with
  # heavy caveats or excluded from the main conclusions.