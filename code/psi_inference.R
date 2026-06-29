library(igraph)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)
library(ggplot2)

# ==============================================================================
# Define incubation periods
# ==============================================================================
# 
# Sources:
#	---
#   COVID-19: Li et al. (2020)
#	  [https://www.nejm.org/doi/full/10.1056/NEJMoa2001316]
#	  Mean 5.2 days, 95% percentile 12.5 days; moment matched SD = 3.9
#	---
#   MERS: Assiri et al. (2013)
#	  [https://www.thelancet.com/journals/laninf/article/
#	    PIIS1473-3099(13)70204-4/fulltext]
#	  Median 5.2, 95% percentile 12.4; 
#	---
#   Measles: Lessler et al. (2009)
#	  [https://www.thelancet.com/journals/laninf/article/
#	     PIIS1473-3099(09)70069-6/fulltext]
#	  Median: 12.5; 95% percentile, 17.7
#	---
#   Ebola: WHO Ebola Response Team (2014)
#	  [https://www.nejm.org/doi/full/10.1056/NEJMoa1411100]
#	  Mean: 9.4, sd: 7.4
#	---
#   Pneumonic plague: Gani & Leach (2004)
#	  [https://wwwnc.cdc.gov/eid/article/10/4/03-0509_article]
#	  Mean: 4.3, SD: 1.8
#	---
#   Norovirus: Lee et al. (2013)
#	  [https://link.springer.com/article/10.1186/1471-2334-13-446]
#	  Median 1.2, 95th percentile 2.6; 
#	---
#   Nipah: Nikolay et al. (2019)
#	  [https://www.nejm.org/doi/full/10.1056/NEJMoa1805376]
#	  [https://www.nejm.org/doi/suppl/10.1056/NEJMoa1805376/
#	    suppl_file/nejmoa1805376_appendix.pdf]
#	  Mean: 9.7, SD: 2.2 
#	---
#   Smallpox: Nishiura & Eichner (2007)
#	  [https://www.jstor.org/stable/4621176?if_data=e30%3D&seq=2]
#	  Mean(log(t)): 2.47; SD(log(t)):0.17
#	---
#   Hepatitis A: CDC (2024)
#	  [https://www.cdc.gov/hepatitis-a/hcp/clinical-overview/index.html]
#	  Mean: 28, Range: 15-50; treat range as 95% interval.
#	---
#   Influenza A: Lessler et al. (2009)
#	  [https://www.thelancet.com/journals/laninf/article/
#	    PIIS1473-3099(09)70069-6/fulltext]
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
# ---
# COVID-19: 
#	[https://www.thelancet.com/journals/lanepe/article/PIIS2666-7762(22)00140-5/
#	  fulltext]
#	Shape mean: 2.39; Scale mean: 2.95
# ---
# MERS: 
#	[https://www.pnas.org/doi/10.1073/pnas.1519235113#st03]
#	Mean: 6.8; SD: 4.1
# ---
# Measles: 
#	[https://www.sciencedirect.com/science/article/abs/pii/
#	  S0022519311003146?via%3Dihub]
#	Mean: 12.2; SD: 3.62
# ---
# Ebola: 
#	[https://www.nejm.org/doi/10.1056/NEJMoa1411100#APPNEJMoa1411100SUP]
#	Mean: 15.3, SD: 9.1
# ---
# Pneumonic plague: 
#	[https://pmc.ncbi.nlm.nih.gov/articles/PMC2566243/]
#	Mean: 5.1, SD: 2.3
# ---
# Norovirus: 
#	[https://pmc.ncbi.nlm.nih.gov/articles/instance/2660689/bin/
#	  08-0299_Techapp1-s1.pdf]
#	alpha: 3.35; beta: 1/1.09
# ---
# Nipah virus: 
#	[https://www.nejm.org/doi/suppl/10.1056/NEJMoa1805376/suppl_file/
#	  nejmoa1805376_appendix.pdf]
#	Gamma mean: 12.7; Gamma SD: 3.0
# ---
# Smallpox: 
#	[https://www.cambridge.org/core/journals/epidemiology-and-infection/article/
#	  infectiousness-of-smallpox-relative-to-disease-age-estimates-based-on-
#	  transmission-network-and-incubation-period/
#	  F8449950FBEFCEED57994B50D18F96FB]
#	Mean: 16.0; SD: 4.0
# ---
# Hepatitis A: 
#	[https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0204201]
#	Mean: 23.9, SD: 20.9
# ---
# Influenza 
#	[https://pmc.ncbi.nlm.nih.gov/articles/PMC11370535/]
#	Mean: 3.2; SD: 2.1

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
		gi_mean   = gi$mean,  gi_sd  = gi$sd,
		inc_mean  = inc$mean, inc_sd = inc$sd,
		alpha_gi  = gi$mean^2  / gi$sd^2,
		beta_gi   = gi$mean    / gi$sd^2,
		a_obs     = inc$mean^2 / inc$sd^2,
		b_obs     = inc$mean   / inc$sd^2,
		gi_source = "literature"
	)
}, generation_params, incubation_params)

# ==============================================================================
# Import and clean OutbreakTrees data
# ==============================================================================

dat <- readRDS("data/data_tibble_trees.RDS") 
has_onset <- grepl("symptom_onset", dat$Attributes)
onset_dat <- dat[has_onset, ]

# Parse onset values to numeric days (relative to earliest onset in tree) 
parse_onset <- function(onset_raw){
	n <- length(onset_raw) 
	result <- rep(NA_real_, n)

	skip <- is.na(onset_raw) | 
		grepl("asymptomatic|unclear|unknown|or|before|after|largely|diagnosed|early|mid|late",
			onset_raw, ignore.case=TRUE)
	working <- onset_raw 
	working[skip] <- NA

	# Try 1: already numeric (relative days) 
	nums <- suppressWarnings(as.numeric(working))
	if(sum(!is.na(nums))>0) {
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
			print("WARNING: early onset dates are cycled to following year")
			early_idx <- which(!is.na(dates) & month_nums <= 3)
			dates[early_idx] <- dates[early_idx] + 365
		}
		ref <- min(dates, na.rm = TRUE)
		result[!is.na(dates)] <- as.numeric(dates[!is.na(dates)] - ref)
		return(result)
	}

	result

}

# Extract clusters of serial intervals from a transmission tree 
extract_clusters <- function(tree){
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

# Process all trees 
cat("\n--- Processing trees ---\n")
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
# from the observed SI pool.
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

		pathogen_params[[disease]]$gi_mean   <- gi_mean_hat
		pathogen_params[[disease]]$gi_sd     <- gi_sd_hat
		pathogen_params[[disease]]$alpha_gi  <- gi_mean_hat^2 / gi_var_hat
		pathogen_params[[disease]]$beta_gi   <- gi_mean_hat  / gi_var_hat
		pathogen_params[[disease]]$gi_source <- "data-derived"
	}
}

# ==============================================================================
# Likelihood machinery
#
# precompute_densities, compute_psi_posterior, and loglik_cluster live in
# code/psi_likelihood.R, shared with code/psi_identifiability.R so that
# empirical estimation and identifiability use identical likelihoods.
# Underlying density helpers dgamma_sum / dgamma_diff live in utils.R.
# ==============================================================================

source("code/psi_likelihood.R")

# ==============================================================================
# Main estimation loop
# ==============================================================================

cat("\n--- Estimating psi ---\n")

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
	cat(sprintf("    GI (%s): mean=%.1f, sd=%.1f => alpha=%.2f, beta=%.3f \n",
	    pp$gi_source, pp$gi_mean, pp$gi_sd, alpha_gi, beta_gi))
	cat(sprintf("    Incubation: mean=%.1f, sd=%.1f\n", pp$inc_mean, pp$inc_sd))

	cat("    Building densities...\n")
	densities <- precompute_densities(psi_grid, alpha_gi, beta_gi, a_obs, b_obs)
	cat("    Computing posterior...\n")

	post <- compute_psi_posterior(clusters, densities, psi_grid)

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

# ==============================================================================
# Figures
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
	# geom_text(data = annot_df,
	#           aes(x = 0.95, y = Inf, label = label),
	#           hjust = 1, vjust = 1.5, size = 3, color = "grey40") +
	facet_wrap(~ disease, scales = "free_y", ncol = 3) +
	coord_cartesian(xlim = c(0, 1)) +
	labs(
		x = expression(psi),
		y = "Posterior density"
		# title = expression("Posterior distribution of" ~ psi),
		# subtitle = "Dashed lines = 95% CI."
	) +
	theme_classic(base_size = 20) +
	theme(strip.text = element_text(face = "bold"),
	      strip.background = element_blank())

n_diseases <- length(disease_order)
fig_height <- ceiling(n_diseases / 3) * 3.5
save_fig(fig1, "psi_empirical_posteriors", width = 14, height = fig_height)


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
	coord_cartesian(xlim = c(0, 1)) +
	labs(
		x = expression(psi ~ "(posterior mean with 95% CI)"),
		y = NULL,
		title = expression("Estimated" ~ psi ~ "across diseases"),
		subtitle = "Bars = 95% credible intervals."
	) +
	theme_classic(base_size = 13) +
	theme(panel.grid.major.y = element_blank())
save_fig(fig2, "psi_empirical_summary", width = 10, height = max(4, n_diseases * 0.6))

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
	         sprintf("psi_empirical_si_%s", gsub("[^A-Za-z0-9]", "_", disease)),
	         width = 8, height = max(3, n_clusters * 0.25))
	cat(sprintf("    Saved SI figure for %s (%d infectors)\n", disease, n_clusters))
}


# --------------------------------------------------------------------------
# Print summary table
# --------------------------------------------------------------------------

cat("\n--- Psi estimation summary ---\n")
cat(sprintf("  %-20s  %5s  %5s  %8s  %8s  %14s  %8s  %s\n",
    "Disease", "Multi", "TotSI", "GI_alpha", "GI_beta", "psi 95% CI", "psi mode", "GI source"))
cat(paste(rep("-", 110), collapse = ""), "\n")

for (d in disease_order) {
	r <- annot_df[annot_df$disease == d, ]
	p <- posterior_results[[d]]
	gi_source <- pathogen_params[[d]]$gi_source
	cat(sprintf("  %-20s  %5d  %5d  %8.2f  %8.3f  [%4.2f, %4.2f]    %8.2f  %s\n",
	    d, r$n_multi, r$n_si,
	    p$alpha_gi[1], p$beta_gi[1],
	    r$ci_lo, r$ci_hi, r$post_mode, gi_source))
}

# Save results
results_out <- posterior_df %>%
	select(disease, psi_grid, posterior, n_clusters, n_multi, n_si,
	       alpha_gi, beta_gi, a_obs, b_obs, 
	       post_mode, post_mean, ci_lo, ci_hi)
write_csv(results_out, file.path("output", "psi_empirical_results.csv"))
cat(sprintf("\n  Saved results to output/psi_empirical_results.csv\n"))

# Save cluster data for downstream scripts (e.g. psi_sensitivity_icc.R)
saveRDS(
	list(all_disease_clusters = all_disease_clusters,
	     pathogen_params      = pathogen_params),
	file.path("output", "psi_cluster_data.RDS")
)
cat("  Saved cluster data to output/psi_cluster_data.RDS\n")

cat("\n=== Empirical psi estimation complete ===\n")

