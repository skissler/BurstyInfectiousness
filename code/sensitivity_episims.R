library(tidyverse)
library(parallel)
source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

# ==============================================================================
# SENSITIVITY ANALYSIS: epidemic simulations under alternative burst models
#
# Replicates the core results of episims.R under two additional burst models:
#   - lognormal: tau_ij = l_i + eps_j, g(tau) log-normal, f_l via
#                Fenton-Wilkinson moment matching
#   - gamma2:    tau_ij = l_i + eps_j, eps_j = sqrt(psi)*tau_tilde (same shape
#                as g(tau) ~ Gamma), f_l via exact compound-Poisson sampler
# These are compared against the original gamma (shared-rate) model.
#
# Cache files use the prefix "sensitivity_episims_" to avoid overwriting the
# main analysis cache.
# ==============================================================================

psivals <- c(0, 0.5, 1)

# Maximum time displayed in the trajectory figures (days). Trajectories are
# filtered to this window before plotting so the 3x3 grid stays readable.
xlim_max <- 100

establishment_threshold <- as.integer(round(establishment_threshold))

n_cores <- getOption("mc.cores",
					 as.integer(Sys.getenv("MC_CORES",
										   max(1L, detectCores() - 1L))))

models <- c("gamma", "lognormal", "gamma2")

model_labels <- c("gamma"     = "Gamma (Type I)",
				  "lognormal" = "Log-normal",
				  "gamma2"    = "Gamma (Type II)")

model_cols <- c("Gamma (Type I)"  = "#1f77b4",
				"Log-normal"      = "#ff7f0e",
				"Gamma (Type II)" = "#2ca02c")

sensitivity_cache_path <- function(pathogen, nsim, popsize) {
	file.path("output",
			  sprintf("sensitivity_episims_%s_n%d_s%d.csv", pathogen, popsize, nsim))
}

sensitivity_cache_path_plot <- function(pathogen, nsim, popsize) {
	file.path("output",
			  sprintf("sensitivity_episims_plot_%s_n%d_s%d.csv", pathogen, popsize, nsim))
}

sensitivity_load_cache <- function(pathogen, nsim, popsize, psivals, models) {
	f_sum  <- sensitivity_cache_path(pathogen, nsim, popsize)
	f_plot <- sensitivity_cache_path_plot(pathogen, nsim, popsize)
	if (!file.exists(f_sum) || !file.exists(f_plot)) return(NULL)
	df <- read_csv(f_sum, show_col_types = FALSE)
	if (!all(psivals %in% unique(df$psi)))   return(NULL)
	if (!all(models  %in% unique(df$model))) return(NULL)
	cat(sprintf("  %s: loading sensitivity cache from %s\n", pathogen, f_sum))
	list(
		summary = df,
		plot    = read_csv(f_plot, show_col_types = FALSE)
	)
}

# ==============================================================================
# Loop over pathogens
# ==============================================================================

all_summaries <- list()

for (pars in parslist[2]) {

	pathogen  <- pars$pathogen
	Tgen      <- pars$Tgen
	Tvar      <- pars$Tvar
	alpha <- pars$alpha
	beta  <- pars$beta
	R0        <- pars$R0

	cat(sprintf("\n===== SENSITIVITY %s: T=%.2f, alpha=%.2f, R0=%.1f =====\n",
				pathogen, Tgen, alpha, R0))

	# ------------------------------------------------------------------
	# Deterministic renewal-equation baseline (model-specific GI kernel)
	# Log-normal parameters moment-matched to Gamma(alpha, beta)
	# ------------------------------------------------------------------
	sigma2_tau_ln <- log(1 + 1 / alpha)
	mu_tau_ln     <- log(Tgen) - sigma2_tau_ln / 2
	sigma_tau_ln  <- sqrt(sigma2_tau_ln)

	ren_baselines <- list(
		"Gamma (Type I)"  = renewal_epidemic(R0, alpha, Tgen, popsize),
		"Log-normal"      = renewal_epidemic(R0, alpha, Tgen, popsize,
								  gi_pdf = function(t) dlnorm(t, mu_tau_ln, sigma_tau_ln)),
		"Gamma (Type II)" = renewal_epidemic(R0, alpha, Tgen, popsize)
	)

	ren_out <- bind_rows(
		mapply(function(df, mod) mutate(df, model = mod),
			   ren_baselines, names(ren_baselines), SIMPLIFY = FALSE)
	) %>%
		mutate(model = factor(model, levels = model_labels[models]))

	ren_daily <- ren_out %>%
		mutate(day = floor(t)) %>%
		group_by(model, day) %>%
		summarise(cuminf = max(cuminf), .groups = "drop") %>%
		arrange(model, day) %>%
		group_by(model) %>%
		mutate(newinf = cuminf - lag(cuminf),
			   newinf = if_else(is.na(newinf), cuminf, newinf)) %>%
		ungroup()

	# ------------------------------------------------------------------
	# Simulations (load from cache or run fresh)
	# ------------------------------------------------------------------
	cached <- sensitivity_load_cache(pathogen, nsim_small, popsize, psivals, models)

	if (!is.null(cached)) {
		sim_df  <- cached$summary
		plot_df <- cached$plot %>%
			mutate(psi   = factor(psi,   levels = c(0, 0.5, 1)),
				   model = factor(model, levels = models, labels = model_labels[models]))
	} else {
		fl_objs <- list()
		for (psi in psivals) {
			if (psi > 0 && psi < 1) {
				fl_objs[[paste0("lognormal_", psi)]] <-
					precompute_fl_lognormal(psi,
						mu_tau    = log(Tgen) - log(1 + 1/alpha) / 2,
						sigma_tau = sqrt(log(1 + 1/alpha)))
			}
		}

		task_grid <- expand_grid(sim   = seq_len(nsim_small),
								 psi   = psivals,
								 model = models)

		RNGkind("L'Ecuyer-CMRG")
		set.seed(sum(utf8ToInt(pathogen)) + 42L)

		run_one <- function(i) {
			sim   <- task_grid$sim[i]
			psi   <- task_grid$psi[i]
			model <- task_grid$model[i]

			fl_key_ln <- paste0("lognormal_", psi)

			gen <- switch(model,
				gamma     = gen_inf_attempts_gamma(Tgen, R0, alpha, psi),
				lognormal = gen_inf_attempts_lognormal(
								Tgen, Tvar, R0, psi,
								fl_obj = if (psi > 0 && psi < 1) fl_objs[[fl_key_ln]] else NULL),
				gamma2    = gen_inf_attempts_gamma2(Tgen, R0, alpha, beta, psi)
			)

			tinf            <- sim_stochastic_fast(n = popsize, gen_inf_attempts = gen)
			infection_times <- sort(tinf[is.finite(tinf)])
			n_inf           <- length(infection_times)
			established     <- as.integer(n_inf >= establishment_threshold)
			establishment_time <- if (established == 1L)
									 infection_times[establishment_threshold]
								  else NA_real_

			summary_row <- tibble(sim = sim, psi = psi, model = model,
								  established = established,
								  establishment_time = establishment_time)

			plot_row <- if (sim <= max_plot_sims && established == 1L) {
				tibble(sim    = sim,
					   psi    = psi,
					   model  = model,
					   tinf   = infection_times,
					   cuminf = seq_along(infection_times))
			} else NULL

			list(summary = summary_row, plot = plot_row)
		}

		results <- mclapply(seq_len(nrow(task_grid)), run_one,
							mc.cores = n_cores, mc.set.seed = TRUE)

		errs <- vapply(results, function(x) inherits(x, "try-error"), logical(1))
		if (any(errs))
			stop(sprintf("%s: %d sensitivity task(s) failed; first error:\n%s",
						 pathogen, sum(errs), as.character(results[[which(errs)[1]]])))

		sim_df  <- bind_rows(lapply(results, `[[`, "summary")) %>%
			mutate(pathogen = pathogen)
		plot_df_raw <- bind_rows(lapply(results, `[[`, "plot"))

		write_csv(sim_df,      sensitivity_cache_path(pathogen, nsim_small, popsize))
		write_csv(plot_df_raw, sensitivity_cache_path_plot(pathogen, nsim_small, popsize))
		cat(sprintf("  %s: sensitivity simulations saved\n", pathogen))

		plot_df <- plot_df_raw %>%
			mutate(psi   = factor(psi,   levels = c(0, 0.5, 1)),
				   model = factor(model, levels = models, labels = model_labels[models]))
	}

	all_summaries[[pathogen]] <- sim_df

	# ------------------------------------------------------------------
	# Aggregate plot subset to daily resolution
	# ------------------------------------------------------------------
	lastday <- ceiling(max(plot_df$tinf))

	dayjoin <- expand_grid(
		psi   = unique(plot_df$psi),
		model = unique(plot_df$model),
		sim   = unique(plot_df$sim),
		day   = 0:lastday)

	dailyinf_df <- plot_df %>%
		mutate(day = floor(tinf)) %>%
		group_by(psi, model, sim, day) %>%
		summarise(ninf = n(), .groups = "drop") %>%
		right_join(dayjoin, by = c("psi", "model", "sim", "day")) %>%
		replace_na(list(ninf = 0)) %>%
		group_by(psi, model, sim) %>%
		arrange(day, .by_group = TRUE) %>%
		mutate(cuminf = cumsum(ninf))

	# ------------------------------------------------------------------
	# Figures: cumulative curves (3x3 grid: rows = burst model, cols = psi)
	# ------------------------------------------------------------------
	fig_cuminf <- plot_df %>%
		filter(tinf <= xlim_max) %>%
		ggplot(aes(x = tinf, y = cuminf, group = sim)) +
			geom_line(alpha = 0.2, colour = "grey") +
			geom_line(data = filter(ren_out, t <= xlim_max),
					  aes(x = t, y = cuminf * popsize, group = model),
					  colour = "black", linewidth = 1, alpha = 0.8,
					  inherit.aes = FALSE) +
			facet_grid(model ~ psi,
					   labeller = labeller(psi = as_labeller(function(x) paste0("psi == ", x), label_parsed))) +
			coord_cartesian(xlim = c(0, xlim_max)) +
			theme_classic(base_size = 13) +
			labs(x = "Time (days)", y = "Cumulative infections"
				 # title = sprintf("Sensitivity: cumulative infections — %s", pathogen)
				 )

	save_fig(fig_cuminf, sprintf("sensitivity_episims_cuminf_%s", pathogen),
			 width = 10, height = 8)

	# ------------------------------------------------------------------
	# Figures: daily incidence (3x3 grid: rows = burst model, cols = psi)
	# ------------------------------------------------------------------
	fig_daily <- dailyinf_df %>%
		filter(day <= xlim_max) %>%
		ggplot(aes(x = day, y = ninf, group = sim)) +
			geom_line(alpha = 0.2, colour = "grey") +
			geom_line(data = filter(ren_daily, day <= xlim_max),
					  aes(x = day, y = newinf * popsize, group = model),
					  colour = "black", linewidth = 1, alpha = 0.8,
					  inherit.aes = FALSE) +
			facet_grid(model ~ psi,
					   labeller = labeller(psi = as_labeller(function(x) paste0("psi == ", x), label_parsed))) +
			coord_cartesian(xlim = c(0, xlim_max)) +
			theme_classic(base_size = 13) +
			labs(x = "Time (days)", y = "Daily new infections"
				 # title = sprintf("Sensitivity: daily incidence — %s", pathogen)
				 )

	save_fig(fig_daily, sprintf("sensitivity_episims_daily_%s", pathogen),
			 width = 10, height = 8)

	cat(sprintf("  %s: trajectory figures saved.\n", pathogen))
}
cat("\nSensitivity epidemic simulations complete.\n")
