library(tidyverse)
library(odin)
library(parallel)
source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

# ==============================================================================
# Set parameters
# ==============================================================================

psivals <- c(0, 0.5, 1)

# Coerce establishment_threshold to integer so it's safe for index-based
# lookups (infection_times[establishment_threshold]) under any choice of
# popsize * establishment_prop.
establishment_threshold <- as.integer(round(establishment_threshold))

# Number of worker cores for the parallel simulation loop. Override by setting
# options(mc.cores=...) before sourcing, or the MC_CORES env var. Defaults to
# all but one core. On Windows, mclapply silently falls back to serial.
n_cores <- getOption("mc.cores",
                     as.integer(Sys.getenv("MC_CORES",
                                           max(1L, detectCores() - 1L))))

# ==============================================================================
# Loop over pathogens
# ==============================================================================

for (pars in parslist) {

pathogen <- pars$pathogen
pathogen_label <- switch(pathogen,
    influenza = "Influenza",
    omicron   = "SARS-CoV-2 Omicron",
    measles   = "Measles",
    pathogen
)
Tgen     <- pars$Tgen
alpha    <- pars$alpha
R0       <- pars$R0

cat(sprintf("\n===== %s: T=%.2f, alpha=%.2f, R0=%.1f =====\n",
            pathogen, Tgen, alpha, R0))

# ==============================================================================
# Generate mean-field solution
# ==============================================================================

# Run the renewal equation model
ren_out <- renewal_epidemic(R0, alpha, Tgen, popsize)

# Aggregate renewal equation output to daily new-infection counts
ren_daily <- ren_out %>%
	mutate(day=floor(t)) %>%
	group_by(day) %>%
	summarise(cuminf=max(cuminf)) %>%
	arrange(day) %>%
	mutate(newinf=cuminf-lag(cuminf)) %>%
	mutate(newinf=case_when(is.na(newinf)~cuminf, TRUE~newinf))

# ==============================================================================
# Generate stochastic simulations (use cache if it exists)
# ==============================================================================

cache <- load_cache(pathogen, nsim, popsize, psivals)

if (!is.null(cache)) {
	sim_summary_df <- cache$summary
	plot_df        <- cache$plot
} else {
	# Run simulations in parallel over the (sim, psi) grid. Each task is fully
	# independent, so we fork the work across n_cores via mclapply.
	#   - Compute per-sim summary immediately, discard raw infection times
	#   - Keep full trajectories only for sim <= max_plot_sims (for plotting)
	# L'Ecuyer-CMRG gives independent, reproducible RNG streams across forks.
	task_grid <- expand_grid(sim = 1:nsim, psi = psivals)

	# Deterministic, pathogen-specific seed so streams are reproducible across
	# runs but distinct across pathogens.
	RNGkind("L'Ecuyer-CMRG")
	set.seed(sum(utf8ToInt(pathogen)) + 1L)

	run_one <- function(i){
	    sim <- task_grid$sim[i]
	    psi <- task_grid$psi[i]

	    tinf <- sim_stochastic_fast(n=popsize,
	                                gen_inf_attempts=gen_inf_attempts_gamma(Tgen, R0, alpha, psi))

	    # Extract sorted finite infection times
	    infection_times <- sort(tinf[is.finite(tinf)])
	    n_inf <- length(infection_times)
	    established <- as.integer(n_inf >= establishment_threshold)

	    # Compute time to establishment threshold
	    establishment_time <- if(established == 1) infection_times[establishment_threshold] else NA_real_

	    summary_row <- tibble(
	        sim = sim,
	        psi = psi,
	        established = established,
	        establishment_time = establishment_time
	    )

	    # Keep full trajectories for a subset of established epidemics for plotting
	    plot_row <- if(sim <= max_plot_sims && established == 1){
	        tibble(
	            sim = sim,
	            psi = psi,
	            tinf = infection_times,
	            cuminf = seq_along(infection_times),
	        )
	    } else NULL

	    list(summary = summary_row, plot = plot_row)
	}

	results <- mclapply(seq_len(nrow(task_grid)), run_one,
	                    mc.cores = n_cores, mc.set.seed = TRUE)

	# Surface worker errors (mclapply returns try-error objects rather than stopping)
	errs <- vapply(results, function(x) inherits(x, "try-error"), logical(1))
	if(any(errs)) stop(sprintf("%s: %d parallel sim task(s) failed; first error:\n%s",
	                           pathogen, sum(errs), as.character(results[[which(errs)[1]]])))

	sim_summary_df <- bind_rows(lapply(results, `[[`, "summary"))
	plot_df        <- bind_rows(lapply(results, `[[`, "plot"))

	# Save to cache
	write_csv(sim_summary_df, cache_path_summary(pathogen, nsim, popsize))
	write_csv(plot_df, cache_path_plot(pathogen, nsim, popsize))
	cat(sprintf("  %s: simulations saved\n", pathogen))

	sim_summary_df <- sim_summary_df %>% mutate(psi = factor(psi))
	plot_df        <- plot_df %>% mutate(psi = factor(psi))
}

# ==============================================================================
# Aggregate plot subset to daily resolution (for trajectory plots)
# ==============================================================================

lastday <- ceiling(max(plot_df$tinf))

dayjoin <- expand_grid(
	psi=unique(plot_df$psi),
	sim=unique(plot_df$sim),
	day=0:lastday)

dailyinf_df <- plot_df %>%
	mutate(day = floor(tinf)) %>%
	group_by(psi, sim, day) %>%
	summarise(ninf = n(), .groups = "drop") %>%
	right_join(dayjoin, by = c("psi", "sim", "day")) %>%
	replace_na(list(ninf = 0)) %>%
	group_by(psi, sim) %>%
	arrange(day, .by_group = TRUE) %>%
	mutate(cuminf = cumsum(ninf))

# ==============================================================================
# Report establishment times 
# ==============================================================================

est_time_table <- sim_summary_df %>%
	filter(established==1, !is.na(establishment_time)) %>%
	group_by(psi) %>%
	summarise(
		et_mean = mean(establishment_time),
		et_sd   = sd(establishment_time),
		et_q05  = quantile(establishment_time, 0.05),
		et_q95  = quantile(establishment_time, 0.95),
		.groups = "drop")

cat(sprintf("  %s: Time to %d cases (established epidemics):\n", pathogen, establishment_threshold))
print(est_time_table)

# ==============================================================================
# Figures — epidemic trajectories (plot subset only)
# ==============================================================================

# Cumulative stochastic epidemic curves (grey) with ODE overlay (black)
fig_cuminf_overlay <- plot_df %>%
	ggplot(aes(x=tinf, y=cuminf, group=sim)) +
		geom_line(alpha=0.2, col="grey") +
		geom_line(data=filter(ren_out, t<=lastday),
			aes(x=t, y=cuminf*popsize, group=1),
			alpha=0.8, linewidth=1, col="black") +
		theme_classic(base_size = 14) +
		theme(strip.background = element_blank()) +
		labs(x="Time (days)", y="Cumulative infections", title = pathogen_label) +
		facet_wrap(~psi, nrow=1,
		           labeller = label_bquote(psi == .(as.numeric(as.character(psi)))))

save_fig(fig_cuminf_overlay, paste0("fig_cuminf_overlay_", pathogen))

# Cumulative curves with time-to-threshold milestone annotations
milestone_df <- plot_df %>%
	filter(cuminf >= establishment_threshold) %>%
	group_by(sim, psi) %>%
	slice(1) %>%
	group_by(psi) %>%
	summarise(
		mean_t = mean(tinf),
		q05_t  = quantile(tinf, 0.05),
		q95_t  = quantile(tinf, 0.95),
		.groups  = "drop")

# cat(sprintf("  %s: Time to reach 5% of the population infected:\n", pathogen))
# print(milestone_df, n=Inf)

fig_cuminf_milestones <- plot_df %>%
	ggplot(aes(x=tinf, y=cuminf, group=sim)) +
		geom_line(alpha=0.2, col="grey") +
		geom_line(data=filter(ren_out, t<=lastday),
			aes(x=t, y=cuminf*popsize, group=1),
			alpha=0.8, linewidth=1, col="black") +
		geom_segment(data=milestone_df,
			aes(x=q05_t, xend=q05_t,
			    y=establishment_threshold*0.85, yend=establishment_threshold*1.15, group=1),
			col="red", linewidth=0.6) +
		geom_segment(data=milestone_df,
			aes(x=q95_t, xend=q95_t,
			    y=establishment_threshold*0.85, yend=establishment_threshold*1.15, group=1),
			col="red", linewidth=0.6) +
		geom_point(data=milestone_df,
			aes(x=mean_t, y=establishment_threshold, group=1),
			col="red", size=2.5) +
		theme_classic(base_size = 14) +
		labs(x="Time (days)", y="Cumulative infections", title = pathogen) +
		facet_wrap(~psi, nrow=1)

save_fig(fig_cuminf_milestones, paste0("fig_cuminf_milestones_", pathogen))

# Daily stochastic incidence curves (grey) with ODE overlay (black)
fig_inf_overlay <- dailyinf_df %>%
	ggplot(aes(x=day, y=ninf, group=sim)) +
		geom_line(alpha=0.2, col="grey") +
		geom_line(data=filter(ren_daily, day<=lastday),
			aes(x=day, y=newinf*popsize, group=1),
			alpha=0.8, linewidth=1, col="black") +
		theme_classic(base_size = 14) +
		theme(strip.background = element_blank()) +
		labs(x="Time (days)", y="Daily new infections", title = pathogen_label) +
		facet_wrap(~psi, nrow=1,
		           labeller = label_bquote(psi == .(as.numeric(as.character(psi)))))

save_fig(fig_inf_overlay, paste0("fig_inf_overlay_", pathogen))

cat(sprintf("  %s: figures saved.\n", pathogen))

} # end pathogen loop
