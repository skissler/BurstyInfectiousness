# ==============================================================================
# Overdispersion: k heatmaps (periodic and Gamma/Poisson contact models)
# ==============================================================================

library(tidyverse)
source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

c_per             <- 7   # Period of the periodic contact model (days)
c_amp_vals_sim    <- seq(0, 1, length.out=25)
psi_vals_sim      <- seq(0, 1, length.out=25)
c_amp_vals_theory <- seq(0, 1, length.out=100)
psi_vals_theory   <- seq(0, 1, length.out=100)

lambda_gp        <- 1    # Poisson switching rate (switches/day)
k_c_vals_sim     <- exp(seq(log(0.1), log(1000), length.out=25))
k_c_vals_theory  <- exp(seq(log(0.1), log(1000), length.out=100))

n_index <- 10000
k_cap   <- 50
od_cap  <- 1

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Start pathogen loop
# ==============================================================================

for (idx_pathogen in seq_along(parslist)) {

pars <- parslist[[idx_pathogen]]
pathogen <- pars$pathogen
pathogen_label <- switch(pathogen,
    influenza = "Influenza",
    omicron   = "SARS-CoV-2 Omicron",
    measles   = "Measles",
    pathogen
)
Tgen     <- pars$Tgen
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

# Per-pathogen seed so results don't depend on the order of earlier pathogens.
set.seed(42 + idx_pathogen)

cat(sprintf("\n===== %s: T=%.2f, alpha=%.2f, R0=%.1f =====\n",
            pathogen, Tgen, alpha, R0))

# ==============================================================================
# Periodic contacts: heatmaps
# ==============================================================================

periodic_sim_cache <- file.path(
	"output",
	sprintf("k_heatmap_periodic_sim_%s_n%d_cper%g.csv", pathogen, n_index, c_per)
)

need_combos <- expand_grid(psi = psi_vals_sim, c_amp = c_amp_vals_sim)
sim_grid_periodic <- NULL
if (file.exists(periodic_sim_cache)) {
	cached <- read_csv(periodic_sim_cache, show_col_types = FALSE)
	missing <- anti_join(need_combos, cached %>% distinct(psi, c_amp),
	                     by = c("psi", "c_amp"))
	if (nrow(missing) == 0) {
		cat(sprintf("  %s: loaded periodic k_sim grid from %s\n", pathogen, periodic_sim_cache))
		sim_grid_periodic <- cached
	}
}

if (is.null(sim_grid_periodic)) {
	sim_grid_periodic <- need_combos %>% mutate(k_sim = NA_real_)

	for(idx in seq_len(nrow(sim_grid_periodic))){
		psi   <- sim_grid_periodic$psi[idx]
		c_amp <- sim_grid_periodic$c_amp[idx]
		z     <- make_contact_fn_periodic(R0, c_amp, c_per)
		z_max <- R0*(1 + c_amp)
		gfun  <- gen_inf_attempts_gamma_contacts(Tgen=Tgen, z=z, z_max=z_max,
		                                         alpha=alpha, psi=psi)
		tinfs      <- c_per * runif(n_index)
		noffspring <- sapply(lapply(tinfs, gfun), length)
		m_s <- mean(noffspring)
		v_s <- var(noffspring)
		sim_grid_periodic$k_sim[idx] <- if(v_s > m_s){ m_s^2 / (v_s - m_s) } else { Inf }
		if(idx %% 50 == 0) cat(sprintf("  [%s] heatmap: %d / %d\n", pathogen, idx, nrow(sim_grid_periodic)))
	}

	write_csv(sim_grid_periodic, periodic_sim_cache)
}

theory_grid_periodic <- expand_grid(
	    psi = psi_vals_theory,
	    c_amp = c_amp_vals_theory) %>%
	mutate(k_theory = (2/c_amp^2)*(1 + (2*pi/c_per)^2 / beta^2)^(psi*alpha))

sim_grid_periodic <- sim_grid_periodic %>%
	mutate(k_capped  = pmin(k_sim, k_cap),
	       od_sim    = 1 / k_sim,
	       od_capped = pmin(od_sim, od_cap))

theory_grid_periodic <- theory_grid_periodic %>%
	mutate(k_capped  = pmin(k_theory, k_cap),
	       od_theory = 1 / k_theory,
	       od_capped = pmin(od_theory, od_cap))

pest_periodic_path <- file.path(
	"output",
	sprintf("pest_periodic_%s_n%d_p%d_cper%g.csv", pathogen, nsim, popsize, c_per)
)
pest_df <- if (file.exists(pest_periodic_path)) read_csv(pest_periodic_path, show_col_types = FALSE) else NULL

fig_heatmap_periodic <- ggplot() +
	geom_tile(data = sim_grid_periodic, aes(x=psi, y=c_amp, fill=k_capped)) +
	geom_contour(data = theory_grid_periodic,
		aes(x=psi, y=c_amp, z=k_theory),
		col = "black", linewidth = 1.0, breaks = c(1, 2, 5, 10, 20, 50)) +
	geom_contour(data = theory_grid_periodic,
		aes(x=psi, y=c_amp, z=k_theory),
		col = "white", linewidth = 0.4, breaks = c(1, 2, 5, 10, 20, 50)) +
	scale_fill_viridis_c(option = "inferno", name = sprintf("k\n(capped\nat %d)", k_cap), limits = c(0, k_cap)) +
	theme_classic(base_size = 14) +
	labs(x = expression(psi),
	     y = expression("Contact amplitude " * italic(A)),
	     title = sprintf("%s (R0 = %g)", pathogen_label, R0))
if (!is.null(pest_df)) fig_heatmap_periodic <- fig_heatmap_periodic +
	geom_point(data = pest_df, aes(x = psi, y = c_amp, size = p_extinct),
	           shape = 21, fill = "white", color = "black", alpha = 0.8) +
	scale_size_area(max_size = 10, name = "Extinction prob.")

save_fig(fig_heatmap_periodic, sprintf("fig_heatmap_periodic_%s", pathogen), width=5, height=4)
cat(sprintf("  Saved fig_heatmap_periodic_%s\n", pathogen))

fig_heatmap_periodic_od <- ggplot() +
	geom_tile(data = sim_grid_periodic, aes(x=psi, y=c_amp, fill=od_capped)) +
	geom_contour(data = theory_grid_periodic,
		aes(x=psi, y=c_amp, z=od_theory),
		col = "black", linewidth = 1.0, breaks = c(0.02, 0.05, 0.1, 0.2, 0.5, 1)) +
	geom_contour(data = theory_grid_periodic,
		aes(x=psi, y=c_amp, z=od_theory),
		col = "white", linewidth = 0.4, breaks = c(0.02, 0.05, 0.1, 0.2, 0.5, 1)) +
	scale_fill_viridis_c(option = "inferno", name = "OD\n(1/k)", limits = c(0, od_cap)) +
	theme_classic(base_size = 14) +
	labs(x = expression(psi),
	     y = expression("Contact amplitude " * italic(A)),
	     title = sprintf("%s (R0 = %g)", pathogen_label, R0))
if (!is.null(pest_df)) fig_heatmap_periodic_od <- fig_heatmap_periodic_od +
	geom_point(data = pest_df, aes(x = psi, y = c_amp, size = p_extinct),
	           shape = 21, fill = "white", color = "black", alpha = 0.8) +
	scale_size_area(max_size = 10, name = "Extinction prob.")

save_fig(fig_heatmap_periodic_od, sprintf("fig_heatmap_periodic_od_%s", pathogen), width=5, height=4)
cat(sprintf("  Saved fig_heatmap_periodic_od_%s\n", pathogen))

# ==============================================================================
# Gamma/Poisson contacts: heatmaps
# ==============================================================================

gp_sim_cache <- file.path(
	"output",
	sprintf("k_heatmap_gp_sim_%s_n%d_lam%g.csv", pathogen, n_index, lambda_gp)
)

need_combos_gp <- expand_grid(psi = psi_vals_sim, k_c = k_c_vals_sim)
sim_grid_gp <- NULL
if (file.exists(gp_sim_cache)) {
	cached <- read_csv(gp_sim_cache, show_col_types = FALSE)
	have <- cached %>% distinct(psi, k_c)
	# Floating-point tolerance on k_c because it's drawn on a log scale.
	missing_psi <- !(psi_vals_sim %in% have$psi)
	missing_kc  <- sapply(k_c_vals_sim, function(v) !any(abs(have$k_c - v) < 1e-8 * max(1, v)))
	if (!any(missing_psi) && !any(missing_kc)) {
		cat(sprintf("  %s: loaded GP k_sim grid from %s\n", pathogen, gp_sim_cache))
		sim_grid_gp <- cached
	}
}

if (is.null(sim_grid_gp)) {
	sim_grid_gp <- need_combos_gp %>% mutate(k_sim = NA_real_)

	for(idx in seq_len(nrow(sim_grid_gp))){
		psi <- sim_grid_gp$psi[idx]
		k_c <- sim_grid_gp$k_c[idx]
		gfun <- gen_inf_attempts_gammapoisson_contacts(Tgen, R0, alpha, psi, k_c, lambda_gp)
		noffspring <- replicate(n_index, length(gfun(0)))
		m_s <- mean(noffspring)
		v_s <- var(noffspring)
		sim_grid_gp$k_sim[idx] <- if(v_s > m_s){ m_s^2 / (v_s - m_s) } else { Inf }
		if(idx %% 50 == 0) cat(sprintf("  [%s] GP heatmap: %d / %d\n", pathogen, idx, nrow(sim_grid_gp)))
	}

	write_csv(sim_grid_gp, gp_sim_cache)
}

theory_grid_gp <- expand_grid(psi=psi_vals_theory, k_c=k_c_vals_theory) %>%
	mutate(k_theory = mapply(k_theory_gammapoisson, k_c, psi, alpha, lambda_gp, beta))

sim_grid_gp <- sim_grid_gp %>%
	mutate(k_capped  = pmin(k_sim, k_cap),
	       od_sim    = 1 / k_sim,
	       od_capped = pmin(od_sim, od_cap))

theory_grid_gp <- theory_grid_gp %>%
	mutate(k_capped  = pmin(k_theory, k_cap),
	       od_theory = 1 / k_theory,
	       od_capped = pmin(od_theory, od_cap))

# Overlay: extinction probabilities from overdispersion_extinction.R, one dot
# per (psi, k_c) pair sampled there.
pest_gp_path <- file.path(
	"output",
	sprintf("pest_gammapoisson_%s_n%d_p%d_lam%g.csv",
	        pathogen, nsim, popsize, lambda_gp)
)
pest_df_gp <- if (file.exists(pest_gp_path)) read_csv(pest_gp_path, show_col_types = FALSE) else NULL

fig_heatmap_gp <- ggplot() +
	geom_tile(data = sim_grid_gp, aes(x=psi, y=k_c, fill=k_capped)) +
	geom_contour(data = theory_grid_gp,
		aes(x=psi, y=k_c, z=k_theory),
		col = "black", linewidth = 1.0, breaks = c(1, 2, 5, 10, 20, 50)) +
	geom_contour(data = theory_grid_gp,
		aes(x=psi, y=k_c, z=k_theory),
		col = "white", linewidth = 0.4, breaks = c(1, 2, 5, 10, 20, 50)) +
	scale_y_log10() +
	scale_fill_viridis_c(option = "inferno", name = sprintf("k\n(capped\nat %d)", k_cap), limits = c(0, k_cap)) +
	theme_classic(base_size = 14) +
	labs(x = expression(psi),
	     y = expression("Contact shape " * sigma),
	     # title = sprintf("%s: Gamma/Poisson contacts (R0 = %g, λ = %g)", pathogen_label, R0, lambda_gp))
	     title = sprintf("%s", pathogen_label))
if (!is.null(pest_df_gp)) fig_heatmap_gp <- fig_heatmap_gp +
	geom_point(data = pest_df_gp, aes(x = psi, y = k_c, size = p_extinct),
	           shape = 21, fill = "white", color = "black", alpha = 0.8) +
	scale_size_area(max_size = 10, name = "Extinction prob.")

save_fig(fig_heatmap_gp, sprintf("fig_heatmap_gammapoisson_%s", pathogen), width=5, height=4)
cat(sprintf("  Saved fig_heatmap_gammapoisson_%s\n", pathogen))

fig_heatmap_gp_od <- ggplot() +
	geom_tile(data = sim_grid_gp, aes(x=psi, y=k_c, fill=od_capped)) +
	geom_contour(data = theory_grid_gp,
		aes(x=psi, y=k_c, z=od_theory),
		col = "black", linewidth = 1.0, breaks = c(0.02, 0.05, 0.1, 0.2, 0.5, 1)) +
	geom_contour(data = theory_grid_gp,
		aes(x=psi, y=k_c, z=od_theory),
		col = "white", linewidth = 0.4, breaks = c(0.02, 0.05, 0.1, 0.2, 0.5, 1)) +
	scale_y_log10() +
	scale_fill_viridis_c(option = "inferno", name = "OD\n(1/k)", limits = c(0, od_cap)) +
	theme_classic(base_size = 14) +
	labs(x = expression(psi),
	     y = expression("Contact shape " * sigma),
	     title = sprintf("%s: Gamma/Poisson contacts (R0 = %g, λ = %g)", pathogen_label, R0, lambda_gp))
if (!is.null(pest_df_gp)) fig_heatmap_gp_od <- fig_heatmap_gp_od +
	geom_point(data = pest_df_gp, aes(x = psi, y = k_c, size = p_extinct),
	           shape = 21, fill = "white", color = "black", alpha = 0.8) +
	scale_size_area(max_size = 10, name = "Extinction prob.")

save_fig(fig_heatmap_gp_od, sprintf("fig_heatmap_gammapoisson_od_%s", pathogen), width=5, height=4)
cat(sprintf("  Saved fig_heatmap_gammapoisson_od_%s\n", pathogen))

} # end pathogen loop
