# =====================================================================================

# ==============================================================================
# COVID-19: 
#   Source: Li et al. (2020) NEJM [https://www.nejm.org/doi/full/10.1056/NEJMoa2001316]
#   Mean: 5.2
#   95% percentile: 12.5

target_mean <- 5.2
target_p95  <- 12.5

f <- function(sigma) {
	mu <- log(target_mean) - sigma^2 / 2
	qlnorm(0.95, mu, sigma) - target_p95
}

sigma_hat <- uniroot(f, interval = c(0.01, 2))$root
mu_hat    <- log(target_mean) - sigma_hat^2 / 2

# Now recover mean and SD of the distribution
mean_hat <- exp(mu_hat + sigma_hat^2 / 2)
sd_hat   <- sqrt((exp(sigma_hat^2) - 1) * exp(2*mu_hat + sigma_hat^2))