# This script analyzes the simulated and fitted models from simulate.R
library(BDAepimodel)

# Plot the simulated disease dynamics
plot_simulated_dynamics <- function(simulated_epimodel, popsize = 500) {
  plot(x = simulated_epimodel$pop_mat[,"time"], y = simulated_epimodel$pop_mat[,"I"], 
       "l", xlab = "Time", ylab = "Prevalence", main = "Simulated disease dynamics", 
       ylim = c(0, popsize))
  points(x = simulated_epimodel$dat[,"time"], y = simulated_epimodel$dat[,"I"])
}

# Plot the fitted disease dynamics
plot_fitted_dynamics <- function(fitted_epimodel) {
  plot_latent_posterior(fitted_epimodel, states = "I", times = fitted_epimodel$obstimes, cm = "mn")
}

# Trace plots for beta, mu, and rho
plot_trace <- function(fitted_epimodel, burn_in = 10000) {
  # Beta
  ts.plot(fitted_epimodel$results$params[ , "beta"], ylab = expression("Sampled" ~ beta), 
          gpars = list(xlab = "Iteration", cex.lab = 5))
  abline(v = burn_in, col = "red", lty = 2)
  
  # Mu
  ts.plot(fitted_epimodel$results$params[ , "mu"], ylab = expression("Sampled" ~ mu), 
          gpars = list(xlab = "Iteration", cex.lab = 5))
  abline(v = burn_in, col = "red", lty = 2)
  
  # Rho
  ts.plot(fitted_epimodel$results$params[ , "rho"], ylab = expression("Sampled" ~ rho), 
          gpars = list(xlab = "Iteration", cex.lab = 5))
  abline(v = burn_in, col = "red", lty = 2)
}

# Obtain posterior median, credible interval, and true value for beta, mu, and rho
get_posmed <- function(fitted_epimodel, burn_in = 10000, credible_lvl = 0.05) {
  burnt_idx <- 1:burn_in
  quants <- c(credible_lvl/2, 1 - credible_lvl/2)
  
  # Beta 
  betas <- fitted_epimodel$results$params[-burnt_idx , "beta"]
  beta_med <- signif(median(betas), 3)
  beta_ci <- signif(quantile(betas, probs = quants), 3)
  
  # Mu
  mus <- fitted_epimodel$results$params[-burnt_idx, "mu"]
  mu_med <- signif(median(mus), 3)
  mu_ci <- signif(quantile(mus, probs = quants), 3)
  
  # Rho
  rhos <- fitted_epimodel$results$params[-burnt_idx, "rho"]
  rho_med <- signif(median(rhos), 3)
  rho_ci <- signif(quantile(rhos, probs = quants), 3)
  
  # Print results
  print(paste0("beta: ", beta_med, " (", beta_ci[1], " - ", beta_ci[2], ")"))
  print(paste0("mu: ", mu_med, " (", mu_ci[1], " - ", mu_ci[2], ")"))
  print(paste0("rho: ", rho_med, " (", rho_ci[1], " - ", rho_ci[2], ")"))
}

# Output the result plots and numbers
result_dirs <- c("results/beta=0.5_mu=0.005_rho=0.8_popsize=500_niter=50000_obsint=7/", 
                 "results/beta=0.005_mu=0.005_rho=0.8_popsize=500_niter=50000_obsint=7/", 
                 "results/beta=0.005_mu=0.5_rho=0.8_popsize=500_niter=50000_obsint=7/")

# Analysis
for (result_dir in result_dirs) {
  simulated_epimodel <- readRDS(paste0(result_dir, "simulated_epimodel.rds"))
  fitted_epimodel <- readRDS(paste0(result_dir, "fitted_epimodel.rds"))
  # Plot and save disease dynamics
  png(filename = paste0(result_dir, "simulated_dynamics.png"))
  plot_simulated_dynamics(simulated_epimodel)
  dev.off()
  png(filename = paste0(result_dir, "fitted_dynamics.png"), width = 600)
  plot_fitted_dynamics(fitted_epimodel)
  dev.off()
  
  # Plot and save trace plots
  png(filename = paste0(result_dir, "traceplots.png"), width = 600, height = 600)
  par(mfrow = c(3, 1))
  plot_trace(fitted_epimodel)
  dev.off()
  
  # Print results
  print(paste0("Result dir: ", result_dir))
  get_posmed(fitted_epimodel)
}
