#!/usr/bin/env Rscript
library(optparse)
library(BDAepimodel)

# Command line arguments ------------------------------
option_list = list(make_option(c("-b", "--beta"), type = "double", default = 0.00035, 
                               help = "Beta infection rate [default = %default]", metavar = "double"), 
                   make_option(c("-m", "--mu"), type = "double", default = 0.15, 
                               help = "Mu recovery rate [default = %default]", metavar = "double"), 
                   make_option(c("-r", "--rho"), type = "double", default = 0.2, 
                               help = "Rho disease detection probability [default = %default]", metavar = "double"), 
                   make_option(c("-N", "--popsize"), type = "integer", default = 20, 
                               help = "Population size [default = %default]", metavar = "integer"), 
                   make_option(c("-t", "--niter"), type = "integer", default = 100, 
                               help = "Number of iterations for MCMC [default = %default]", metavar = "integer"),
                   make_option(c("-o", "--obsint"), type = "integer", default = 7, 
                               help = "Length of observation interval [default = %default] (total period = [0, 105])", metavar = "integer"), 
                   make_option(c("-v", "--verbose"), action = "store_true", default = FALSE, 
                               help = "Flag for monitoring MCMC progress")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Read in the command line arguments
beta <- opt$beta
mu <- opt$mu
rho <- opt$rho
popsize <- opt$popsize
niter <- opt$niter
obsint <- opt$obsint
verbose <- opt$verbose

# NOTE: Code from the vignette The BDAepimodel package for fitting stochastic epidemic models to 
# disease prevalence data via Bayesian data augmentaiton
# by Jon Fintzi, Xiang Cui, Jon Wakefield, Vladimir Minin

# Simulate data ------------------------------
# Randomly sample for the measurement process
r_meas_process <- function(state, meas_vars, params){
  # in our example, rho will be the name of the binomial sampling probability parameter.
  # this function returns a matrix of observed counts
  rbinom(n = nrow(state), 
         size = state[,meas_vars], # binomial sample of the unobserved prevalenc
         prob = params["rho"])     # sampling probability
}

# Density for the measurement process
d_meas_process <- function(state, meas_vars, params, log = TRUE) {
  # note that the names of the measurement variables are endowed with suffixes "_observed" and "_augmented". This is required.
  # we will declare the names of the measurement variables shortly.
  dbinom(x = state[, "I_observed"], 
         size = state[, "I_augmented"], 
         prob = params["rho"], log = log)
}

# Initialize the stochastic epidemic dynamic object
epimodel <- init_epimodel(obstimes = seq(0, 105, by = obsint),                        # vector of observation times
                          popsize = popsize,                                          # population size
                          states = c("S", "I", "R"),                                  # compartment names
                          params = c(beta = beta,                                     # infectivity parameter
                                     mu = mu,                                         # recovery rate
                                     rho = rho,                                       # binomial sampling probability
                                     S0 = 0.9, I0 = 0.03, R0 = 0.07),                 # initial state probabilities
                          rates = c("beta * I", "mu"),                                # unlumped transition rates
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T),  # flow matrix
                          meas_vars = "I",                                            # name of measurement variable
                          r_meas_process = r_meas_process,                            # measurement process functions
                          d_meas_process = d_meas_process)

# Simulate the epidemic and observed data
simulated_epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)

# Rcpp helper function for computing sufficient statistics ------------------------------
Rcpp::cppFunction("Rcpp::NumericVector getSuffStats_SIR(const Rcpp::NumericMatrix& pop_mat, const int ind_final_config) {
                  
                  // initialize sufficient statistics
                  int num_inf = 0;       // number of infection events
                  int num_rec = 0;       // number of recovery events
                  double beta_suff = 0;  // integrated hazard for the infectivity
                  double mu_suff = 0;    // integrated hazard for the recovery
                  // initialize times
                  double cur_time = 0;              // current time
                  double next_time = pop_mat(0,0);  // time of the first event
                  double dt = 0;                    // time increment
                  
                  // compute the sufficient statistics - loop through the pop_mat matrix until
                  // reaching the row for the final observation time
                  for(int j = 0; j < ind_final_config - 1; ++j) {
                  
                  cur_time = next_time;         
                  next_time = pop_mat(j+1, 0); // grab the time of the next event
                  dt = next_time - cur_time;   // compute the time increment
                  
                  beta_suff += pop_mat(j, 3) * pop_mat(j, 4) * dt; // add S*I*(t_{j+1} - t_j) to beta_suff
                  mu_suff += pop_mat(j, 4) * dt;                   // add I*(t_{j+1} - t_j) to mu_suff
                  
                  // increment the count for the next event
                  if(pop_mat(j + 1, 2) == 1) {  
                  num_inf += 1;
                  } else if(pop_mat(j + 1, 2) == 2) {
                  num_rec += 1;
                  }
                  }
                  
                  // return the vector of sufficient statistics for the rate parameters
                  return Rcpp::NumericVector::create(num_inf, beta_suff, num_rec, mu_suff);
                  }")

# Gibbs sampler ------------------------------
gibbs_kernel_SIR <- function(epimodel) {
  
  # get sufficient statistics using the previously compiled getSuffStats_SIR function (above)
  suff_stats <- getSuffStats_SIR(epimodel$pop_mat, epimodel$ind_final_config)
  
  # update parameters from their univariate full conditional distributions
  # Priors: beta ~ gamma(0.3, 1000)
  #         mu   ~ gamma(1, 8)
  #         rho  ~ beta(2, 7)
  proposal          <- epimodel$params # params is the vector of ALL model parameters
  proposal["beta"]  <- rgamma(1, 0.3 + suff_stats[1], 1000 + suff_stats[2])
  proposal["mu"]    <- rgamma(1, 1 + suff_stats[3], 8 + suff_stats[4])
  proposal["rho"]   <- rbeta(1, shape1 = 2 + sum(epimodel$obs_mat[, "I_observed"]),
                             shape2 = 7 + sum(epimodel$obs_mat[, "I_augmented"] - epimodel$obs_mat[, "I_observed"]))
  
  # update array of rate matrices
  epimodel <- build_new_irms(epimodel, proposal)
  
  # update the eigen decompositions (This function is built in and computes eigen decompositions analytically)
  buildEigenArray_SIR(real_eigenvals = epimodel$real_eigen_values,
                      imag_eigenvals = epimodel$imag_eigen_values,
                      eigenvecs      = epimodel$eigen_vectors, 
                      inversevecs    = epimodel$inv_eigen_vectors, 
                      irm_array      = epimodel$irm, 
                      n_real_eigs    = epimodel$n_real_eigs, 
                      initial_calc   = FALSE)
  
  # get log-likelihood of the observations under the new parameters
  obs_likelihood_new  <- calc_obs_likelihood(epimodel, params = proposal, log = TRUE) #### NOTE - log = TRUE
  
  # get the new population level CTMC log-likelihood
  pop_likelihood_new  <- epimodel$likelihoods$pop_likelihood_cur +
    suff_stats[1] * (log(proposal["beta"]) - log(epimodel$params["beta"])) +
    suff_stats[3] * (log(proposal["mu"]) - log(epimodel$params["mu"])) -
    suff_stats[2] * (proposal["beta"] - epimodel$params["beta"]) - 
    suff_stats[4] * (proposal["mu"] - epimodel$params["mu"])
  
  # update parameters, likelihood objects, and eigen decompositions
  epimodel  <-
    update_params(
      epimodel,
      params = proposal,
      pop_likelihood = pop_likelihood_new,
      obs_likelihood = obs_likelihood_new
    )
  
  return(epimodel)
}

# Fit the stochastic epidemic model ------------------------------
# Grab the data that was simulated previously. No need to redefine the measurement process functions, they remain unchanged.
dat <- simulated_epimodel$dat

# Initial values for initial state parameters
init_dist <- MCMCpack::rdirichlet(1, c(9,0.5,0.1))
epimodel <- init_epimodel(popsize = popsize,                                                   # population size
                          states = c("S", "I", "R"),                                           # compartment names
                          params = c(beta = abs(rnorm(1, beta, beta*0.1)),                     # per-contact infectivity rate
                                     mu = abs(rnorm(1, mu, mu*0.002)),                         # recovery rate
                                     rho = rbeta(1, 21, 75),                                   # binomial sampling probability
                                     S0 = init_dist[1], I0 = init_dist[2], R0 = init_dist[3]), # initial state probabilities
                          rates = c("beta * I", "mu"),                                         # unlumped transition rates
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T),           # flow matrix
                          dat = dat,                                                           # dataset
                          time_var = "time",                                                   # name of time variable in the dataset
                          meas_vars = "I",                                                     # name of measurement variable
                          initdist_prior = c(90, 2, 5), ### Parameters for the dirichlet prior distribution for the initial state probs
                          r_meas_process = r_meas_process,
                          d_meas_process = d_meas_process)

# Set up MCMC settings
epimodel <- init_settings(epimodel,
                          niter = niter,  # this was set to 100,000 in the paper
                          save_params_every = 1, 
                          save_configs_every = niter/10, # this was set to 250 in the paper 
                          kernel = list(gibbs_kernel_SIR),
                          configs_to_redraw = 75, # this was set to 75 in the paper
                          analytic_eigen = "SIR", # compute eigen decompositions and matrix inverses analytically
                          ecctmc_method = "unif")   # sample subject paths in interevent intervals via modified rejection sampling
tic <- proc.time()
fitted_epimodel <- fit_epimodel(epimodel, monitor = verbose)
toc <- proc.time()

# Save files ------------------------------
file_suffix <- paste0("beta=", beta, "_", 
                      "mu=", mu, "_", 
                      "rho=", rho, "_", 
                      "popsize=", popsize, "_", 
                      "niter=", niter, "_", 
                      "obsint=", obsint)

file_dir <- paste0("results/", file_suffix, "/")
if (!dir.exists(file_dir)) {
  dir.create(file_dir)
}

saveRDS(simulated_epimodel, file = paste0(file_dir, "simulated_epimodel.rds"))
saveRDS(fitted_epimodel, file = paste0(file_dir, "fitted_epimodel.rds"))

profile_file <- paste0(file_dir, "profile.txt")
cat("MCMC runtime", "\n", file = profile_file)
cat("User:", toc[1] - tic[1], "\n", file = profile_file, append = TRUE)
cat("System:", toc[2] - tic[2], "\n", file = profile_file, append = TRUE)
cat("Elapsed:", toc[3] - tic[3], file = profile_file, append = TRUE)
