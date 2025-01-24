# Inverse Leverage Effect and Direct Effect Analysis

This repository provides the code from our Bayesian data analysis used to assess differences in the leverage effect for stochastic volatility models with discontinuous jumps for assets such as cryptocurrencies and meme stocks when modelled univariately or bivariately with a market index, such as the S&P 500.  The results of this research can be found in our preprint "Inverse Leverage Effect for Cryptocurrencies and Meme Stocks: A Comprehensive Framework".  

The code is broken down into the following files:
- pgas_2d.cpp: The C++ code for simulating a random draw from the particular full conditional posterior distribution of interest either using a Gibbs sampler for a model parameter or using a Particle Gibbs with Ancestor Sampling algorithm for the latent time series variables such as volatility and jump sizes
- run_mcmc_2c.R: The R code for running our MCMC algorithm for a prespecified burnin period, number of desired samples, and thinning parameter (if specified) at a given set of starting values specified in starting_values_2d.R
- simulation_study_2d.R: The R code for conducting our thorough simulation study
- 
