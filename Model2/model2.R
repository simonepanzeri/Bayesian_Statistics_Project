setwd("C:/Users/beacr/Politecnico di Milano/Simone Panzeri - Bayesian Project/model2")

library(rstan)
library(coda)
library(rjags)

library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)
require(gplots)
require(ggpubr)

load ('N_intervals.dat')
load ('mean_intensity.dat')
load ('geom_cut.dat')
load ('sigma2.dat')
load ('mub.dat')
load ('B_matrix.dat')
load ('g.dat')

#### AUTOREGRESSIVE MODEL ####

# p <- 170  # median value of geom and g
# p <- 147 # min value of geom and max value of g 
# p <- 188 # max value of geom
# p <- 406 # min value of g

  regression_data <- list(Y = mean_intensity[1:N_intervals[p],p],
                           G = c(geom[p],g[p]),
                           sigma = S[p],
                           mub = mub,
                           B = B,
                           N = N_intervals[p])
  
  fit <- stan(file='model2.stan', 
               data = regression_data,
               chains = 2,
               iter = 10000,
               warmup = 100, 
               thin = 10,
               algorithm = 'NUTS',
               verbose = TRUE,
               seed = 42)
  
  x11()
  rstan::traceplot(fit, pars = "beta", inc_warmup = TRUE)

  coda_chain <- As.mcmc.list(fit, pars = "beta")
  summary(coda_chain)

  gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
  geweke.diag(coda_chain, frac1=0.1, frac2=0.5)

  x11()
  acfplot(coda_chain, lag.max = 30)

  x11()
  plot_post <- fit%>%
    rstan::extract("beta") %>%
    as.data.frame() %>%
    map_df(as_data_frame, .id = 'param')
  plot_post %>%
    ggplot(aes(value, fill = param)) +
    geom_density() +
    facet_wrap(~param, scales = 'free') +
    scale_fill_locuszoom() +
    theme_minimal() +
    theme(legend.position="none")

  b = fit%>%
    rstan::extract("beta")
  c = exp(b[[1]][,1]*geom[p]+b[[1]][,2]*g[p])/(1+exp(b[[1]][,1]*geom[p]+b[[1]][,2]*g[p]));
  a = (exp(2*c)-1)/(exp(2*c)+1);

  x11()
  plot_post <- a %>%
    as.data.frame() %>%
    map_df(as_data_frame, .id = 'param')
  plot_post %>%
    ggplot(aes(value, fill = param)) +
    geom_density() +
    facet_wrap(~param, scales = 'free') +
    scale_fill_locuszoom() +
    theme_minimal() +
    theme(legend.position="none")