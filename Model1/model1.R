setwd("C:/Users/beacr/Politecnico di Milano/Simone Panzeri - Bayesian Project/model1")

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

load ('mean_intensity.dat')
load ('N_on.dat')
load ('N_pass.dat')
load ('T_raff.dat')
load ('matrice_B.dat')
Bb <- B
load ('sigma0.dat')
load ('geom_cut.dat')
load ('g.dat')
load ('N_intervals.dat')
load ('mub.dat')
mug <- mub
load ('B_matrix.dat')
Bg <- B
load ('mu.dat')


for(p in 1:420)
{
  for (f in (N_intervals[p]+1):12)
    if (f < 13)
  {
    N_on[f,p] <- 0
    N_pass[f,p] <- 0
    T_raff[f+1,p] <- 0
    mean_intensity[f+1,p] <- 0
  }
}

D<-diag(12)
Uno <- rep(1,10)

regression_data <- list(K = max(N_intervals),
                        N = 420,
                        N_int = N_intervals,
                        x1 = t(N_on),
                        x2 = t(N_pass),
                        x3 = t(T_raff[1:12,]),
                        Y = t(mean_intensity[2:13,]),
                        Sigma = s2_err,
                        D = D,
                        G = cbind(rep(1,420),geom,g),
                        mug = mug,
                        Bg = Bg,
                        mub = mub,
                        Bb = Bb)


model <- stan("model1.stan",
              data = regression_data,
              chains = 1,
              iter = 1000,
              algorithm = 'NUTS',
              verbose = TRUE,
              seed = 42)

x11()
rstan::traceplot(model, pars = "beta")

coda_chain <- As.mcmc.list(model, pars = "beta")
summary(coda_chain)

x11()
acfplot(coda_chain, lag.max = 30)

x11()
plot_post <- model%>%
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