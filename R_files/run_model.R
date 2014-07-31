library(MCMCpack)
library(parallel)
library(rstan)

setwd("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/Stan/My Code/farhang_katznelson_rep")


# load prepared data
loadData.path <- "Data/prepared_data.RData"
load(loadData.path)


# to run in parallel
PARS <- c("b_0","b_URB","b_AA","b_UNION", "b_DEM", "b_LABORCOMM")
rng.seed <- 123456
stanFile.path <- "Stan_files/Attempt1_fix.stan"
fit <- stan(file = stanFile.path, data = data.list, pars = PARS, chains = 0)
sflist <- mclapply(1:4, mc.cores = 4, 
                   function(k) stan(fit = fit, 
                                    seed = rng.seed, 
                                    data = data.list,
                                    pars = PARS,
                                    chains = 1,
                                    iter = 1000,
                                    chain_id = k,
                                    refresh = -1))
fits_Attempt1_fix <- sflist2stanfit(sflist)
saveFits.path <- "Stan_fits/stanfit_Attempt1_fix.RData"
save(fits_Attempt1_fix, file = saveFits.path)





# to run the model one chain at a time
n.chains <- 4
n.iter <- 3000
fit1 <- stan(fit = fit, data = data.list, iter = n.iter, chains = n.chains, init="random", verbose = TRUE)

saveStan.path <- "Stan_fits/stanfit_no_tau.RData"
save(fit1, file = saveStan.path)



