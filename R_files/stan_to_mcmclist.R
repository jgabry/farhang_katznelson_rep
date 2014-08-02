library(MCMCpack)
library(rstan)
library(xtable)


# load the stan output ----------------------------------------------------
setwd("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/Rsync")

# load("sflist1.RData")
# load("sflist2.RData")
# stanfit1 <- sflist2stanfit(sflist1)
# stanfit2 <- sflist2stanfit(sflist2)
# stanfit <- sflist2stanfit(list(stanfit1, stanfit2))
# save(stanfit, file = "all_chains_stanfit_local.RData")

# load("all_chains_stanfit_local.RData")

MCMCstanfit <- mcmc.list(lapply(1:ncol(stanfit), function(x) mcmc(as.array(stanfit)[,x,])))



# plot autocorrelations for each chain
pdf(file = "autocorr_plots.pdf")
for (i in 1:8){
  autocorr.plot(MCMCstanfit[i], ask = FALSE, sub = list(paste0("Chain_",i),col = "blue") )
}
dev.off()


# effective sample sizes
esample <- effectiveSize(MCMCstanfit)
esample.tab <- matrix(esample[1:96],24,4)
esample.tab <- as.data.frame(esample.tab)
colnames(esample.tab) <- c("b0", "URB", "AA", "UN")

NS <- seq(1, 22, 3) 
BS <- seq(2, 23, 3) 
DS <- seq(3, 24, 3)
rownames(esample.tab)[as.numeric(rownames(esample.tab)) %in% NS] <- paste0("NS",1:8)
rownames(esample.tab)[as.numeric(rownames(esample.tab)) %in% BS] <- paste0("BS",1:8)
rownames(esample.tab)[as.numeric(rownames(esample.tab)) %in% DS] <- paste0("DS",1:8)

print(xtable(esample.tab, digits = 0, include.rownames = T), 
      sanitize.text.function = function(x){x}, booktabs = T)

  

# gelman-rubin Rhat
rhat <- gelman.diag(MCMCstanfit)
rhat.tab <- as.data.frame(rhat[1])
colnames(rhat.tab) <- c("Point.Est", "Upper95")
rhat.tab_b0 <- rhat.tab[1:24,]
rhat.tab_URB <- rhat.tab[25:48,]
rhat.tab_AA <- rhat.tab[49:72,]
rhat.tab_UN <- rhat.tab[73:96,]

rownames(rhat.tab_b0)[NS] <- paste0("NS",1:8)
rownames(rhat.tab_b0)[BS] <- paste0("BS",1:8)
rownames(rhat.tab_b0)[DS] <- paste0("DS",1:8)
rownames(rhat.tab_URB)[NS] <- paste0("NS",1:8)
rownames(rhat.tab_URB)[BS] <- paste0("BS",1:8)
rownames(rhat.tab_URB)[DS] <- paste0("DS",1:8)
rownames(rhat.tab_AA)[NS] <- paste0("NS",1:8)
rownames(rhat.tab_AA)[BS] <- paste0("BS",1:8)
rownames(rhat.tab_AA)[DS] <- paste0("DS",1:8)
rownames(rhat.tab_UN)[NS] <- paste0("NS",1:8)
rownames(rhat.tab_UN)[BS] <- paste0("BS",1:8)
rownames(rhat.tab_UN)[DS] <- paste0("DS",1:8)

print(xtable(rhat.tab_b0, digits = 4, include.rownames = T), 
      sanitize.text.function = function(x){x}, booktabs = T)

print(xtable(rhat.tab_URB, digits = 4, include.rownames = T), 
      sanitize.text.function = function(x){x}, booktabs = T)

print(xtable(rhat.tab_AA, digits = 4, include.rownames = T), 
      sanitize.text.function = function(x){x}, booktabs = T)

print(xtable(rhat.tab_UN, digits = 4, include.rownames = T), 
      sanitize.text.function = function(x){x}, booktabs = T)




library(mcmcplots)
png(file = "b_0_traceplots.png", w = 8, h = 10.5, units = "in", res = 300)
traplot(MCMCstanfit, parms = "b_0")
dev.off()

png(file = "b_URB_traceplots.png", w = 8, h = 10.5, units = "in", res = 300)
traplot(MCMCstanfit, parms = "b_URB")
dev.off()

png(file = "b_AA_traceplots.png", w = 8, h = 10.5, units = "in", res = 300)
traplot(MCMCstanfit, parms = "b_AA")
dev.off()

png(file = "b_UNION_traceplots.png", w = 8, h = 10.5, units = "in", res = 300)
traplot(MCMCstanfit, parms = "b_UNION")
dev.off()