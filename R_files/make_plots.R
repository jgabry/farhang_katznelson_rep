library(MCMCpack)
library(rstan)


# load the stan output ----------------------------------------------------
setwd("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/Rsync")

# load("sflist1.RData")
# load("sflist2.RData")
# stanfit1 <- sflist2stanfit(sflist1)
# stanfit2 <- sflist2stanfit(sflist2)
# stanfit <- sflist2stanfit(list(stanfit1, stanfit2))
# save(stanfit, file = "all_chains_stanfit_local.RData")

load("all_chains_stanfit_local.RData")
print(stanfit)
monitor(stanfit)
# convergence diagnostics -------------------------------------------------
PARS <- c("b_0","b_URB","b_AA","b_UNION", "b_DEM", "b_LABORCOMM")

# summary statistics
monitor(stanfit)
library(xtable)
print(xtable(monitor(stanfit), digits = 1), 
      sanitize.text.funtion = function(x){x}, 
      tabular.environment = 'longtable', 
      booktabs = TRUE,
      format.args = list(digits = 3))



# traceplots
traceplot(stanfit, pars = PARS)


# dependence plot 
source("dependence.R")

svg(filename = "plots/dependence_plot.svg")
invisible(dependence(stanfit))
dev.off()

# extract the posterior samples -------------------------------------------
extractStan <- extract(stanfit, pars = PARS)
b_DEM <- extractStan$b_DEM
b_LABORCOMM <- extractStan$b_LABORCOMM
b_0 <- extractStan$b_0
b_UNION <- extractStan$b_UNION
b_AA <- extractStan$b_AA
b_URB <- extractStan$b_URB


# define regression.plots function ----------------------------------------
regression.plots <- function (ests, varname, regionname, congnum = seq(73,80,1), conf.level = 0.90) {
  lb <- (1 - conf.level)/2
  ub <- 1 - lb
  lbub <- apply(ests, 2, quantile, (probs = c(lb, ub)))
  b.range <- range(lbub[1, ], lbub[2, ]) + range(lbub[1, ], lbub[2, ])*.01
  
  medians <- apply(ests, 2, median)
  plot(congnum, medians,
       ylim = b.range,
       pch = 20,
       cex.axis = 0.8,
       xlab = "Congress",
       ylab = "Estimates",
       main = list(varname, cex = 0.9),
       sub = regionname)
  
  segments(congnum, lbub[1, ], congnum, lbub[2, ])
  
  ll <- length(congnum)
  lines (c(congnum[1] - 1, congnum[ll] + 1), c(0, 0), lty = 3, lwd = 1.5, col = "turquoise4")
  
  cat(regionname, varname, "\n")
  print(cbind(congnum, medians, lbub[1, ], lbub[2, ]))
}



# make the plots ----------------------------------------------------------
NS <- seq(1, 22, 3) 
BS <- seq(2, 23, 3) 
DS <- seq(3, 24, 3)
r.names <- c("Deep South", "Border South", "Non-South")
names(r.names) <- c("DS", "BS", "NS")
v.names <- c("Region-Period Effect", "Union Pct", "African-American Pct", "Urban Pct")
names(v.names) <- c("b_0", "b_UNION", "b_AA", "b_URB")


pdf(file = "plots/votes_regression_figs_Stan_probit.pdf")
par(mfcol=c(4,3))
for (X in names(r.names)) {
  for (Y in names(v.names)) {
    args <- list(ests = get(Y)[,get(X)], varname = v.names[Y], regionname = r.names[X])
    do.call(regression.plots, args)
  }
}
dev.off()
par(mfcol = c(1,1))






library(MASS)

dens_dem <- density(b_DEM)
dens_dem$x[which.max(dens_dem$y)]

dens_lab <- density(b_LABORCOMM)
dens_lab$x[which.max(dens_lab$y)]

leg.text_dem <- c(paste("Mean:", round(mean(b_DEM), 3)),
                  paste("Median:", round(median(b_DEM), 3)),
                  paste("MAP:", round(dens_dem$x[which.max(dens_dem$y)], 3)))


leg.text_lab <- c(paste("Mean:", round(mean(b_LABORCOMM), 3)),
                  paste("Median:", round(median(b_LABORCOMM), 3)),
                  paste("MAP:", round(dens_lab$x[which.max(dens_lab$y)], 3)))



pdf(file = "plots/b_DEM_b_LABORCOMM_stan_new.pdf")
par(mfrow = c(1,2))
truehist(b_DEM, col = "skyblue", border = "skyblue", yaxt = "n")
# abline(v = dens_dem$x[which.max(dens_dem$y)])
legend("topright", leg.text_dem, bty = "n",title = "Estimates:", cex = 0.7)


truehist(b_LABORCOMM, col = "seagreen", border = "seagreen", yaxt = "n")
legend("topleft", leg.text_lab, bty = "n", title = "Estimates:", cex = 0.7)
# abline(v = dens_lab$x[which.max(dens_lab$y)])
dev.off()
par(mfrow = c(1,1))

