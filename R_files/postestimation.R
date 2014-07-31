setwd("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/Stan/My Code/farhang_katznelson_rep")

fit.path <- "stanfit_Attempt1_fix.RData"
fit.path <- paste0("Stan_fits/",fit.path)

fit.name <- "fits_yeti"

plotStan.path <- "Stan_fits/Plots/"

dependence.plot.name <- paste0("dependence_",fit.name,".svg")
estimates.plot.name <- paste0("votes_regression_figs_",fit.name,".pdf") 
trace.plot.name <- paste0("traceplots_",fit.name,".pdf") 

pars.list <- c("b_0","b_URB","b_AA","b_UNION", "b_DEM", "b_LABORCOMM")




# Summary of MCMC simulation samples  -------------------------------------

monitor(get(fit.name), pars = pars.list)
# plot(get(fit.name), pars = pars.list)
# pairs(get(fit.name), pars = pars.list)



# Make plots --------------------------------------------------------------

load(fit.path)


# dependence plot
source("R_files/dependence.R")
svg(file=paste0(plotStan.path,dependence.plot.name))
invisible(dependence(get(fit.name)))
dev.off()




# traceplots
pdf(file=paste0(plotStan.path, trace.plot.name))
traceplot(get(fit.name), pars = pars.list, inc_warmup=FALSE)
dev.off()



# extract posterior samples
extractStan <- extract(get(fit.name), pars = pars.list)

b_DEM <- extractStan$b_DEM
b_LABORCOMM <- extractStan$b_LABORCOMM
b_0 <- extractStan$b_0
b_UNION <- extractStan$b_UNION
b_AA <- extractStan$b_AA
b_URB <- extractStan$b_URB

NS <- seq(1,22, 3); BS <- seq(2,23,3); DS <- seq(3,24,3)

b_0_all <- cbind(b_0[,NS], b_0[,BS], b_0[,DS])
b_UNION_all <- cbind(b_UNION[,NS], b_UNION[,BS], b_UNION[,DS])
b_AA_all <- cbind(b_AA[,NS], b_AA[,BS], b_AA[,DS])
b_URB_all <- cbind(b_URB[,NS], b_URB[,BS], b_URB[,DS])


congnum <- seq(73,80,1)
lb <- .05
ub <- .95


source("R_files/regression_plots_function.R")
pdf(file = paste0(plotStan.path, estimates.plot.name))
par(mfrow=c(4,3))

regression.plots(b_0_all[,17:24],"Region-Period Effect","Deep South")
regression.plots(b_0_all[,9:16],"Region-Period Effect","Border South")
regression.plots(b_0_all[,1:8],"Region-Period Effect","Non-South")

regression.plots(b_UNION_all[,17:24],"Union Pct","Deep South")
regression.plots(b_UNION_all[,9:16],"Union Pct","Border South")
regression.plots(b_UNION_all[,1:8],"Union Pct","Non-South")

regression.plots(b_AA_all[,17:24],"African-American Pct","Deep South")
regression.plots(b_AA_all[,9:16],"African-American Pct","Border South")
regression.plots(b_AA_all[,1:8],"African-American Pct","Non-South")

regression.plots(b_URB_all[,17:24],"Urban Pct","Deep South")
regression.plots(b_URB_all[,9:16],"Urban Pct","Border South")
regression.plots(b_URB_all[,1:8],"Urban Pct","Non-South")

dev.off()
par(mfrow = c(1,1))





library(RColorBrewer)
colors <- c(brewer.pal(8, "Accent"), brewer.pal(8, "Pastel1"), brewer.pal(8, "Dark2"))
par(mfrow = c(6,4))
for(i in 1:24){
  hist(b[, i], col = colors[i],
       xlab = "", ylab = "", yaxt = "n", main = bquote(alpha[.(i)]))
}
par(mfrow = c(1,1))
