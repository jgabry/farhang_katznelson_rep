library(MCMCpack)
library(rstan)

# check out library(ggmcmc)


setwd("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/Rsync/The_Real_Deal")


# functions to set everything up ------------------------------------------
# _________________________________________________________________________

# function to make new subdirectory in plots folder for the current date
make_todays_plots_dir <- function(plot.dir, setwd = TRUE) {
  date <- Sys.Date()
  date <- gsub("-","_", date)
  new.dir <- paste(plot.dir, date, sep = "/")
  dir.create(new.dir)
  dir.msg <- paste0("New directory created: ", getwd(),"/",new.dir)
  if (setwd == TRUE) {
    setwd(new.dir)
    print(dir.msg)
    print("Has been set as working directory")
  }
  else {
    print(dir.msg)
  }
}

# function to load stan object, rename to stanfit (converting from sflist to stanfit if needed)
prepare_stan <- function(stan.object.name, path, convert = TRUE, monitor = FALSE) {
  # sflist.name = name of sflist object
  # rel.path = relative path to sflist.name if not in wd
  
  require(rstan)
  
  if (missing(path)) {
    to_load <- paste0(sflist.name, ".RData")
  }
  else {
    to_load <- paste0(path,"/",sflist.name, ".RData")
  }
  do.call("load", list(file = to_load))
  stanfit <- get(stan.object.name)
  
  if (convert == TRUE) {
    stanfit <- sflist2stanfit(stanfit)
    convert.msg <- paste(stan.object.name, "converted to class stanfit")
    print(convert.msg)
  }
  
  if (monitor == TRUE) {
    monitor(stanfit)
  }
  
  return(stanfit)
}


# load stan object, and create today's plot directory
sflist.name <- "mvn_prec_fit_spatial_Aug7"
stanfit <- prepare_stan(sflist.name, convert = TRUE)

plot.dir <- "plots"
make_todays_plots_dir(plot.dir, setwd = TRUE)





# make trace plots --------------------------------------------------------
# _________________________________________________________________________
make_trace_plots <- function(dir = "wd"){
  require(rstan)
  
  if (dir != "wd") {
    wd <- getwd()
    setwd(dir)
  }
  
  print("Making trace plots. Please wait ...")
  
  # beta_urb traceplots
  png(file = "traceplots_URB.png", 
      w = 8, h = 10.5, units = "in", res = 300)
  traceplot(stanfit, pars = paste0("beta[1,",1:24,"]"), nr = 4, nc = 6)
  dev.off()
  
  # beta_aa traceplots
  png(file = "traceplots_AA.png", 
      w = 8, h = 10.5, units = "in", res = 300)
  traceplot(stanfit, pars = paste0("beta[2,",1:24,"]"), nr = 4, nc = 6)
  dev.off()
  
  # beta_union traceplots
  png(file = "traceplots_UNION.png", 
      w = 8, h = 10.5, units = "in", res = 300)
  traceplot(stanfit, pars = paste0("beta[3,",1:24,"]"), nr = 4, nc = 6)
  dev.off()
  
  # beta_0 traceplots
  png(file = "traceplots_b0.png", 
      w = 8, h = 10.5, units = "in", res = 300)
  traceplot(stanfit, pars = paste0("beta[4,",1:24,"]"), nr = 4, nc = 6)
  dev.off()
  
  # delta & delta_noise traceplots
  png(file = "traceplots_delta.png", 
      w = 8, h = 10.5, units = "in", res = 300)
  traceplot(stanfit, pars = c("delta","delta_noise"), nr=4,nc=1)
  dev.off()
  
  # tau & tau_mean traceplots
  png(file = "traceplots_tau.png", 
      w = 8, h = 10.5, units = "in", res = 300)
  traceplot(stanfit, pars = c("tau","tau_mean"), nr=5,nc=1)
  dev.off()
  
  # tau_unit traceplots
  png(file = "traceplots_tau_unit.png", 
      w = 8, h = 10.5, units = "in", res = 300)
  traceplot(stanfit, pars = c("tau_unit"), nr=4,nc=1)
  dev.off()
  
  # rho & rho_mean traceplots
  png(file = "traceplots_rho.png", 
      w = 8, h = 10.5, units = "in", res = 300)
  traceplot(stanfit, pars = c("rho", "rho_mean"), nr=5,nc=1)
  dev.off()
  
  # alpha0 & lp__ traceplot
#   png(file = "traceplots_alpha0_and_lp.png", 
#       w = 8, h = 10.5, units = "in", res = 300)
#   traceplot(stanfit, pars = c("alpha0", "lp__"), nr=2)
#   dev.off()

  # lp__ traceplot
  png(file = "traceplots_lp.png", 
      w = 8, h = 10.5, units = "in", res = 300)
  traceplot(stanfit, pars = "lp__")
  dev.off()
  
  if (dir != "wd") {
    setwd(wd)
  }
}

make_trace_plots()

# summary statistics
monitor(stanfit, digits = 4)
library(xtable)
print(xtable(summary(stanfit)$summary)), booktabs = TRUE), 
#       sanitize.text.funtion = function(x){x}, 
      tabular.environment = 'longtable', 
      booktabs = TRUE,
      format.args = list(digits = 3))



# Ben's dependence plot 
source("../dependence_new.R")
svg(filename = paste0("dependence.svg"))
invisible(dependence(stanfit))
dev.off()


# Ben's convergence test
source("../convergence_new.R")
convergence_test(stanfit, 
                 thin = stop("check the dependence plot"),
                 R = stop("see ?disco"))




# extract posterior samples -----------------------------------------------
# _________________________________________________________________________
extractStan <- extract(stanfit)
b_DEM <- extractStan$delta[,1]
b_LABORCOMM <- extractStan$delta[,2]
b_URB <- extractStan$beta[,1,]
b_AA <- extractStan$beta[,2,]
b_UNION <- extractStan$beta[,3,]
b_0 <- extractStan$beta[,4,]


# define regression.plots function ----------------------------------------
# _________________________________________________________________________
regression.plots <- function (ests, varname, regionname, 
                              congnum = seq(73,80,1), conf.level = 0.90) {
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


# make regression figs plots ----------------------------------------------
# _________________________________________________________________________
NS <- seq(1, 22, 3) 
BS <- seq(2, 23, 3) 
DS <- seq(3, 24, 3)
r.names <- c("Deep South", "Border South", "Non-South")
names(r.names) <- c("DS", "BS", "NS")
v.names <- c("Region-Period Effect", "Union Pct", "African-American Pct", "Urban Pct")
names(v.names) <- c("b_0", "b_UNION", "b_AA", "b_URB")


pdf(file = "votes_regression_figs.pdf")
par(mfcol=c(4,3))
for (X in names(r.names)) {
  for (Y in names(v.names)) {
    args <- list(ests = get(Y)[,get(X)], varname = v.names[Y], regionname = r.names[X])
    do.call(regression.plots, args)
  }
}
dev.off()
par(mfcol = c(1,1))




# plot dem and laborcomm posterior densities
library(ggplot2)
library(reshape2)

pdf(file = "gg_DEM_LABORCOMM.pdf",w = 10.5, h = 7.5, paper = "USr")

molten.df <- melt(data.frame(DemocraticParty = b_DEM, LaborCommittee = b_LABORCOMM))
(ggplot(molten.df, aes(x = value, group = variable, color = variable, fill = variable)) 
 + geom_density() 
 + xlab("") + ggtitle("Poseterior Kernel Density Estimates"))
#  + theme(legend.position = "bottom"))
dev.off()



# make gg density plots for betas

source("../make_gg_densities.R")

pdf(file = "gg_beta_densities.pdf",
    w = 10.5, h = 7.5, paper = "USr")
make_gg_densities(ests = b_0, plot_title = "Region-Period Effect")
make_gg_densities(b_URB, "Urban Pct")
make_gg_densities(b_AA, "African-American Pct")
make_gg_densities(b_UNION, "Union Pct")
dev.off()


# make gg line and ribbon plots for betas

source("../make_gg_lines.R")

Bnames <- c("b_0","b_URB", "b_AA", "b_UNION")
titles <- c("Region-period","Urban", "African-American", "Union")
my.pal <- c("purple4", "#D55E00", "#56B4E9")

gg_lines <- make_gg_lines(Bnames = Bnames, titles = titles, palette = my.pal)

pdf(file = "gg_beta_lines.pdf", w = 10.5, h = 7.5, paper = "USr")
gg_all <- grid.arrange(arrangeGrob(  
  gg_lines[[1]], gg_lines[[2]], gg_lines[[3]], gg_lines[[4]],
  nrow = 2, main = "Posterior Medians & 90% Credible Intervals"))
dev.off()


