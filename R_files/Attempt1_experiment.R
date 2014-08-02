# Load data ---------------------------------------------------------------
# _________________________________________________________________________

setwd("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/farhang_katznelson_rep")

loadData.path <- "Data/rollcall.logit.inter.all.names.txt"
Data <- read.table(loadData.path, header = TRUE)


# Construct map matrix ----------------------------------------------------
# _________________________________________________________________________

nr <- with(Data, length(unique(region)))
nc <- with(Data, length(unique(t)))

map.matrix <- mat.or.vec(nr,nc)
rownames(map.matrix) <- c("NS", "BS", "DS")
colnames(map.matrix) <- paste0("t",1:length(unique(Data$t)))

map.matrix[1, ] <- c(3, rep(5, nc - 2), 3)
map.matrix[2, ] <- c(5, rep(8, nc - 2), 5)
map.matrix[3, ] <- c(3, rep(5, nc - 2), 3)
map.matrix

# R = number of unique regions (i.e. number cells in map.matrix)
R <- prod(dim(map.matrix))

# N = total number of neighbors
N <- 3*sum(map.matrix == 3) + 5*sum(map.matrix == 5) + 8*sum(map.matrix == 8)

# just as a guide/reference: number the map cells from 1 to 24
map.index <- matrix(1:R, nr, nc)
colnames(map.index) <- colnames(map.matrix)
rownames(map.index) <- rownames(map.matrix)
map.index




# Fill in map vector to pass to Stan --------------------------------------
# _________________________________________________________________________

map <- c(2, 4, 5)               # neighbors for cell 1
map <- c(map, c(1, 3, 4, 5, 6)) # neighbors for cell 2
map <- c(map, c(2, 5, 6))       # neighbors for cell 3

# for filling in neighbors for cells 4 through 21
NS <- seq(4, 19, 3); top <- c(-2, -1, 2, 4, 5)
BS <- seq(5, 20, 3); mid <- c(-2:0, 1, 3, 4:6)
DS <- seq(6, 21, 3); bot <- c(-1, 0, 2, 5, 6)

i <- 4
while(i <= 21) {
  if (i %in% NS) map <- c(map, top + 3*which(NS == i)) 
  if (i %in% BS) map <- c(map, mid + 3*which(BS == i)) 
  if (i %in% DS) map <- c(map, bot + 3*which(DS == i))
  i <- i + 1
}

map <- c(map, c(19, 20, 23))          # neighbors for cell 22
map <- c(map, c(19, 20, 21, 22, 24))  # neighbors for cell 23
map <- c(map, c(20, 21, 23))          # neighbors for cell 24



# Prepare rest of data for Stan -------------------------------------------
# _________________________________________________________________________

# offsets vector
off <- c(0, cumsum(map.matrix))

# I = number of observations
I <- nrow(Data)

# Create region vector with possible values 1 to R (1 to 24 in this case)
# Note: region codes in data: DS = 1, BS = 2, NS, = 3

library(plyr)
Data$region.name <- mapvalues(Data$region, from = c(1:3), to = c("DS", "BS", "NS"))

region <- rep(NA, 24)
for(i in 1:8){
  region[with(Data, region.name == "NS" & period == i)] <- seq(1, 22, 3)[i]
  region[with(Data, region.name == "BS" & period == i)] <- seq(2, 23, 3)[i]
  region[with(Data, region.name == "DS" & period == i)] <- seq(3, 24, 3)[i]
}




# Create Nneighs vector containing number of neighbors for each region
Nneighs <- sapply(1:R, function(r) off[r+1] - off[r]) # or equivalently Nneighs <- c(map.matrix)




# Make data list to pass to Stan ------------------------------------------
# _________________________________________________________________________

data.list <- list(I = I, # number of obs
                  R = R, # number of regions
                  N = N, # length of map vector (total number of neighbors)
                  K = 4, # number of parameters to get instrinsic autoregressive priors
                  Nneighs = Nneighs,
                  map = map, 
                  off = off,
                  region = region,
                  LABORCOMM = Data$laborcomm,
                  DEM = Data$dem,
                  URB = Data$urbanpct,
                  AA = Data$aapct,
                  UNION = Data$unionpop,
                  VOTE = Data$prolabor)






# Run the model -----------------------------------------------------------
# _________________________________________________________________________

library(MCMCpack)
library(parallel)
library(rstan)

setwd("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/Github/farhang_katznelson_rep/Stan_files")

# compile C++ code
stanFile.path <- "Attempt1_experiment.stan"
PARS <- c("B", "b_DEM", "b_LABORCOMM")
fit.test <- stan(file = stanFile.path, data = data.list, chains = 1, iter = 50, pars = PARS)


# to run the model one chain at a time
n.chains <- 1
n.iter <- 500
fit1 <- stan(fit = fit, data = data.list, iter = n.iter, chains = n.chains, init="random", verbose = TRUE)


# to run in parallel
library(parallel)
rng.seed <- 1234
fit <- stan(file = stanFile.path, data = data.list, pars = PARS, chains = 0)
sflist <- mclapply(1:4, mc.cores = 4, 
                   function(k) stan(fit=fit, 
                                    seed = rng.seed, 
                                    data= data.list,
                                    pars=PARS,
                                    chains = 1,
                                    iter = 1000,
                                    chain_id = k,
                                    verbose = TRUE,
                                    refresh = -1))
fits <- sflist2stanfit(sflist)
saveFits.path <- "Stan_fits/stanfit_tau.RData"
save(fits, file = saveFits.path)


monitor(fits)
pars.list <- c("B", "b_DEM", "b_LABORCOMM")
plot(fits, pars = pars.list)
pairs(fits, pars = pars.list)

plotStan.path <- "Stan_fits/Plots/"
pdf(file=paste0(plotStan.path,"traceplots.pdf"))
traceplot(fits, pars = pars.list, inc_warmup=FALSE)
dev.off()

source("R_files/dependence.R")
svg(file=paste0(plotStan.path,"dependence.svg"))
invisible(dependence(fits))
dev.off()



# extracting the posterior samples
extractStan <- extract(fits, pars = pars.list)

b_DEM <- extractStan$b_DEM
b <- extractStan$B[,1,]
b_UNION <- extractStan$B[,2,]
b_AA <- extractStan$B[,3,]
b_URB <- extractStan$B[,4,]

NS <- seq(1,22, 3); BS <- seq(2,23,3); DS <- seq(3,24,3)

b_all <- cbind(b[,NS], b[,BS], b[,DS])
b_UNION_all <- cbind(b_UNION[,NS], b_UNION[,BS], b_UNION[,DS])
b_AA_all <- cbind(b_AA[,NS], b_AA[,BS], b_AA[,DS])
b_URB_all <- cbind(b_URB[,NS], b_URB[,BS], b_URB[,DS])


congnum <- seq(73,80,1)
lb <- .05
ub <- .95
regression.plots <- function (ests,varname,regionname) {
  lbub <- apply(ests,2,quantile,(probs=c(lb,ub)))
  b.range <- range(lbub[1,],lbub[2,])+range(lbub[1,],lbub[2,])*.01
  meds <- apply(ests,2,median)
  plot(congnum,meds,ylim=b.range,pch=20, cex.axis=.8,xlab="Congress",ylab="Estimates",main=list(varname,cex=.9),sub=regionname)
  segments(congnum,lbub[1,],congnum,lbub[2,])
  lines (c(72,81), c(0,0),lty=3, lwd = 1.5, col = "turquoise4")
  cat(regionname,varname,"\n")
  print(cbind(congnum,apply(ests,2,median),lbub[1,],lbub[2,]))
  
}


pdf(paste0(plots.path,"votes_regression_figs3.pdf"))
par(mfrow=c(4,3))

regression.plots(b_all[,17:24],"Region-Period Effect","Deep South")
regression.plots(b_all[,9:16],"Region-Period Effect","Border South")
regression.plots(b_all[,1:8],"Region-Period Effect","Non-South")

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
