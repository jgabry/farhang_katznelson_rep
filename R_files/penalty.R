### see equation 18 (pp. 14) http://www.stat.uni-muenchen.de/~bayesx/manual/methodology_manual.pdf

setwd("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/Rsync")

loadData.path <- "rollcall.logit.inter.all.names.txt"
Data <- read.table(loadData.path, header = TRUE)


# Construct map matrix ----------------------------------------------------
# _________________________________________________________________________

nr <- with(Data, length(unique(region)))
nc <- with(Data, length(unique(period)))

map.matrix <- mat.or.vec(nr,nc)
rownames(map.matrix) <- c("NS", "BS", "DS")
colnames(map.matrix) <- paste0("t", 1:nc)

# just as a guide/reference: number the map cells from 1 to 24
map.matrix[1, ] <- c(3, rep(5, nc - 2), 3)
map.matrix[2, ] <- c(5, rep(8, nc - 2), 5)
map.matrix[3, ] <- c(3, rep(5, nc - 2), 3)
map.matrix

# R = number of unique region/period combinations (i.e. number cells in map.matrix)
R <- prod(dim(map.matrix))

map.index <- matrix(1:R, nr, nc)
colnames(map.index) <- colnames(map.matrix)
rownames(map.index) <- rownames(map.matrix)
map.index




# Fill in map vector ------------------------------------------------------
# _________________________________________________________________________


map <- c(2, 4, 5)               # neighbors for cell 1
map <- c(map, c(1, 3, 4, 5, 6)) # add neighbors for cell 2
map <- c(map, c(2, 5, 6))       # add neighbors for cell 3

# for filling in neighbors for cells 4 through 21
NS <- seq(4, 19, 3); top <- c(-2, -1, 2, 4, 5)
BS <- seq(5, 20, 3); mid <- c(-2, -1, 0, 1, 3, 4, 5, 6)
DS <- seq(6, 21, 3); bot <- c(-1, 0, 2, 5, 6)

i <- 4
while(i <= 21) {
  if (i %in% NS) map <- c(map, top + 3*which(NS == i)) 
  if (i %in% BS) map <- c(map, mid + 3*which(BS == i)) 
  if (i %in% DS) map <- c(map, bot + 3*which(DS == i))
  i <- i + 1
}

map <- c(map, c(19, 20, 23))          # add neighbors for cell 22
map <- c(map, c(19, 20, 21, 22, 24))  # add neighbors for cell 23
map <- c(map, c(20, 21, 23))          # add neighbors for cell 24

map


# offsets vector
off <- c(0, cumsum(map.matrix))

# N = number of observations
N <- nrow(Data)


# Create region vector with possible values 1 to R (1 to 24 in this case). Regions are 
# assigned to correspond to the cell numbers in map.index

library(plyr) 
Data$region.name <- mapvalues(Data$region, 
                              from = 1:3, # Note: region codes in data: DS = 1, BS = 2, NS, = 3
                              to = c("DS", "BS", "NS"))

region <- rep(NA, 24)
for(i in 1:8){
  region[with(Data, region.name == "NS" & period == i)] <- seq(1, 22, 3)[i]
  region[with(Data, region.name == "BS" & period == i)] <- seq(2, 23, 3)[i]
  region[with(Data, region.name == "DS" & period == i)] <- seq(3, 24, 3)[i]
}



# design matrix
XX <- mat.or.vec(N,R)
for(n in 1:N){
  for(r in 1:R){
    XX[n,r] <- ifelse(region[n] == r, 1, 0)
  }
}


XX_urb <- XX*Data$urbanpct
XX_aa <- XX*Data$aapct
XX_union <- XX*Data$union 


# penalty/area-adjancency matrix
KK <- mat.or.vec(R, R)
for (i in 1:R) {
  for (j in 1:R) {
    KK[i, map[ (off[i] + 1) : off[i+1]] ] <- 1
  }
}

D <- matrix(NA, R, R)
for (i in 1:R) {
  for (j in 1:R) {
    D[i, j] <- ifelse(i == j, sum(KK[i, ]), 0.0)
  }
}



library(Matrix)
rankMatrix(KK)

data.list <- list(XX = XX, 
                  XX_urb = XX_urb, 
                  XX_aa = XX_aa, 
                  XX_union = XX_union,
                  KK = KK,
                  D = D,
                  R = R,
                  N = N,
                  DEM = Data$dem,
                  LAB = Data$laborcomm,
                  VOTE = Data$prolabor)

library(rstan)
init <- function(){
  list(beta = matrix(0.0, 4, R), 
       tau_sq = rep(1.0, 4), 
       delta_noise = rnorm(2), 
       rho = runif(4))
}



library(rstan)
fit.compile <- stan(file = "MRF_attempt.stan", data = data.list, chains = 0)
fit.test <- stan(fit = fit.compile, data = data.list, chains = 1, iter = 50)
fit.test.init <- stan(fit = fit.compile, data = data.list, chains = 1, iter = 50, init = init)

PARS <- c("beta", "delta", "tau_sq", "rho")
rng.seed <- 123456

library(parallel)
sflist_MRF <- mclapply(1:6, mc.cores = 6, 
                       function(k) stan(fit = fit.compile, 
                                        seed = rng.seed, 
                                        data = data.list,
                                        pars = PARS,
                                        chains = 1,
                                        iter = 2000,
                                        chain_id = k,
                                        verbose = TRUE,
                                        refresh = -1))
save(sflist_MRF, file = "sflist_MRF.RData")


stanfit_MRF <- sflist2stanfit(sflist_MRF)


#########






fit.compile <- stan(file = "MRF_attempt_rank.stan", data = data.list, chains = 0)
fit.test <- stan(fit = fit.compile, data = data.list, chains = 1, iter = 50, init = init)


KK2 <- KK
for(i in 1:24){
  KK2[i,i] <- 1
}

rank <- as.numeric(rankMatrix(KK2))

library(MASS)
M <- ginv(KK2)
M_block <- M[1:rank, 1:rank]

data.list <- list(XX = XX, 
                  XX_urb = XX_urb, 
                  XX_aa = XX_aa, 
                  XX_union = XX_union,
                  M = M_block,
                  rank = rank,
                  R = R,
                  N = I,
                  DEM = Data$dem,
                  LAB = Data$laborcomm,
                  VOTE = Data$prolabor)

