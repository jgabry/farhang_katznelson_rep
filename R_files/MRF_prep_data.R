### see equation 18 (pp. 14) http://www.stat.uni-muenchen.de/~bayesx/manual/methodology_manual.pdf

setwd("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/Rsync/The_Real_Deal")

loadData.path <- "rollcall.logit.inter.all.names.txt"
Data <- read.table(loadData.path, header = TRUE)


#### Construct map matrix ####
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

map.index <- matrix(1:nr*nc, nr, nc)
colnames(map.index) <- colnames(map.matrix)
rownames(map.index) <- rownames(map.matrix)
map.index


#### Fill in map vector ####
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


#### offsets vector ####
off <- c(0, cumsum(map.matrix))

#### constants to pass to stan ####
N <- nrow(Data) # number of observations
R <- prod(dim(map.matrix)) # number of unique 'regions' (i.e. region/period combinations)


#### Create region vector #### 
  # Regions (from 1 to R) are assigned to correspond to the cell numbers in map.index

library(plyr) 
Data$region.name <- mapvalues(Data$region, 
                                # Note: region codes in data: DS = 1, BS = 2, NS, = 3
                              from = 1:3, 
                              to = c("DS", "BS", "NS"))

region <- rep(NA, R)
for(i in 1:8){
  region[with(Data, region.name == "NS" & period == i)] <- seq(1, 22, 3)[i]
  region[with(Data, region.name == "BS" & period == i)] <- seq(2, 23, 3)[i]
  region[with(Data, region.name == "DS" & period == i)] <- seq(3, 24, 3)[i]
}



#### incidence and design matrices ####
M <- mat.or.vec(N,R)
for(n in 1:N){
  for(r in 1:R){
    M[n,r] <- ifelse(region[n] == r, 1, 0)
  }
}

X_urb <- M*Data$urbanpct
X_aa <- M*Data$aapct
X_union <- M*Data$unionpop



X <- array(NA, dim = c(3,N,R))
X[1,,] <- X_urb
X[2,,] <- X_aa
X[3,,] <- X_union


#### area-adjancency matrix ####
A <- mat.or.vec(R, R)
for (i in 1:R) {
  for (j in 1:R) {
    A[i, map[ (off[i] + 1) : off[i+1]] ] <- 1
  }
}

#### degree matrix ####
D <- mat.or.vec(R, R)
for (i in 1:R) {
  for (j in 1:R) {
    D[i, j] <- ifelse(i == j, sum(A[i, ]), 0.0)
  }
}

map.matrix


#### other covariates ####
Z <- with(Data, cbind(dem, laborcomm))

#### outcome variable ####
Y <- Data$prolabor


# STAN --------------------------------------------------------------------
# _________________________________________________________________________



data.list <- list(Z = Z,
                  X = X,
                  Y = Y,
                  M = M,
                  A = A, 
#                   D = D,
                  R = R,
                  N = N)

PARS <- c("beta", "delta","tau", "rho")

library(rstan)
fit.compile <- stan(file = "MRF_new.stan", 
                    data = data.list, pars = PARS, chains = 0)
fit.test <- stan(fit = fit.compile, 
                 data = data.list, pars = PARS, chains = 1, iter = 30)
# fit.test.init <- stan(fit = fit.compile, data = data.list, chains = 1, iter = 50, init = init)



library(parallel)
rng.seed <- 123456
sflist_MRF_new <- mclapply(X = 1:6, 
                           mc.cores = 6, 
                           FUN = function(k) {
                             stan(fit = fit.compile, 
                                  seed = rng.seed, 
                                  data = data.list,
                                  pars = PARS,
                                  chains = 1,
                                  iter = 2000,
                                  chain_id = k,
                                  verbose = TRUE,
                                  refresh = -1) })


save(sflist_MRF_new, file = "sflist_MRF_new.RData")



load("sflist_MRF_new.RData")
stanfit_MRF_new <- sflist2stanfit(sflist_MRF_new)
monitor(stanfit_MRF_new)