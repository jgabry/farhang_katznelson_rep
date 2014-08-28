## CHANGES: ###

# made urbanpct, aapct AND unionpop into proportions 

# setwd("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/Rsync/The_Real_Deal")

rm(list = ls())
setwd("/home/jsg2201/wawro/The_Real_Deal/Aug27")

# setwd("/home/jsg2201/wawro/The_Real_Deal/Aug25")
# rm(list = ls())

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

map.index <- matrix(1:(nr*nc), nr, nc)
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
NS <- map.index[1,]
BS <- map.index[2,]
DS <- map.index[3,]

for(i in 1:8){
  region[with(Data, region.name == "NS" & period == i)] <- NS[i]
  region[with(Data, region.name == "BS" & period == i)] <- BS[i]
  region[with(Data, region.name == "DS" & period == i)] <- DS[i]
}


#### incidence and design matrices ####
M <- matrix(NA, nrow = N, ncol = R)
for(n in 1:N){
  for(r in 1:R){
    M[n,r] <- ifelse(region[n] == r, 1, 0)
  }
}

urbanprop <- Data$urbanpct / 100 # make into proportion
unionprop <- Data$unionpop / 100 # make into proportion
X_urb <- M*urbanprop
X_union <- M*unionprop

X <- array(NA, dim = c(3, N, R))
X[1,,] <- X_urb
X[2,,] <- X_union
X[3,,] <- M



region.names <- c("NS", "BS", "DS")
M_aa <- matrix(NA, nrow = N, ncol = length(region.names))
for(n in 1:N){
  for(r in 1:length(region.names)){
    M_aa[n,r] <- ifelse(Data$region.name[n] == region.names[r], 1, 0)
  }
}

aaprop <- Data$aapct / 100 # make into proportion
X_aa <- M_aa*aaprop


#### area-adjancency matrix ####
A <- mat.or.vec(R, R)
for (i in 1:R) {
  for (j in 1:R) {
    A[i, map[ (off[i] + 1) : off[i+1]] ] <- 1
  }
}

#### degree matrix ####
D <- diag(rowSums(A))


#### other covariates ####
Z <- cbind(Data$dem, Data$laborcomm)


#### outcome variable ####
Y <- Data$prolabor


### sparse matrix stuff ###
# A_sparse <- which(A == 1, arr.ind=TRUE)
# # remove duplicates (because matrix is symmetric)
# A_sparse <- A_sparse[A_sparse[,1] < A_sparse[,2],]
# A_N <- dim(A_sparse)[1]
# A1 <- A_sparse[,1]
# A2 <- A_sparse[,2]
# 

### data list for stan ###
mvn_prec_data <- list(Y = Y,
                      X = X,
                      X_aa = X_aa,
                      Z = Z,
                      A = A, 
                      d = diag(D),
                      R = R,
                      N = nrow(M),
                      C = ncol(Z),
                      J = dim(X)[1],
                      P = ncol(M_aa))


library(rstan)
set_cppo("fast")


fit.compile <- stan(file = "the_slow_way_08_27.stan", data = mvn_prec_data, chains = 0)

fit.test <- stan(fit = fit.compile, data = mvn_prec_data, chains = 1, iter = 30, refresh = 5)

chain_1 <- stan(fit = fit.compile, data = mvn_prec_data, chains = 1, chain_id = 1, iter = 3000, refresh = 25)
save(chain_1, file = "Aug27_chain_1.RData")
chain_2 <- stan(fit = fit.compile, data = mvn_prec_data, chains = 1, chain_id = 2, iter = 3000, refresh = 25)
save(chain_2, file = "Aug27_chain_2.RData")


fit.compile <- stan(file = "the_slow_way_08_26.stan", data = mvn_prec_data, chains = 0)
chain_3 <- stan(fit = fit.compile, data = mvn_prec_data, chains = 1, chain_id = 3, iter = 3000, refresh = 50)
save(chain_3, file = "Aug27_chain_3.RData")
chain_4 <- stan(fit = fit.compile, data = mvn_prec_data, chains = 1, chain_id = 4, iter = 3000, refresh = 50)
save(chain_4, file = "Aug27_chain_4.RData")



fit.compile <- stan(file = "the_slow_way_08_26.stan", data = mvn_prec_data, chains = 0)
chain_5 <- stan(fit = fit.compile, data = mvn_prec_data, chains = 1, chain_id = 5, iter = 3000, refresh = 50)
save(chain_5, file = "Aug27_chain_5.RData")
chain_6 <- stan(fit = fit.compile, data = mvn_prec_data, chains = 1, chain_id = 6, iter = 3000, refresh = 50)
save(chain_6, file = "Aug27_chain_6.RData")


# 
# fit.compile <- stan(file = "the_slow_way_08_26.stan", data = mvn_prec_data, chains = 0)
# chain_7 <- stan(fit = fit.compile, data = mvn_prec_data, chains = 1, chain_id = 7, iter = 3000, refresh = 50)
# save(chain_7, file = "Aug27_chain_7.RData")
# chain_8 <- stan(fit = fit.compile, data = mvn_prec_data, chains = 1, chain_id = 8, iter = 3000, refresh = 50)
# save(chain_8, file = "Aug27_chain_8.RData")
# 
# 
# 
# fit.compile <- stan(file = "the_slow_way_08_26.stan", data = mvn_prec_data, chains = 0)
# chain_9 <- stan(fit = fit.compile, data = mvn_prec_data, chains = 1, chain_id = 9, iter = 3000, refresh = 50)
# save(chain_9, file = "Aug27_chain_9.RData")
# chain_10 <- stan(fit = fit.compile, data = mvn_prec_data, chains = 1, chain_id = 10, iter = 3000, refresh = 50)
# save(chain_10, file = "Aug27_chain_10.RData")
# 
# 
# 
# fit.compile <- stan(file = "the_slow_way_08_26.stan", data = mvn_prec_data, chains = 0)
# chain_11 <- stan(fit = fit.compile, data = mvn_prec_data, chains = 1, chain_id = 11, iter = 3000, refresh = 50)
# save(chain_11, file = "Aug27_chain_11.RData")
# chain_12 <- stan(fit = fit.compile, data = mvn_prec_data, chains = 1, chain_id = 12, iter = 3000, refresh = 50)
# save(chain_12, file = "Aug27_chain_12.RData")



# load finished chains
sflist <- list()
dir.path <- "Stan_tests/Aug_27th/all_chains/"

# for (i in 1:12) {
for (i in 1:6) {
  file <- paste0(dir.path,"Aug27_chain_",i,".RData")  
  load(file)
  sflist[i] <- get(paste0("chain_",i))
}


library(rstan)
stanfit <- sflist2stanfit(sflist)



mean_n_divergent <- sapply(get_sampler_params(stanfit), FUN = function(x) {
  mean(x[,"n_divergent__"])
})

round(mean_n_divergent, 3)



source("plots/pairs_new.R")

graphics.off()



for(i in 1:2){
  rho <- extract(stanfit, pars = paste0("rho[",i,"]"), permuted = F)
  pars <- c(paste0("beta[", i , ",", seq(1,24,4), "]"),
            paste0("rho[",i,"]"), "lp__")
  
  png(file = paste0("Stan_tests/Aug_27th/pairs_",i,".png"), 
      w = 10, h = 10, un = "in", res = 150) 
  pairs(stanfit, pars = pars, condition = rho > median(rho)) 
  dev.off()
}

rho <- extract(stanfit, pars = "rho[3]", permuted = F)
pars <- c(paste0("alpha[", seq(1,24,4), "]"),
          "rho[3]", "lp__")

png(file = paste0("Stan_tests/Aug_25th/pairs_",3,".png"), 
    w = 10, h = 10, un = "in", res = 150) 
pairs(stanfit, pars = pars, condition = rho > median(rho)) 
dev.off()


png(file = paste0("plots/2014_08_26/pairs_n_divergent.png"), 
    w = 10, h = 10, un = "in", res = 150) 
pairs(stanfit, pars = c(paste0("rho[", 1:3, "]"), paste0("tau[", 1:3, "]"), "lp__"),
      condition = "n_divergent__")
dev.off()

