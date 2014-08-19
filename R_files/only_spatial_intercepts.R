## CHANGES: ###

# made urbanpct and aapct into proportions and unionpop into z-score

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
urbanprop <- Data$urbanpct / 100 # make into proportion
aaprop <- Data$aapct / 100 # make into proportion
unionpop_z <- as.vector(scale(Data$unionpop)) # make into z-score

Z <- cbind(Data$dem, Data$laborcomm, urbanprop, aaprop, unionpop_z)


#### outcome variable ####
Y <- Data$prolabor


### sparse matrix stuff ###
A_sparse <- which(A == 1, arr.ind=TRUE)
# remove duplicates (because matrix is symmetric)
A_sparse <- A_sparse[A_sparse[,1] < A_sparse[,2],]
A_N <- dim(A_sparse)[1]
A1 <- A_sparse[,1]
A2 <- A_sparse[,2]


### data list for stan ###
mvn_prec_data <- list(Y = Y,
                      M = M,
                      Z = Z,
                      A = A, 
                      A_N = A_N, 
                      A1 = A1, 
                      A2 = A2,
                      d = diag(D),
                      R = R,
                      N = nrow(M),
                      C = ncol(Z))


library(rstan)
set_cppo("fast")

fit.compile <- stan(file = "only_spatial_intercepts.stan", data = mvn_prec_data, chains = 0)


fit.test <- stan(fit = fit.compile, data = mvn_prec_data, chains = 1, iter = 30, refresh = 5)

for(j in 11:12){
  stan_out <- stan(fit = fit.compile, data = mvn_prec_data, chains = 1, chain_id = j, iter = 2500, refresh = 50)
  chn <- paste0("chain_",j)
  assign(chn, stan_out)
}

save(chain_11, file = "Aug19chain_11.RData")
save(chain_12, file = "Aug19chain_12.RData")


