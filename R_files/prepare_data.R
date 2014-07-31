# Load data ---------------------------------------------------------------
# _________________________________________________________________________

setwd("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/Stan/My Code/farhang_katznelson_rep")

loadData.path <- "Data/rollcall.logit.inter.all.names.txt"
Data <- read.table(loadData.path, header = TRUE)


# Construct map matrix ----------------------------------------------------
# _________________________________________________________________________

nr <- with(Data, length(unique(region)))
nc <- with(Data, length(unique(period)))

map.matrix <- mat.or.vec(nr,nc)
rownames(map.matrix) <- c("NS", "BS", "DS")
colnames(map.matrix) <- paste0("t", 1:nc)

map.matrix[1, ] <- c(3, rep(5, nc - 2), 3)
map.matrix[2, ] <- c(5, rep(8, nc - 2), 5)
map.matrix[3, ] <- c(3, rep(5, nc - 2), 3)
map.matrix

# R = number of unique region/period combinations (i.e. number cells in map.matrix)
R <- prod(dim(map.matrix))

# N = total number of neighbors (the length of the map vector to be created)
N <- 3*sum(map.matrix == 3) + 5*sum(map.matrix == 5) + 8*sum(map.matrix == 8)

# just as a guide/reference: number the map cells from 1 to 24
map.index <- matrix(1:R, nr, nc)
colnames(map.index) <- colnames(map.matrix)
rownames(map.index) <- rownames(map.matrix)
map.index




# Fill in map vector to pass to Stan --------------------------------------
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


# Prepare rest of data for Stan -------------------------------------------
# _________________________________________________________________________

# offsets vector
off <- c(0, cumsum(map.matrix))

# I = number of observations
I <- nrow(Data)



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


# Create Nneighs vector containing number of neighbors for each region
Nneighs <- sapply(1:R, function(r) off[r+1] - off[r]) # or equivalently Nneighs <- c(map.matrix)




# Make data list to pass to Stan ------------------------------------------
# _________________________________________________________________________

data.list <- list(I = I, 
                  R = R, 
                  N = N,
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

saveData.path <- "Data/prepared_data.RData"
save(data.list, file = saveData.path)
