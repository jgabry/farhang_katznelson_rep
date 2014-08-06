setwd("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/Rsync/The_Real_Deal")

loadData.path <- "rollcall.logit.inter.all.names.txt"
Data <- read.table(loadData.path, header = TRUE)

#### Create region vector #### 
# Regions (from 1 to 24) are assigned to correspond to the cell numbers in map.index

library(plyr) 
Data$region.name <- mapvalues(Data$region, 
                              # Note: region codes in data: DS = 1, BS = 2, NS, = 3
                              from = 1:3, 
                              to = c("DS", "BS", "NS"))

Region <- rep(NA, 24)
for(i in 1:8){
  Region[with(Data, region.name == "NS" & period == i)] <- seq(1, 22, 3)[i]
  Region[with(Data, region.name == "BS" & period == i)] <- seq(2, 23, 3)[i]
  Region[with(Data, region.name == "DS" & period == i)] <- seq(3, 24, 3)[i]
}



FNs <- c("find_neighbors","make_adjacency_matrix","make_degree_matrix",
         "make_design_matrix_array","make_incidence_matrix","make_sparse",
         "make_stan_data_list")

for(fn in FNs){
  source(paste0("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/functions/",fn,".R"))  
}

data.list <- do.call("make_stan_data_list", 
                     args = list(
                       df = Data, 
                       map.dims = c(3,8), 
                       region.vector = Region, 
                       Xvarlist = c("urbanpct","aapct","unionpop"),
                       Zvarlist = c("dem", "laborcomm"),
                       Yname = "prolabor")
                     )




