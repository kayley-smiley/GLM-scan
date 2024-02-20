#comparing results
library(sf)
data(nysf)

#quasi-poisson, no significant clusters, but this is the MLC 
#this was with the LLR subtraction order switched (null - alt)
cluster1 <- c(91, 88, 86, 92, 87, 89, 93, 90, 277)

cluster_col = rep(NA, nrow(nysf))
cluster_col[cluster1] <- "cadetblue"


plot(st_geometry(nysf), col = cluster_col, border = "lightgrey")

#quasi-poisson, one significant cluster
#this was with the LLR subtraction order fixed (alt - null)
cluster1 <- c(42, 43, 45, 44, 41, 46, 40, 39, 38, 53, 37, 54, 52)

cluster_col = rep(NA, nrow(nysf))
cluster_col[cluster1] <- "cadetblue"


plot(st_geometry(nysf), col = cluster_col, border = "lightgrey")


#poisson glm
cluster1 <- c(8, 86, 87, 89, 92, 91, 85, 93, 84, 90, 259)
cluster2 <- c(33, 32, 18, 9, 8, 7, 17, 11, 10, 6, 34, 5, 12, 4, 14, 13, 16, 3, 27, 15, 2,  
              1, 35, 50, 49, 31, 25, 47, 48, 52, 24, 51, 26, 36, 28, 55, 37, 38, 39, 53, 40, 44, 43, 46)

cluster_col = rep(NA, nrow(nysf))
cluster_col[cluster1] <- "cadetblue"
cluster_col[cluster2] <- "pink"

plot(st_geometry(nysf), col = cluster_col, border = "lightgrey")



#what happens if there are no cases in a region and we use the glm function,
#specifically, what is beta? the closed form solutions give -infinity for beta

#candidate zone 363 has no cases, it only consists of region 8, which has no cases
ind <- indmat[,363]

test <- glm(formula = cases ~ ind, offset = log(pop), data =  nysf, family = poisson)










