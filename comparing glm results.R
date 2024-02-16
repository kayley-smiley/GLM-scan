#comparing results
library(sf)
data(nysf)

#quasi-poisson, no significant clusters, but this is the MLC
cluster1 <- c(91, 88, 86, 92, 87, 89, 93, 90, 277)

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
