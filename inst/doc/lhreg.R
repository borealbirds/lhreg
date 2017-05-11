## ------------------------------------------------------------------------
data(lhreg_data)
str(lhreg_data)

## ----eval=FALSE----------------------------------------------------------
#  library(ape)
#  mph <- read.nexus("11960.tre") # 1000 trees with Ericson backbone
#  cph <- consensus(mph)
#  table(sapply(mph, function(z) length(z$tip.label)))
#  CORR <- TRUE
#  vv <- list()
#  vv[[1]] <- vcv(mph[[1]], corr=CORR)
#  for (i in 2:length(mph)) {
#      v <- vcv(mph[[i]], corr=CORR)
#      v <- v[rownames(vv[[1]]), colnames(vv[[1]])]
#      vv[[i]] <- v
#  }
#  vvv <- v
#  for (i in 1:length(v)) {
#      vvv[i] <- mean(sapply(vv, function(z) z[i]))
#  }
#  spp <- intersect(rownames(lhreg_data), rownames(vvv))
#  vvv <- vvv[spp,spp]
#  cor_matrix <- as.matrix(nearPD(vvv, corr=TRUE)$mat)

## ------------------------------------------------------------------------
data(cor_matrix)
str(cor_matrix)
heatmap(cor_matrix)

