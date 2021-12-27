# Install ManiNetCluster
# If error delete r-miniconda in AppData/Local
#library(reticulate)
#use_python("C:\\Users\\nck\\AppData\\Local\\r-miniconda\\envs\\r-reticulate\\python.exe", required=TRUE)
#conda_install("r-reticulate", c("pandas", "scipy", "matplotlib"))
#conda_install("r-reticulate", c("sklearn"), pip=TRUE)
#devtools::install_github("namtk/ManiNetCluster", INSTALL_opts=c("--no-multiarch"))
library("ManiNetCluster")

# Load data
load("bulk1.RData")
load("bulk2.RData")
write.csv(mat1, "data/mat1.csv")
write.csv(mat2, "data/mat2.csv")
write.csv(knn_ini, "data/corr.csv")
write.csv(sel.meta1, "data/meta1.csv")
write.csv(sel.meta2, "data/meta2.csv")
mat1 = as.matrix(read.csv("data/mat1.csv", row.names=1))
mat2 = as.matrix(read.csv("data/mat2.csv", row.names=1))
knn_ini = as.matrix(read.csv("data/corr.csv", row.names=1))

# Use KNN as correspondence
dim(mat1)
dim(mat2)
cor(mat1, mat2)
XY_corr=Correspondence(matrix=knn_ini)
# Run NLMA
ManiNetCluster(
  mat1,mat2,
  nameX='sample1',nameY='sample2',
  corr=XY_corr,
  d=3L,
  method='nonlinear manifold aln',
  k_NN=3L,
  k_medoids=6L
)
