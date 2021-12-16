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
mat1 = as.matrix(read.csv("mat1.csv", row.names=1))
mat2 = as.matrix(read.csv("mat2.csv", row.names=1))
knn_ini = as.matrix(read.csv("corr.csv", row.names=1))
write.csv(mat1, "mat1_2.csv")
write.csv(mat2, "mat2_2.csv")
write.csv(knn_ini, "corr_2.csv")

# Use KNN as correspondence
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

