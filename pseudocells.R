### To create Pseudocells from multiple GEX datasets...</h3>
### Filter the datasets and generate clusters
for (data in datasets) {
	#Generate a Seurat object
	obj <- CreateSeuratObject(~)

	#Normalize the data and perform PCA on most variable features
	obj <- FindVariableFeatures(obj, nfeatures=6000)
	obj <- ScaleData(obj, features=row.names(obj))
	obj <- RunPCA(obj, features = VariableFeatures(object=obj))

	#Compute the KNN and cluster
	obj <- FindNeighbors(obj, dims=1:10)
	obj <- FindClusters(obj, resolution=0.5)

	#Scale the counts and produce the final output
	obj <- RelativeCounts(data=obj[['RNA']]@data, scale.factor=1e6)
	output_data <- input_data[row.names(input_data) %in% VariableFeatures(obj),]
	output_meta <- input_meta[row.names(input_meta) %in% colnames(output_data),]
	output_meta$seurat_clusters = obj@meta.data$seurat_clusters
}


### Combine the datasets
gex <- cbind(~)
meta <- cbind(~)


### Generate the pseudocells
## Perform hierarchical clustering on PCA data or use existing time metadata
# Use automatic
n = dim(gex)[1]
pca.res = prcomp(gex)
pcs = gex %*% pca.res$rotation
tpc = pcs[, 1:20]
d = rdist(tpc, metric="euclidean")
psm = cutree(hclust(d),k=min(floor(.2 * n), n))
# Average the expression of each cluster
avg.gex = matrix(0, nrow=length(unique(psm)), ncol=dim(gex)[2])
for (id in unique(psm)) {avg.gex[id,] = colMeans(gex[psm==id,])}
return avg.gex

# Use time
for (time in unique(meta$time)) {
	for (clust in unique(meta$seurat_clusters)) {
		# Find and assign average gex[in time/clust,] to avg.gex
		# Record tag, clust to avg.meta
	}
}
return List(avg.gex, avg.meta)