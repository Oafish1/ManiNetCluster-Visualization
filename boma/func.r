library(Seurat)
library(rdist)


#functions
get_form_seurat<-function(rdata,meta,type='COUNTS',hvg=6000,resolution=0.5) {
	form.data = c()
	if (type=='COUNTS') {
		obj <- CreateSeuratObject(counts = rdata, min.cells = 30, min.features = 200)
		obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
		obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
		obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = hvg)
		obj <- ScaleData(obj, features = row.names(obj))
		obj <- RunPCA(obj, features = VariableFeatures(object = obj))
		obj <- FindNeighbors(obj, dims = 1:10)
		obj <- FindClusters(obj, resolution = resolution)
		
		form.data <- RelativeCounts(data = obj[['RNA']]@data,scale.factor=1e6)
	} else {
		stop('for now, can only accept raw counts')
	}
	form.data = form.data[row.names(form.data) %in% VariableFeatures(obj),]
	#print(rownames(obj[['RNA']]@data))
	#print(row.names(meta))
	form.meta = meta[row.names(meta) %in% colnames(form.data),]
	form.meta$seurat_clusters = obj@meta.data$seurat_clusters

	return(list(form.data,form.meta))
}


get_form_data<-function(rdata,meta,type='NORM',hvg=6000,add=NULL) {#prefiltering and formalize the expression data,type='NORM' or 'COUNTS'
	form.data = c()
	if (type=='COUNTS') {
		obj <- CreateSeuratObject(counts = rdata, min.cells = 30, min.features = 200)
		obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
		obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
		obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = hvg)
		form.data <- RelativeCounts(data = obj[['RNA']]@data,scale.factor=1e6)
	} else if (type=='NORM') {
		form.data <- rdata
	}
	if (is.null(add)) {
		form.data = form.data[row.names(form.data) %in% VariableFeatures(obj),]
		form.meta = meta[row.names(meta) %in% colnames(form.data),]
	} else {
		form.data = form.data[row.names(form.data) %in% unique(unlist(c(VariableFeatures(obj),add))),]
		form.meta = meta[row.names(meta) %in% colnames(form.data),]
	}
	return(list(form.data,form.meta))
}

eu.dist<-function(p1,p2) { #given 2 vectors of positions, calculate  distance[e.g.eucleadian]
	return(sqrt(sum((p1-p2)^2)))
}

inn.dist<-function(p1,p2) { #normalized inner product
	p1=p1/norm(p1,type="2")
	p2=p2/norm(p2,type="2")
	return(p1 %*% p2)
}


get_pseudo_cells<-function(gex,keep.perc=.2,tag=null,pcn=20) {#form/average cells into pseudo cells---original sc exps too stochastic,input is the log-transformed matrix; keep.perc: percentage of pseudo cells VS original cells(%)
	n = dim(gex)[1]
	pca.res = prcomp(gex)
	pcs = gex %*% pca.res$rotation
	tpc = pcs[,1:pcn] #top PCs only 
	d = rdist(tpc,metric="euclidean")
	psm = cutree(hclust(d),k=min(floor(n*keep.perc),n))  #membership of cell clusters
	
	#average expression of each cell cluster
	avg.gex = matrix(0,nrow=length(unique(psm)),ncol=dim(gex)[2])
	for (id in unique(psm)) {
		tmp = gex[psm==id,]
		if (!is.null(dim(tmp))) {
			avg.gex[id,] = colMeans(tmp)
		} else {
			avg.gex[id,] = tmp #need or not? throw singletons?
		}
	}
	tmp.tag = rep(tag,dim(avg.gex)[1])
	ps.names = paste(c(1:dim(avg.gex)[1]),tmp.tag,sep='-')
	row.names(avg.gex) = ps.names
	colnames(avg.gex) = colnames(gex)
	
	return(list(avg.gex,tmp.tag))
}


get_pseudo_cells2<-function(gex,keep.perc=.2,tag=null,pcn=20,ctp=null,ctp.collapse.thr=.8) {
#form/average cells into pseudo cells---original sc exps too stochastic,input is the log-transformed matrix; keep.perc: percentage of pseudo cells VS original cells(%)
#if cell types list provided, also collapse on the cell type of pseudocell
	n = dim(gex)[1]
	pca.res = prcomp(gex)
	ctp = as.character(ctp);ctp[is.na(ctp)]='NA'
	pcs = gex %*% pca.res$rotation
	tpc = pcs[,1:pcn] #top PCs only 
	d = rdist(tpc)
	psm = cutree(hclust(d),k=min(floor(n*keep.perc),n))  #membership of cell clusters
	
	#average expression of each cell cluster
	avg.gex = matrix(0,nrow=length(unique(psm)),ncol=dim(gex)[2])
	ps.ctp = rep('NA',length(unique(psm)))
	for (id in unique(psm)) {
		tmp = gex[psm==id,]
		tmp.ctp = ctp[psm==id]
		if (!is.null(dim(tmp))) {
			avg.gex[id,] = colMeans(tmp)
			sorts = sort(table(tmp.ctp)/sum(table(tmp.ctp)),decreasing=T)
			if (sorts[1] >= ctp.collapse.thr) {
				ps.ctp[id] = names(sorts)[1]
			}
		} else {#need or not? throw singletons?
			avg.gex[id,] = tmp 
			ps.ctp[id] = tmp.ctp
		}
	}
	tmp.tag = rep(tag,dim(avg.gex)[1])
	ps.names = paste(c(1:dim(avg.gex)[1]),tmp.tag,sep='-')
	row.names(avg.gex) = ps.names
	colnames(avg.gex) = colnames(gex)
	tmp.tag = data.frame('tag'=rep(tag,dim(avg.gex)[1]),'ctp'=ps.ctp)
	
	return(list(avg.gex,tmp.tag))
}


get_pseudo_cells3<-function(gex,keep.perc=.2,tag=null,pcn=20,ctp=null,ctp.collapse.thr=.8) {
#form/average cells into pseudo cells---original sc exps too stochastic,input is the log-transformed matrix; keep.perc: percentage of pseudo cells VS original cells(%)
#if cell types list provided, also collapse on the cell type of pseudocell``
#record which cell belongs to which pseudocell
	n = dim(gex)[1]
	pca.res = prcomp(gex)
	ctp = as.character(ctp);ctp[is.na(ctp)]='NA'
	pcs = gex %*% pca.res$rotation
	tpc = pcs[,1:pcn] #top PCs only 
	d = rdist(tpc)
	psm = cutree(hclust(d),k=min(floor(n*keep.perc),n))  #membership of cell clusters
	
	#average expression of each cell cluster
	avg.gex = matrix(0,nrow=length(unique(psm)),ncol=dim(gex)[2])
	ps.ctp = rep('NA',length(unique(psm)))
	for (id in unique(psm)) {
		tmp = gex[psm==id,]
		tmp.ctp = ctp[psm==id]
		if (!is.null(dim(tmp))) {
			avg.gex[id,] = colMeans(tmp)
			sorts = sort(table(tmp.ctp)/sum(table(tmp.ctp)),decreasing=T)
			if (sorts[1] >= ctp.collapse.thr) {
				ps.ctp[id] = names(sorts)[1]
			}
		} else {#need or not? throw singletons?
			avg.gex[id,] = tmp 
			ps.ctp[id] = tmp.ctp
		}
	}
	tmp.tag = rep(tag,dim(avg.gex)[1])
	ps.names = paste(c(1:dim(avg.gex)[1]),tmp.tag,sep='-')
	row.names(avg.gex) = ps.names
	colnames(avg.gex) = colnames(gex)
	tmp.tag = data.frame('tag'=rep(tag,dim(avg.gex)[1]),'ctp'=ps.ctp)
	
	tmp.belong = paste(psm,tag,sep='-')
	
	return(list(avg.gex,tmp.tag,tmp.belong,gex))
}

##generate pseudo cells based on 'tag' info in the meta table;
pcell.from.seurat<-function(gex,meta) { #function to generate pseudo-cells
	tgex = t(gex)
	tmeta = meta

	ps.mat = c()
	ps.tag = c()
	ps.belong = data.frame(meta$tag2)
	names = c()
	ps.gex = gex

	for(t in unique(tmeta$tag)) {
		tmp.gex = tgex[tmeta$tag==t,]
		tmp.meta = tmeta[tmeta$tag==t,]
		if (!is.null(dim(tmp.gex))) {
			#for (t in unique(tmeta$tag)) {
			for(cls in unique(tmeta$seurat_clusters)) {
	        		print(paste(t,cls))
				tmp.gex2 = tmp.gex[tmp.meta$seurat_clusters==cls,]
				tmp.meta2 = tmp.meta[tmp.meta$seurat_clusters==cls,]
				if (is.null(dim(tmp.gex2))) {
					#ps.mat = rbind(ps.mat,tmp.gex2)
					#ps.tag = rbind(ps.tag,c(cls,t,tmp.meta2[1,]$time,names(sort(table(tmp.meta2$ctp),decreasing=T))[1]))
					#single cell throw out
				} else if (dim(tmp.gex2)[1] <= 1) { 
				        #no data
				}else {
					ps.mat = rbind(ps.mat,colMeans(tmp.gex2))
					ps.tag = rbind(ps.tag,c(cls,t,tmp.meta2[1,]$time,names(sort(table(tmp.meta2$ctp),decreasing=T))[1]))
					names = c(names,tmp.meta2$tag2[1])
				}
				#print(paste(t,tmp.meta2[1,]$Age))
			}
			#ps.belong = c(ps.belong,ps.res[[3]])
			#ps.gex = rbind(ps.gex,ps.res[[4]])
		}
	}
	row.names(ps.tag) = names
	return(list(ps.mat,ps.tag,ps.belong,ps.gex))
}


##generate pseudo cells based on 'tag' info in the meta table;[tag includes sample source & celltype]
pcell.from.fctp_seurat<-function(gex,meta) { #function to generate pseudo-cells
	tgex = t(gex)
	tmeta = meta

	ps.mat = c()
	ps.tag = c()
	ps.belong = data.frame(meta$tag2)
	names = c()
	ps.gex = gex


	for(t in unique(tmeta$tag)) {
		tmp.gex = tgex[tmeta$tag==t,]
		tmp.meta = tmeta[tmeta$tag==t,]
		if (!is.null(dim(tmp.gex))) {
			#for (t in unique(tmeta$tag)) {
			for(cls in unique(tmeta$seurat_clusters)) {
	        		print(paste(t,cls))
				tmp.gex2 = tmp.gex[tmp.meta$seurat_clusters==cls,]
				tmp.meta2 = tmp.meta[tmp.meta$seurat_clusters==cls,]
				if (is.null(dim(tmp.gex2))) {
					#do nothing
				} else if (dim(tmp.gex2)[1] <= 1) { 
				        #do nothing
				}else {
					for (ctp in unique(tmp.meta2$fctp)) {
						print(paste(t,cls,ctp))
						tmp.gex3 = tmp.gex2[tmp.meta2$fctp==ctp,]
						tmp.meta3 = tmp.meta2[tmp.meta2$fctp==ctp,]
						if (is.null(dim(tmp.gex3))) {
							#do nothing
						} else if (dim(tmp.gex3)[1] == 1) {
							##only one cell,do nothing
							ps.mat = rbind(ps.mat,colMeans(tmp.gex3))
							ps.tag = rbind(ps.tag,c(cls,t,tmp.meta3[1,]$time,tmp.meta3[1,]$fctp,1,tmp.meta3[1,]$time_weeks))

						} else {
							ps.mat = rbind(ps.mat,colMeans(tmp.gex3))
							ps.tag = rbind(ps.tag,c(cls,t,tmp.meta3[1,]$time,names(sort(table(tmp.meta3$fctp),decreasing=T))[1],dim(tmp.gex3)[1],tmp.meta3[1,]$time_weeks))
							names = c(names,tmp.meta3$tag2[1])
						}

					}
				}
			}
		}
	}
	row.names(ps.tag) = names
	return(list(ps.mat,ps.tag,ps.belong,ps.gex))
}



runMSMA_cor<-function(mat1,mat2,k=5) {
#mat has genes on the column !!!!!!!!!!!!!!!!!
	sim_mat = c() #matrix to store the similarity results after 1st step
	knn_ini = c() #matrix to store the correspondence matrix after 1st step align
	df1 = NULL     #for validation

	sim_mat = cor(t(mat1),t(mat2))
	#knn
	knn_ini=t(matrix(0,nrow=dim(sim_mat)[2],ncol=dim(sim_mat)[1]))
	cm = apply(sim_mat,2,function(x) order(x,decreasing=T)[1:k])
	if (k > 1) {
		for (i in c(1:dim(cm)[2])) {
			knn_ini[cbind(cm[,i]),i] = 1
		}
	} else if (k==1) {
		for (i in c(1:length(cm))) {
			knn_ini[cm[i],i] = 1
		}
	}
	
	#after get the knn_ini, use it as correspondence matrix for 2nd align
	XY_corr=Correspondence(matrix=knn_ini)
	df2=ManiNetCluster(mat1,mat2,nameX='sample1',nameY='sample2',corr=XY_corr,d=3L,method='nonlinear manifold aln',k_NN=3L,k_medoids=6L) #NMA
	
	return(list(sim_mat,knn_ini,df2,df1))

}


runMSMA_mw<-function(mat1,mat2,k=5) {
#mat has genes on the column !!!!!!!!!!!!!!!!!
	sim_mat = c() #matrix to store the similarity results after 1st step
	knn_ini = c() #matrix to store the correspondence matrix after 1st step align
	df1 = NULL     #for validation

	df=ManiNetCluster(mat1,mat2,nameX='sample1',nameY='sample2',d=10L,method='manifold warping',k_NN=3L,k_medoids=3L)
	pair_dist = apply(df[df$data=='sample1',c(3:12)],1,function(x) {
			d = apply(df[df$data=='sample2',c(3:12)],1,function(y) eu.dist(x,y))
	})
	row.names(pair_dist)=row.names(mat2)
	colnames(pair_dist)=row.names(mat1)
	range01 <- function(x){(x-min(x))/(max(x)-min(x))}
	pair_dist = range01(pair_dist)
	sim_mat = 1/(1+pair_dist)
	df1 = df  #range01 OR NOT?
	#knn
	k1=k;k2=k
	knn_mat1 = matrix(0,nrow=dim(pair_dist)[1],ncol=dim(pair_dist)[2]) #for row
	knn_mat2 = matrix(0,nrow=dim(pair_dist)[1],ncol=dim(pair_dist)[2]) #for col
	for(i in 1:dim(knn_mat1)[1]) {#for row
		knn_mat1[i,order(pair_dist[i,])[1:k1]]=1
	}
	for(i in 1:dim(knn_mat2)[2]) {#for col
			knn_mat2[order(pair_dist[,i])[1:k2],i]=1
	}
	knn_ini = t(1*(knn_mat1 & knn_mat2))
	
	#after get the knn_ini, use it as correspondence matrix for 2nd align
	XY_corr=Correspondence(matrix=knn_ini)
	df2=ManiNetCluster(mat1,mat2,nameX='sample1',nameY='sample2',corr=XY_corr,d=3L,method='nonlinear manifold aln',k_NN=3L,k_medoids=6L) #NMA
	
	return(list(sim_mat,knn_ini,df2,df1))

}


runMSMA_dtw<-function(mat1,mat2) { #method can be 'cor','mw' or 'dtw'
#mat has genes on the column !!!!!!!!!!!!!!!!!
#1st step has 3 options:1)correlation;2)Manifod warping;3)Dynamic Time Warping
	sim_mat = c() #matrix to store the similarity results after 1st step
	knn_ini = c() #matrix to store the correspondence matrix after 1st step align
	df1 = NULL     #for validation

	library(dtw)
	cor_mat = cor(t(mat1),t(mat2))
	dist_mat = 1/(1+cor_mat)
	dtw_res=dtw(dist_mat,keep=TRUE,step=asymmetric,open.end=T,open.begin=T)  #open beginning and ending for flexibility
	knn_ini = matrix(0,nrow=nrow(dist_mat),ncol=ncol(dist_mat))
	knn_ini[cbind(dtw_res$index1,dtw_res$index2)]=1
	sim_mat=NULL #position holder
	
	#after get the knn_ini, use it as correspondence matrix for 2nd align
	XY_corr=Correspondence(matrix=knn_ini)
	df2=ManiNetCluster(mat1,mat2,nameX='sample1',nameY='sample2',corr=XY_corr,d=3L,method='nonlinear manifold aln',k_NN=3L,k_medoids=6L) #NMA
	
	return(list(sim_mat,knn_ini,df2,df1))

}



runMSMA_diag<-function(mat1,mat2) { #method can be 'cor','mw' or 'dtw'
#mat has genes on the column !!!!!!!!!!!!!!!!!
#1st step has 3 options:1)correlation;2)Manifod warping;3)Dynamic Time Warping
	sim_mat = c() #matrix to store the similarity results after 1st step
	knn_ini = c() #matrix to store the correspondence matrix after 1st step align
	df1 = NULL     #for validation

	if (nrow(mat1) != nrow(mat2)) {
		stop('for diag case---the numerb of samples need to be same')
	}else {
		knn_ini = diag(nrow(mat1))
	}
	
	#after get the knn_ini, use it as correspondence matrix for 2nd align
	XY_corr=Correspondence(matrix=knn_ini)
	df2=ManiNetCluster(mat1,mat2,nameX='sample1',nameY='sample2',corr=XY_corr,d=3L,method='nonlinear manifold aln',k_NN=3L,k_medoids=6L) #NMA
	
	return(list(sim_mat,knn_ini,df2,df1))

}








runMSMA<-function(mat1,mat2,method='cor',k=5) { #method can be 'cor','mw' or 'dtw'
#mat has genes on the column !!!!!!!!!!!!!!!!!
#1st step has 3 options:1)correlation;2)Manifod warping;3)Dynamic Time Warping
	sim_mat = c() #matrix to store the similarity results after 1st step
	knn_ini = c() #matrix to store the correspondence matrix after 1st step align
	df1 = NULL     #for validation
	if (method=='cor') {
		sim_mat = cor(t(mat1),t(mat2))
		
		#knn
		knn_ini=t(matrix(0,nrow=dim(sim_mat)[2],ncol=dim(sim_mat)[1]))
		cm = apply(sim_mat,2,function(x) order(x,decreasing=T)[1:k])
		for (i in c(1:dim(cm)[2])) {
			knn_ini[cbind(cm[,i]),i] = 1
		}
	}else if (method=='mw') {
		df=ManiNetCluster(mat1,mat2,nameX='sample1',nameY='sample2',d=10L,method='manifold warping',k_NN=3L,k_medoids=3L)
		pair_dist = apply(df[df$data=='sample1',c(3:12)],1,function(x) {
				d = apply(df[df$data=='sample2',c(3:12)],1,function(y) eu.dist(x,y))
		})
		row.names(pair_dist)=row.names(mat2)
		colnames(pair_dist)=row.names(mat1)
		range01 <- function(x){(x-min(x))/(max(x)-min(x))}
		pair_dist = range01(pair_dist)
		sim_mat = 1/(1+pair_dist)
		df1 = df  #range01 OR NOT?
		#knn
		k1=k;k2=k
		knn_mat1 = matrix(0,nrow=dim(pair_dist)[1],ncol=dim(pair_dist)[2]) #for row
		knn_mat2 = matrix(0,nrow=dim(pair_dist)[1],ncol=dim(pair_dist)[2]) #for col
		for(i in 1:dim(knn_mat1)[1]) {#for row
			knn_mat1[i,order(pair_dist[i,])[1:k1]]=1
		}
		for(i in 1:dim(knn_mat2)[2]) {#for col
				knn_mat2[order(pair_dist[,i])[1:k2],i]=1
		}
		knn_ini = t(1*(knn_mat1 & knn_mat2))

	}else if (method=='dtw') {
		library(dtw)
		cor_mat = cor(t(mat1),t(mat2))
		dist_mat = 1/(1+cor_mat)
		dtw_res=dtw(dist_mat,keep=TRUE,step=asymmetric,open.end=T,open.begin=T)  #open beginning and ending for flexibility
		knn_ini = matrix(0,nrow=nrow(dist_mat),ncol=ncol(dist_mat))
		knn_ini[cbind(dtw_res$index1,dtw_res$index2)]=1
		sim_mat=NULL #position holder
	}else if (method=='diag') {#use diagnal matrix as correspondence matrix---for 1-1 match
		if (nrow(mat1) != nrow(mat2)) {
			stop('for diag case---the numerb of samples need to be same')
		}else {
			knn_ini = diag(nrow(mat1))
		}
	}else {
		stop('not recognized method')
	}
	
	#after get the knn_ini, use it as correspondence matrix for 2nd align
	XY_corr=Correspondence(matrix=knn_ini)
	df2=ManiNetCluster(mat1,mat2,nameX='sample1',nameY='sample2',corr=XY_corr,d=3L,method='nonlinear manifold aln',k_NN=3L,k_medoids=6L) #NMA
	
	return(list(sim_mat,knn_ini,df2,df1))

}








##########FUNCTION that identifies DEGs across two separate experiments#########










