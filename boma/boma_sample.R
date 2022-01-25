#given two matrices[normalized], perform alignment
#Steps:
#1. Seurat scale both matrices & clean the matrices [organize meta info, select HVGs]
#2. generate pseudo-cells & correspondence based on seurat clusters
#3. perform MA alignment on the matrix
#4. output figures>>distance heatmap;cell-type enrichment

library(pheatmap)
library(Matrix)
library(stringr)
library(reticulate)
library(ManiNetCluster)
source('./boma/func.r')


##############################start#################
load('all.data.RData')

nTopGenes = 15000
##Seurat to normalize & scale data matrix 1
tmp.form=get_form_seurat(org.data1,org.meta1,type='COUNTS',hvg=nTopGenes,resolution=10)
form.data1 = tmp.form[[1]]
form.meta1 = tmp.form[[2]]

##Seurat to normalize & scale data matrix 2
tmp.form=get_form_seurat(human.data1,human.meta1,type='COUNTS',hvg=nTopGenes,resolution=10)
form.data2 = tmp.form[[1]]
form.meta2 = tmp.form[[2]]

# CHANGED UP TO HERE

########################################################
##overlap with human brain developmental genes##
sel.genes = intersect(intersect(row.names(form.data1),row.names(form.data2)),unique(all.rec$gene))
##normalize 
sel.data1 = form.data1[row.names(form.data1) %in% sel.genes,]; sel.data1 = sel.data1[!duplicated(row.names(sel.data1)),]; sel.meta1=form.meta1
sel.data2 = form.data2[row.names(form.data2) %in% sel.genes,]; sel.data2 = sel.data2[!duplicated(row.names(sel.data2)),]; sel.meta2=form.meta2 
#log transform
exp1 = log(sel.data1+1)
exp2 = log(sel.data2+1)
#reorder
exp1 = exp1[order(row.names(exp1)),order(sel.meta1$time)];sel.meta1=sel.meta1[order(sel.meta1$time),]
exp2 = exp2[order(row.names(exp2)),order(sel.meta2$time)];sel.meta2=sel.meta2[order(sel.meta2$time),]

print('generating pseudo cells')
source('src/func.r')
###########################################
##pcells for sample1
sel.meta1$tag2 = paste(sel.meta1$tag,sel.meta1$seurat_clusters,sel.meta1$fctp,sep='-')
pcell1.info = pcell.from.fctp_seurat(exp1,sel.meta1)

##pcells for sample2
sel.meta2$tag2 = paste(sel.meta2$tag,sel.meta2$seurat_clusters,sel.meta2$fctp,sep='-')
sel.meta2$fctp[is.na(sel.meta2$fctp)] = 'Unknown' 
pcell2.info = pcell.from.fctp_seurat(exp2,sel.meta2)

ps.mat1=pcell1.info[[1]]
ps.tag1 = pcell1.info[[2]];
ps.belong1 = pcell1.info[[3]]
ps.gex1 = pcell1.info[[4]]
#ps.sc.time1 = tmeta$time
#ps.sc.ctp1 = data.frame('orig_ctp'=tmeta$WGCNAcluster, 'ctp'=tmeta$ctp)
ps.mat2=pcell2.info[[1]]
ps.tag2 = pcell2.info[[2]]
ps.belong2 = pcell2.info[[3]]
ps.gex2 = pcell2.info[[4]]


print('performing alignment')
###########perform alignment#############


##focus on overlapped cell types between sample sources only!!!!!!!!!!!!!!!##############
ps.mat01=ps.mat1;ps.mat02=ps.mat2;ps.tag01=ps.tag1;ps.tag02=ps.tag2
ps.mat1 = ps.mat1[(ps.tag1[,4] %in% c('Astro','EN','IN','IPC','OPC','RG')),]
ps.tag1 = ps.tag1[(ps.tag1[,4] %in% c('Astro','EN','IN','IPC','OPC','RG')),]
ps.mat2 = ps.mat2[(ps.tag2[,4] %in% c('Astro','EN','IN','IPC','OPC','RG')),]
ps.tag2 = ps.tag2[(ps.tag2[,4] %in% c('Astro','EN','IN','IPC','OPC','RG')),]

########
algn_res=runMSMA_cor(ps.mat1,ps.mat2,k=1)