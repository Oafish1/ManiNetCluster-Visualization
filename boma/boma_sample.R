#given two matrices[normalized], perform alignment
#Steps:
#1. Seurat scale both matrices & clean the matrices [organize meta info, select HVGs]
#2. generate pseudo-cells & correspondence based on seurat clusters
#3. perform MA alignment on the matrix
#4. output figures>>distance heatmap;cell-type enrichment

library(pheatmap)
library(Matrix)
library(stringr)
Sys.setenv(RETICULATE_PYTHON = "/home/hechenfon/anaconda3/bin/python")
library(reticulate)
path_to_python<-'/home/hechenfon/anaconda3/bin/python'
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

#relabel the time manually---for human brain development stages
sts=as.numeric(as.character(form.meta2$time_weeks))
for (i in c(1:length(sts))) {
  t = sts[i]
  if (t<8) {
    sts[i]=1
  } else if (t>=8 & t<10) {
    sts[i]=2
  } else if (t>=10 & t<13) {
    sts[i]=3
  } else if (t>=13 & t<16) {
    sts[i] = 4
  } else if (t>=16 & t<19) {
    sts[i] = 5
  } else if (t>=19 & t<24) {
    sts[i] = 6
  } else if (t>=24 & t<38) {
    sts[i] = 7
  } else {
    sts[i] = 8
  }
}
form.meta2$time = sts

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

print('analyzing algnment res')
#########analyze aligned res.#############
library(RColorBrewer)
df2=algn_res[[3]]
df2$time = as.numeric(as.character(c(ps.tag1[,3],ps.tag2[,3])))
df2$tag = as.character(c(ps.tag1[,4],ps.tag2[,4]))
df2$time_weeks = as.numeric(as.character(c(ps.tag1[,6],ps.tag2[,6])))
row.names(df2) = c(row.names(ps.tag1),row.names(ps.tag2))



pair_dist = apply(df2[df2$data=='sample1',c(3:5)],1,function(x) {
  d = apply(df2[df2$data=='sample2',c(3:5)],1,function(y) eu.dist(x,y))
})
row.names(pair_dist)=row.names(ps.tag2)
colnames(pair_dist)=row.names(ps.tag1)
sim_dist = 1/(1+pair_dist)
#############FOR cells+cells similarity bi-clustering#######
sim_mat = t(sim_dist)
write.table(sim_mat,file='dist.csv',row.names=T,col.names=T,sep=',',quote=F)

####################analyzing the alignment results###############################
#pheatmap(t(sim_dist),show_rownames=F,show_colnames=F,cluster_rows=F,cluster_cols=F)

#plot 3D on single source 
library(plot3D)
library(plyr)

#plot neural cells only
int.ctp = c('EN','IN','RG','IPC','OPC')
#pdf('sample2.algn.neural.pdf')
#res = data.frame(df2[df2$data=='sample1' & df2$tag %in% c('EN','IN','IPC'),])1res = data.frame(df2[df2$data=='sample2',])
#res = data.frame(df2[df2$data=='sample2' & df2$cellPerc>10e-4,])
res = data.frame(df2[df2$data=='sample2' & df2$tag %in% int.ctp,])
#tag.cols = colorRampPalette(brewer.pal(9, "Set1"))(10)[c(1,3,5:10)] #color for smaple1
tag.cols = colorRampPalette(brewer.pal(6, "Set1"))(6)[c(4,2,3,5,1)] #colors
s3d<-scatter3D(x=res[,3],y=res[,4],z=res[,5],col='grey',bg=tag.cols[as.numeric(mapvalues(res$tag,names(table(res$tag)),c(1:length(unique(res$tag)))))],lwd=1,pch=21,
               colkey=F,theta = 90, phi = 0,cex=1.5,
               xlim=c(min(df2$Val0),max(df2$Val0)),
               ylim=c(min(df2$Val1),max(df2$Val1)),
               zlim=c(min(df2$Val2),max(df2$Val2)))
legend("top", legend = levels(as.factor(res$data)), pch = c(16, 17),inset = -0.1, xpd = TRUE, horiz = TRUE)
legend("left", legend = levels(as.factor(res$tag)), col = tag.cols,pch=16,inset = -0.1, xpd = TRUE, horiz = FALSE)
#dev.off()





