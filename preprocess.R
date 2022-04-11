# Load data
load("bulk1.RData")
load("bulk2.RData")
meta1 = sel.meta1
meta2 = sel.meta2
# mat1, meta1, mat2, meta2

load("sc2sc.small.RData")
meta1 = ps.meta1
meta2 = ps.meta2
# mat1, meta1, mat2, meta2

load("all.data.RData")
# human.data1, human.data2, human.data3
# human.meta1, human.meta2, human.meta3
# org.data1, org.data2
# org.meta1, org.meta2


### Misc
write.csv(knn_ini, "corr.csv")


### Works
# Vars
prefix = "data/NowakowskiScience2017/"
mat = t(human.data1)
meta = human.meta1

# Write Data
paste("Mat Dim:", dim(mat)[1], "x", dim(mat)[2], "  Meta Dim:", dim(meta)[1], "x", dim(meta)[2])
dir.create(prefix)
write.csv(mat, paste(prefix, "mat.csv", sep=""))
write.csv(meta, paste(prefix, "meta.csv", sep=""))


# Vars
prefix = "data/BireyNature2017/"
mat = t(org.data1)
meta = org.meta1

# Write Data
paste("Mat Dim:", dim(mat)[1], "x", dim(mat)[2], "  Meta Dim:", dim(meta)[1], "x", dim(meta)[2])
dir.create(prefix)
write.csv(mat, paste(prefix, "mat.csv", sep=""))
write.csv(meta, paste(prefix, "meta.csv", sep=""))



### Doesn't work
# Vars
prefix = "data/TrevinoCell2021/"
mat = t(as.matrix(human.data2)) # This is problematic
meta = human.meta2

# Write Data
paste("Mat Dim:", dim(mat)[1], "x", dim(mat)[2], "  Meta Dim:", dim(meta)[1], "x", dim(meta)[2])
dir.create(prefix)
write.csv(mat, paste(prefix, "mat.csv", sep=""))
write.csv(meta, paste(prefix, "meta.csv", sep=""))


# Vars
prefix = "data/BhaduriNature2020Brain/"
mat = human.data3
meta = human.meta3

# Write Data
paste("Mat Dim:", dim(mat)[1], "x", dim(mat)[2], "  Meta Dim:", dim(meta)[1], "x", dim(meta)[2])
dir.create(prefix)
write.csv(mat, paste(prefix, "mat.csv", sep=""))
write.csv(meta, paste(prefix, "meta.csv", sep=""))


# Vars
prefix = "data/BhaduriNature2020Org/"
mat = org.data2
meta = org.meta2

# Write Data
paste("Mat Dim:", dim(mat)[1], "x", dim(mat)[2], "  Meta Dim:", dim(meta)[1], "x", dim(meta)[2])
dir.create(prefix)
write.csv(mat, paste(prefix, "mat.csv", sep=""))
write.csv(meta, paste(prefix, "meta.csv", sep=""))

mat
