library(readr)
require(ggplot2)
require(reshape2)
require(gplots)
require(dendextend)
library(RColorBrewer)

DF = read.csv('Results/Heatmap/htmp_data_gt1Tmt_D5.csv')

colnames(DF)
new_colnames <- c( "id", "BDQ_5_0625","Cla_5_4",
                   "INH_5_125","Rif_5_125", "BDQ_5_125",
                   "Cla_5_1", "Emb_5_1","Cla_5_2",
                   "Rif_5_25","INH_5_25","Emb_5_25",
                   "Levo_5_125", "Strep_5_1","Linez_5_25",
                   "INH_5_5","Linez_5_0625", "BDQ_5_25",
                   "Levo_5_5", "Vanc_5_25", "Vanc_5_125",
                   "Strep_5_25", "Vanc_5_0625", "Emb_5_5",
                   "Levo_5_25", "Strep_5_5","Rif_5_0625",
                   "Linez_5_125")

colnames(DF) <- new_colnames
rownames(DF) <- DF[,c('id')]
DF <- DF[, names(DF) != "id"]
head(DF)

distance.row = dist(as.matrix(DF), method = "euclidean") # NAs introduced by coercion - change to zeros
cluster.row = hclust(distance.row, method = "ward.D") 

distance.col = as.dist((1-cor(DF))/2)
cluster.col = hclust(distance.col, method = "ward.D") 

desired_order_5 = c(
  "Cla_5_1", "Cla_5_2", "Cla_5_4",
  "Levo_5_125", "Levo_5_25", "Levo_5_5",
  "INH_5_125", "INH_5_25", "INH_5_5",
  "BDQ_5_0625", "BDQ_5_125", "BDQ_5_25",
  "Emb_5_25", "Emb_5_5", "Emb_5_1",
  "Linez_5_0625", "Linez_5_125", "Linez_5_25",
  "Strep_5_25", "Strep_5_5", "Strep_5_1",
  "Vanc_5_0625", "Vanc_5_125", "Vanc_5_25",
  "Rif_5_0625", "Rif_5_125", "Rif_5_25")


############ MAKE SURE USING RIGHT LIST ####################
ordered.cluster.col = rotate(cluster.col, desired_order_5)


############# MAKE SURE LABELLED CORRECTLY ################

path = 'Results/Heatmap/'
output = "htmp_Ward_D_gt1Tmt_D5.png"


png(paste(path, output, sep=''), width=9, height=10, units="in", res=2000)

h <- heatmap.2(as.matrix(DF),
          trace="none",
          margins = c(8, 8),
          srtCol=90,
          col=colorRampPalette(c('royalblue3', 'white', 'tomato3'))(100),
          Rowv=as.dendrogram(cluster.row),
          Colv=as.dendrogram(ordered.cluster.col),
          cexCol = 0.5,
          cexRow = 0.03,
          breaks=seq(-2, 2, length.out=101))

dev.off()
