library(Mfuzz)
total = read.table('Huebner201609_XENLAtxV2_DMZ.mean_p_rpkm.txt.humanized_with_XEN.txt',sep = '\t', quote = "")
total_XEN = total[,c(-1)]
require(Biobase)
tf.tbl.temp = new('ExpressionSet', exprs = as.matrix(total_XEN))
tf.tbl.s = standardise(tf.tbl.temp)
tf.mfuzz = mfuzz(tf.tbl.s,c = 16, m = 1.25)
mfuzzColorBar()
mfuzz.plot(tf.tbl.s,cl = tf.mfuzz, mfrow = c(4,4))
dev.copy(pdf,'Total_data_cluster.pdf')
dev.off()
cluster.output = as.data.frame(tf.mfuzz$cluster)
tf.tbl.out = cbind(total, cluster.output)
#tf.tbl.out = cbind(total_human, tf.tbl.out)
colnames(tf.tbl.out)[12] = 'group'

library(gplots)
library(RColorBrewer)
tbl = read.table('Cluster1_transcription_factor.txt', sep = '\t')
tbl = tbl[,c(-1)]
tbl_cor = cor(t(tbl))
heatmap(x=tbl_cor,sym = TRUE)
Mem.data = tf.mfuzz$membership
std.data = exprs(tf.tbl.s)
mean.data = as.data.frame(tf.mfuzz$centers)
mean.data = mean.data[order(-mean.data$DMZ110),]
mean.big2small = as.vector(rownames(mean.data))
confident_list = c()
cluster_color = c()
for (i in mean.big2small){
  count = length(rownames(Mem.data[Mem.data[,i] > 0.995,]))
  confident_list = c(confident_list,rownames(Mem.data[Mem.data[,i] > 0.995,]))
  print(count)
  cluster_color = c(cluster_color, rep(setcol[as.numeric(i)],count))  
}
std.ready = std.data[confident_list,]
std.ready = cbind(std.ready,data.frame(cluster_list))
std.ready.sort = std.ready[with(std.ready, order(cluster_list,-DMZ110)),]
setcol <- colorRampPalette(brewer.pal(9,'Set1'))(16)
heat.temp = heatmap.2(as.matrix(std.ready),main = "cluster_output", Colv = NA, Rowv = NA , dendrogram = 'none',
                      margins = c(6,10), trace = 'none', key = TRUE,lwid = c(1,5),lhei = c(1,6), scale = 'none',col = rev(brewer.pal(11, "RdBu")),
                      RowSideColors = cluster_color)
