library(data.table)
library(fgsea)
library(ggplot2)
ranks = read.table("DAY8_vs_TPP-GA_DAY8.entrez.rnk",header= TRUE, colClasses = c("character","numeric"))
ranks = ranks[!duplicated(ranks$geneID),]
ranks = setNames(ranks$log2FC,ranks$geneID)
pathways = gmtPathways("BP_mouse_gsea.c5.gmt")
set.seed(123)
fgseaRes <- fgsea(pathways = pathways, stats = ranks,minSize=15, maxSize=500,nperm=10000)
GSEA_pl = plotEnrichment(pathways[["GO_WHITE_FAT_CELL_DIFFERENTIATION"]], ranks) + labs(title="GO_GO_WHITE_FAT_CELL_DIFFERENTIATION") + 
  theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5), axis.title.x = element_text(size = 16), 
          axis.title.y = element_text(size = 16))
pdf(file = "GSEA_white_fat_cell_differentiation", height = 8, width = 11, family = "Arial")
op = par(mar = c(5, 4, 0.05, 0.05) + 0.1)
plot(GSEA_pl)
par(op)
dev.off()
#fwrite(fgseaRes[fgseaRes$padj < 0.05,], file="DAY2_vs_DAY8.BP_fgsea.txt", sep = '\t', sep2 = c(""," ",""))
