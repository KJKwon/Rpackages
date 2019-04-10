library(ggplot2)
library(extrafont)
se_list <- read.csv("selected_GO_total_gene_list_entrez.clean.csv",sep = '\t', head = FALSE, row.names = 1)
FC_tbl <- read.table("DAY8_vs_TPP-GA_DAY8.edgeR_out",sep = '\t', head = TRUE,row.names = NULL)
FC_tbl <- FC_tbl[FC_tbl$FDR < 0.05, ]
names(FC_tbl)[1] = "gene_name"
FC_tbl = FC_tbl[,1:2]
uniq_test <- !duplicated(regmatches(FC_tbl[,1],regexpr(('[A-Za-z0-9]+'),FC_tbl[,1])))
FC_tbl = FC_tbl[uniq_test,]
FC_tbl[,1] = regmatches(FC_tbl[,1],regexpr(('[A-Za-z0-9]+'),FC_tbl[,1]))
GO_term = rownames(se_list[12,])
gene_list = as.vector(se_list[12,][se_list[12,] != ""])
gene_location = c()
for (gene_nm in gene_list){
  gene_location = c(gene_location, which(tolower(FC_tbl$gene_name) %in% tolower(gene_nm)))
  #For family or subtype
  gene_location = c(gene_location, which(startsWith(tolower(FC_tbl$gene_name), tolower(gene_nm))))
}
gene_location = unique(gene_location)
FC_tbl_ready = FC_tbl[gene_location,]
FC_tbl_ready = FC_tbl_ready[abs(FC_tbl_ready$logFC) >= 1, ]
#Even though name is similar, totally different type of gene may exist
FC_tbl_ready = FC_tbl_ready[!(FC_tbl_ready$gene_name %in% c("Aspn","Aspm","Il13ra1","Il17rd","Il11ra1","Il15ra")),]
FC_tbl_ready$value = ifelse(FC_tbl_ready$logFC >= 0, "positive", "negative")
FC_tbl_ready$gene_name = factor(FC_tbl_ready$gene_name, levels = FC_tbl_ready$gene_name[order(FC_tbl_ready$logFC)])
##Two color mode
bar_pl = ggplot(FC_tbl_ready, aes(gene_name,logFC, label = gene_name, fill = factor(value))) + geom_text(aes(y = ifelse(logFC > 0, logFC+0.2, logFC-0.2)), size = 5, family = "Arial")
bar_pl = bar_pl + geom_col() + theme_bw() + coord_flip() + labs(x = "") + theme(axis.text.y = element_blank(), axis.text = element_text(size = 14),
                                                                      axis.title = element_text(size = 16, face = "bold"), legend.position = 'none') + 
  labs(title = GO_term, y = 'log2 fold change')

##One color mode
#bar_pl = ggplot(FC_tbl_ready, aes(gene_name,logFC, label = gene_name )) + geom_text(aes(y = ifelse(logFC > 0, logFC+0.1, logFC-0.1)), size = 5, family = "Arial")
#bar_pl = bar_pl + geom_col(fill = "#00BFC4") + theme_bw() + coord_flip() + labs(x = "") + theme(axis.text.y = element_blank(), axis.text = element_text(size = 14),
#                                                                      axis.title = element_text(size = 16, face = "bold"), legend.position = 'none') + 
#  labs(title = GO_term, y = 'log2 fold change')

pdf(GO_term, height = 8, width = 11, family = "Arial")
op = par(mar = c(5, 4, 0.05, 0.05) + 0.1)
plot(bar_pl)
par(op)
dev.off()
