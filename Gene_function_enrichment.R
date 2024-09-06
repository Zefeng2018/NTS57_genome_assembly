library(clusterProfiler)
library(readxl)


## 导入GO 注释信息
setwd("/DATA4T/Brassica_napus_assembly/assembly/result_202301300547/Annotation/06_release/04_Gene_Function_Annotation/GO/")
data <-read.delim("youcai.go.classification.xls")
write.table(data,file = "youcai.go.classification.tab.txt",quote = FALSE,sep="\t",col.names = T,row.names = F)
# awk 'OFS="\t"{FS="\t";u=split($5,m,",");for(x=1;x<=length(m);x++)print $1,$2,$3,m[x]}' youcai.go.classification.tab.txt > GO_reform.txt

# 导入基因列表
genes <- read.csv("/DATA4T/Brassica_napus_RNA_Seq/3.5tmp_deg/t2_t0.deg.csv")
genes <- genes$gene_id

# 导入注释文件
go_info <- read.table("GO_reform.txt",header = T,sep="\t",quote = "")
CC <- go_info[go_info$namespace=="cellular_component",]
BP <- go_info[go_info$namespace=="biological_process",]
MF <- go_info[go_info$namespace=="molecular_function",]

term2gene <- CC[,c(1,4)]
term2name <- CC[,c(1,2)]
# 富集分析
CC_enrich <- enricher(genes,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)

term2gene <- MF[,c(1,4)]
term2name <- MF[,c(1,2)]
# 富集分析
MF_enrich <- enricher(genes,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)

term2gene <- BP[,c(1,4)]
term2name <- BP[,c(1,2)]
# 富集分析
BP_enrich <- enricher(genes,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)



library(ggpubr)
BP <- dotplot(BP_enrich, title = "BP",label_format=100)
CC <- dotplot(CC_enrich, title = "CC",label_format=100)
MF <- dotplot(MF_enrich, title = "MF",label_format=100)
ggarrange(BP, CC, MF, ncol = 1, nrow = 3, align = "hv")


## elegant  plot

CC_enrich@result$Category <- "CC"
MF_enrich@result$Category <- "MF"
BP_enrich@result$Category <- "BP"

GO_enrich <- rbind(CC_enrich@result[1:7,],MF_enrich@result[1:10,],BP_enrich@result[1:10,])


## plot
#先自定义主题：
mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  plot.title = element_text(size = 14,
                            hjust = 0.5,
                            face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5,
                       r = 10,
                       l = 5.5,
                       b = 5.5)
)
#常规富集条形图绘图：
p <- ggplot(data = GO_enrich, aes(x = Count, y = Description, fill = -log10(pvalue))) +
  scale_fill_distiller(palette = "RdPu",direction = 1) + #更改配色
  geom_bar(stat = "identity", width = 0.8) + #绘制条形图
  labs(x = "Number of Gene", y = "", title = "KEGG enrichment barplot") + #修改/添加各标题
  theme_bw() + mytheme+ #主题更改+
  facet_wrap(.~Category,nrow = 3,scales = "free")


mytheme2 <- mytheme + theme(axis.text.y = element_blank())


p1 <- ggplot(data = GO_enrich, aes(y = reorder(Description,-p.adjust),x=Count,fill = -log10(p.adjust))) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "Number of genes", y = "GO terms") +
  geom_text(aes(x = 0.03, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  mytheme2+
  facet_wrap(Category~.,scales = "free_y",nrow = 3,strip.position = "left")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12))
ggsave(filename = "/DATA4T/Brassica_napus_RNA_Seq/3.5tmp_deg/T2toT0_DEG.go.pdf",plot = p1,width=12, height=10)


#### kegg  analysis


