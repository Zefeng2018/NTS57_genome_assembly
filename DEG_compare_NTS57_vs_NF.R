setwd("/DATA4T/Brassica_napus_RNA_Seq_published/")
load("6tmp_deg/EdgeR_t1_vs_t0.Rdata")
NTS57_T1 <- tTag
rm(tTag) # 46040

load("7NF_tmp_deg/EdgeR_t1_vs_t0.Rdata")
NF_T1 <- tTag
rm(tTag) # 44274


load("6tmp_deg/EdgeR_t2_vs_t0.Rdata")
NTS57_T2 <- tTag2 # 46040
rm(tTag2)

load("7NF_tmp_deg/EdgeR_t2_vs_t0.Rdata")
NF_T2 <- tTag2  # 44274


#cm <- intersect(rownames(NTS57_T1),rownames(NF_T1))
cm <- unique(c(rownames(NTS57_T1),rownames(NF_T1))) # 48584
## compare 
df <- data.frame(CK=cm,stringsAsFactors = F)
df$NTS57_T1 <- as.character(NTS57_T1[df$CK,]$change)
df$NTS57_T1[is.na(df$NTS57_T1)] <- "NOT"

df$NTS57_T2 <- as.character(NTS57_T2[df$CK,]$change)
df$NTS57_T2[is.na(df$NTS57_T2)] <- "NOT"

df$NF_T1 <- as.character(NF_T1[df$CK,]$change)
df$NF_T1[is.na(df$NF_T1)] <- "NOT"

df$NF_T2 <- as.character(NF_T2[df$CK,]$change)
df$NF_T2[is.na(df$NF_T2)] <- "NOT"

df$CK <-"All"


library(ggsankey)
library(ggplot2)
df <- df %>% make_long(NTS57_T2,NTS57_T1,CK,NF_T1,NF_T2)

ggplot(df, aes(x = x, next_x = next_x, node = node, 
               next_node = next_node, 
               fill = factor(node), label = node)) +
  geom_alluvial(flow.alpha = .8,width = 0.2) +
  geom_alluvial_text(size = 4, color = "black") +
  scale_fill_manual(values = c("steelblue", "tomato", "gray", "orange"))+
  theme_alluvial(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5),
  ) +
  ggtitle("Expression change")


### specific genes
t1_cm <- intersect(rownames(NTS57_T1),rownames(NF_T1)) # 41730
plot(NTS57_T1[t1_cm,]$logFC,NF_T1[t1_cm,]$logFC)
df1 <- data.frame(GN=t1_cm,NTS57=NTS57_T1[t1_cm,]$logFC,NF=NF_T1[t1_cm,]$logFC)
library(ggpubr)
p1 <-ggscatter(data = df1,
          y = "NTS57",
          x = "NF",
          rug = F,
          ellipse = TRUE,
          ellipse.border.remove = F,
          alpha=0.5,color  = "steelblue",
          add = "reg.line",
          add.params = list(color = "brown", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"))+
  
  geom_vline(xintercept = c(-1,1),linetype = "dashed",color="orange")+
  geom_hline(yintercept = c(-1,1),linetype = "dashed",color="orange")+
  xlab("log(Fold change) in NF")+
  ylab("log(Fold change) in NTS57")


t2_cm <- intersect(rownames(NTS57_T2),rownames(NF_T2)) # 41730
plot(NTS57_T2[t2_cm,]$logFC,NF_T2[t2_cm,]$logFC)
df2 <- data.frame(GN=t2_cm,NTS57=NTS57_T2[t2_cm,]$logFC,NF=NF_T2[t2_cm,]$logFC)
library(ggpubr)
p2 <- ggscatter(data = df2,
          y = "NTS57",
          x = "NF",
          rug = F,
          ellipse = TRUE,
          ellipse.border.remove = F,
          alpha=0.5,color  = "steelblue",
          add = "reg.line",
          add.params = list(color = "brown", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"))+
  
  geom_vline(xintercept = c(-1,1),linetype = "dashed",color="orange")+
  geom_hline(yintercept = c(-1,1),linetype = "dashed",color="orange")+
  xlab("log(Fold change) in NF")+
  ylab("log(Fold change) in NTS57")

ggarrange(p1,p2,nrow = 1)

### 
t1_NTS57_spcific <- df1$GN[abs(df1$NTS57)>1&abs(df1$NF)<1]
write.table(t1_NTS57_spcific,file = "8NTS57_specific/NTS57_spec_t1.txt",quote = FALSE,row.names = F,col.names = F)
t2_NTS57_spcific <- df2$GN[abs(df2$NTS57)>1&abs(df2$NF)<1]
write.table(t2_NTS57_spcific,file = "8NTS57_specific/NTS57_spec_t2.txt",quote = FALSE,row.names = F,col.names = F)



### GO enrichment
genes_nts57_t1 <- df1$GN[abs(df1$NTS57)>1&abs(df1$NF)<1]
go_enrich_function <- function(genes){
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
  
  CC_enrich@result$Category <- "CC"
  MF_enrich@result$Category <- "MF"
  BP_enrich@result$Category <- "BP"
  
  GO_enrich <- rbind(CC_enrich@result[1:nrow(CC_enrich),],MF_enrich@result[1:nrow(MF_enrich),],BP_enrich@result[1:nrow(BP_enrich),])
  
  
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
  return(GO_enrich)
}
GO_enrich1 <-go_enrich_function(genes_nts57_t1)
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


genes_nts57_t2 <- df2$GN[df2$NTS57>1&abs(df2$NF)<1]
GO_enrich2 <-go_enrich_function(genes_nts57_t2)
p2 <- ggplot(data = GO_enrich2, aes(y = reorder(Description,-p.adjust),x=Count,fill = -log10(p.adjust))) +
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

##################
genes_NF_t1 <- read.csv("/DATA4T/Brassica_napus_RNA_Seq_published/7NF_tmp_deg/t1_t0.deg.csv")
genes_NF_t1 <- genes_NF_t1$gene_id
GO_enrich <-go_enrich_function(genes_NF_t1)
p3 <- ggplot(data = GO_enrich, aes(y = reorder(Description,-p.adjust),x=Count,fill = -log10(p.adjust))) +
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

genes_NF_t2 <- read.csv("/DATA4T/Brassica_napus_RNA_Seq_published/7NF_tmp_deg/t2_t0.deg.csv")
genes_NF_t2 <- genes_NF_t2$gene_id
GO_enrich <-go_enrich_function(genes_NF_t2)
p4 <- ggplot(data = GO_enrich, aes(y = reorder(Description,-p.adjust),x=Count,fill = -log10(p.adjust))) +
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

ggarrange(p3,p4)


## gene expression for tube ralated genes
T1_tube_genes <- GO_enrich1$geneID[GO_enrich1$Description%in%c("microtubule","microtubule-based process")]
T1_tube_genes <- unlist(str_split(T1_tube_genes,pattern = "/"))

T2_tube_genes <- GO_enrich2$geneID[GO_enrich2$Description%in%c("microtubule","microtubule-based process")]
T2_tube_genes <- unlist(str_split(T2_tube_genes,pattern = "/"))

tube_genes <- unique(c(T1_tube_genes,T2_tube_genes))

df <- data.frame(genes = tube_genes)
df$NF_T1 <- NF_T1[df$genes,]$logFC
df$NF_T2 <- NF_T2[df$genes,]$logFC
df$NTS57_T1 <- NTS57_T1[df$genes,]$logFC
df$NTS57_T2 <- NTS57_T2[df$genes,]$logFC

anno_data<-read.delim("../Brassica_napus_assembly/assembly/result_202301300547/Annotation/06_release/04_Gene_Function_Annotation/Total/Annotation_Summary.xls",header = T,sep="\t")

df$GN <- anno_data$Swissprot_Annotation[match(df$genes,anno_data$gene_id)]
df$GN2 <- apply(df,1,function(x)rev(unlist(str_split(x[6],pattern = " ")))[3])
df$GN2 <- str_replace(df$GN2,pattern = "GN=",replacement = "")
df <- df[,c(1,2,3,4,5,7)]
df[25,6] <-"dlc6"
df$genes <- str_split_fixed(df$genes,pattern = "\\.",n = 2)[,1]
rownames(df) <- paste(df$genes,df$GN2,sep = "(")
rownames(df) <- paste(rownames(df),")",sep="")
df <- df[,c(2,3,4,5)]
pheatmap::pheatmap(df,display_numbers = T)


### veendigram
library(ggVennDiagram)
NTS57_T1 <- read.csv("3.5tmp_deg/t1_t0.deg.csv")
NTS57_T2 <- read.csv("3.5tmp_deg/t2_t0.deg.csv")
NF_T1 <- read.csv("7NF_tmp_deg/t1_t0.deg.csv")
NF_T2 <- read.csv("7NF_tmp_deg/t2_t0.deg.csv")

ggVennDiagram(list(NTS57_T1$gene_id,NTS57_T2$gene_id,NF_T1$gene_id,NF_T2$gene_id))
