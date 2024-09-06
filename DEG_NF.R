# modified by 2021.8.12
setwd("/DATA4T/Brassica_napus_RNA_Seq_published/")
library(tidyverse)
library(ggrain)
library(ggpubr)
library(EnvStats)


## reads count load
fc_count_dir <-"3featurecount/"
fc_count_files <- grep(dir(fc_count_dir,full.names = TRUE),pattern = ".summary",invert = TRUE,value = TRUE)

fc<-c()
for (file in fc_count_files){
  temp_file <- read.table(file,header = TRUE,stringsAsFactors = FALSE)[]
  fc<-cbind(fc,temp_file[,7])
  }

rownames(fc)<-temp_file$Geneid
colnames(fc) <- basename(fc_count_files)
sample_info  <- read.csv("/DATA4T/Brassica_napus_RNA_Seq_published/sample_info.csv",header = F)

colnames(fc)<-sample_info$V2[match(colnames(fc),sample_info$V1)]
fc <- fc[,grep(pattern = "17",colnames(fc))]
## gene length

## pca
library(mixOmics)
library(edgeR)
library(stringr)
MyResult.pca <- pca(t(log(fc[rowSums(cpm(fc)>=1)>=3,]+1,2)))
#MyResult.pca <- pca(sub_data,scale = TRUE)
group <- str_split_fixed(colnames(fc),pattern = "_",n = 2)[,1]
plotIndiv(MyResult.pca,group = group,ellipse = TRUE,ellipse.level = 0.6,title = "Samples",size.axis = 20,size.xlabel = 20,size.ylabel = 20)


##DEG analysis
library(DESeq2)
library(edgeR)
library(stringr)

fc <- fc[,grep(pattern = "NF17",colnames(fc))]
group_list<-factor(str_split_fixed(colnames(fc),pattern = "_",n = 3)[,2]) # make group list
coldata<-data.frame(row.names = colnames(fc),condition=group_list)  # make coldata

expdat<-round(fc,0) # make round data 30879
expdat = expdat[rowSums(cpm(expdat)>=1) >= 3,]

dds<-DESeqDataSetFromMatrix(countData = expdat,colData = coldata,design = ~ condition) # make DEseq dataset
 # filter
dds<-DESeq(dds) # normalization

res<-results(dds,contrast = c("condition","t1","t0"))   ## tumor/normal
res2 <-results(dds,contrast = c("condition","t2","t0"))
res3 <- results(dds,contrast = c("condition","t2","t1"))

resOrdered <- res[order(res$pvalue),]
DEG <- as.data.frame(resOrdered)
DEG <- na.omit(DEG)

logFC_cutoff<-1      ##  mean + 2sd 
DEG$change<-as.factor(ifelse(DEG$padj<0.01&abs(DEG$log2FoldChange)>logFC_cutoff,ifelse(DEG$log2FoldChange>logFC_cutoff,"UP","DOWN"),"NOT")) # label the change
save(DEG,file = "7NF_tmp_deg//DEGseq2_t1vs_t0.Rdata")


resOrdered2 <- res2[order(res2$pvalue),]
DEG2 <- as.data.frame(resOrdered2)
DEG2 <- na.omit(DEG2)

logFC_cutoff<-1      ##  mean + 2sd 
DEG2$change<-as.factor(ifelse(DEG2$padj<0.01&abs(DEG2$log2FoldChange)>logFC_cutoff,ifelse(DEG2$log2FoldChange>logFC_cutoff,"UP","DOWN"),"NOT")) # label the change
save(DEG2,file = "7NF_tmp_deg//DEGseq2_t2vs_t0.Rdata")

########## t2 to t1
resOrdered3 <- res3[order(res3$pvalue),]
DEG3 <- as.data.frame(resOrdered3)
DEG3 <- na.omit(DEG3)

logFC_cutoff<-1      ##  mean + 2sd 
DEG3$change<-as.factor(ifelse(DEG3$padj<0.01&abs(DEG3$log2FoldChange)>logFC_cutoff,ifelse(DEG3$log2FoldChange>logFC_cutoff,"UP","DOWN"),"NOT")) # label the change
save(DEG3,file = "7NF_tmp_deg//DEGseq2_t2vs_t1.Rdata")



#load("3.5tmp_deg/DEGseq2.Rdata")
## plot vacano plot
library(ggplot2)
require("ggrepel")


df <- rbind(DEG[DEG$change!="NOT",][1:10,],DEG[DEG$change=="DOWN",][1:10,])

this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up-regulated genes is ',nrow(DEG[DEG$change =='UP',]) ,
                     '\nThe number of down-regulated genes is ',nrow(DEG[DEG$change =='DOWN',]))
p1 <- ggplot(data=DEG,aes(x=log2FoldChange,y=-log10(padj),color=change))+
  geom_point(alpha=0.8,size=1.75)+
  geom_hline(yintercept = -log10(0.01),color="#990000",linetype="dashed")+
  geom_vline(xintercept = -1,color="#990000",linetype="dashed")+
  geom_vline(xintercept = 1,color="#990000",linetype="dashed")+
  labs(x="log2 fold change")+ 
  ylab("-log10 pvalue")+
  ggtitle(this_title)+
  theme_bw(base_size = 15)+
  theme(plot.title = element_text(size=15,hjust=0.5))+
  scale_color_manual(values=c('orange','gray','steelblue'))+
  
  geom_label_repel(
    data = df,
    aes(label = rownames(df),
        fill=as.factor(change)),
    color="white",
    show.legend = FALSE ,
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = 'grey50',
    max.overlaps = 15
  )+
  scale_fill_manual(values=c("orange","steelblue"))

## get DEG list from degseq2
library(tidyverse)
DEG_out <- DEG%>%arrange(desc(log2FoldChange),padj,change)%>%filter(abs(log2FoldChange)>1,padj<0.01)
DEG_out2 <- DEG2%>%arrange(desc(log2FoldChange),padj,change)%>%filter(abs(log2FoldChange)>1,padj<0.01)
DEG_out3 <- DEG3%>%arrange(desc(log2FoldChange),padj,change)%>%filter(abs(log2FoldChange)>1,padj<0.01)

#### 
## edgeR analysis
library(edgeR)
## make deglist
pheatmap::pheatmap(cor(log(cpm(fc)+1,2)))
expdat = fc[rowSums(cpm(fc)>=1) >= 3,]  # filter
exprSet <- DGEList(counts = expdat, group = group_list) # make dgelist
exprSet <- calcNormFactors(exprSet)  # normlaized factor

exprSet <- estimateCommonDisp(exprSet) # 
exprSet <- estimateTagwiseDisp(exprSet)

et <- exactTest(exprSet,pair = c("t0","t1"))
et2 <- exactTest(exprSet,pair = c("t0","t2"))
et3 <- exactTest(exprSet,pair = c("t1","t2"))

tTag <- topTags(et, n=nrow(exprSet))
tTag <- as.data.frame(tTag)

tTag2 <- topTags(et2, n=nrow(exprSet))
tTag2 <- as.data.frame(tTag2)

tTag3 <- topTags(et3, n=nrow(exprSet))
tTag3 <- as.data.frame(tTag3)


#logFC_cutoff<-with(tTag,mean(abs(logFC))+2*sd(abs(logFC)))       ##  mean + 2sd 
logFC_cutoff<-1
tTag$change<-as.factor(ifelse(tTag$FDR<0.01&abs(tTag$logFC)>logFC_cutoff,
                              ifelse(tTag$logFC>logFC_cutoff,"UP","DOWN"),"NOT"))

tTag2$change<-as.factor(ifelse(tTag2$FDR<0.01&abs(tTag2$logFC)>logFC_cutoff,
                              ifelse(tTag2$logFC>logFC_cutoff,"UP","DOWN"),"NOT"))

tTag3$change<-as.factor(ifelse(tTag3$FDR<0.01&abs(tTag3$logFC)>logFC_cutoff,
                               ifelse(tTag3$logFC>logFC_cutoff,"UP","DOWN"),"NOT"))

save(tTag,file = "7NF_tmp_deg//EdgeR_t1_vs_t0.Rdata")
save(tTag2,file = "7NF_tmp_deg//EdgeR_t2_vs_t0.Rdata")
save(tTag3,file = "7NF_tmp_deg//EdgeR_t2_vs_t1.Rdata")



df <- rbind(tTag[tTag$change!="NOT",][1:10,],tTag[tTag$change=="DOWN",][1:10,])

this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up-regulated gene is ',nrow(tTag[tTag$change =='UP',]) ,
                     '\nThe number of down-regulated gene is ',nrow(tTag[tTag$change =='DOWN',]))

p2<-ggplot(data=tTag,aes(x=logFC,y=-log10(FDR),color=change))+
  geom_point(alpha=0.8,size=1.75)+
  geom_hline(yintercept = -log10(0.01),color="#990000",linetype="dashed")+
  geom_vline(xintercept = -1,color="#990000",linetype="dashed")+
  geom_vline(xintercept = 1,color="#990000",linetype="dashed")+
  labs(x="log2 fold change")+ 
  ylab("-log10 pvalue")+
  theme_bw(base_size = 25)+
  theme(plot.title = element_text(size=15,hjust=0.5),legend.position = c(0.2,0.8))+
  scale_color_manual(values=c('orange','gray','steelblue'))+
  
  geom_label_repel(
    data = df,
    aes(label = rownames(df),
        fill=as.factor(change)),
    color="white",
    show.legend = FALSE ,
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = 'grey50',
    max.overlaps = 15
  )+
  scale_fill_manual(values=c("orange","steelblue"))

## plot vacano for both methods
library(ggpubr)
ggarrange(p1,p2,labels = "AUTO",common.legend = TRUE)

### get differentially expressed genes  from edgeR
library(tidyverse)
tag_out <- tTag%>%arrange(desc(logFC),FDR,change)%>%filter(abs(logFC)>1,FDR<0.01)
tag_out2 <- tTag2%>%arrange(desc(logFC),FDR,change)%>%filter(abs(logFC)>1,FDR<0.01)
tag_out3 <- tTag3%>%arrange(desc(logFC),FDR,change)%>%filter(abs(logFC)>1,FDR<0.01)


## ggveen plot
edgeR_up  <- rownames(tTag)[tTag$logFC>1&tTag$FDR<0.01]
edgeR_down <- rownames(tTag)[tTag$logFC< -1&tTag$FDR<0.01]
  
DEGseq2_up <- rownames(DEG)[DEG$log2FoldChange>1&DEG$padj<0.01]
DEGseq2_down <- rownames(DEG)[DEG$log2FoldChange< -1&DEG$padj<0.01]

library(ggVennDiagram)
p1<- ggVennDiagram(list(edgeR_up=edgeR_up,DEGseq2_up=DEGseq2_up),
                   label = "count",
                   edge_size = 0)+
  scale_fill_gradient(low="steelblue",high = "orange")+
  ggtitle("T1 vs. T0")+
  theme(plot.title = element_text(size=12,hjust=0.5))
p2<-ggVennDiagram(list(edgeR_down=edgeR_down,DEGseq2_down=DEGseq2_down),
                  label = "count",
                  edge_size = 0)+
  scale_fill_gradient(low="steelblue",high = "orange")+
  ggtitle("T1 vs. T0")+
  theme(plot.title = element_text(size=12,hjust=0.5))



############ t2to t0
edgeR_up2  <- rownames(tTag2)[tTag2$logFC>1&tTag2$FDR<0.01]
edgeR_down2 <- rownames(tTag2)[tTag2$logFC< -1&tTag$FDR<0.01]

DEGseq2_up2 <- rownames(DEG2)[DEG2$log2FoldChange>1&DEG2$padj<0.01]
DEGseq2_down2 <- rownames(DEG2)[DEG2$log2FoldChange< -1&DEG2$padj<0.01]

p3<- ggVennDiagram(list(edgeR_up=edgeR_up2,DEGseq2_up=DEGseq2_up2),
                   label = "count",
                   edge_size = 0)+
  scale_fill_gradient(low="steelblue",high = "orange")+
  ggtitle("T2 vs. T0")+
  theme(plot.title = element_text(size=12,hjust=0.5))
p4<-ggVennDiagram(list(edgeR_down=edgeR_down2,DEGseq2_down=DEGseq2_down2),
                  label = "count",
                  edge_size = 0)+
  scale_fill_gradient(low="steelblue",high = "orange")+
  ggtitle("T2 vs. T0")+
  theme(plot.title = element_text(size=12,hjust=0.5))



############ t2to t1
edgeR_up3  <- rownames(tTag3)[tTag3$logFC>1&tTag3$FDR<0.01]
edgeR_down3 <- rownames(tTag3)[tTag3$logFC< -1&tTag3$FDR<0.01]

DEGseq2_up3 <- rownames(DEG3)[DEG3$log2FoldChange>1&DEG3$padj<0.01]
DEGseq2_down3 <- rownames(DEG3)[DEG3$log2FoldChange< -1&DEG3$padj<0.01]

p5<- ggVennDiagram(list(edgeR_up=edgeR_up3,DEGseq2_up=DEGseq2_up3),
                   label = "count",
                   edge_size = 0)+
  scale_fill_gradient(low="steelblue",high = "orange")+
  ggtitle("T2 vs. T1")+
  theme(plot.title = element_text(size=12,hjust=0.5))
p6<-ggVennDiagram(list(edgeR_down=edgeR_down3,DEGseq2_down=DEGseq2_down3),
                  label = "count",
                  edge_size = 0)+
  scale_fill_gradient(low="steelblue",high = "orange")+
  ggtitle("T2 vs. T1")+
  theme(plot.title = element_text(size=12,hjust=0.5))

ggarrange(p1,p2,p3,p4,p5,p6,nrow = 3,ncol = 2)


##########
#common DEG
##########
up_reg <- tTag[intersect(edgeR_up,DEGseq2_up),]%>%tibble::rownames_to_column()%>%arrange(desc(logFC),FDR,change)%>%filter(abs(logFC)>1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
down_reg <- tTag[intersect(edgeR_down,DEGseq2_down),]%>%tibble::rownames_to_column()%>%arrange(logFC,FDR,change)%>%filter(abs(logFC)>1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
degs_out <- rbind(up_reg,down_reg) # 2415 DEGs
anno_data<-read.delim("../Brassica_napus_assembly/assembly/result_202301300547/Annotation/06_release/04_Gene_Function_Annotation/Total/Annotation_Summary.xls",header = T,sep="\t")
anno_data_deg <- anno_data[match(rownames(degs_out),anno_data$gene_id),] 
deg_out_df <- cbind(degs_out,anno_data_deg)

write_excel_csv(deg_out_df,"7NF_tmp_deg/t1_t0.deg.csv")


## T2 to T0
up_reg2 <- tTag2[intersect(edgeR_up2,DEGseq2_up2),]%>%tibble::rownames_to_column()%>%arrange(desc(logFC),FDR,change)%>%filter(abs(logFC)>1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
down_reg2 <- tTag2[intersect(edgeR_down2,DEGseq2_down2),]%>%tibble::rownames_to_column()%>%arrange(logFC,FDR,change)%>%filter(abs(logFC)>1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
degs_out2 <- rbind(up_reg2,down_reg2) # 2415 DEGs
anno_data<-read.delim("../Brassica_napus_assembly/assembly/result_202301300547/Annotation/06_release/04_Gene_Function_Annotation/Total/Annotation_Summary.xls",header = T,sep="\t")
anno_data_deg2 <- anno_data[match(rownames(degs_out2),anno_data$gene_id),] 
deg_out_df2 <- cbind(degs_out2,anno_data_deg2)
write_excel_csv(deg_out_df2,"7NF_tmp_deg/t2_t0.deg.csv")


## T2 to T1
up_reg3 <- tTag3[intersect(edgeR_up3,DEGseq2_up3),]%>%tibble::rownames_to_column()%>%arrange(desc(logFC),FDR,change)%>%filter(abs(logFC)>1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
down_reg3 <- tTag3[intersect(edgeR_down3,DEGseq2_down3),]%>%tibble::rownames_to_column()%>%arrange(logFC,FDR,change)%>%filter(abs(logFC)>1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
degs_out3 <- rbind(up_reg3,down_reg3) # 2415 DEGs
anno_data<-read.delim("../Brassica_napus_assembly/assembly/result_202301300547/Annotation/06_release/04_Gene_Function_Annotation/Total/Annotation_Summary.xls",header = T,sep="\t")
anno_data_deg3 <- anno_data[match(rownames(degs_out3),anno_data$gene_id),] 
deg_out_df3 <- cbind(degs_out3,anno_data_deg3)
write_excel_csv(deg_out_df3,"7NF_tmp_deg/t2_t1.deg.csv")


## flow chart
library(ggsankey)

df <- data.frame(CK=row.names(expdat),stringsAsFactors = F)
df$T1 <- as.character(deg_out_df[df$CK,]$change)
df$T1[is.na(df$T1)]<-"Unchanged" # UO

df$T2 <- as.character(deg_out_df2[df$CK,]$change)
df$T2[is.na(df$T2)]<-"Unchanged"
df$CK <-"All"

df <- df %>% make_long(CK,T1,T2)


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


## Compare AA and CC gene expression
library(stringr)
gene_exp <- read.csv("../Brassica_napus_RNA_Seq/all_FPKM.csv",row.names = 1)
group <- str_split_fixed(colnames(gene_exp),pattern = "_",n = 3)[,c(1,2)]
group <-apply(group,1,function(x)paste(x[1],x[2],sep = "_"))

for (m in group){
  gene_exp[,m] <- rowMeans(gene_exp[,grep(pattern = m,colnames(gene_exp))])
}

gene_exp <- gene_exp[!startsWith(rownames(gene_exp),prefix = "Contig"),]
AA_exp <- subset(gene_exp,startsWith(rownames(gene_exp),prefix = "A"))
CC_exp <- subset(gene_exp,startsWith(rownames(gene_exp),prefix = "C"))

df1 <- data.frame(expression=log(AA_exp$NS17_t0+1,2),group="AA")
df2 <- data.frame(expression=log(CC_exp$NS17_t0+1,2),group="CC")
data = rbind(df1,df2)

p2 <- ggplot(data, aes(group, expression, fill = group)) +
  geom_rain(alpha = .5,) +
  theme_pubr(base_size = 18) +
  scale_fill_brewer(palette = 'Dark2') +
  guides(fill = 'none', color = 'none')+
  ggsignif::geom_signif(
    comparisons = list(c("AA", "CC")),
    map_signif_level = T)+
  stat_n_text()



p<-ggplot(data,aes(x=group,y=expression,col=group))+
  geom_violin(trim = FALSE,aes(fill=group),alpha=0.5,col="white")+
  geom_boxplot(width=0.05)+
  scale_fill_aaas()+
  scale_color_aaas()+
  geom_signif(comparisons = list(c("AA","CC")),
              test="wilcox.test", test.args=list(alternative="greater"),
              step_increase = 0.05,tip_length = 0.01)+
  theme_bw(base_size = 20)+
  scale_x_discrete(labels=c("AA","CC"))+
  theme(legend.position="none",axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5,size=20,face = "bold"))+ # title posistion
  ylab("Gene expression (log2FPKM)")+
  stat_n_text()
