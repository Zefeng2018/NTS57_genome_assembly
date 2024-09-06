library(rtracklayer)
library(Biostrings)
library(BSgenome)
library(stringr)
library(GenomicFeatures)
library(methylKit)
library(ggpubr)

setwd("/DATA4T/Brassica_napus_methylation_published/")
# load genomic features
gtf <- makeTxDbFromGFF("0reference_genome/youcai.genome.gff")
genes <- genes(gtf)

genes<-genes[width(genes)>=100]
genes$gl <- width(genes)
genes <- genes[!startsWith(x = names(genes),prefix = "Contig")]
################################################################################
####################### 0. Methylkit data load (multiple samples analysis)===== 
file.names<- dir("4methylkit/",pattern = "CpG.txt",full.names = T)[5]
sample_id <- str_replace_all(basename(file.names),pattern = "_CpG.txt",replacement = "")
file.list <-as.list(file.names)

myobj=methRead(file.list,
               sample.id=as.list(sample_id),
               assembly="N.napus",
               treatment=c(rep(1,length(file.list))))
#############################################################################
###################### 1. dna methylation profiles without groups###############################
for (m in seq(1:length(myobj))){
  sample_name = sample_id[m]
  message(sample_name)
  chromsome_meth <- as(myobj[[m]],"GRanges")
  
  flank= 2000
  bins = 100
  date()
  
  ## classify genes into positive and negative
  ## 1 gene body 
  genes_pos <- genes[strand(genes)=="+"]
  genes_pos_tile <- tile(genes_pos,n = bins) 
  genes_pos_tile_unlist <- unlist(genes_pos_tile)
  genes_pos_tile_unlist$id <- rep(seq(1,bins),length(genes_pos))
  
  genes_neg <- genes[strand(genes)=="-"]
  genes_neg_tile <- tile(genes_neg,n = bins) 
  genes_neg_tile_unlist <- unlist(genes_neg_tile)
  genes_neg_tile_unlist$id <- rep(seq(bins,1),length(genes_neg)) # reverse bins
  
  genes_body_tile <- c(genes_pos_tile_unlist,genes_neg_tile_unlist)
  
  overlaps <- findOverlaps(chromsome_meth, genes_body_tile)           # find the overlapped methylated site in each bin 
  signal <- chromsome_meth[queryHits(overlaps)]                    
  sumSignal_Cs <- aggregate(signal$numCs, list(subjectHits(overlaps)), sum) ## very important: get methylated C
  sumSignal_Ts <- aggregate(signal$numTs, list(subjectHits(overlaps)), sum) ## very important: get unmethylated T
  
  genes_body_tile$coverage_Cs <- NA
  genes_body_tile$coverage_Ts <- NA
  
  genes_body_tile$coverage_Cs[sumSignal_Cs$Group.1] <- sumSignal_Cs$x
  genes_body_tile$coverage_Ts[sumSignal_Ts$Group.1] <- sumSignal_Ts$x
  
  genes_body_tile_coverage<- as.data.frame(mcols(genes_body_tile))
  genes_body_tile_coverage$gid <- str_split_fixed(rownames(genes_body_tile_coverage),"\\.",2)[,1]
  genes_body_tile_coverage<- genes_body_tile_coverage[!is.na(genes_body_tile_coverage$coverage_Ts),] # remove bins that not have no covered reads
  genes_body_tile_coverage$meth_level <- genes_body_tile_coverage$coverage_Cs/(genes_body_tile_coverage$coverage_Cs+genes_body_tile_coverage$coverage_Ts)
  
  genes_body_tile_density <-aggregate(genes_body_tile_coverage$meth_level,by=list(genes_body_tile_coverage$id),mean) # get the mean methylation of each bin
  message("Gene body methylation calculation is finished!")
  ## 2 upstream
  genes_pos_up <- flank(genes_pos,width = flank,both = FALSE)
  genes_pos_up_tile <- tile(genes_pos_up,n = bins) 
  genes_pos_up_tile_unlist <- unlist(genes_pos_up_tile)
  genes_pos_up_tile_unlist$id <- rep(seq(1,bins),length(genes_pos_up))
  
  genes_neg_up <- flank(genes_neg,width = flank,both = FALSE) #  
  genes_neg_up_tile <- tile(genes_neg_up,n =bins) 
  genes_neg_up_tile_unlist <- unlist(genes_neg_up_tile)
  genes_neg_up_tile_unlist$id <- rep(seq(bins,1),length(genes_neg_up))
  
  genes_up_tile <- c(genes_pos_up_tile_unlist,genes_neg_up_tile_unlist)
  
  
  overlaps <- findOverlaps(chromsome_meth, genes_up_tile)           # find the overlapped methylated site in each bin 
  signal <- chromsome_meth[queryHits(overlaps)]                    
  sumSignal_Cs <- aggregate(signal$numCs, list(subjectHits(overlaps)), sum) ## very important: get methylated C
  sumSignal_Ts <- aggregate(signal$numTs, list(subjectHits(overlaps)), sum) ## very important: get unmethylated T
  
  genes_up_tile$coverage_Cs <- NA
  genes_up_tile$coverage_Ts <- NA
  
  genes_up_tile$coverage_Cs[sumSignal_Cs$Group.1] <- sumSignal_Cs$x
  genes_up_tile$coverage_Ts[sumSignal_Ts$Group.1] <- sumSignal_Ts$x
  
  genes_up_tile_coverage<- as.data.frame(mcols(genes_up_tile))
  genes_up_tile_coverage$gid <- str_split_fixed(rownames(genes_up_tile_coverage),"\\.",2)[,1]
  genes_up_tile_coverage <- genes_up_tile_coverage[!is.na(genes_up_tile_coverage$coverage_Ts),] # remove bins that not have no covered reads
  genes_up_tile_coverage$meth_level <- genes_up_tile_coverage$coverage_Cs/(genes_up_tile_coverage$coverage_Cs+genes_up_tile_coverage$coverage_Ts)
  
  genes_up_tile_density <-aggregate(genes_up_tile_coverage$meth_level,by=list(genes_up_tile_coverage$id),mean)
  
  message("Gene upstream methylation calculation is finished!")
  ## 3. downstream 1kb
  
  genes_pos_down <- flank(genes_pos,width = flank,both = FALSE,start = FALSE)
  genes_pos_down_tile <- tile(genes_pos_down,n = bins) 
  genes_pos_down_tile_unlist <- unlist(genes_pos_down_tile)
  genes_pos_down_tile_unlist$id <- rep(seq(1,bins),length(genes_pos_down))
  
  genes_neg_down <- flank(genes_neg,width = flank,both = FALSE,start = FALSE)
  genes_neg_down_tile <- tile(genes_neg_down,n = bins) 
  genes_neg_down_tile_unlist <- unlist(genes_neg_down_tile)
  genes_neg_down_tile_unlist$id <- rep(seq(bins,1),length(genes_neg_down))
  
  genes_down_tile <- c(genes_pos_down_tile_unlist,genes_neg_down_tile_unlist)
  
  overlaps <- findOverlaps(chromsome_meth, genes_down_tile)           # find the overlapped methylated site in each bin 
  signal <- chromsome_meth[queryHits(overlaps)]                    
  sumSignal_Cs <- aggregate(signal$numCs, list(subjectHits(overlaps)), sum) ## very important: get methylated C
  sumSignal_Ts <- aggregate(signal$numTs, list(subjectHits(overlaps)), sum) ## very important: get unmethylated T
  
  genes_down_tile$coverage_Cs <- NA
  genes_down_tile$coverage_Ts <- NA
  
  genes_down_tile$coverage_Cs[sumSignal_Cs$Group.1] <- sumSignal_Cs$x
  genes_down_tile$coverage_Ts[sumSignal_Ts$Group.1] <- sumSignal_Ts$x
  
  genes_down_tile_coverage<- as.data.frame(mcols(genes_down_tile))
  genes_down_tile_coverage$gid <- str_split_fixed(rownames(genes_down_tile_coverage),"\\.",2)[,1]
  genes_down_tile_coverage <- genes_down_tile_coverage[!is.na(genes_down_tile_coverage$coverage_Ts),] # remove bins that not have no covered reads
  genes_down_tile_coverage$meth_level <- genes_down_tile_coverage$coverage_Cs/(genes_down_tile_coverage$coverage_Cs+genes_down_tile_coverage$coverage_Ts)
  
  genes_down_tile_density <-aggregate(genes_down_tile_coverage$meth_level,by=list(genes_down_tile_coverage$id),mean)
  message("Gene downstream methylation calculation is finished!")
  
  dd<-rbind(genes_up_tile_density,genes_body_tile_density,genes_down_tile_density)
  dd$n <- seq(1:300)
  dd$sample <- sample_name
  print(ggplot(dd,aes(x=n,y=x))+
    geom_smooth(span=0.1)+
    theme_pubr(base_size = 20)+
    xlab("")+
    ylab("CpG methylation level")+
    scale_x_continuous(limits=c(0,bins*3),
                       breaks = seq(0,bins*3,bins),
                       labels=c("-2 kb","TSS","TTS","2 kb")))
  write.table(dd,paste("results/gene_metaplot_chh/",sample_name,sep = ""),row.names = F,col.names = T,quote = F,sep="\t")  
  rm(chromsome_meth)
  gc()
}
################################################################################
####################### 2.class gene into two classes:A and C-sub genomes #####
for (m in seq(1:length(myobj))){
  message("Sample: ",m)
  sample_name = sample_id[m]
  message(sample_name)
  chromsome_meth <- as(myobj[[m]],"GRanges")
  
  flank= 2000
  bins = 100
  date()
  
  
  for (x in c("A","C")){
    message("Subgenome:",x)
    ## classify genes into positive and negative
    ## 1 gene body 
    genes_class <- genes[grepl(genes$gene_id,pattern = x)]
    genes_pos <- genes_class[strand(genes_class)=="+"]
    genes_pos_tile <- tile(genes_pos,n = bins) 
    genes_pos_tile_unlist <- unlist(genes_pos_tile)
    genes_pos_tile_unlist$id <- rep(seq(1,bins),length(genes_pos))
    
    genes_neg <- genes_class[strand(genes_class)=="-"]
    genes_neg_tile <- tile(genes_neg,n = bins) 
    genes_neg_tile_unlist <- unlist(genes_neg_tile)
    genes_neg_tile_unlist$id <- rep(seq(bins,1),length(genes_neg)) # reverse bins
    
    genes_body_tile <- c(genes_pos_tile_unlist,genes_neg_tile_unlist)
    
    overlaps <- findOverlaps(chromsome_meth, genes_body_tile)           # find the overlapped methylated site in each bin 
    signal <- chromsome_meth[queryHits(overlaps)]                    
    sumSignal_Cs <- aggregate(signal$numCs, list(subjectHits(overlaps)), sum) ## very important: get methylated C
    sumSignal_Ts <- aggregate(signal$numTs, list(subjectHits(overlaps)), sum) ## very important: get unmethylated T
    
    genes_body_tile$coverage_Cs <- NA
    genes_body_tile$coverage_Ts <- NA
    
    genes_body_tile$coverage_Cs[sumSignal_Cs$Group.1] <- sumSignal_Cs$x
    genes_body_tile$coverage_Ts[sumSignal_Ts$Group.1] <- sumSignal_Ts$x
    
    genes_body_tile_coverage<- as.data.frame(mcols(genes_body_tile))
    genes_body_tile_coverage$gid <- str_split_fixed(rownames(genes_body_tile_coverage),"\\.",2)[,1]
    genes_body_tile_coverage<- genes_body_tile_coverage[!is.na(genes_body_tile_coverage$coverage_Ts),] # remove bins that not have no covered reads
    genes_body_tile_coverage$meth_level <- genes_body_tile_coverage$coverage_Cs/(genes_body_tile_coverage$coverage_Cs+genes_body_tile_coverage$coverage_Ts)
    
    genes_body_tile_density <-aggregate(genes_body_tile_coverage$meth_level,by=list(genes_body_tile_coverage$id),mean) # get the mean methylation of each bin
    message("Gene body methylation calculation is finished!")
    ## 2 upstream
    genes_pos_up <- flank(genes_pos,width = flank,both = FALSE)
    genes_pos_up_tile <- tile(genes_pos_up,n = bins) 
    genes_pos_up_tile_unlist <- unlist(genes_pos_up_tile)
    genes_pos_up_tile_unlist$id <- rep(seq(1,bins),length(genes_pos_up))
    
    genes_neg_up <- flank(genes_neg,width = flank,both = FALSE) #  
    genes_neg_up_tile <- tile(genes_neg_up,n =bins) 
    genes_neg_up_tile_unlist <- unlist(genes_neg_up_tile)
    genes_neg_up_tile_unlist$id <- rep(seq(bins,1),length(genes_neg_up))
    
    genes_up_tile <- c(genes_pos_up_tile_unlist,genes_neg_up_tile_unlist)
    
    
    overlaps <- findOverlaps(chromsome_meth, genes_up_tile)           # find the overlapped methylated site in each bin 
    signal <- chromsome_meth[queryHits(overlaps)]                    
    sumSignal_Cs <- aggregate(signal$numCs, list(subjectHits(overlaps)), sum) ## very important: get methylated C
    sumSignal_Ts <- aggregate(signal$numTs, list(subjectHits(overlaps)), sum) ## very important: get unmethylated T
    
    genes_up_tile$coverage_Cs <- NA
    genes_up_tile$coverage_Ts <- NA
    
    genes_up_tile$coverage_Cs[sumSignal_Cs$Group.1] <- sumSignal_Cs$x
    genes_up_tile$coverage_Ts[sumSignal_Ts$Group.1] <- sumSignal_Ts$x
    
    genes_up_tile_coverage<- as.data.frame(mcols(genes_up_tile))
    genes_up_tile_coverage$gid <- str_split_fixed(rownames(genes_up_tile_coverage),"\\.",2)[,1]
    genes_up_tile_coverage <- genes_up_tile_coverage[!is.na(genes_up_tile_coverage$coverage_Ts),] # remove bins that not have no covered reads
    genes_up_tile_coverage$meth_level <- genes_up_tile_coverage$coverage_Cs/(genes_up_tile_coverage$coverage_Cs+genes_up_tile_coverage$coverage_Ts)
    
    genes_up_tile_density <-aggregate(genes_up_tile_coverage$meth_level,by=list(genes_up_tile_coverage$id),mean)
    
    message("Gene upstream methylation calculation is finished!")
    ## 3. downstream 1kb
    
    genes_pos_down <- flank(genes_pos,width = flank,both = FALSE,start = FALSE)
    genes_pos_down_tile <- tile(genes_pos_down,n = bins) 
    genes_pos_down_tile_unlist <- unlist(genes_pos_down_tile)
    genes_pos_down_tile_unlist$id <- rep(seq(1,bins),length(genes_pos_down))
    
    genes_neg_down <- flank(genes_neg,width = flank,both = FALSE,start = FALSE)
    genes_neg_down_tile <- tile(genes_neg_down,n = bins) 
    genes_neg_down_tile_unlist <- unlist(genes_neg_down_tile)
    genes_neg_down_tile_unlist$id <- rep(seq(bins,1),length(genes_neg_down))
    
    genes_down_tile <- c(genes_pos_down_tile_unlist,genes_neg_down_tile_unlist)
    
    overlaps <- findOverlaps(chromsome_meth, genes_down_tile)           # find the overlapped methylated site in each bin 
    signal <- chromsome_meth[queryHits(overlaps)]                    
    sumSignal_Cs <- aggregate(signal$numCs, list(subjectHits(overlaps)), sum) ## very important: get methylated C
    sumSignal_Ts <- aggregate(signal$numTs, list(subjectHits(overlaps)), sum) ## very important: get unmethylated T
    
    genes_down_tile$coverage_Cs <- NA
    genes_down_tile$coverage_Ts <- NA
    
    genes_down_tile$coverage_Cs[sumSignal_Cs$Group.1] <- sumSignal_Cs$x
    genes_down_tile$coverage_Ts[sumSignal_Ts$Group.1] <- sumSignal_Ts$x
    
    genes_down_tile_coverage<- as.data.frame(mcols(genes_down_tile))
    genes_down_tile_coverage$gid <- str_split_fixed(rownames(genes_down_tile_coverage),"\\.",2)[,1]
    genes_down_tile_coverage <- genes_down_tile_coverage[!is.na(genes_down_tile_coverage$coverage_Ts),] # remove bins that not have no covered reads
    genes_down_tile_coverage$meth_level <- genes_down_tile_coverage$coverage_Cs/(genes_down_tile_coverage$coverage_Cs+genes_down_tile_coverage$coverage_Ts)
    
    genes_down_tile_density <-aggregate(genes_down_tile_coverage$meth_level,by=list(genes_down_tile_coverage$id),mean)
    message("Gene downstream methylation calculation is finished!")
    
    dd<-rbind(genes_up_tile_density,genes_body_tile_density,genes_down_tile_density)
    dd$n <- seq(1:300)
    dd$sample <- sample_name
    print(ggplot(dd,aes(x=n,y=x))+
            geom_smooth(span=0.1)+
            theme_pubr(base_size = 20)+
            xlab("")+
            ylab("CpG methylation level")+
            scale_x_continuous(limits=c(0,bins*3),
                               breaks = seq(0,bins*3,bins),
                               labels=c("-2 kb","TSS","TTS","2 kb")))
    write.table(dd,paste("results/gene_metaplot/CHG/",sample_name,x, sep = "_"),row.names = F,col.names = T,quote = F,sep="\t")  
  }
  #rm(meth)
  #rm(chromsome_meth)
  gc()
}
################################################################################
### 3class gene into three classes and methylation profile
gene_exp <- read.csv("../Brassica_napus_RNA_Seq/all_FPKM.csv",row.names = 1)
group <- str_split_fixed(colnames(gene_exp),pattern = "_",n = 3)[,c(1,2)]
group <-apply(group,1,function(x)paste(x[1],x[2],sep = "_"))

for (m in group){
  gene_exp[,m] <- rowMeans(gene_exp[,grep(pattern = m,colnames(gene_exp))])
}
average_exp <- gene_exp[,tail(c(1:ncol(gene_exp)),5)]
average_exp<-average_exp %>% mutate(new_bin_NS17_to = ntile(NS17_t0, n=3))
#ggplot(average_exp,aes(x=as.factor(new_bin_NS17_to)),y=NS17_t0)+geom_boxplot() 
genes_low <- rownames(average_exp[average_exp$new_bin_NS17_to==1,])
genes_median <- rownames(average_exp[average_exp$new_bin_NS17_to==2,])
genes_high <- rownames(average_exp[average_exp$new_bin_NS17_to==3,])

genes_list<- list()
genes_list[[1]]<-genes_low
genes_list[[2]]<-genes_median
genes_list[[3]]<-genes_high

for (m in seq(1:length(myobj))){
  message("Sample: ",m)
  sample_name = sample_id[m]
  message(sample_name)
  chromsome_meth <- as(myobj[[m]],"GRanges")
  
  flank= 2000
  bins = 100
  date()
  
  
  for (x in seq(1,3)){
    message("Gene list:",x)
    ## classify genes into positive and negative
    ## 1 gene body 
    genes_class <- genes[genes$gene_id%in%genes_list[[x]]]
    genes_pos <- genes_class[strand(genes_class)=="+"]
    genes_pos_tile <- tile(genes_pos,n = bins) 
    genes_pos_tile_unlist <- unlist(genes_pos_tile)
    genes_pos_tile_unlist$id <- rep(seq(1,bins),length(genes_pos))
    
    genes_neg <- genes_class[strand(genes_class)=="-"]
    genes_neg_tile <- tile(genes_neg,n = bins) 
    genes_neg_tile_unlist <- unlist(genes_neg_tile)
    genes_neg_tile_unlist$id <- rep(seq(bins,1),length(genes_neg)) # reverse bins
    
    genes_body_tile <- c(genes_pos_tile_unlist,genes_neg_tile_unlist)
    
    overlaps <- findOverlaps(chromsome_meth, genes_body_tile)           # find the overlapped methylated site in each bin 
    signal <- chromsome_meth[queryHits(overlaps)]                    
    sumSignal_Cs <- aggregate(signal$numCs, list(subjectHits(overlaps)), sum) ## very important: get methylated C
    sumSignal_Ts <- aggregate(signal$numTs, list(subjectHits(overlaps)), sum) ## very important: get unmethylated T
    
    genes_body_tile$coverage_Cs <- NA
    genes_body_tile$coverage_Ts <- NA
    
    genes_body_tile$coverage_Cs[sumSignal_Cs$Group.1] <- sumSignal_Cs$x
    genes_body_tile$coverage_Ts[sumSignal_Ts$Group.1] <- sumSignal_Ts$x
    
    genes_body_tile_coverage<- as.data.frame(mcols(genes_body_tile))
    genes_body_tile_coverage$gid <- str_split_fixed(rownames(genes_body_tile_coverage),"\\.",2)[,1]
    genes_body_tile_coverage<- genes_body_tile_coverage[!is.na(genes_body_tile_coverage$coverage_Ts),] # remove bins that not have no covered reads
    genes_body_tile_coverage$meth_level <- genes_body_tile_coverage$coverage_Cs/(genes_body_tile_coverage$coverage_Cs+genes_body_tile_coverage$coverage_Ts)
    
    genes_body_tile_density <-aggregate(genes_body_tile_coverage$meth_level,by=list(genes_body_tile_coverage$id),mean) # get the mean methylation of each bin
    message("Gene body methylation calculation is finished!")
    ## 2 upstream
    genes_pos_up <- flank(genes_pos,width = flank,both = FALSE)
    genes_pos_up_tile <- tile(genes_pos_up,n = bins) 
    genes_pos_up_tile_unlist <- unlist(genes_pos_up_tile)
    genes_pos_up_tile_unlist$id <- rep(seq(1,bins),length(genes_pos_up))
    
    genes_neg_up <- flank(genes_neg,width = flank,both = FALSE) #  
    genes_neg_up_tile <- tile(genes_neg_up,n =bins) 
    genes_neg_up_tile_unlist <- unlist(genes_neg_up_tile)
    genes_neg_up_tile_unlist$id <- rep(seq(bins,1),length(genes_neg_up))
    
    genes_up_tile <- c(genes_pos_up_tile_unlist,genes_neg_up_tile_unlist)
    
    
    overlaps <- findOverlaps(chromsome_meth, genes_up_tile)           # find the overlapped methylated site in each bin 
    signal <- chromsome_meth[queryHits(overlaps)]                    
    sumSignal_Cs <- aggregate(signal$numCs, list(subjectHits(overlaps)), sum) ## very important: get methylated C
    sumSignal_Ts <- aggregate(signal$numTs, list(subjectHits(overlaps)), sum) ## very important: get unmethylated T
    
    genes_up_tile$coverage_Cs <- NA
    genes_up_tile$coverage_Ts <- NA
    
    genes_up_tile$coverage_Cs[sumSignal_Cs$Group.1] <- sumSignal_Cs$x
    genes_up_tile$coverage_Ts[sumSignal_Ts$Group.1] <- sumSignal_Ts$x
    
    genes_up_tile_coverage<- as.data.frame(mcols(genes_up_tile))
    genes_up_tile_coverage$gid <- str_split_fixed(rownames(genes_up_tile_coverage),"\\.",2)[,1]
    genes_up_tile_coverage <- genes_up_tile_coverage[!is.na(genes_up_tile_coverage$coverage_Ts),] # remove bins that not have no covered reads
    genes_up_tile_coverage$meth_level <- genes_up_tile_coverage$coverage_Cs/(genes_up_tile_coverage$coverage_Cs+genes_up_tile_coverage$coverage_Ts)
    
    genes_up_tile_density <-aggregate(genes_up_tile_coverage$meth_level,by=list(genes_up_tile_coverage$id),mean)
    
    message("Gene upstream methylation calculation is finished!")
    ## 3. downstream 1kb
    
    genes_pos_down <- flank(genes_pos,width = flank,both = FALSE,start = FALSE)
    genes_pos_down_tile <- tile(genes_pos_down,n = bins) 
    genes_pos_down_tile_unlist <- unlist(genes_pos_down_tile)
    genes_pos_down_tile_unlist$id <- rep(seq(1,bins),length(genes_pos_down))
    
    genes_neg_down <- flank(genes_neg,width = flank,both = FALSE,start = FALSE)
    genes_neg_down_tile <- tile(genes_neg_down,n = bins) 
    genes_neg_down_tile_unlist <- unlist(genes_neg_down_tile)
    genes_neg_down_tile_unlist$id <- rep(seq(bins,1),length(genes_neg_down))
    
    genes_down_tile <- c(genes_pos_down_tile_unlist,genes_neg_down_tile_unlist)
    
    overlaps <- findOverlaps(chromsome_meth, genes_down_tile)           # find the overlapped methylated site in each bin 
    signal <- chromsome_meth[queryHits(overlaps)]                    
    sumSignal_Cs <- aggregate(signal$numCs, list(subjectHits(overlaps)), sum) ## very important: get methylated C
    sumSignal_Ts <- aggregate(signal$numTs, list(subjectHits(overlaps)), sum) ## very important: get unmethylated T
    
    genes_down_tile$coverage_Cs <- NA
    genes_down_tile$coverage_Ts <- NA
    
    genes_down_tile$coverage_Cs[sumSignal_Cs$Group.1] <- sumSignal_Cs$x
    genes_down_tile$coverage_Ts[sumSignal_Ts$Group.1] <- sumSignal_Ts$x
    
    genes_down_tile_coverage<- as.data.frame(mcols(genes_down_tile))
    genes_down_tile_coverage$gid <- str_split_fixed(rownames(genes_down_tile_coverage),"\\.",2)[,1]
    genes_down_tile_coverage <- genes_down_tile_coverage[!is.na(genes_down_tile_coverage$coverage_Ts),] # remove bins that not have no covered reads
    genes_down_tile_coverage$meth_level <- genes_down_tile_coverage$coverage_Cs/(genes_down_tile_coverage$coverage_Cs+genes_down_tile_coverage$coverage_Ts)
    
    genes_down_tile_density <-aggregate(genes_down_tile_coverage$meth_level,by=list(genes_down_tile_coverage$id),mean)
    message("Gene downstream methylation calculation is finished!")
    
    dd<-rbind(genes_up_tile_density,genes_body_tile_density,genes_down_tile_density)
    dd$n <- seq(1:300)
    dd$sample <- sample_name
    print(ggplot(dd,aes(x=n,y=x))+
            geom_smooth(span=0.1)+
            theme_pubr(base_size = 20)+
            xlab("")+
            ylab("CpG methylation level")+
            scale_x_continuous(limits=c(0,bins*3),
                               breaks = seq(0,bins*3,bins),
                               labels=c("-2 kb","TSS","TTS","2 kb")))
    write.table(dd,paste("results/gene_metaplot/CHH/",sample_name,x, sep = ""),row.names = F,col.names = T,quote = F,sep="\t")  
  }
  #rm(meth)
  #rm(chromsome_meth)
  gc()
}



###############################################################################
##################################4 all plots
dd_low <- read.table("results/gene_metaplot/CpG/D1806119_L3_J90351",header = T)
dd_low$class <- "low"
dd_mediam <- read.table("results/gene_metaplot/CpG/D1806119_L3_J90352",header = T)
dd_mediam$class <-"medium"
dd_high <- read.table("results/gene_metaplot/CpG/D1806119_L3_J90353",header = T)
dd_high$class <- "high"

dd1 <- rbind(dd_low,dd_mediam,dd_high)
dd1$type <- "CpG"


dd_low <- read.table("results/gene_metaplot/CHG/D18061191",header = T)
dd_low$class <- "low"
dd_mediam <- read.table("results/gene_metaplot/CHG/D18061192",header = T)
dd_mediam$class <-"medium"
dd_high <- read.table("results/gene_metaplot/CHG/D18061193",header = T)
dd_high$class <- "high"

dd2 <- rbind(dd_low,dd_mediam,dd_high)
dd2$type <- "CHG"


dd_low <- read.table("results/gene_metaplot/CHH//D18061191",header = T)
dd_low$class <- "low"
dd_mediam <- read.table("results/gene_metaplot/CHH/D18061192",header = T)
dd_mediam$class <-"medium"
dd_high <- read.table("results/gene_metaplot/CHH/D18061193",header = T)
dd_high$class <- "high"

dd3 <- rbind(dd_low,dd_mediam,dd_high)
dd3$type <- "CHH"

dd <- rbind(dd1,dd2,dd3)

p<-ggplot(dd,aes(x=n,y=x,col=class))+
  geom_smooth(span=0.1,se = F)+
  theme_pubr(base_size = 12)+
  xlab("")+
  ylab("Methylation level")+
  scale_x_continuous(limits=c(0,bins*3),
                     breaks = seq(0,bins*3,bins),
                     labels=c("-2 kb","TSS","TTS","2 kb"))+
  facet_wrap(.~type,scales = "free")

ggsave(p,filename = "results/gene_expression_methylation.pdf",width = 7,height = 3)



### AA and CC  

dd_A1 <- read.table("results/gene_metaplot/CpG/D1806119_L3_J9035A",header = T)
dd_A1$class <- "AA"
dd_C1 <- read.table("results/gene_metaplot/CpG/D1806119_L3_J9035C",header = T)
dd_C1$class <-"CC"


dd_A2 <- read.table("results/gene_metaplot_chg_AACC/_D1806119_A",header = T)
dd_A2$class <- "AA"
dd_C2 <- read.table("results/gene_metaplot_chg_AACC/_D1806119_C",header = T)
dd_C2$class <-"CC"

dd_A3 <- read.table("results/gene_metaplot_chh_AACC/_D1806119_A",header = T)
dd_A3$class <- "AA"
dd_C3 <- read.table("results/gene_metaplot_chh_AACC/_D1806119_C",header = T)
dd_C3$class <-"CC"


dd1 <- rbind(dd_A1,dd_C1)
dd1$type <- "CpG"

dd2 <- rbind(dd_A2,dd_C2)
dd2$type <- "CHG"
dd3 <- rbind(dd_A3,dd_C3)
dd3$type <- "CHH"

dd <- rbind(dd1,dd2,dd3)

p2<-ggplot(dd,aes(x=n,y=x,col=class))+
  geom_smooth(span=0.1,se = F)+
  theme_pubr(base_size = 12)+
  xlab("")+
  ylab("Methylation level")+
  scale_x_continuous(limits=c(0,bins*3),
                     breaks = seq(0,bins*3,bins),
                     labels=c("-2 kb","TSS","TTS","2 kb"))+
facet_wrap(.~type,scales = "free")
ggsave(p2,filename = "results/methylation_AACC.pdf",width = 7,height = 3)
