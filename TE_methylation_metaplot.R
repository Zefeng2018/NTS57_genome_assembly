library(rtracklayer)
library(Biostrings)
library(BSgenome)
library(stringr)
library(GenomicFeatures)
library(methylKit)
library(ggpubr)

setwd("/DATA4T/Brassica_napus_methylation/")
# load genomic features
Repeats<- import.gff3("../Brassica_napus_assembly/assembly/result_202301300547/Annotation/06_release/01_Genome_Repeats/youcai_repeats.gff")

################################################################################
####################### 0. Methylkit data load (multiple samples analysis)===== 
file.names<- dir("4methylkit/",pattern = "CHH.txt",full.names = T)[3:10]
sample_id <- str_replace_all(basename(file.names),pattern = "_CHH.txt",replacement = "")
file.list <-as.list(file.names)

myobj=methRead(file.list,
               sample.id=as.list(sample_id),
               assembly="N.napus",
               treatment=c(rep(1,8)),
               context = "CHH")

### class TEs into different types

TE_list<- list()
TE_list[[1]]<- Repeats[grep(pattern = "LTR",Repeats$type)] # 400301
TE_list[[2]]<- Repeats[grep(pattern = "DNA",Repeats$type)] # 490882
TE_list[[3]]<- Repeats[grep(pattern = "tandem_repeat",Repeats$type)] # 109213
strand(TE_list[[3]]) <-"+"
names(TE_list) <- c("LTR","DNA","tandem_repeat")

for (m in seq(1,length(myobj))){
  message("Sample: ",m)
  sample_name = sample_id[m]
  message(sample_name)
  chromsome_meth <- as(myobj[[m]],"GRanges")
  
  flank= 2000
  bins = 10
  date()
  
  
  for (x in names(TE_list)){
    message("TE list:",x)
    ## classify TE into positive and negative
    ## 1 gene body 
    TE_class <- TE_list[[x]]
    #TE_class <- TE_class[width(TE_class)>=100]
    TE_pos <- TE_class[strand(TE_class)=="+"]
    TE_pos_tile <- tile(TE_pos,n = bins) 
    TE_pos_tile_unlist <- unlist(TE_pos_tile)
    TE_pos_tile_unlist$id <- rep(seq(1,bins),length(TE_pos))
    
    TE_neg <- TE_class[strand(TE_class)=="-"]
    TE_neg_tile <- tile(TE_neg,n = bins) 
    TE_neg_tile_unlist <- unlist(TE_neg_tile)
    TE_neg_tile_unlist$id <- rep(seq(bins,1),length(TE_neg)) # reverse bins
    
    TE_body_tile <- c(TE_pos_tile_unlist,TE_neg_tile_unlist)
    
    overlaps <- findOverlaps(chromsome_meth, TE_body_tile)           # find the overlapped methylated site in each bin 
    signal <- chromsome_meth[queryHits(overlaps)]                    
    sumSignal_Cs <- aggregate(signal$numCs, list(subjectHits(overlaps)), sum) ## very important: get methylated C
    sumSignal_Ts <- aggregate(signal$numTs, list(subjectHits(overlaps)), sum) ## very important: get unmethylated T
    
    TE_body_tile$coverage_Cs <- NA
    TE_body_tile$coverage_Ts <- NA
    
    TE_body_tile$coverage_Cs[sumSignal_Cs$Group.1] <- sumSignal_Cs$x
    TE_body_tile$coverage_Ts[sumSignal_Ts$Group.1] <- sumSignal_Ts$x
    
    TE_body_tile_coverage<- as.data.frame(mcols(TE_body_tile))
    TE_body_tile_coverage$gid <- str_split_fixed(rownames(TE_body_tile_coverage),"\\.",2)[,1]
    TE_body_tile_coverage<- TE_body_tile_coverage[!is.na(TE_body_tile_coverage$coverage_Ts),] # remove bins that not have no covered reads
    TE_body_tile_coverage$meth_level <- TE_body_tile_coverage$coverage_Cs/(TE_body_tile_coverage$coverage_Cs+TE_body_tile_coverage$coverage_Ts)
    
    TE_body_tile_density <-aggregate(TE_body_tile_coverage$meth_level,by=list(TE_body_tile_coverage$id),mean) # get the mean methylation of each bin
    message("Gene body methylation calculation is finished!")
    ## 2 upstream
    TE_pos_up <- flank(TE_pos,width = flank,both = FALSE)
    TE_pos_up_tile <- tile(TE_pos_up,n = bins) 
    TE_pos_up_tile_unlist <- unlist(TE_pos_up_tile)
    TE_pos_up_tile_unlist$id <- rep(seq(1,bins),length(TE_pos_up))
    
    TE_neg_up <- flank(TE_neg,width = flank,both = FALSE) #  
    TE_neg_up_tile <- tile(TE_neg_up,n =bins) 
    TE_neg_up_tile_unlist <- unlist(TE_neg_up_tile)
    TE_neg_up_tile_unlist$id <- rep(seq(bins,1),length(TE_neg_up))
    
    TE_up_tile <- c(TE_pos_up_tile_unlist,TE_neg_up_tile_unlist)
    
    
    overlaps <- findOverlaps(chromsome_meth, TE_up_tile)           # find the overlapped methylated site in each bin 
    signal <- chromsome_meth[queryHits(overlaps)]                    
    sumSignal_Cs <- aggregate(signal$numCs, list(subjectHits(overlaps)), sum) ## very important: get methylated C
    sumSignal_Ts <- aggregate(signal$numTs, list(subjectHits(overlaps)), sum) ## very important: get unmethylated T
    
    TE_up_tile$coverage_Cs <- NA
    TE_up_tile$coverage_Ts <- NA
    
    TE_up_tile$coverage_Cs[sumSignal_Cs$Group.1] <- sumSignal_Cs$x
    TE_up_tile$coverage_Ts[sumSignal_Ts$Group.1] <- sumSignal_Ts$x
    
    TE_up_tile_coverage<- as.data.frame(mcols(TE_up_tile))
    TE_up_tile_coverage$gid <- str_split_fixed(rownames(TE_up_tile_coverage),"\\.",2)[,1]
    TE_up_tile_coverage <- TE_up_tile_coverage[!is.na(TE_up_tile_coverage$coverage_Ts),] # remove bins that not have no covered reads
    TE_up_tile_coverage$meth_level <- TE_up_tile_coverage$coverage_Cs/(TE_up_tile_coverage$coverage_Cs+TE_up_tile_coverage$coverage_Ts)
    
    TE_up_tile_density <-aggregate(TE_up_tile_coverage$meth_level,by=list(TE_up_tile_coverage$id),mean)
    
    message("Gene upstream methylation calculation is finished!")
    ## 3. downstream 1kb
    
    TE_pos_down <- flank(TE_pos,width = flank,both = FALSE,start = FALSE)
    TE_pos_down_tile <- tile(TE_pos_down,n = bins) 
    TE_pos_down_tile_unlist <- unlist(TE_pos_down_tile)
    TE_pos_down_tile_unlist$id <- rep(seq(1,bins),length(TE_pos_down))
    
    TE_neg_down <- flank(TE_neg,width = flank,both = FALSE,start = FALSE)
    TE_neg_down_tile <- tile(TE_neg_down,n = bins) 
    TE_neg_down_tile_unlist <- unlist(TE_neg_down_tile)
    TE_neg_down_tile_unlist$id <- rep(seq(bins,1),length(TE_neg_down))
    
    TE_down_tile <- c(TE_pos_down_tile_unlist,TE_neg_down_tile_unlist)
    
    overlaps <- findOverlaps(chromsome_meth, TE_down_tile)           # find the overlapped methylated site in each bin 
    signal <- chromsome_meth[queryHits(overlaps)]                    
    sumSignal_Cs <- aggregate(signal$numCs, list(subjectHits(overlaps)), sum) ## very important: get methylated C
    sumSignal_Ts <- aggregate(signal$numTs, list(subjectHits(overlaps)), sum) ## very important: get unmethylated T
    
    TE_down_tile$coverage_Cs <- NA
    TE_down_tile$coverage_Ts <- NA
    
    TE_down_tile$coverage_Cs[sumSignal_Cs$Group.1] <- sumSignal_Cs$x
    TE_down_tile$coverage_Ts[sumSignal_Ts$Group.1] <- sumSignal_Ts$x
    
    TE_down_tile_coverage<- as.data.frame(mcols(TE_down_tile))
    TE_down_tile_coverage$gid <- str_split_fixed(rownames(TE_down_tile_coverage),"\\.",2)[,1]
    TE_down_tile_coverage <- TE_down_tile_coverage[!is.na(TE_down_tile_coverage$coverage_Ts),] # remove bins that not have no covered reads
    TE_down_tile_coverage$meth_level <- TE_down_tile_coverage$coverage_Cs/(TE_down_tile_coverage$coverage_Cs+TE_down_tile_coverage$coverage_Ts)
    
    TE_down_tile_density <-aggregate(TE_down_tile_coverage$meth_level,by=list(TE_down_tile_coverage$id),mean)
    message("Gene downstream methylation calculation is finished!")
    
    dd<-rbind(TE_up_tile_density,TE_body_tile_density,TE_down_tile_density)
    dd$n <- seq(1:(bins*3))
    dd$sample <- sample_name
    print(ggplot(dd,aes(x=n,y=x))+
            #geom_smooth(span=0.1)+
            geom_line()+
            theme_pubr(base_size = 20)+
            xlab("")+
            ylab("CpG methylation level")+
            scale_x_continuous(limits=c(0,bins*3),
                               breaks = seq(0,bins*3,bins),
                               labels=c("-2 kb","TSS","TTS","2 kb")))
    write.table(dd,paste("results/TE_metaplot_chh/",sample_name,x, sep = ""),row.names = F,col.names = T,quote = F,sep="\t")  
  }
  #rm(meth)
  #rm(chromsome_meth)
  gc()
}
