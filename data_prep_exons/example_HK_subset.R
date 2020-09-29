library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
human <- readRDS("/home/kunqi/conservation/result/human.rds")
mcols(human) <- NULL

Methylation_level <- readRDS("/home/kunqi/meta_analysis/human_meRIP_Level.rds")
mcols(human) <- Methylation_level[,c(8,9,14,15,17,18,29,30,36,37,57,58,83,84,97,98)]
YGasd <- list()
for(join in 1:10){
  po <- (join-1)*5133+1
  mo <-(join)*5133
  asda <- readRDS("/home/kunqi/meta_analysis/HK_exons.rds")[po:mo]
  human_exon <- subsetByOverlaps(human,asda) 
  
  human_exon <- human_exon[-unique(c(which(is.na(human_exon$V1)),which(is.na(human_exon$V2)),
                                   which(is.na(human_exon$V3)),which(is.na(human_exon$V4)),
                                   which(is.na(human_exon$V5)),which(is.na(human_exon$V6)),
                                   which(is.na(human_exon$V7)),which(is.na(human_exon$V8)),
                                   which(is.na(human_exon$V9)),which(is.na(human_exon$V10)),
                                   which(is.na(human_exon$V11)),which(is.na(human_exon$V12)),
                                   which(is.na(human_exon$V13)),which(is.na(human_exon$V14)),
                                   which(is.na(human_exon$V15)),which(is.na(human_exon$V16))))] 
  ##human_exon  can be consider as experiment group
  
  # length(which(c(human_exon$V1-human_exon$V2)>0))
  
  motif <- m6ALogisticModel::sample_sequence("DRACH",subsetByOverlaps(asda,human),Hsapiens)
  motif <- motif-2
  motif <- motif[-which(motif%in%human_exon)]
  set.seed(213)
  background_motif <- motif[sample(1:length(motif))][1:length(human_exon)]
  
  XXX <- mcols(human_exon)
  ix <- c("GSE46705/HeLa_Ctrl","GSE55572/HEK293T_Ctrl","GSE55572/A549_Ctrl","GSE94808/GSC_Ctrl",
          "GSE94613/MOLM13_Ctrl","GSE110320/HepG2_Ctrl","GSE132306/EndoC-bH1_Ctrl","GSE93911/HEC-1-A_Ctrl")
  ix2 <- c("GSE46705/HeLa_M3KO","GSE55572/HEK293T_M3KO","GSE55572/A549_M3KO","GSE94808/GSC_M3KO",
           "GSE94613/MOLM13_METTL3KD","GSE110320/HepG2_M3KO","GSE132306/EndoC-bH1_M3KO","GSE93911/HEC-1-A_METTL3KD")
  for(i in 1:8){
    Ctrl_peak <- read.csv(paste0("/home/kunqi/human_peakCalling/",ix[i],"/Mod.csv"))
    Ctrl_peak_gr <- GRanges(seqnames = Ctrl_peak$chr,IRanges(start = Ctrl_peak$chromStart,end=Ctrl_peak$chromEnd),
                            strand = Ctrl_peak$strand) 
    mcols(Ctrl_peak_gr) <- Ctrl_peak$log2FoldChange
    M3KO_peak <- read.csv(paste0("/home/kunqi/human_peakCalling/",ix2[i],"/Mod.csv"))
    M3KO_peak_gr <- GRanges(seqnames = M3KO_peak$chr,IRanges(start = M3KO_peak$chromStart,end=M3KO_peak$chromEnd),
                            strand = M3KO_peak$strand) 
    mcols(M3KO_peak_gr) <- M3KO_peak$log2FoldChange
    
    background_motif$Ctrl <- 0 
    background_motif$M3KO <- 0 
    background_motif$Ctrl[queryHits(findOverlaps(background_motif,Ctrl_peak_gr))] <- (Ctrl_peak_gr$X[subjectHits(findOverlaps(background_motif,Ctrl_peak_gr))])
    background_motif$M3KO[queryHits(findOverlaps(background_motif,M3KO_peak_gr))] <- (M3KO_peak_gr$X[subjectHits(findOverlaps(background_motif,M3KO_peak_gr))])
    XXX[,c((i*2-1),i*2)] <- mcols(background_motif)
  }
  
  XXX1 <- matrix(NA,ncol = 4,nrow = 8) 
  
  for(i in 1:8){
    XXX1[i,1] <- length(which((mcols(human_exon)[,(i*2-1)]-mcols(human_exon)[,(i*2)])>0))
    XXX1[i,2] <- length(which((mcols(human_exon)[,(i*2-1)]-mcols(human_exon)[,(i*2)])<=0))
    XXX1[i,3] <- length(which((XXX[,(i*2-1)]-XXX[,(i*2)])>0))
    XXX1[i,4] <- length(which((XXX[,(i*2-1)]-XXX[,(i*2)])<=0))
  }
  YGasd[[join]] <- XXX1
}

YGasd1 <- YGasd[[1]]
for(kj in 2:10){
  YGasd1 <- rbind(YGasd1,YGasd[[kj]])
}

colnames(YGasd1) <- c( "exp_pos","exp_neg","ctrl_pos", "ctrl_neg")
res = rma(ai=ctrl_pos, bi=ctrl_neg, ci=exp_pos, di=exp_neg, data=YGasd1[c(2,10,18,26,50,7,15,23,31,55),], measure="RR")
print(res$I2)
forest(res)
#typ1
asdiuh <- YGasd1[c(2,10,18,26,50,7,15,23,31,55),]
asdiuh <- as.data.frame(asdiuh)
asdiuh$gene <- c(1,2,3,4,5,1,2,3,4,5)
asdiuh$class <- c(1,1,1,1,1,2,2,2,2,2)
reg = rma(ai=ctrl_pos, bi=ctrl_neg, ci=exp_pos, di=exp_neg,data=asdiuh,mods=~gene+class, measure="RR")
reg$I2
library(metafor)
#type2
asdiuh <- YGasd1[c(2,10,50,7,15,55,8,16,56),]
asdiuh <- as.data.frame(asdiuh)
asdiuh$gene <- c(1,2,3,1,2,3,1,2,3)
asdiuh$class <- c(1,1,1,2,2,2,3,3,3)
reg = rma(ai=ctrl_pos, bi=ctrl_neg, ci=exp_pos, di=exp_neg,data=asdiuh,mods=~class, measure="RR")
reg$I2
forest(reg)

library(meta)
colnames(YGasd1) <- c( "exp_pos","exp_neg","ctrl_pos", "ctrl_neg")
res = rma(ai=ctrl_pos, bi=ctrl_neg, ci=exp_pos, di=exp_neg, data=YGasd1[c(2,10,50,7,15,55,8,16,56),], measure="RR")
res = rma(ci=ctrl_pos, di=ctrl_neg, ai=exp_pos, bi=exp_neg,mods=~class, data=asdiuh, measure="RR")
forest(res)
print(res$I2)

forest(res)
