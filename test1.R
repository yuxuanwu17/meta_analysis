library(metafor)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
human <- readRDS("/home/kunqi/conservation/result/human.rds")
mcols(human) <- NULL

Methylation_level <- readRDS("/home/kunqi/meta_analysis/human_meRIP_Level.rds")
mcols(human) <- Methylation_level[,c(8,10,17,19,29,31,48,49,57,59,83,85)]
YGasd <- list()

  hk_exon <- readRDS("/home/kunqi/meta_analysis/HK_exons.rds")
  human_exon <- subsetByOverlaps(human,hk_exon) 
  
  human_exon <- human_exon[-unique(c(which(is.na(human_exon$V1)),which(is.na(human_exon$V2)),
                                     which(is.na(human_exon$V3)),which(is.na(human_exon$V4)),
                                     which(is.na(human_exon$V5)),which(is.na(human_exon$V6)),
                                     which(is.na(human_exon$V7)),which(is.na(human_exon$V8)),
                                     which(is.na(human_exon$V9)),which(is.na(human_exon$V10)),
                                     which(is.na(human_exon$V11)),which(is.na(human_exon$V12))))] 
  ##human_exon  can be consider as experiment group
  
  # length(which(c(human_exon$V1-human_exon$V2)>0))
  
  motif <- m6ALogisticModel::sample_sequence("DRACH",subsetByOverlaps(hk_exon,human),Hsapiens)
  motif <- motif-2
  motif <- motif[-which(motif%in%human_exon)]
  set.seed(213)
  background_motif <- motif[sample(1:length(motif))][1:length(human_exon)]
  
  XXX <- mcols(human_exon)
  ix <- c("GSE46705/HeLa_Ctrl","GSE55572/A549_Ctrl","GSE94808/GSC_Ctrl",
          "GSE90642/HepG2_Ctrl","GSE110320/HepG2_Ctrl","GSE132306/EndoC-bH1_Ctrl")
  ix2 <- c("GSE46705/HeLa_M14KO","GSE55572/A549_M14KO","GSE94808/GSC_M14KO",
           "GSE90642/HepG2_METTL14KD","GSE110320/HepG2_M14KO","GSE132306/EndoC-bH1_M14KO")
  for(i in 1:6){
    Ctrl_peak <- read.csv(paste0("/home/kunqi/human_peakCalling/",ix[i],"/Mod.csv"))
    Ctrl_peak_gr <- GRanges(seqnames = Ctrl_peak$chr,IRanges(start = Ctrl_peak$chromStart,end=Ctrl_peak$chromEnd),
                            strand = Ctrl_peak$strand) 
    mcols(Ctrl_peak_gr) <- Ctrl_peak$log2FoldChange
    M14KO_peak <- read.csv(paste0("/home/kunqi/human_peakCalling/",ix2[i],"/Mod.csv"))
    M14KO_peak_gr <- GRanges(seqnames = M14KO_peak$chr,IRanges(start = M14KO_peak$chromStart,end=M14KO_peak$chromEnd),
                             strand = M14KO_peak$strand) 
    mcols(M14KO_peak_gr) <- M14KO_peak$log2FoldChange
    
    background_motif$Ctrl <- 0 
    background_motif$M14KO <- 0 
    background_motif$Ctrl[queryHits(findOverlaps(background_motif,Ctrl_peak_gr))] <- (Ctrl_peak_gr$X[subjectHits(findOverlaps(background_motif,Ctrl_peak_gr))])
    background_motif$M14KO[queryHits(findOverlaps(background_motif,M14KO_peak_gr))] <- (M14KO_peak_gr$X[subjectHits(findOverlaps(background_motif,M14KO_peak_gr))])
    XXX[,c((i*2-1),i*2)] <- mcols(background_motif)
  }
  
  XXX1 <- matrix(NA,ncol = 4,nrow = 6) 
  
  for(i in 1:6){
    XXX1[i,1] <- length(which((mcols(human_exon)[,(i*2-1)]-mcols(human_exon)[,(i*2)])>0))
    XXX1[i,2] <- length(which((mcols(human_exon)[,(i*2-1)]-mcols(human_exon)[,(i*2)])<=0))
    XXX1[i,3] <- length(which((XXX[,(i*2-1)]-XXX[,(i*2)])>0))
    XXX1[i,4] <- length(which((XXX[,(i*2-1)]-XXX[,(i*2)])<=0))
  }
  YGasd[[join]] <- XXX1


YGasd1 <- YGasd[[1]]
for(kj in 2:10){
  YGasd1 <- rbind(YGasd1,YGasd[[kj]])
}

colnames(YGasd1) <- c( "exp_pos","exp_neg","ctrl_pos", "ctrl_neg")

## find the suitable combination (3,6)
colnames(XXX1) <- c( "exp_pos","exp_neg","ctrl_pos", "ctrl_neg")
XXX1

