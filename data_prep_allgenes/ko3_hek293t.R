library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
human <- readRDS("/home/kunqi/conservation/result/human.rds")
mcols(human) <- NULL
## huamn: m6a site in experiment 

Methylation_level <- readRDS("/home/kunqi/meta_analysis/human_meRIP_Level.rds")

## add the GSE55572 Mettle3 knockdown and control comparison 
mcols(human) <- Methylation_level[,c(14:15)]

human_exon <- subsetByOverlaps(human,exons(txdb)) #有表达的exon

#eliminate the NAs in V1 and V2
human_exon <- human_exon[-unique(which(is.na(human_exon$V1)),which(is.na(human_exon$V2)))] 

##human_exon  can be consider as experiment group
exp_pos <- length(which(c(human_exon$V1-human_exon$V2)>0))
exp_neg <- length(which(c(human_exon$V1-human_exon$V2)<=0))

# find the DRACH sequence motif 
motif <- m6ALogisticModel::sample_sequence("DRACH",subsetByOverlaps(exons(txdb),human),Hsapiens)
motif <- motif-2
motif <- motif[-which(motif%in%human_exon)]   ##此时里面应该全是negative，因为是减掉true

set.seed(213)
# set the negative control, which restrict to the same size as human_exon
background_motif <- motif[sample(1:length(motif))][1:length(human_exon)]

Ctrl_peak <- read.csv("/home/kunqi/human_peakCalling/GSE55572/HEK293T_Ctrl/Mod.csv")
Ctrl_peak_gr <- GRanges(seqnames = Ctrl_peak$chr,IRanges(start = Ctrl_peak$chromStart,end=Ctrl_peak$chromEnd),
                        strand = Ctrl_peak$strand) 

mcols(Ctrl_peak_gr) <- Ctrl_peak$log2FoldChange
M3KO_peak <- read.csv("/home/kunqi/human_peakCalling/GSE55572/HEK293T_M3KO/Mod.csv")
M3KO_peak_gr <- GRanges(seqnames = M3KO_peak$chr,IRanges(start = M3KO_peak$chromStart,end=M3KO_peak$chromEnd),
                        strand = M3KO_peak$strand) 
mcols(M3KO_peak_gr) <- M3KO_peak$log2FoldChange

background_motif$Ctrl <- 0 
background_motif$M3KO <- 0 
background_motif$Ctrl[queryHits(findOverlaps(background_motif,Ctrl_peak_gr))] <- (Ctrl_peak_gr$X[subjectHits(findOverlaps(background_motif,Ctrl_peak_gr))])
background_motif$M3KO[queryHits(findOverlaps(background_motif,M3KO_peak_gr))] <- (M3KO_peak_gr$X[subjectHits(findOverlaps(background_motif,M3KO_peak_gr))])
##background_motif  can be consider as control group 

ctrl_pos <- length(which(c(background_motif$Ctrl-background_motif$M3KO)>0))
ctrl_neg <- length(which(c(background_motif$M3KO-background_motif$Ctrl)>=0))

exp_pos 
exp_neg 
ctrl_pos
ctrl_neg

