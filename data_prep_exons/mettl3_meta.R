library(metafor)

# Mettl3 <- read.csv("/home/yuxuan.wu/meta_analysis/mettl3_exon_divided.csv")
# Mettl3 <- read.csv("/home/yuxuan.wu/meta_analysis/Mettl3.csv")
# Mettl3 <- read.csv("/home/yuxuan.wu/meta_analysis/mettl3_exon_ctrlChange.csv")
# Mettl3 <- read.csv("/home/yuxuan.wu/meta_analysis/mettl3_exon_divided_ctrl.csv")
Mettl3 <- read.csv("/home/yuxuan.wu/meta_analysis/Mettl3_exon.csv")
# Mettl3 <- read.csv("/home/yuxuan.wu/meta_analysis/mettl3_total.csv")

# Mettl3 <- Mettl3[c(1,5,2),]
dat <- escalc(measure = "RR", ai=exp_pos, bi=exp_neg, ci=ctrl_pos, di=ctrl_neg, data = Mettl3, append = TRUE)
res_reml <- rma.uni(yi, vi, data=dat, slab=paste(dat$X,digits=2, method="REML"))
res_reml$I2
baujat(res_reml)

# mettl14_exon
Mettl14_exon <- readRDS("~/meta_analysis/Mettl14_exon.rds")
res = rma(ci=ctrl_pos, di=ctrl_neg, ai=exp_pos, bi=exp_neg,data=Mettl14_exon, measure="RR", 
          method = "REML")
baujat(res)


dat <- escalc(measure = "RR", ai=exp_pos, bi=exp_neg, ci=ctrl_pos, di=ctrl_neg, data = Mettl3, append = TRUE)
res_reml <- rma.uni(yi, vi, data=dat, slab=paste(dat$X,digits=2, method="REML"))
res_reml$I2
forest(res_reml)
baujat(res_reml)

Mettl3$X <- c("GSE46705_HeLa","GSE55572_HEK293T","GSE55572_A549","GSE94808_GSC",
              "GSE94613_MOLM13","GSE110320_HepG2","GSE132306_EndoC-bH1","GSE93911_HEC-1-A")
res = rma(ci=ctrl_pos, di=ctrl_neg, ai=exp_pos, bi=exp_neg,data=Mettl3, measure="RR",slab = X,
          method = "REML")
# forest(res)

forest(res, xlim=c(-16, 6), at=log(c(1, 4, 16)), atransf=exp,
       ilab=cbind(Mettl3$exp_pos,Mettl3$exp_neg, Mettl3$ctrl_pos,Mettl3$ctrl_neg),
       ilab.xpos=c(-8.5,-6.5,-3.5,-1.5), cex=0.75, ylim=c(-1, 12),
       xlab="Log Risk Ratio", mlab="", psize=1, header = "dsafdsafjdashflaksd")

text(-16, -1, pos=4, cex=0.75, bquote(paste("RE Model for All Studies (Q = ",
                                            .(formatC(res$QE, digits=2, format="f")), ", df = ", .(res$k - res$p),
                                            ", p = ", .(formatC(res$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                            .(formatC(res$I2, digits=1, format="f")), "%)")))


# y = metabin(exp_pos, exp_total,ctrl_pos, ctrl_total, data = Mettl3, sm = "RR", studlab = X)
# forest(y)
# # plot(res,qqplot=TRUE)
# 
# res_DL <- rma.uni(yi, vi, data=dat, slab=paste(dat$X,digits=2, method="DL"))
# forest(res_DL)
# res_DL$I2
# 
# 
# res_HE <- rma.uni(yi, vi, data=dat, slab=paste(dat$X,digits=2, method="HE"))
# forest(res_HE)
# res_HE$I2
# 
# 
# res_HS <- rma.uni(yi, vi, data=dat, slab=paste(dat$X,digits=2, method="HS"))
# forest(res_HS)
# res_HS$I2
