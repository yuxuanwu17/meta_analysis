library(metafor)
Mettl14 <- read.csv("/Users/hopezhu/Desktop/meta-analysis/Mettl14_exon_group.csv")
Mettl14$X <- c("GSE94808_GSC_class1","GSE94808_GSC_class2","GSE94808_GSC_class3",
               "GSE94808_GSC_class4","GSE94808_GSC_class5",
               "GSE132306_EndoC−bH1_class1",
               "GSE132306_EndoC−bH1_class2",
               "GSE132306_EndoC−bH1_class3",
               "GSE132306_EndoC−bH1_class4",
               "GSE132306_EndoC−bH1_class5")
res = rma(ci=ctrl_pos, di=ctrl_neg, ai=exp_pos, bi=exp_neg,data=Mettl14, measure="RR",slab = X,
          method = "REML")
# forest(res)
par(cex=0.75, font=2)
par(mar=c(4,4,1,2))
forest(res, xlim=c(-16, 6), at=log(c(1, 4, 16)), atransf=exp,
       ilab=cbind(Mettl14$exp_pos,Mettl14$exp_neg, Mettl14$ctrl_pos,Mettl14$ctrl_neg),
       ilab.xpos=c(-8.5,-6.5,-3.5,-1.5), cex=0.75, ylim=c(-1, 20),
       rows = c(11:15,2:6),
       xlab="Log Risk Ratio", mlab="", psize=1)
par(font=2)
text(-16, -1, pos=4, cex=0.85, bquote(paste("RE Model for All Studies (Q = ",
                                            .(formatC(res$QE, digits=2, format="f")), ", df = ", .(res$k - res$p),
                                            ", p = ", .(formatC(res$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                            .(formatC(res$I2, digits=1, format="f")), "%)")))

op <- par(cex=0.75, font=2)
# 
text(-12.75,19, "GEO datasets ")
text( 4.00,19, "Correlation [95% CI]",cex = 0.8)
text(c(-8.5,-6.5,-3.25,-1.25), 19, c("KO+", "KO-", "KO+", "KO-"))
text(c(-7.35,-2.25),20, c("Experiment", "Control"))
#
op <- par(cex=0.75, font=4)
text(-16, c(16.5,7.5), pos=4, c("GSE55572_HEK293T","GSE132306_EndoC−bH1"))


res.s <- rma(ci=ctrl_pos, di=ctrl_neg, ai=exp_pos, bi=exp_neg,data=Mettl14, measure="RR",
             subset=(class=="1"), method="REML")
res.r <- rma(ci=ctrl_pos, di=ctrl_neg, ai=exp_pos, bi=exp_neg,data=Mettl14, measure="RR",
             subset=(class=="2"), method="REML")
# 
# 
addpoly(res.s, row=10, cex=0.75, atransf=exp, mlab="")
addpoly(res.r, row= 1, cex=0.75, atransf=exp, mlab="")


text(-16, 10, pos=4, cex=0.75, bquote(paste("RE Model for Subgroup (Q = ",
                                            .(formatC(res.s$QE, digits=2, format="f")), ", df = ", .(res.s$k - res.s$p),
                                            ", p = ", .(formatC(res.s$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                            .(formatC(res.s$I2, digits=1, format="f")), "%)")))
text(-16, 1, pos=4, cex=0.75, bquote(paste("RE Model for Subgroup (Q = ",
                                           .(formatC(res.r$QE, digits=2, format="f")), ", df = ", .(res.r$k - res.r$p),
                                           ", p = ", .(formatC(res.r$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                           .(formatC(res.r$I2, digits=1, format="f")), "%)")))

