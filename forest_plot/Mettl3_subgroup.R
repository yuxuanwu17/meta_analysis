library(metafor)
Mettl3 <- read.csv("/Users/hopezhu/Desktop/meta-analysis/Mettl3_exon_group.csv")
Mettl3$X <- c("GSE55572_HEK293T_class1","GSE55572_HEK293T_class2","GSE55572_HEK293T_class3",
              "GSE132306_EndoC−bH1_class1","GSE132306_EndoC−bH1_class2","GSE132306_EndoC−bH1_class3",
              "GSE93911_HEC−1−A_class1","GSE93911_HEC−1−A_class2","GSE93911_HEC−1−A_class3")
res = rma(ci=ctrl_pos, di=ctrl_neg, ai=exp_pos, bi=exp_neg,data=Mettl3, measure="RR",slab = X,
          method = "REML")
# forest(res)
par(cex=0.75, font=2)
par(mar=c(4,4,1,2))
forest(res, xlim=c(-16, 6), at=log(c(1, 4, 16)), atransf=exp,
       ilab=cbind(Mettl3$exp_pos,Mettl3$exp_neg, Mettl3$ctrl_pos,Mettl3$ctrl_neg),
       ilab.xpos=c(-8.5,-6.5,-3.5,-1.5), cex=0.75, ylim=c(-1, 20),
       rows = c(14:16,8:10,2:4),
       xlab="Log Risk Ratio", mlab="", psize=1)
par(font=2)
text(-16, -1, pos=4, cex=0.85, bquote(paste("RE Model for All Studies (Q = ",
                                            .(formatC(res$QE, digits=2, format="f")), ", df = ", .(res$k - res$p),
                                            ", p = ", .(formatC(res$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                            .(formatC(res$I2, digits=1, format="f")), "%)")))

op <- par(cex=0.75, font=2)

text(-12.75,19, "GEO datasets ")
text( 4.00,19, "Correlation [95% CI]",cex = 0.8)
text(c(-8.5,-6.5,-3.25,-1.25), 19, c("KO+", "KO-", "KO+", "KO-"))
text(c(-7.35,-2.25),20, c("Experiment", "Control"))
# 
op <- par(cex=0.75, font=4)
text(-16, c(17,11,5), pos=4, c("GSE55572_HEK293T",
                               "GSE132306_EndoC−bH1",
                               "GSE93911_HEC−1−A"))


res.s <- rma(ci=ctrl_pos, di=ctrl_neg, ai=exp_pos, bi=exp_neg,data=Mettl3, measure="RR",
             subset=(class=="1"), method="REML")
res.r <- rma(ci=ctrl_pos, di=ctrl_neg, ai=exp_pos, bi=exp_neg,data=Mettl3, measure="RR",
             subset=(class=="2"), method="REML")
res.a <- rma(ci=ctrl_pos, di=ctrl_neg, ai=exp_pos, bi=exp_neg,data=Mettl3, measure="RR",
             subset=(class=="3"), method="REML")


addpoly(res.s, row=13, cex=0.75, atransf=exp, mlab="")
addpoly(res.r, row= 7, cex=0.75, atransf=exp, mlab="")
addpoly(res.a, row= 0.75, cex=0.75, atransf=exp, mlab="")

text(-16, 13, pos=4, cex=0.75, bquote(paste("RE Model for Subgroup (Q = ",
                                              .(formatC(res.s$QE, digits=2, format="f")), ", df = ", .(res.s$k - res.s$p),
                                              ", p = ", .(formatC(res.s$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                              .(formatC(res.s$I2, digits=1, format="f")), "%)")))
text(-16, 7, pos=4, cex=0.75, bquote(paste("RE Model for Subgroup (Q = ",
                                             .(formatC(res.r$QE, digits=2, format="f")), ", df = ", .(res.r$k - res.r$p),
                                             ", p = ", .(formatC(res.r$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                             .(formatC(res.r$I2, digits=1, format="f")), "%)")))
text(-16, .75, pos=4, cex=0.75, bquote(paste("RE Model for Subgroup (Q = ",
                                             .(formatC(res.a$QE, digits=2, format="f")), ", df = ", .(res.a$k - res.a$p),
                                             ", p = ", .(formatC(res.a$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                             .(formatC(res.a$I2, digits=1, format="f")), "%)")))
