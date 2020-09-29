# Meta-analysis for knockdown of METTL3 or METTL14 affecting N6-methyladenosine methylation level

## Abstract
N6-methyladenosine (m6A) is the most prevalent internal chemical modification present in multiple eukaryotic mRNAs, which is dynamically installed by methyltransferases (“writers”) and removed by demethylases (“erasers”). METTL3 and METTL14 are two typical m6A writers and several studies have found that knockdown or depletion of METTL3 and METTL14 genes influence the methylation level. A meta-analysis study was performed to assess the strength and quality of current evidence regarding to it and assess the how the knockdown of METTL3 or METTL14 affect the mRNA methylation level. Data was collected and classified into groups. Metafor, a R package was used to generate the Random-effects meta-analysis model and calculate the risk ratios (RRs) and 95% CIs. Since the heterogeneity of the overall data is extremely high (I2 > 99%), specific samples were excluded to minimize the downgrade the heterogeneity to acceptable level. Three and two published trails for METTL3 and METTLE14 samples, respectively, were eligible for review and analysis due to the acceptable low heterogeneity of their combination (METTL3 (Q = 12.85, I2 = 37.8%) and METTL14 (Q = 45.48, I2 = 82.9%) ). The positive relationship between knockdown of METTL3/METTL14 and methylation level was successfully verified and the difference among cell lines and gene types was identified.

## Usage:
The full script was in final file. But due to the size limitation, I could not upload the data set, please email me if you need the data set

## Results

![Alt text](https://github.com/yuxuanwu17/meta_analysis/blob/master/figure/WechatIMG26.jpeg)

Figure 1. Forest plot of the basic model. This is a Random Effect model which accommodates differences in study sample. Q statistic is based on the chi-square distribution, most commonly testing heterogeneity. I2 index is the total variability in a set of effect sizes due to true heterogeneity. Log Risk Ratio is a measurement of effect size. 

![Alt text](https://github.com/yuxuanwu17/meta_analysis/blob/master/figure/WechatIMG30.png)

Figure 2. The Baujat plot to identify studies contributing to heterogeneity. The study ID numbers represent eight studies, respectively.


![Alt text](https://github.com/yuxuanwu17/meta_analysis/blob/master/figure/WechatIMG31.jpeg)

Figure 3. Forest plot of moderator analysis.

![Alt text](https://github.com/yuxuanwu17/meta_analysis/blob/master/figure/WechatIMG32.png)

Figure 4. The Funnel plot.
