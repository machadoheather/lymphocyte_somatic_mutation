## Nov 2020
# running gam of S9 mutations and genomic properties
# Using only the Gm blood line for replication timing, and importantly keeping only sites 
#with sufficient data (sum density signal > 95)


setwd("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/mutsig_byregion/regression_analysis_hdp_denovo_1Mb/results_S9_majority_reptimingGm")

location="local"
if (location=="farm"){
  root="/lustre/scratch116/casm/cgp/users/hm8"
} else if (location=="local"){
  root="/Users/hm8/volumes/hm8_network"
} else {
  warning("location not set- using current directory")
  root=getwd()
}

library(BSgenome.Hsapiens.UCSC.hg19)
library(selectiveInference)
library(gamsel)
library(mgcv)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)

load(file=paste("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/mutsig_byregion/regression_analysis_hdp_denovo_1Mb/results/resultsGR10Kb_hdp_denovo_1Mb.Rdata", sep=""))
s9_results_filtered = data.frame(s9_resultsGR)
# remove unnecessary columns
s9_results_noquant = s9_results_filtered[ ,grep(colnames(s9_results_filtered), pattern="q_", invert=TRUE)]
#todrop = c("H3K9ac", "H3K4me2", "H3K4me3") ## don't drop
#s9_results_noquant_histdrop = s9_results_noquant[,-which(colnames(s9_results_noquant) %in% todrop)]

myrep = read.table("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/mutsig_byregion/regression_analysis_hdp_denovo_1Mb/wgEncodeUwRepliSeqGm12878WaveSignalRep1_averagebed.txt", header=F, stringsAsFactors=FALSE)
mysum = read.table("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/mutsig_byregion/regression_analysis_hdp_denovo_1Mb/wgEncodeUwRepliSeqGm12878SumSignalRep1_averagebed.txt", header=F, stringsAsFactors=FALSE)
colnames(myrep) =  c("name","size","covered","sum","mean0","mean")
colnames(mysum) = c("name","size","covered","sumsum","mean0sum","meansum")
myrepsum = merge(myrep, mysum[,c("name","sumsum","mean0sum","meansum")], by="name")
coordinates_resultsGR = read.table(file="/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/mutsig_byregion/regression_analysis_hdp_denovo_1Mb/coordinates_resultsGR.bed", sep="\t", header=F, stringsAsFactors = T)
colnames(coordinates_resultsGR) = c("chr","start","end","name")
myrep_name = data.frame(merge(myrepsum, coordinates_resultsGR, by="name"))
myrep_name$rep_timing_Gm = myrep_name$mean0
myrep_name$rep_timing_Gm[myrep_name$mean0sum < 95] = NA 
myrep_name$chr_start = paste(myrep_name$chr, as.numeric(as.character(myrep_name$start))+1, sep="_")


# combining all data with new rep timing info
s9_results_noquant$chr_start = paste(s9_results_noquant$seqnames, s9_results_noquant$start, sep="_")
s9_results_noquant_newrep = na.omit(merge(s9_results_noquant, myrep_name[,c("chr_start","rep_timing_Gm")], by="chr_start"))

# looking at strange feature in replication timing data
# s9_results_noquant_newrep_maj = subset(s9_results_noquant_newrep, meanS9 > 0.2)
# ggplot(s9_results_noquant_newrep_maj, aes(sumS9, rep_timing_Gm))+
#   geom_point(alpha=0.2)+
#   geom_smooth()+
#   theme_light()#+
#   #xlim(c(0,5))
# subset(s9_results_noquant_newrep_maj, sumS9>1.2 & sumS9<2 & rep_timing_Gm>25)
## caused by the boundary of high meanS9 for 1 mutation (sumS9 just below 1 and low replication timing value) and
#    low meanS9 for 2 mutations (sumS9 just above 1 and high replication timing value)

# Histograms of each variable (only sites with mutations)
s9_majority_props_all = s9_results_noquant_newrep[,-which(colnames(s9_results_noquant_newrep) %in% c("seqnames", "start","end", "strand","width", "Nsamples", "Ig", "size","hist","chr","meanX1", "sumX1" ,"meanSIg" ,"sumSIg","chr_start","rep_timing_value_value"))]
s9_majority_props_all$reptimingNA = "include"
s9_majority_props_all_reptimingNA = na.omit(s9_majority_props_all)
s9_majority_props_all_reptimingNA$reptimingNA = "omit"
s9_majority_props_all_both = rbind(s9_majority_props_all, s9_majority_props_all_reptimingNA)
s9_majority_props_all_bothDF = melt(s9_majority_props_all_both, id.vars=c("sumS9", "reptimingNA", "Nmut", "meanS9"))
ggplot(subset(s9_majority_props_all_bothDF, Nmut>0 & meanS9 >= 0.5), aes(x=value) ) +
  geom_histogram()+
  facet_wrap(reptimingNA~variable, nrow=2, scales="free") + 
  theme_light()
ggsave("rawdata_histogram_S9maj_byreptimingNA.pdf", width=49, height=4)

# Subset 
s9_results_noquant_zero = subset(s9_results_noquant_newrep, Nmut==0)
s9_results_noquant_nonzero = subset(s9_results_noquant_newrep, Nmut>0 & meanS9>0.5) # 28929
s9_results_noquant_zero10K = s9_results_noquant_zero[sample(1:nrow(s9_results_noquant_zero), size=round(nrow(s9_results_noquant_nonzero)/9) ), ]
s9_results_noquant_nonzero_zero10K = rbind(s9_results_noquant_zero10K, s9_results_noquant_nonzero)
s9_majority_props = s9_results_noquant_nonzero_zero10K[,-which(colnames(s9_results_noquant_nonzero_zero10K) %in% c("chr_start","seqnames", "start","end", "strand","width", "Nmut", "Nsamples", "meanS9","Ig", "size","hist","chr","meanX1", "sumX1" ,"meanSIg" ,"sumSIg","rep_timing_value_value"))]
mean(s9_majority_props$sumS9==0) # 0.09999067

# X1
s1_results_noquant_nonzero = subset(s9_results_noquant_newrep, Nmut>0 & meanX1>0.5) # 28929
s1_results_noquant_zero10K = s9_results_noquant_zero[sample(1:nrow(s9_results_noquant_zero), size=round(nrow(s1_results_noquant_nonzero)/9) ), ]
s1_results_noquant_nonzero_zero10K = rbind(s1_results_noquant_zero10K, s1_results_noquant_nonzero)
s1_majority_props = s1_results_noquant_nonzero_zero10K[,-which(colnames(s1_results_noquant_nonzero_zero10K) %in% c("chr_start","seqnames", "start","end", "strand","width", "Nmut", "Nsamples", "meanS9","Ig", "size","hist","chr","meanX1", "sumSIg" ,"meanSIg" ,"sumS9","rep_timing_value_value"))]
mean(s1_majority_props$sumX1==0) # 0.1

# SIg
sIg_results_noquant_nonzero = subset(s9_results_noquant_newrep, Nmut>0 & meanSIg>0.5) # 28929
sIg_results_noquant_zero10K = s9_results_noquant_zero[sample(1:nrow(s9_results_noquant_zero), size=round(nrow(sIg_results_noquant_nonzero)/9) ), ]
sIg_results_noquant_nonzero_zero10K = rbind(sIg_results_noquant_zero10K, sIg_results_noquant_nonzero)
sIg_majority_props = sIg_results_noquant_nonzero_zero10K[,-which(colnames(sIg_results_noquant_nonzero_zero10K) %in% c("chr_start","seqnames", "start","end", "strand","width", "Nmut", "Nsamples", "meanS9","Ig", "size","hist","chr","meanX1", "sumX1" ,"meanSIg" ,"sumS9","rep_timing_value_value"))]
mean(sIg_majority_props$sumSIg==0) #  0.1000063

# write.table(s9_majority_props, file="s9_majority_props_reptimingGm.txt", col.names=T, row.names=F, quote=F)
# write.table(s1_majority_props, file="s1_majority_props_reptimingGm.txt", col.names=T, row.names=F, quote=F)
# write.table(sIg_majority_props, file="sIg_majority_props_reptimingGm.txt", col.names=T, row.names=F, quote=F) # 354

# HSC
hscNmut = read.table("../hscNmut10KB.txt", header=T, stringsAsFactors = F)
colnames(hscNmut) = c("chr" ,  "start", "end"  , "HSC_Nmut",  "ID" )
hscNmut_info = merge(hscNmut[,c("HSC_Nmut","ID")], s9_results_noquant_newrep, by.x="ID", by.y="chr_start")

hsc_results_noquant_zero = subset(hscNmut_info, HSC_Nmut==0)
hsc_results_noquant_nonzero = subset(hscNmut_info, HSC_Nmut>0) # 28929
hsc_results_noquant_zero10K = hsc_results_noquant_zero[sample(1:nrow(hsc_results_noquant_zero), size=round(nrow(hsc_results_noquant_nonzero)/9) ), ]
hsc_results_noquant_nonzero_zero10K = rbind(hsc_results_noquant_zero10K, hsc_results_noquant_nonzero)
hsc_majority_props = hsc_results_noquant_nonzero_zero10K[,-which(colnames(hsc_results_noquant_nonzero_zero10K) %in% c("ID","chr_start","seqnames", "start","end", "strand","width", "Nmut", "Nsamples", "meanS9","Ig", "size","hist","chr","meanX1", "sumX1" ,"meanSIg" ,"sumSIg","rep_timing_value_value","sumS9"))]
#write.table(hsc_majority_props, file="hsc_props_reptimingGm.txt", col.names=T, row.names=F, quote=F) # 354


######################################################################
############## Colinearity analysis (produce table of colinearity for supplement)
######################################################################
s9_majority_props_complete = s9_majority_props[complete.cases(s9_majority_props[,1]),]
s9_majority_props_scale = scale(s9_majority_props_complete) # center and scale 
# 8687/ 6215540 bins have an NA value for gc_content_value, so make sure to remove NA's
s9_majority_props_scale = na.omit(s9_majority_props_scale)
cor_mat = matrix(ncol=ncol(s9_majority_props_scale), nrow=ncol(s9_majority_props_scale))
colnames(cor_mat) = colnames(s9_majority_props)
rownames(cor_mat) = colnames(s9_majority_props)
for (i in 1:ncol(s9_majority_props_scale)){
  for (j in 1:ncol(s9_majority_props_scale)){
    cor_mat[i,j] = cor(s9_majority_props_scale[,i], s9_majority_props_scale[,j])
  }
}
write.table(cor_mat, file="correlation_coef.S9_majority_props_reptimingGm_Jan2021.txt", col.names=T, row.names=T, quote=F)

# plotting
cor_mat_var = cor_mat[-1,-1]
diag(cor_mat_var) = NA
cor_mat_melt = melt(cor_mat_var)
cor_mat_melt$value2 = cor_mat_melt$value
cor_mat_melt$value2[is.na(cor_mat_melt$value2)] = 1
ggplot(cor_mat_melt, aes(x=Var1, y=Var2, fill=value) )+
  geom_tile() +
  geom_text(aes(label=sprintf("%0.1f", round(value2, digits = 1)) ), size=2)+
  #scale_fill_gradient(low="red", high="blue", na.value = "white", limits = c(-1, 1) )
  scale_fill_distiller(expression(italic(r)), palette = "Spectral", na.value = "white", limits = c(-1, 1))+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(file="correlation_coef.heatmap.S9_majority_props_reptimingGm_Jan2021_Jan2021.pdf", width=12, height=10)  

#cor_mat = read.table(file="correlation_coef.S9glm_Nov2020.txt", header=T, stringsAsFactors = F)
cor_matDF =data.frame(data.frame(cor_mat)$rep_timing_Gm, colnames(cor_mat))
cor_matDF[order(abs(cor_matDF[,1]), decreasing = T),]
# data.frame.cor_mat..rep_timing_Gm               colnames.cor_mat.
# 37                       1.000000000                   rep_timing_Gm
# 4                       -0.532588562          cpg_islands_dist_log10
# 33                      -0.421163899                TAD_b_dist_log10
# 1                       -0.411971230                           sumS9
# 12                       0.379583367                   gene_dens_1e6
# 11                       0.339824440                gc_content_value
# 26                      -0.275029753                    LAD_dens_1e6
# 2                       -0.166158500              ALU_rep_dist_log10
# 10                      -0.150832367                   g4_dist_log10
# 27                       0.136790733              LTR_rep_dist_log10
# 28                      -0.129039920              MIR_rep_dist_log10


######################################################################
### Individual regressions- no scaling
######################################################################
s9_majority_props = read.table(file="s9_majority_props_reptimingGm.txt", header=T, stringsAsFactors=F)
pvalue_cat = coef_cat = r_cat = c()
for (i in 1:(ncol(s9_majority_props)-1)){
  myaov = lm(s9_majority_props[,1] ~ s9_majority_props[,i+1])
  pvalue_cat[i] = summary(myaov)$coefficients[2,4]
  coef_cat[i] = summary(myaov)$coefficients[2,1]
  r_cat[i] = summary(myaov)$r.squared
}
indlm_S9 = data.frame(category=colnames(s9_majority_props)[2:ncol(s9_majority_props)], coefficient=coef_cat, pvalue=pvalue_cat, R2=r_cat)
write.table(indlm_S9, file="individual_lm.S9_majority_props_reptimingGm.txt", quote=F, col.names=T, row.names=F)
head(indlm_S9[order(indlm_S9$R2, decreasing=T),], n=10)
#                           category coefficient        pvalue          R2
# 36                   rep_timing_Gm -0.01837326  0.000000e+00 0.169720295
# 3           cpg_islands_dist_log10  0.38463678  0.000000e+00 0.084561381
# 10                gc_content_value -3.05304271  0.000000e+00 0.045104884
# 32                TAD_b_dist_log10  0.36310834  0.000000e+00 0.042056766
# 11                   gene_dens_1e6 -0.55292731  0.000000e+00 0.037668148
# 25                    LAD_dens_1e6  0.22263286 2.546015e-123 0.011809294
# 9                    g4_dist_log10  0.19383501 4.883498e-122 0.011685172

### spline regressions- no scaling
pvalue_cat_p = pvalue_cat_s = coef_cat_p = r_cat = c()
for (i in 1:(ncol(s9_majority_props)-1)){
  myaov = gam(s9_majority_props[,1] ~ s(s9_majority_props[,i+1]) )
  pvalue_cat_p[i] = summary(myaov)$p.pv
  pvalue_cat_s[i] = summary(myaov)$s.pv
  coef_cat_p[i] = summary(myaov)$p.coeff
  r_cat[i] = summary(myaov)$r.sq
}
indgam_S9 = data.frame(category=colnames(s9_majority_props)[2:ncol(s9_majority_props)], coefficient_p=coef_cat_p, pvalue_p=pvalue_cat_p, pvalue_s=pvalue_cat_s, R2=r_cat)
write.table(indgam_S9, file="individual_gam.S9_majority_props_reptimingGm.txt", quote=F, col.names=T, row.names=F)
head(indgam_S9[order(indgam_S9$R2, decreasing=T),], n=10)
#                 category coefficient_p pvalue_p      pvalue_s          R2
# 36          rep_timing_Gm      1.153279        0  0.000000e+00 0.173599520
# 3  cpg_islands_dist_log10      1.153279        0  0.000000e+00 0.086911571
# 11          gene_dens_1e6      1.153279        0  0.000000e+00 0.066381474
# 32       TAD_b_dist_log10      1.153279        0  0.000000e+00 0.053319993
# 10       gc_content_value      1.153279        0  0.000000e+00 0.045500134
# 25           LAD_dens_1e6      1.153279        0 4.280706e-240 0.024417950
# 9           g4_dist_log10      1.153279        0 4.301762e-120 0.012340867
# 33    telomere_dist_log10      1.153279        0 3.424888e-105 0.010752635

## these are colinear: H3K4me2, H3K9ac, H3K4me3
# category coefficient_p pvalue_p  pvalue_s            R2
# 12    H2A.Z      1.153279        0 0.4971891  6.477923e-05
# 17  H3K4me2      1.153279        0 0.4604000  2.708475e-05
# 20   H3K9ac      1.153279        0 0.8410607 -2.044094e-05
# 18  H3K4me3      1.153279        0 0.9229986 -2.109846e-05



#############################################
# Full GAM model
########################### 
s9_majority_props = read.table(file="s9_majority_props_reptimingGm.txt", header=T, stringsAsFactors=F)
#todrop = c("H3K9ac", "H3K4me2", "H3K4me3") 
#s9_majority_props = s9_majority_props[,-which(colnames(s9_majority_props) %in% todrop)]

focal_ind = which(colnames(s9_majority_props)=="sumS9")
predictors = s9_majority_props[,-focal_ind]
response = s9_majority_props[,focal_ind]

# response vector and covariate matrix
response = response[complete.cases(predictors)]
covs = predictors[complete.cases(predictors),]

# model matrix WITHOUT intercept, scaled and centered
xx = model.matrix(response~.,covs)[,-1]
xx = scale(covs) # center and scale 
yy <- response

# train on 2/3, test on 1/3
train <- sample(1:nrow(xx), ceiling(nrow(xx)/3*2))
test <- (-train)

    # fit GAM with lasso penalties, do 5-fold cross-validation on training data
    # to find lamba (penalisation strength) that minimises classification error
    gam_cv_train <- cv.gamsel(xx[train,], yy[train], family='gaussian',
                        nfolds=5, type.measure='mse')
    # Had to include custom lambdas to test, as the model maxed out at the last one tested
    mylambdas = c(100,50,20,10,9:2,1.5,1.25,1,seq(from=0.9, to=0.2, by=-0.1))
    gam_cv_train <- cv.gamsel(xx[train,], yy[train], family='gaussian',
                          nfolds=5, type.measure='mse', lambda = mylambdas )
    gam_bestind <- gam_cv_train$index.min
    gam_bestlam <- gam_cv_train$lambda[gam_bestind] # same as  gam_cv_train$lambda.min
    
    # fit GAM on training data and use CV calculated best lamba index to predict
    # response on the held-out test data
    gam_fit_train <- gamsel(xx[train,], yy[train], family='gaussian', lambda = mylambdas)
    par(mfrow=c(1,2)); summary(gam_fit_train)
    if(gam_bestind != which.min(abs(gam_fit_train$lambdas - gam_bestlam))) stop('problem with gam_bestind')
    gam_predcat <- predict(gam_fit_train, index=gam_bestind, newdata=xx[test,], type='response')

    # fit GAM to full dataset, choosing bestlam from CV on 2/3 data
    gam_fit_full <- gamsel(xx, yy, family='gaussian', lambda = mylambdas)
    pdf(file="components_S9majority_gam_reptimingGm_fits.pdf", width=8, height=4)
    par(mfrow = c(1, 2))
    summary(gam_fit_full, label = TRUE)
    dev.off()
    gam_full_bestind <- which.min(abs(gam_fit_full$lambdas-gam_bestlam))
    gam_full_deviance = gam_fit_full$dev.ratio[gam_full_bestind] #  0.2016055
    variables_keep = colnames(predictors)[abs(gam_fit_full$alphas[,gam_full_bestind])>0]
    data.frame(variable=colnames(predictors), alpha=gam_fit_full$alphas[,gam_full_bestind], betas=gam_fit_full$betas[1:ncol(predictors),gam_full_bestind])
    # outputs
    ans <- list(realcat=yy[test], # real categories of test set
                gam_predcat = gam_predcat,
                gam_mod = gam_fit_full,
                gam_ind = gam_full_bestind,
                gam_dev = gam_full_deviance,
                variables_keep = variables_keep
                )
save(gam_fit_full, ans, file = 'S9majority_gam_reptimingGm_fits.RData')
#load(file = 'S9majority_gam_reptimingGm_fits.RData')


# $gam_dev
# [1] 0.1998495

# variables to keep in model
# 1] "centromere_dist_log10"           "cpg_islands_dist_log10"          "cruciform_inverted_rep_dens_3e3"
# [4] "DNA_rep_dist_log10"              "DNase"                           "g4_dist_log10"                  
# [7] "gc_content_value"                "gene_dens_1e6"                   "H3K27ac"                        
# [10] "H3K4me1"                         "H3K9me3"                         "L1_rep_dist_log10"              
# [13] "LAD_dens_1e6"                    "LTR_rep_dist_log10"              "recomb_rate_nearest_value"      
# [16] "SIMPLE_REPEAT_rep_dist_log10"    "telomere_dist_log10"             "rep_timing_Gm"

# individual gam results for significant factors in gam
mysig = subset(indgam_S9, category %in%  ans$variables_keep)
mysig[order(abs(mysig$R2), decreasing=T),]
# category coefficient_p pvalue_p      pvalue_s            R2
# 36                   rep_timing_Gm      1.153279        0  0.000000e+00  1.735995e-01
# 3           cpg_islands_dist_log10      1.153279        0  0.000000e+00  8.691157e-02
# 11                   gene_dens_1e6      1.153279        0  0.000000e+00  6.638147e-02
# 10                gc_content_value      1.153279        0  0.000000e+00  4.550013e-02
# 25                    LAD_dens_1e6      1.153279        0 4.280706e-240  2.441795e-02
# 9                    g4_dist_log10      1.153279        0 4.301762e-120  1.234087e-02
# 33             telomere_dist_log10      1.153279        0 3.424888e-105  1.075263e-02
# 26              LTR_rep_dist_log10      1.153279        0  6.517013e-67  6.892148e-03
# 4  cruciform_inverted_rep_dens_3e3      1.153279        0  1.369302e-64  6.738131e-03
# 28       recomb_rate_nearest_value      1.153279        0  1.838545e-35  3.932500e-03
# 2            centromere_dist_log10      1.153279        0  3.840816e-27  2.994625e-03
# 6               DNA_rep_dist_log10      1.153279        0  1.215045e-23  2.607900e-03
# 31    SIMPLE_REPEAT_rep_dist_log10      1.153279        0  2.216447e-13  1.451479e-03
# 23               L1_rep_dist_log10      1.153279        0  4.659429e-12  1.271714e-03
# 16                         H3K4me1      1.153279        0  6.308456e-03  3.158146e-04
# 8                            DNase      1.153279        0  3.668353e-01  6.869611e-05
# 21                         H3K9me3      1.153279        0  4.447974e-01  6.719430e-05
# 13                         H3K27ac      1.153279        0  7.178004e-01 -5.830739e-06

mysig[order(abs(mysig$R2), decreasing=T),]$category


########## Model comparisons
##### checking most appropriate regression for replication timing
RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}
## Linear
replication_lm = lm(sumS9 ~ rep_timing_Gm, data=s9_majority_props)
plot(replication_lm)
ggplot(focal_dataset, aes(sumS9, rep_timing_Gm) ) +
  geom_point() +
  stat_smooth()
model <- lm(sumS9 ~ rep_timing_Gm, data = s9_majority_props[train,])
predictions <- model %>% predict(s9_majority_props[test,])
RMSE(predictions, s9_majority_props$sumS9[test])
# 0.7709544

### polynomial
model <- lm(sumS9 ~ poly(rep_timing_Gm, 2, raw = TRUE), data = s9_majority_props[train,])
predictions <- model %>% predict(s9_majority_props[test,])
RMSE(predictions, s9_majority_props$sumS9[test])
# 0.7697171

model <- lm(sumS9 ~ poly(rep_timing_Gm, 3, raw = TRUE), data = s9_majority_props[train,])
predictions <- model %>% predict(s9_majority_props[test,])
RMSE(predictions, s9_majority_props$sumS9[test])
# 0.7695249



############################# X1: SBSblood + SBS1
######################################################################
### Individual regressions- no scaling
######################################################################
s1_majority_props = read.table(file="s1_majority_props_reptimingGm.txt", header=T, stringsAsFactors=F)
focal_ind = which(colnames(s1_majority_props)=="sumX1")
s1_majority_props = data.frame(sumX1=s1_majority_props[,focal_ind], s1_majority_props[,-focal_ind])
pvalue_cat = coef_cat = r_cat = c()
for (i in 1:(ncol(s1_majority_props)-1)){
  myaov = lm(s1_majority_props[,1] ~ s1_majority_props[,i+1])
  pvalue_cat[i] = summary(myaov)$coefficients[2,4]
  coef_cat[i] = summary(myaov)$coefficients[2,1]
  r_cat[i] = summary(myaov)$r.squared
}
indlm_s1 = data.frame(category=colnames(s1_majority_props)[2:ncol(s1_majority_props)], coefficient=coef_cat, pvalue=pvalue_cat, R2=r_cat)
write.table(indlm_s1, file="individual_lm.s1_majority_props_reptimingGm.txt", quote=F, col.names=T, row.names=F)
head(indlm_s1[order(indlm_s1$R2, decreasing=T),], n=10)
# category   coefficient       pvalue           R2
# 25                  LAD_dens_1e6  0.0460651243 3.795827e-16 0.0013903949
# 36                 rep_timing_Gm -0.0006358814 1.664131e-09 0.0007617134
# 26            LTR_rep_dist_log10 -0.0267296834 1.500040e-07 0.0005784740
# 27            MIR_rep_dist_log10 -0.0266087647 6.820721e-07 0.0005171772
# 34 triplex_mirror_rep_dist_log10  0.0258234551 1.628940e-06 0.0004820498
# 33           telomere_dist_log10 -0.0234734316 5.200613e-06 0.0004353533

### spline regressions- no scaling
pvalue_cat_p = pvalue_cat_s = coef_cat_p = r_cat = c()
for (i in 1:(ncol(s1_majority_props)-1)){
  myaov = gam(s1_majority_props[,1] ~ s(s1_majority_props[,i+1]) )
  pvalue_cat_p[i] = summary(myaov)$p.pv
  pvalue_cat_s[i] = summary(myaov)$s.pv
  coef_cat_p[i] = summary(myaov)$p.coeff
  r_cat[i] = summary(myaov)$r.sq
}
indgam_s1 = data.frame(category=colnames(s1_majority_props)[2:ncol(s1_majority_props)], coefficient_p=coef_cat_p, pvalue_p=pvalue_cat_p, pvalue_s=pvalue_cat_s, R2=r_cat)
write.table(indgam_s1, file="individual_gam.s1_majority_props_reptimingGm.txt", quote=F, col.names=T, row.names=F)
head(indgam_s1[order(indgam_s1$R2, decreasing=T),], n=10)
# category coefficient_p pvalue_p     pvalue_s           R2
# 36          rep_timing_Gm     0.7961655        0 3.957187e-61 0.0064472665
# 3  cpg_islands_dist_log10     0.7961655        0 1.141629e-32 0.0035593928
# 32       TAD_b_dist_log10     0.7961655        0 9.369144e-14 0.0016129370
# 25           LAD_dens_1e6     0.7961655        0 3.793548e-16 0.0013694509
# 11          gene_dens_1e6     0.7961655        0 6.972106e-10 0.0011643906
# 2   centromere_dist_log10     0.7961655        0 5.937946e-08 0.0009650810
# 33    telomere_dist_log10     0.7961655        0 3.474460e-07 0.0008772822



#############
# Full GAM model
########################### 
focal_ind = which(colnames(s1_majority_props)=="sumX1")
predictors = s1_majority_props[,-focal_ind]
response = s1_majority_props[,focal_ind]

# response vector and covariate matrix
response = response[complete.cases(predictors)]
covs = predictors[complete.cases(predictors),]

# model matrix WITHOUT intercept, scaled and centered
xx = model.matrix(response~.,covs)[,-1]
xx = scale(covs) # center and scale 
yy <- response

# train on 2/3, test on 1/3
train <- sample(1:nrow(xx), ceiling(nrow(xx)/3*2))
test <- (-train)

# fit GAM with lasso penalties, do 5-fold cross-validation on training data
# to find lamba (penalisation strength) that minimises classification error
gam_cv_train <- cv.gamsel(xx[train,], yy[train], family='gaussian',
                          nfolds=5, type.measure='mse')
# Had to include custom lambdas to test, as the model maxed out at the last one tested
mylambdas = c(100,50,20,10,9:2,1.5,1.25,1,seq(from=0.9, to=0.2, by=-0.1))
gam_cv_train <- cv.gamsel(xx[train,], yy[train], family='gaussian',
                          nfolds=5, type.measure='mse', lambda = mylambdas )
gam_bestind <- gam_cv_train$index.min
gam_bestlam <- gam_cv_train$lambda[gam_bestind] # same as  gam_cv_train$lambda.min

# fit GAM on training data and use CV calculated best lamba index to predict
# response on the held-out test data
gam_fit_train <- gamsel(xx[train,], yy[train], family='gaussian', lambda = mylambdas)
par(mfrow=c(1,2)); summary(gam_fit_train)
if(gam_bestind != which.min(abs(gam_fit_train$lambdas - gam_bestlam))) stop('problem with gam_bestind')
gam_predcat <- predict(gam_fit_train, index=gam_bestind, newdata=xx[test,], type='response')

# fit GAM to full dataset, choosing bestlam from CV on 2/3 data
gam_fit_full <- gamsel(xx, yy, family='gaussian', lambda = mylambdas)
pdf(file="components_s1majority_gam_reptimingGm_fits.pdf", width=8, height=4)
par(mfrow = c(1, 2))
summary(gam_fit_full, label = TRUE)
dev.off()
gam_full_bestind <- which.min(abs(gam_fit_full$lambdas-gam_bestlam))
gam_full_deviance = gam_fit_full$dev.ratio[gam_full_bestind] #  0.2016055
variables_keep = colnames(predictors)[abs(gam_fit_full$alphas[,gam_full_bestind])>0]
# [1] "ALU_rep_dist_log10"              "cruciform_inverted_rep_dens_3e3" "g4_dist_log10"                  
# [4] "gc_content_value"                "gene_dens_1e6"                   "H3K27me3"                       
# [7] "H3K4me1"                         "H3K79me2"                        "H3K9me3"                        
# [10] "L1_rep_dist_log10"               "LAD_dens_1e6"                    "LTR_rep_dist_log10"             
# [13] "MIR_rep_dist_log10"              "TAD_b_dist_log10"                "telomere_dist_log10"            
# [16] "triplex_mirror_rep_dist_log10" 

# outputs
ans <- list(realcat=yy[test], # real categories of test set
            gam_predcat = gam_predcat,
            gam_mod = gam_fit_full,
            gam_ind = gam_full_bestind,
            gam_dev = gam_full_deviance,
            variables_keep = variables_keep
)
save(gam_fit_full, ans, file = 's1majority_gam_reptimingGm_fits.RData')


############################# HSC
######################################################################
### Individual regressions- no scaling
######################################################################
hsc_majority_props = read.table(file="hsc_props_reptimingGm.txt", header=T, stringsAsFactors=F)
focal_ind = which(colnames(hsc_majority_props)=="HSC_Nmut")
hsc_majority_props = data.frame(HSC_Nmut=hsc_majority_props[,focal_ind], hsc_majority_props[,-focal_ind])
pvalue_cat = coef_cat = r_cat = c()
for (i in 1:(ncol(hsc_majority_props)-1)){
  myaov = lm(hsc_majority_props[,1] ~ hsc_majority_props[,i+1])
  pvalue_cat[i] = summary(myaov)$coefficients[2,4]
  coef_cat[i] = summary(myaov)$coefficients[2,1]
  r_cat[i] = summary(myaov)$r.squared
}
indlm_hsc = data.frame(category=colnames(hsc_majority_props)[2:ncol(hsc_majority_props)], coefficient=coef_cat, pvalue=pvalue_cat, R2=r_cat)
write.table(indlm_hsc, file="individual_lm.hsc_majority_props_reptimingGm.txt", quote=F, col.names=T, row.names=F)
head(indlm_hsc[order(indlm_hsc$R2, decreasing=T),], n=10)
# category  coefficient       pvalue           R2
# 36          rep_timing_Gm -0.001020102 6.962323e-14 0.0016449478
# 11          gene_dens_1e6 -0.084240738 2.886261e-13 0.0015630366
# 3  cpg_islands_dist_log10  0.027059522 4.736221e-10 0.0011379037
# 32       TAD_b_dist_log10  0.027399134 9.397974e-06 0.0005761119

### spline regressions- no scaling
pvalue_cat_p = pvalue_cat_s = coef_cat_p = r_cat = c()
for (i in 1:(ncol(hsc_majority_props)-1)){
  myaov = gam(hsc_majority_props[,1] ~ s(hsc_majority_props[,i+1]) )
  pvalue_cat_p[i] = summary(myaov)$p.pv
  pvalue_cat_s[i] = summary(myaov)$s.pv
  coef_cat_p[i] = summary(myaov)$p.coeff
  r_cat[i] = summary(myaov)$r.sq
}
indgam_hsc = data.frame(category=colnames(hsc_majority_props)[2:ncol(hsc_majority_props)], coefficient_p=coef_cat_p, pvalue_p=pvalue_cat_p, pvalue_s=pvalue_cat_s, R2=r_cat)
write.table(indgam_hsc, file="individual_gam.hsc_majority_props_reptimingGm.txt", quote=F, col.names=T, row.names=F)
head(indgam_hsc[order(indgam_hsc$R2, decreasing=T),], n=10)
# category coefficient_p pvalue_p     pvalue_s           R2
# 11                 gene_dens_1e6     0.9881407        0 5.477277e-12 0.0016459464
# 36                 rep_timing_Gm     0.9881407        0 7.010452e-14 0.0016156387
# 33           telomere_dist_log10     0.9881407        0 1.080684e-07 0.0012700652
# 3         cpg_islands_dist_log10     0.9881407        0 4.734887e-10 0.0011085806
# 25                  LAD_dens_1e6     0.9881407        0 2.608195e-06 0.0008376951



#############
# Full GAM model
########################### 
focal_ind = which(colnames(hsc_majority_props)=="HSC_Nmut")
predictors = hsc_majority_props[,-focal_ind]
response = hsc_majority_props[,focal_ind]

# response vector and covariate matrix
response = response[complete.cases(predictors)]
covs = predictors[complete.cases(predictors),]

# model matrix WITHOUT intercept, scaled and centered
xx = model.matrix(response~.,covs)[,-1]
xx = scale(covs) # center and scale 
yy <- response

# train on 2/3, test on 1/3
train <- sample(1:nrow(xx), ceiling(nrow(xx)/3*2))
test <- (-train)

# fit GAM with lasso penalties, do 5-fold cross-validation on training data
# to find lamba (penalisation strength) that minimises classification error
gam_cv_train <- cv.gamsel(xx[train,], yy[train], family='gaussian',
                          nfolds=5, type.measure='mse')
# Had to include custom lambdas to test, as the model maxed out at the last one tested
mylambdas = c(100,50,20,10,9:2,1.5,1.25,1,seq(from=0.9, to=0.2, by=-0.1))
gam_cv_train <- cv.gamsel(xx[train,], yy[train], family='gaussian',
                          nfolds=5, type.measure='mse', lambda = mylambdas )
gam_bestind <- gam_cv_train$index.min
gam_bestlam <- gam_cv_train$lambda[gam_bestind] # same as  gam_cv_train$lambda.min

# fit GAM on training data and use CV calculated best lamba index to predict
# response on the held-out test data
gam_fit_train <- gamsel(xx[train,], yy[train], family='gaussian', lambda = mylambdas)
par(mfrow=c(1,2)); summary(gam_fit_train)
if(gam_bestind != which.min(abs(gam_fit_train$lambdas - gam_bestlam))) stop('problem with gam_bestind')
gam_predcat <- predict(gam_fit_train, index=gam_bestind, newdata=xx[test,], type='response')

# fit GAM to full dataset, choosing bestlam from CV on 2/3 data
gam_fit_full <- gamsel(xx, yy, family='gaussian', lambda = mylambdas)
pdf(file="components_hscmajority_gam_reptimingGm_fits.pdf", width=8, height=4)
par(mfrow = c(1, 2))
summary(gam_fit_full, label = TRUE)
dev.off()
gam_full_bestind <- which.min(abs(gam_fit_full$lambdas-gam_bestlam))
gam_full_deviance = gam_fit_full$dev.ratio[gam_full_bestind] #  0.2016055
variables_keep = colnames(predictors)[abs(gam_fit_full$alphas[,gam_full_bestind])>0]
# [1] "centromere_dist_log10"         "cpg_islands_dist_log10"        "DNase"                         "gene_dens_1e6"                
# [5] "H4K20me1"                      "L2_rep_dist_log10"             "LAD_dens_1e6"                  "LTR_rep_dist_log10"           
# [9] "RNAseq"                        "SIMPLE_REPEAT_rep_dist_log10"  "telomere_dist_log10"           "triplex_mirror_rep_dist_log10"
# [13] "z_dna_motif_dist_log10"        "rep_timing_Gm" 

# outputs
ans <- list(realcat=yy[test], # real categories of test set
            gam_predcat = gam_predcat,
            gam_mod = gam_fit_full,
            gam_ind = gam_full_bestind,
            gam_dev = gam_full_deviance,
            variables_keep = variables_keep
)
save(gam_fit_full, ans, file = 'hscmajority_gam_reptimingGm_fits.RData')




########### Plotting regressions
### functions
mean_by_quantile = function(focalvalues, secondaryvalues, Nquantiles, fulldataset){
  #Nquantiles=20
  #focalvalues = per_bin_sig_prop[,focalcolumn]
  #secondaryvalues = per_bin_sig_prop[,my_categories]
  if (length(focalvalues) != nrow(secondaryvalues)){warning("Length of focalvalues do not equal number of rows of secondary values.")}
  allvalues = data.frame(focalvalues, secondaryvalues)
  allquantiles = seq(from=0, to=1, by=1/Nquantiles)
  myquantiles = quantile(focalvalues, probs=allquantiles)
  outlist = outlist_percentile = list()
  for (i in 1:Nquantiles){
    if (myquantiles[i] ==myquantiles[i+1]){
      focalquantile = subset(allvalues, focalvalues==myquantiles[i])
    } else {
      focalquantile = subset(allvalues, focalvalues>=myquantiles[i] & focalvalues<myquantiles[i+1]) 
    }
    my_means = apply(focalquantile[,2:ncol(focalquantile)], MARGIN=2, FUN=mean, na.rm=TRUE)
    my_percentiles = vector()
    for (j in 1:length(my_means)){
      percentile = ecdf(fulldataset[,colnames(secondaryvalues)[j]])
      my_percentiles[j] = percentile(my_means[j])
    }
    outlist[[i]] = c(allquantiles[i+1],myquantiles[i+1], my_means)
    outlist_percentile[[i]] = c(allquantiles[i+1],myquantiles[i+1], my_percentiles)
  }
  outDF = data.frame(do.call(rbind, outlist))
  out_percentileDF = data.frame(do.call(rbind, outlist_percentile))
  colnames(outDF) = c("quantile","focal_value",colnames(secondaryvalues))
  colnames(out_percentileDF) = c("quantile","focal_value",colnames(secondaryvalues))
  outDF$valuetype = "mean"
  out_percentileDF$valuetype = "percentile"
  outDF2 = data.frame(rbind(outDF, out_percentileDF))
  return(outDF2)
}

format_regression_analyses = function(per_bin_sig_prop, focalcolumn, ind_regression_results, N, fulldataset, loess_span=0.5){
#format_regression_analyses = function(per_bin_sig_prop, focalcolumn, mult_regression_results, ind_regression_results, N, fulldataset){
    
  #per_bin_sig_prop=s9_majority_props
  #mult_regression_results=stepAIC_res
  #ind_regression_results=indgam_S9
  #fulldataset=s9_majority_props
  #focalcolumn="sumS9"
  #N=20
  require("plyr")
  
  my_categories = c("ALU_rep_dist_log10","centromere_dist_log10" , "cpg_islands_dist_log10" , "cruciform_inverted_rep_dens_3e3","direct_rep_dist_log10" , "DNA_rep_dist_log10",  "DNAMethylSBS" , "DNase" ,"g4_dist_log10","gc_content_value" , "gene_dens_1e6" ,  "H2A.Z" ,"H3K27ac" , "H3K27me3" , "H3K36me3", "H3K4me1"  ,"H3K79me2", "H3K9me3", "H4K20me1" , "L1_rep_dist_log10" ,"L2_rep_dist_log10" , "LAD_dens_1e6","LTR_rep_dist_log10" ,"MIR_rep_dist_log10"  ,"recomb_rate_nearest_value" , "rep_timing_Gm","RNAseq", "short_tandem_rep_dens_3e3" ,"SIMPLE_REPEAT_rep_dist_log10" ,   "TAD_b_dist_log10" , "telomere_dist_log10", "triplex_mirror_rep_dist_log10" , "z_dna_motif_dist_log10")
 
  plotting_quantiles = seq(from=1/N, to=1, by=1/N)
  #quantileN_focal = mean_by_quantile(focalvalues = non_zero[,focalcolumn], secondaryvalues = non_zero[,my_categories], Nquantiles=N, fulldataset=fulldataset)
  quantileN_focal_all = mean_by_quantile(focalvalues = per_bin_sig_prop[,focalcolumn], secondaryvalues = per_bin_sig_prop[,my_categories], Nquantiles=N, fulldataset=fulldataset)
  quantileN_focal_melt = melt(quantileN_focal_all, id.vars=c("quantile", "focal_value","valuetype") )
  quantileN_focal_melt_lm = merge(quantileN_focal_melt, ind_regression_results, by.x="variable", by.y="category")
  quantileN_focal_melt_lm$variable = factor(quantileN_focal_melt_lm$variable, levels=unique(quantileN_focal_melt_lm$variable[order(quantileN_focal_melt_lm$R2, decreasing = T)]) )
  quantileN_focal_melt_lm$variable_new = revalue(quantileN_focal_melt_lm$variable, c("rep_timing_Gm"="Replication timing","cpg_islands_dist_log10"="CpG island", "gc_content_value" = "GC content", "LAD_dens_1e6" = "LAD density", "gene_dens_1e6" = "Gene density", "TAD_b_dist_log10" = "TAD boundaries","telomere_dist_log10" = "Telomere", "centromere_dist_log10" = "Centromere", "cruciform_inverted_rep_dens_3e3" = "Cruciform repeats", "L1_rep_dist_log10" = "L1 repeats", "SIMPLE_REPEAT_rep_dist_log10" = "Simple repeat", "DNA_rep_dist_log10"= "DNA replication", "ALU_rep_dist_log10" = "ALU repeats", "MIR_rep_dist_log10" = "MIR repeats", "g4_dist_log10" = "G-quadruplex", "LTR_rep_dist_log10" = "LTR repeats", "triplex_mirror_rep_dist_log10" = "Triplex mirror") )
  
  
  ## loess smoothing
  myvars = unique(quantileN_focal_melt_lm$variable)
  means_list = percent_list = list()
  for (i in 1:length(myvars)){
    myvar_focal = subset(quantileN_focal_melt_lm, variable==myvars[i])
    myvar_focal_percent = subset(myvar_focal, valuetype=="percentile")
    myvar_focal_mean = subset(myvar_focal, valuetype=="mean")
    myvar_focal_percent$value_loess = predict(loess(value~quantile, data=myvar_focal_percent, span=loess_span))
    myvar_focal_mean$value_loess = predict(loess(value~quantile, data=myvar_focal_mean, span=loess_span))
    means_list[[i]] = myvar_focal_mean
    percent_list[[i]] = myvar_focal_percent
  }
  quantileN_focal_melt_loess = data.frame(rbind(do.call(rbind, means_list), do.call(rbind, percent_list)))
  
  
  ### polygon plot prep
  myminquantile = min(quantileN_focal_melt_loess$quantile)
  quantileN_focal_melt_lm2 = quantileN_focal_melt_loess[order(quantileN_focal_melt_loess$quantile),]
  quantileN_focal_melt_lm2$plotting_only = FALSE
  tmp1 = quantileN_focal_melt_lm2[nrow(quantileN_focal_melt_lm2):1,]
  tmp1$value = 0
  tmp1$value_loess = 0
  quantileN_focal_melt_lm2$plotting_only = TRUE
  quantileN_focal_melt_lm3 = rbind(quantileN_focal_melt_lm2, tmp1)
  quantileN_focal_melt_lm3$zero_means[quantileN_focal_melt_lm3$quantile != myminquantile | quantileN_focal_melt_lm3$plotting_only == TRUE] = 0
  quantileN_focal_melt_lm3$zero_percentiles[quantileN_focal_melt_lm3$quantile != myminquantile | quantileN_focal_melt_lm3$plotting_only == TRUE] = 0
  quantileN_focal_melt_lm3$zero = 0
  return(quantileN_focal_melt_lm3)
}


## read in data S9
indgam_S9 = read.table(file="individual_gam.S9_majority_props_reptimingGm.txt", header=T, stringsAsFactors = F)
s9_majority_props = read.table(file="s9_majority_props_reptimingGm.txt", header=T, stringsAsFactors=F)
# stepAIC_res
load(file="S9majority_gam_reptimingGm_fits.RData")
s9_ans = ans
## X1 (S1 + Sblood)
# stepAIC_res_Sblood
s1_majority_props = read.table(file="s1_majority_props_reptimingGm.txt", header=T, stringsAsFactors = F)
indgam_s1 = read.table(file="individual_gam.s1_majority_props_reptimingGm.txt", header=T, stringsAsFactors = F)
load(file = 's1majority_gam_reptimingGm_fits.RData')
s1_ans = ans

# Calculate summaries per quantile and loess smoothing (span=0.3)
quantile20_S9 = format_regression_analyses(per_bin_sig_prop=s9_majority_props, focalcolumn="sumS9", ind_regression_results=indgam_S9, N=100, fulldataset=s9_majority_props)
quantile20_S9$dataset = "SBS9"
quantile20_S9$glm_sig = "no"
quantile20_S9$glm_sig[quantile20_S9$variable %in% s9_ans[[6]] ] = "yes"

quantile20_S1 = format_regression_analyses(per_bin_sig_prop=s1_majority_props, focalcolumn="sumX1", ind_regression_results=indgam_s1, N=100, fulldataset=s1_majority_props)
quantile20_S1$dataset = "SBS1/blood"
quantile20_S1$glm_sig = "no"
quantile20_S1$glm_sig[quantile20_S1$variable %in% s1_ans[[6]] ] = "yes"


### Plotting both
focalres = rbind(quantile20_S9, quantile20_S1)
focalres$sig_text = ""
focalres$sig_text[focalres$glm_sig=="yes"] = "*"
focalres$dataset = factor(focalres$dataset, levels=c("SBS9"  ,     "SBS1/blood"))
focalres$variable_new2 = focalres$variable_new
focalres$variable_new2 = revalue(focalres$variable_new2, c("Replication timing" = "Replication\n timing", "Cruciform repeats"= "Cruciform\n repeats", "recomb_rate_nearest_value" = "Recomb.\n rate", "G-quadruplex"="G-\n 4plex", "Gene density"="Gene\n density", "GC content"="GC\n content", "LAD density"="LAD\n density", "LTR repeats"="LTR\n repeats", "CpG island"="CpG\n island", "DNA replication"="DNA\n replication", "Centromere" = "Centromere") )
#save(focalres, file="S9_S1_table_forplotting_loess5.Rdata")

#library(extrafont)
#font_import()
ggplot(subset(focalres, variable_new %in% unique(focalres$variable_new[focalres$glm_sig=="yes" & focalres$R2>0.0025] ) & valuetype=="percentile" ), aes(x=quantile, y=value_loess, fill=R2))+
  geom_polygon(mapping=aes(x=quantile, y=value_loess), color="black")+
  #geom_col(mapping=aes(x=zero, y=zero_percentiles), width=0.05, fill="black") +
  geom_text(mapping=aes(x=0.5, y=0.1, label=sig_text)) +
  scale_fill_gradient2("R2",  low = "#0571b0",
                       mid = "white",
                       high = "#ca0020" )+
  facet_grid(dataset~variable_new2, switch="y")+
  theme_cowplot()+
  ylab("Genomic attribute quantile")+
  xlab("Mutation burden quantile")+
  theme(
    axis.ticks=element_blank(),
    axis.text=element_blank(), 
    #axis.text=element_text(size=8), 
    panel.border=element_blank(),
    panel.background=element_blank(),
    #strip.text.x = element_text(angle=90, hjust = 0),
    strip.text = element_text(size=8),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    axis.line=element_blank(),
    #axis.text=element_text(size=8),
    axis.title=element_text(size=10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8), 
    legend.key.size = unit(0.35, 'cm'),
    panel.spacing = unit(0.2, "lines")
    #strip.text = element_text(margin = margin(0, 0,0,0, "cm"))
    #plot.margin = unit(c(0, 0, 0, 0), "cm") )
  )
ggsave(file="S9S1bloodreg_polygonperattribute_NmutQuant_attributeQuant_reptimingGm_loess5.pdf", width=8.25, height=2)
#ggsave(file="S9S1bloodreg_polygonperattribute_NmutQuant_attributeQuant_reptimingGm_loess5.png", width=8, height=3.5)







###############################################################
####### Alternative plots:
## One large plot for SBS9 vs replication timing
## Two sets of barplots, height being R2, per significant factor in SBS9 and SBSblood
library("tidymv")

## read in data S9
indgam_S9 = read.table(file="individual_gam.S9_majority_props_reptimingGm.txt", header=T, stringsAsFactors = F)
s9_majority_props = read.table(file="s9_majority_props_reptimingGm.txt", header=T, stringsAsFactors=F)
# stepAIC_res
load(file="S9majority_gam_reptimingGm_fits.RData")
s9_ans = ans
## X1 (S1 + Sblood)
# stepAIC_res_Sblood
s1_majority_props = read.table(file="s1_majority_props_reptimingGm.txt", header=T, stringsAsFactors = F)
indgam_s1 = read.table(file="individual_gam.s1_majority_props_reptimingGm.txt", header=T, stringsAsFactors = F)
load(file = 's1majority_gam_reptimingGm_fits.RData')
s1_ans = ans

hsc_majority_props = read.table(file="hsc_props_reptimingGm.txt", header=T, stringsAsFactors = F)
indgam_hsc = read.table(file="individual_gam.hsc_majority_props_reptimingGm.txt", header=T, stringsAsFactors = F)
load(file = 'hscmajority_gam_reptimingGm_fits.RData')
hsc_ans = ans

mygam = gam(sumS9~s(rep_timing_Gm), data=s9_majority_props)
model_p = predict_gam(mygam)
model_p %>%
  ggplot(aes(rep_timing_Gm,fit)) +
  geom_jitter(data=s9_majority_props, aes(y=sumS9, x=rep_timing_Gm), width=0, height=0.5, alpha=0.03, size=0.2)+
  geom_smooth_ci(col="#7a0177", ci_alpha=0, lwd=1, lty=1)+
  xlab("Replication timing")+
  ylab("# SBS9 mutations")+
  theme_bw()+
  ylim(c(0,5))+
  coord_flip()+
  border(color = "black")+
  theme(panel.grid = element_blank(),
        axis.text=element_text(size=8, color="black"),
        axis.title=element_text(size=10))
#ggsave("reptiming_SBS9_scatter_jitter5_gam.pdf", width=3.5, height=2.5)
ggsave("reptiming_SBS9_scatter_jitter5_gam_small.pdf", width=3, height=2.5)

## Non-jitter- not used
model_p %>%
  ggplot(aes(rep_timing_Gm,fit)) +
  geom_jitter(data=s9_majority_props, aes(y=sumS9, x=rep_timing_Gm), width=0, height=0,  alpha=0.03, size=0.2)+
  geom_smooth_ci(col="#7a0177", ci_alpha=0, lwd=1, lty=1)+
  xlab("Replication timing")+
  ylab("# SBS9 mutations")+
  theme_bw()+
  ylim(c(0,5))+
  coord_flip()+
  theme(panel.grid = element_blank(),
        axis.text=element_text(size=8, color="black"),
        axis.title=element_text(size=10))
ggsave("reptiming_SBS9_scatter_gam.pdf",  width=3.5, height=2.5)


## barplots of R2 for factors significant per full gam
mylabels = c(
  "rep_timing_Gm"="Replication timing",
  "cpg_islands_dist_log10"="CpG island", 
  "gc_content_value" = "GC content", 
  "LAD_dens_1e6" = "LAD density", 
  "gene_dens_1e6" = "Gene density", 
  "TAD_b_dist_log10" = "TAD boundaries",
  "telomere_dist_log10" = "Telomere", 
  "centromere_dist_log10" = "Centromere", 
  "cruciform_inverted_rep_dens_3e3" = "Cruciform repeats", 
  "L1_rep_dist_log10" = "L1 repeats", 
  "SIMPLE_REPEAT_rep_dist_log10" = "Simple repeat", 
  "DNA_rep_dist_log10"= "DNA replication", 
  "ALU_rep_dist_log10" = "ALU repeats", 
  "MIR_rep_dist_log10" = "MIR repeats", 
  "g4_dist_log10" = "G-quadruplex", 
  "LTR_rep_dist_log10" = "LTR repeats", 
  "triplex_mirror_rep_dist_log10" = "Triplex mirror",
  "DNase" = "DNase" ,
  "H3K27ac" = "H3K27ac",
  "H3K4me1" = "H3K4me1",
  "H3K9me3" = "H3K9me3" ,
  "recomb_rate_nearest_value" = "Recomb. rate",
  "z_dna_motif_dist_log10" = "Z-DNA motifs",
  "L2_rep_dist_log10" = "L2 repeats"
  )


indgam_S9$fullgam_sig = FALSE
indgam_S9$fullgam_sig[indgam_S9$category %in% s9_ans$variables_keep] = TRUE
indgam_S9$category = factor(indgam_S9$category, levels=indgam_S9[order(indgam_S9$R2, decreasing=TRUE),]$category )
indgam_S9$category = revalue(indgam_S9$category, mylabels )

#gplot_sbs9 = 
ggplot( subset(indgam_S9, fullgam_sig==TRUE), aes(y=R2, x=category) )+  ## 18 for SBS9
  geom_col(width=0.7)+
  xlab("")+
  ylab(expression(paste("R"^2, sep="")) )+
  theme_bw()+
  #expand_limits(y = 0)+
  scale_y_continuous(limits = c(0,max(indgam_S9$R2)) )+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1, size=6, color="black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x=element_text(size=6, color="black"),
        axis.text.y=element_text(size=8, color="black"),
        axis.title=element_text(size=10),
        panel.border = element_blank(), 
        axis.line.x = element_line(), 
        axis.line.y = element_line())
#ggsave("ind_gamR2_fullgamsig_SBS9_barplot.pdf", width=4, height=1.5)
ggsave("ind_gamR2_fullgamsig_SBS9_barplot_short.pdf", width=3.1, height=1.5)

gplot_sbs9 = 
ggplot( subset(indgam_S9, fullgam_sig==TRUE), aes(y=R2, x=category) )+  ## 18 for SBS9
  geom_col(width=0.8)+
  xlab("")+
  ylab(expression(paste("R"^2, sep="")) )+
  theme_bw()+
  #expand_limits(y = 0)+
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(indgam_S9$R2)) )+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1, size=6, color="black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x=element_text(size=6, color="black"),
        axis.text.y=element_text(size=8, color="black"),
        axis.title=element_text(size=10),
        panel.border = element_blank(), 
        axis.line.x = element_line(), 
        axis.line.y = element_line())
#ggsave("ind_gamR2_fullgamsig_SBS9_barplot.pdf", width=4, height=1.5)

  
# add 2 rows of dummy variables to sbsblood, as it has 5 fewer features significant in the model but I want the plot to be the same size
indgam_s1$fullgam_sig = FALSE  ## 16 for sbsBlood/S1
indgam_s1$fullgam_sig[indgam_s1$category %in% s1_ans$variables_keep] = TRUE
indgam_s1$category = factor(indgam_s1$category, levels=indgam_s1[order(indgam_s1$R2, decreasing=TRUE),]$category )

mydummy = data.frame(paste("dummy", 1:2, sep=""), rep(0,times=2), rep(0,times=2), rep(0,times=2), rep(0,times=2), rep(TRUE,times=2))
colnames(mydummy) = colnames(indgam_s1)
indgam_s12 = rbind(indgam_s1[order(indgam_s1$R2, decreasing=TRUE),], mydummy)
indgam_s12$category = factor(indgam_s12$category, levels=indgam_s12$category)
#indgam_s12$category = factor(indgam_s12$category, levels=indgam_s12[order(indgam_s12$R2, decreasing=TRUE),]$category )
indgam_s12$category = revalue(indgam_s12$category, mylabels)

#gplot_sbs1 = 
ggplot( subset(indgam_s12, fullgam_sig==TRUE), aes(y=R2, x=category) )+ 
      geom_col()+ 
      xlab("")+ 
      ylab(expression(paste("R"^2, sep="")) )+ 
  theme_bw()+ 
  #scale_y_continuous(expand = c(0, 0)  )+ 
  #scale_y_continuous(limits = c(0,max(indgam_S9$R2)), breaks=c(0.01,0.05,0.1, 0.15) )+ 
  scale_y_continuous(limits = c(0,max(indgam_S9$R2)))+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1, size=6, color="black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=8, color="black"), 
        axis.title=element_text(size=10), 
        panel.border = element_blank(),  
        axis.line.x = element_line(),  
        axis.line.y = element_line()) 
#ggsave("ind_gamR2_fullgamsig_SBSblood_barplot.pdf", width=4, height=1.5)
ggsave("ind_gamR2_fullgamsig_SBSblood_barplot_short.pdf", width=3.1, height=1.5)
  
  
  
# add 4 rows of dummy variables to hsc, as it has 4 fewer features significant in the model but I want the plot to be the same size
indgam_hsc$fullgam_sig = FALSE  ## 16 for sbsBlood/hsc
indgam_hsc$fullgam_sig[indgam_hsc$category %in% hsc_ans$variables_keep] = TRUE
indgam_hsc$category = factor(indgam_hsc$category, levels=indgam_hsc[order(indgam_hsc$R2, decreasing=TRUE),]$category )
  
  mydummy = data.frame(paste("dummy", 1:4, sep=""), rep(0,times=4), rep(0,times=4), rep(0,times=4), rep(0,times=4), rep(TRUE,times=4))
  colnames(mydummy) = colnames(indgam_hsc)
  indgam_hsc2 = rbind(indgam_hsc[order(indgam_hsc$R2, decreasing=TRUE),], mydummy)
  indgam_hsc2$category = factor(indgam_hsc2$category, levels=indgam_hsc2$category)
  #indgam_hsc2$category = factor(indgam_hsc2$category, levels=indgam_hsc2[order(indgam_hsc2$R2, decreasing=TRUE),]$category )
  indgam_hsc2$category = revalue(indgam_hsc2$category, mylabels)
  
gplot_hsc = 
  ggplot( subset(indgam_hsc2, fullgam_sig==TRUE), aes(y=R2, x=category) )+ 
    geom_col()+ 
    xlab("")+ 
    ylab(expression(paste("R"^2, sep="")) )+ 
    theme_bw()+ 
    scale_y_continuous(limits = c(0,max(indgam_S9$R2)))+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1, size=6, color="black"),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y=element_text(size=8, color="black"), 
          axis.title=element_text(size=10), 
          panel.border = element_blank(),  
          axis.line.x = element_line(),  
          axis.line.y = element_line()) 
#ggsave("ind_gamR2_fullgamsig_HSC_barplot.pdf", width=4, height=1.5)
  
plot_grid(gplot_sbs9+ggtitle("SBS9")+theme(plot.title = element_text(size = 10)), gplot_sbs1+ggtitle("SBSblood")+theme(plot.title = element_text(size = 10)), gplot_hsc+ggtitle("HSPC")+theme(plot.title = element_text(size = 10)), nrow=3)
ggsave("ind_gamR2_fullgamsig_SBS9_SBSblood_HSC_barplot.pdf", width=4, height=6)
