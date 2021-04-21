############################################################################################
## File: mutsig_byregion_hdp_denovo_Oct2020_clean.R
## Project: lymphocyte_somatic_mutation
## Description: hdp denovo extraction of the SHM signature
##
## Date: April 2021
## Author: Heather Machado
############################################################################################


#setwd("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/mutsig_byregion/mutsig_byregion_hdp_denovo")
library("GenomicRanges")
library("Rsamtools")
library("MASS")
library(hdp)
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)


################################################################################################
################ Run hdp (with aid priors (and cosmic?) priors)
################################################################################################
groupname = "by_1MB"

## Mutation matrix of memory B cells per 1MB bin
load(file="../data/mutmats_byregion_bin1000000.Rdata")

## Choose dataset to use
my_prior="noprior"
name_sigs = "memB"
mutmax = mutmax
mutmaxinfo = mutmaxinfo
samp_type = data.frame(mutmaxinfo, cellType="Memory B")
samp_type$ppindex = 1
samp_type$cpindex = 1+1

####### Run the program
### All samples
hdp_mut <- hdp_init(ppindex = c(0, samp_type$ppindex),  # index of parental nodes
                    cpindex = c(1, samp_type$cpindex), # index of concentration param
                    hh=rep(1, 96), # prior is uniform over the 96 mutation categories
                    alphaa=rep(1, max(samp_type$cpindex)), # shape hyperparams for different CPs
                    alphab=rep(1, max(samp_type$cpindex))) # rate hyperparams for different CPs

hdp_mut <- hdp_setdata(hdp_mut,
                       dpindex = 2:numdp(hdp_mut),   # index of nodes to add data to
                       t(mutmax) )   ## flip so samples are in rows

samp_type$dpindex = 2:numdp(hdp_mut)



# Run multiple posterior sampling chains
chlist <- vector("list", 4)

for (i in 1:4) {
  hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10, seed=i*200)
  chlist[[i]] <- hdp_posterior(hdp_activated,
                               burnin=20000,
                               n=50,
                               space=200,
                               cpiter=3,
                               seed=i*1e3)
}


luad_multi <- hdp_multi_chain(chlist)
save(luad_multi, samp_type, file=paste("luad_multi_first.", name_sigs, ".", my_prior,".",groupname, ".Rdata", sep="") )


par(mfrow=c(2,2)) 
p1 <- lapply(chains(luad_multi), plot_lik, bty='L', start=10) # check for convergence (if not, increase burnin)- 10K is good
p2 <- lapply(chains(luad_multi), plot_numcluster, bty='L')  # flat- why?
p3 <- lapply(chains(luad_multi), plot_data_assigned, bty='L') # error: 1:xmax : NA/NaN argument- why?

# Extract components (mutational signatures)
mut_colours <- c(RColorBrewer::brewer.pal(10, 'Paired')[seq(1,10,2)], 'grey70')

luad_multi <- hdp_extract_components(luad_multi)
dp_distn = comp_dp_distn(luad_multi)
comp_trinuc = comp_categ_distn(luad_multi)
Nsig = ncol(dp_distn$mean)

par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
plot_comp_size(luad_multi, bty="L")

trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))

Ncolumns = ceiling(sqrt(Nsig))
Nrows = ceiling(Nsig/Ncolumns)

pdf(file=paste("figures/signatures",name_sigs,"_prior", my_prior, "_", groupname, "_componentsSpectra.pdf", sep=""), width=Ncolumns*6, height=Nrows*3)
  par(mfrow=c(Nrows,Ncolumns))
  par(mar=c(3,3,3,1))
  plot_comp_distn(luad_multi, cat_names=trinuc_context,
                  grouping=group_factor, col=mut_colours,
                  col_nonsig="grey80", show_group_labels=TRUE)
dev.off()
  
## Plot exposures
pdf(file=paste("figures/signatures",name_sigs,"_prior", my_prior, "_", groupname, "_exposures.pdf", sep=""), width=7, height=5)
plot_dp_comp_exposure(luad_multi, samp_type$dpindex, incl_nonsig = F,incl_numdata_plot=F,
                      col=c('black', RColorBrewer::brewer.pal(8, "Set1")))
dev.off()


################# Extracting components for plotting (sig_prop per sample, CI's)
inds = samp_type$dpindex
nsites = apply(mutmax, MARGIN=2, FUN=sum)
samp_type$Nmut = nsites
  
sig_prop = data.frame(dp_distn$mean[samp_type$dpindex, ])
sig_prop$samples = paste(samp_type[,1], samp_type[,2], sep="_")
Nsamp = length(sig_prop$samples)

dp_distn_lowerCI = data.frame(matrix(ncol=ncol(dp_distn$cred.int[[1]]), nrow=Nsamp))
dp_distn_upperCI = data.frame(matrix(ncol=ncol(dp_distn$cred.int[[1]]), nrow=Nsamp))
for (i in inds){
  focal = dp_distn$cred.int[[i]]
  dp_distn_lowerCI[i-(min(inds)-1), ] = focal[1,] 
  dp_distn_upperCI[i-(min(inds)-1), ] = focal[2,] 
}
dp_distn_lowerCI$samples = sig_prop$samples
dp_distn_upperCI$samples = sig_prop$samples
colnames(dp_distn_lowerCI) = colnames(sig_prop)
colnames(dp_distn_upperCI) = colnames(sig_prop)
s1 = melt(sig_prop[,2:ncol(sig_prop)], id=c("samples") )
s2 = melt(dp_distn_lowerCI[,2:ncol(sig_prop)], id=c("samples") )
s3 = melt(dp_distn_upperCI[,2:ncol(sig_prop)], id=c("samples") )
sig_propCI = cbind(s1, s2$value, s3$value)
colnames(sig_propCI) = c("samples","variable","freq","lower","upper")
sig_propCI$samples = factor(sig_propCI$samples, levels(as.factor(sig_propCI$samples)))
sig_prop$Nmut = samp_type$Nmut  # 230 bins had no mutations



############## Calculate cosine similarity with known signatures
pcawg.sigs <- read.csv(paste('../data/sigProfiler_SBS_signatures_2019_05_22.csv', sep=""), header=TRUE)
#  sort by Substitution Type and Trinucleotide
pcawg.sigs <- pcawg.sigs[order(pcawg.sigs$Type, pcawg.sigs$SubType),]
prior_sigsP1 <- as.matrix(pcawg.sigs[,grep('SB', colnames(pcawg.sigs))])
# number of prior signatures to condition on (65)
# each col (signatures) must sum to 1 (over all mutation types)
colsum1 = apply(prior_sigsP1, MARGIN=2, FUN=sum)
prior_sigsP = prior_sigsP1
for (i in 1:length(colsum1)){
  prior_sigsP[,i] = prior_sigsP1[,i]/colsum1[i]
}

#### adding signatures of IGH/L region from CLL samples to check similarity
francAID = read.table("../data/c-AID_maura.txt", stringsAsFactors = F, header=T)

## blood sig
bloodsig = read.table("../data/sigfit_cosmic3_bloodsig_Aug2020.txt" , stringsAsFactors = F, header=T)
selectsigs = data.frame(prior_sigsP, Signature.Ig=francAID$x, SBSblood = bloodsig$Signature.Blood)


library("lsa")
cosmic_cosine = matrix(ncol=Nsig, nrow = ncol(selectsigs))
for (i in 1:Nsig){
  cosmic_cosine[,i] = apply(selectsigs, MARGIN=2, FUN=function(X) cosine(X, comp_trinuc[[1]][i,]) )  
}
rownames(cosmic_cosine) = colnames(selectsigs)


pdf(file=paste("figures/signatures",name_sigs,"_prior", my_prior, "_", groupname, "_cosine_sim_pcawg_cAID.pdf", sep=""), width=Ncolumns*1.8, height=Nrows*1.5)
par(mfrow=c(Nrows,Ncolumns))
par(mar=c(3,3,3,1))
for (i in 1:Nsig){
  hist(cosmic_cosine[,i], main=paste("C:",row.names(comp_trinuc[[1]])[i], " cos=",signif(max(cosmic_cosine[,i]), 2), " ", colnames(selectsigs)[which(cosmic_cosine[,i]==max(cosmic_cosine[,i]))], sep=""), xlab="cosine sim.", col="grey", breaks=15)
} 
dev.off()


####################### Write signatures
sigIDs = c("X0", "X1","S9_cos88","S9_cos94", "X4","SIg")
for (i in 1:Nsig){
  write.table(comp_trinuc[[1]][i,], file=paste("comp_trinuc_", sigIDs[i], ".txt", sep=""),  quote=F, col.names=F, row.names=F)  
}
write.table(t(comp_trinuc[[1]]), file=paste("comp_trinuc_allsigs.txt", sep=""),  quote=F, col.names=T, row.names=F)  

######################## Plotting proportions
sig_prop_type = data.frame(samples=sig_prop$samples, Nmut=sig_prop$Nmut, Cell.type2 = "B Memory", sig_prop[,c(1:Nsig)] )
colnames(sig_prop_type) = c("samples", "Nmut", "Cell.type2", "X0", "X1","S9_cos88","S9_cos94", "X4","SIg")
sig_prop_type$S9 = sig_prop_type$S9_cos88 + sig_prop_type$S9_cos94
save(sig_prop_type, sig_propCI, file=paste("sig_prop.", name_sigs, ".", my_prior,".",groupname, ".Rdata", sep=""))


