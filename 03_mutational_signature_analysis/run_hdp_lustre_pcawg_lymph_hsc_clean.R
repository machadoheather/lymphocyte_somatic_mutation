############################################################################################
## File: run_hdp_lustre_pcawg_lymph_hsc_clean.R
## Project: lymphocyte_somatic_mutation
## Description: hdp signature extraction
##
## Date: April 2021
## Author: Heather Machado
############################################################################################

library(hdp)
library(ggplot2)
library(foreach)
library(doParallel)
library(reshape2)


location="farm"
#location="local"

if (location=="farm"){
  root="/lustre/scratch116/casm/cgp/users/hm8"
} else if (location=="local"){
  root="/Users/hm8/volumes/hm8_network"
} else {
  warning("location not set- using current directory")
  root=getwd()
}




###### load in data
#setwd("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg")
groupname="pcawg_lymph_hsc"

## We know that the colony_info file and the mutcounts_matrix samples are in the same order.
# Otherwise, must correctly order
mutcounts_matrix = read.table(file="../data/mutcounts_matrix_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
colonyinfo_all = read.table(file="../data/colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")

############# Make samp_type object
Nsamp = ncol(mutcounts_matrix)
samp_type = colonyinfo_all
samp_type$ppindex = 1
samp_type$cpindex = 1+1
samp_type$dpindex = 2:(Nsamp+1)
samp_type$samp_names = samp_type$colony
samp_typeA = samp_type[,c("samp_names","Cell.type2","Donor","Nmut_adj_as")]
colnames(samp_typeA) = c("samp_names","celltype","group","Nmut")
  
############## Read in pcawg focal7 June2020 matrix
mutpcawg = read.table(paste("../data/mutcounts_matrix_pcawgfocal7_mm_Aug2020.txt",sep=""), stringsAsFactors = T, header=T)
samppcawg = read.table(paste(root, "../data/samplesgroups_mutcounts_matrix_pcawgfocal7_mm_Aug2020.txt", sep=""), stringsAsFactors = T, header=F)
colnames(samppcawg) = c("samp_names","celltype")
samppcawg$group = samppcawg$celltype
samppcawg$Nmut = apply(mutpcawg, MARGIN=2, FUN=sum)

# exclude those with more than 100K mutations:
keepsamples = which(samppcawg$Nmut < 60000 & samppcawg$celltype %in% c("Lymph-BNHL","Lymph-CLL","mm","Myeloid-AML") )
#samp_names         celltype            group
#102 f9dc999f-6dde-448d-9cf1-2897ddcf7b0b       Lymph-BNHL       Lymph-BNHL
#216 04aa6b77-8074-480c-872e-a1a47afa5314    Skin-Melanoma    Skin-Melanoma
#219 0ab4d782-9a50-48b9-96e4-6ce42b2ea034    Skin-Melanoma    Skin-Melanoma
#220 0dd0718d-5ddf-4c59-8c47-0f51303daeb5    Skin-Melanoma    Skin-Melanoma
#221 108749d2-5c62-4ef1-92df-aec6941ba53b    Skin-Melanoma    Skin-Melanoma
#236 00aa769d-622c-433e-8a8a-63fb5c41ea42 ColoRect-AdenoCA ColoRect-AdenoCA
#241 0980e7fd-051d-45e9-9ca6-2baf073da4e8 ColoRect-AdenoCA ColoRect-AdenoCA
#244 14c5b81d-da49-4db1-9834-77711c2b1d38 ColoRect-AdenoCA ColoRect-AdenoCA
#245 154f80bd-984c-4792-bb89-20c4da0c08e0 ColoRect-AdenoCA ColoRect-AdenoCA
#Nmut
#102  112529
#216  744553
#219  214810
#220  159349
#221  105876
#236  234231
#241  849282
#244 2426470
#245  259867

mutpcawg2 = mutpcawg[,keepsamples]
samppcawg2 = samppcawg[keepsamples,]


####### Combine two datasets
mutcounts_matrix_both = cbind(mutcounts_matrix, mutpcawg2)
samp_type = rbind(samp_typeA, samppcawg2)
samp_names = samp_type$samp_names


####### testing
# mutcounts_matrix_both = mutcounts_matrix_both[,1:10]
# samp_type = samp_type[1:10,]
# samp_names = samp_names[1:10]



######### Running hdp without priors
####### Set the indicies for the samples
groupname="pcawg_lymph_hsc"

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
                       t(mutcounts_matrix_both[,1:ncol(mutcounts_matrix_both)]) )   ## flip so samples are in rows

samp_type$dpindex = 2:numdp(hdp_mut)



# Run multiple posterior sampling chains
cl <- makeCluster(4) 
registerDoParallel(cl)


#chlist <- vector("list", 4)

chlist = foreach(i=1:4) %dopar% {
  require(hdp)
      ###### run in parallel, remember to request enough cores when running (plus one extra?)
  hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10, seed=i*200)
  chlist_part <- hdp_posterior(hdp_activated,
                               burnin=20000,
                               n=50,
                               space=200,
                               cpiter=3,
                               seed=i*1e3)
  chlist_part
}


luad_multi <- hdp_multi_chain(chlist)
name_sigs = "noprior"
my_prior=""
save(luad_multi, samp_type, file=paste("luad_multi_first.", name_sigs, ".", groupname, ".Rdata", sep="") )



############# More extraction
luad_multi <- hdp_extract_components(luad_multi)
dp_distn = comp_dp_distn(luad_multi)
comp_trinuc = comp_categ_distn(luad_multi)
Nsig = ncol(dp_distn$mean)
trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))
inds = samp_type$dpindex
sig_prop = data.frame(dp_distn$mean[samp_type$dpindex, ])
samp_names = samp_type$samp_names
sig_prop$samples = samp_names
Nsamp = length(samp_names)
dp_distn_lowerCI = data.frame(matrix(ncol=ncol(dp_distn$cred.int[[1]]), nrow=Nsamp))
dp_distn_upperCI = data.frame(matrix(ncol=ncol(dp_distn$cred.int[[1]]), nrow=Nsamp))
for (i in inds){
  focal = dp_distn$cred.int[[i]]
  dp_distn_lowerCI[i-(min(inds)-1), ] = focal[1,] 
  dp_distn_upperCI[i-(min(inds)-1), ] = focal[2,] 
}
dp_distn_lowerCI$samples = samp_names
dp_distn_upperCI$samples = samp_names
colnames(dp_distn_lowerCI) = colnames(sig_prop)
colnames(dp_distn_upperCI) = colnames(sig_prop)

s1 = melt(sig_prop[,2:ncol(sig_prop)], id=c("samples") )
s2 = melt(dp_distn_lowerCI[,2:ncol(sig_prop)], id=c("samples") )
s3 = melt(dp_distn_upperCI[,2:ncol(sig_prop)], id=c("samples") )
sig_propCI = cbind(s1, s2$value, s3$value)
colnames(sig_propCI) = c("samples","variable","freq","lower","upper")
sig_propCI$samples = factor(sig_propCI$samples, levels(as.factor(sig_propCI$samples)))

#### Looking at the contributions of each signature to the different cell types
#sig_prop_type = merge(m3, sig_prop, by.x="colony", by.y="samples")
save(luad_multi, samp_type, sig_propCI, dp_distn, comp_trinuc, Nsig, file=paste("luad_multi.", name_sigs, ".", groupname, ".Rdata", sep="") )
save(sig_propCI,samp_type,sig_prop, file=paste("misc_samp_prop_objects.", groupname, ".Rdata", sep="") )

# ### Checking burn-in, etc
# par(mfrow=c(2,2))
# p1 <- lapply(chains(luad_multi), plot_lik, bty='L', start=10) # check for convergence (if not, increase burnin)- needs 20K
# p2 <- lapply(chains(luad_multi), plot_numcluster, bty='L')
# p3 <- lapply(chains(luad_multi), plot_data_assigned, bty='L')
# # Only needs 10K, for future reference
# 
# # Extract components (mutational signatures)
# mut_colours <- c(RColorBrewer::brewer.pal(10, 'Paired')[seq(1,10,2)], 'grey70')
# par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
# plot_comp_size(luad_multi, bty="L")
# trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
# group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
#                               each=16))
# Ncolumns = ceiling(sqrt(Nsig))
# Nrows = ceiling(Nsig/Ncolumns)
# pdf(file=paste("signatures",name_sigs,"_prior", my_prior, "_", groupname, "_componentsSpectra.pdf", sep=""), width=Ncolumns*3, height=Nrows*2)
# par(mfrow=c(Nrows,Ncolumns))
# par(mar=c(3,3,3,1))
# plot_comp_distn(luad_multi, cat_names=trinuc_context,
#                   grouping=group_factor, col=mut_colours,
#                   col_nonsig="grey80", show_group_labels=TRUE)
# dev.off()
# 
# # Save mutational spectrum for components
# comp_distn <- comp_categ_distn(luad_multi)
# write.table(t(comp_distn$mean), file=paste("mutspectrum_components.comp_distn_mean.", name_sigs, ".", groupname, ".Rdata", sep=""), col.names = T, row.names = F, quote=F)

