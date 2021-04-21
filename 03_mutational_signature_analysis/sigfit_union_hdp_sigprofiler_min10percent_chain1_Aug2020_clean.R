############################################################################################
## File: sigfit_union_hdp_sigprofiler_min10percent_chain1_Aug2020_clean.R
## Project: lymphocyte_somatic_mutation
## Description: sigfit per-genome signature attribution: signatures min10%
##
## Date: April 2021
## Author: Heather Machado
############################################################################################

library(sigfit)
data("cosmic_signatures_v3")

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
groupname="pcawg_mm_lymph_hsc"

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
mutpcawg = read.table("../data/mutcounts_matrix_pcawgfocal7_mm_Aug2020.txt", stringsAsFactors = T, header=T)
samppcawg = read.table("../data/samplesgroups_mutcounts_matrix_pcawgfocal7_mm_Aug2020.txt", stringsAsFactors = T, header=F)
colnames(samppcawg) = c("samp_names","celltype")
samppcawg$group = samppcawg$celltype
samppcawg$Nmut = apply(mutpcawg, MARGIN=2, FUN=sum)

# exclude those with more than 100K mutations:
keepsamples = which(samppcawg$Nmut < 60000 & samppcawg$celltype %in% c("Lymph-BNHL","Lymph-CLL","mm","Myeloid-AML") )
mutpcawg2 = mutpcawg[,keepsamples]
samppcawg2 = samppcawg[keepsamples,]


####### Combine two datasets
mutcounts_matrix_both = cbind(mutcounts_matrix, mutpcawg2)
samp_type = rbind(samp_typeA, samppcawg2)
samp_names = samp_type$samp_names


##### adding blood signature to cosmic signatures (they are in the same order)
bloodsig = read.table("../data/sigfit_cosmic3_bloodsig_Aug2020.txt", stringsAsFactors = F, header=T)

cosmic_signatures_v3_blood = rbind(cosmic_signatures_v3, bloodsig$Signature.Blood)
rownames(cosmic_signatures_v3_blood)[68] = "Signature.blood"


##### selecting the signatures that are the union of hdp and sigprofiler results
hdp_sigprofiler = read.table(file="../data/signatures_union_hdp_sigprofiler_min10percent.txt", header=F, stringsAsFactors = F)
selectsigs = cosmic_signatures_v3_blood[rownames(cosmic_signatures_v3_blood) %in% hdp_sigprofiler[,1],]


###### fit cosmic + blood signatures
mcmc_samples_fit = fit_signatures(counts = t(mutcounts_matrix_both),
    signatures = selectsigs,
    iter = 20000,
    warmup = 2000,
    chains = 1,
    seed = 1113)


## refit to abundant signatures
exposures <- retrieve_pars(mcmc_samples_fit, par = "exposures", hpd_prob = 0.90)
save(samp_type, exposures, file=paste("exposures_",groupname, ".sigfit_union_hdp_sigprofiler_min10percent.chain1.Rdata", sep="") )
