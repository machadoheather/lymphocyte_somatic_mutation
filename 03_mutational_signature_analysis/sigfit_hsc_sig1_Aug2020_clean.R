############################################################################################
## File: sigfit_hsc_sig1_Aug2020_clean.R
## Project: lymphocyte_somatic_mutation
## Description: SBSblood signature extraction using sigfit
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



####### Select HSCs (exclude the stemcell ones that have too few mutations to be helpful)
myinds = which(colonyinfo_all$Cell.type2 == "HSC" & colonyinfo_all$Donor %in% c("AX001","KX001","KX002","KX003") )  # 73 HSCs
mutcounts_matrix_both = mutcounts_matrix[,myinds]
samp_type = samp_typeA[myinds, ]
samp_names = samp_type$samp_names


##### adding blood signature to cosmic signatures
cosmic_signatures_v3_blood = cosmic_signatures_v3[1,]

mcmc_samples_fit_1_extra <- fit_extract_signatures(t(mutcounts_matrix_both),
signatures = cosmic_signatures_v3_blood,
iter = 50000,
warmup = 10000,
chains = 1,
seed = 1757,
num_extra_sigs = 1)

myexposures <- retrieve_pars(mcmc_samples_fit_1_extra, par = "exposures", hpd_prob = 0.90)
#plot_all(mcmc_samples = mcmc_samples_fit_2_extra, out_path = getwd(), prefix = "test_mcmc_samples_fit_2_extra_cosmic3")
mysigs <- retrieve_pars(mcmc_samples_fit_1_extra, "signatures")

save(samp_type, myexposures, mysigs, file="exposures_hsc_S1_extra1.sigfit.Rdata")
save(samp_type, mcmc_samples_fit_1_extra, file="hsc_S1_extra1.sigfit.Rdata")
