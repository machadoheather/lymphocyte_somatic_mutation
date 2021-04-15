# ---
#   title: "run_hdp_lustre_pcawg_lymph_hsc"
# ---
# Feb 2019

#setwd('/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_pcawg_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg_June2020')

library(sigfit)
#data("cosmic_signatures_v2")
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




###### load in ARG380 data
#setwd("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg")
groupname="pcawg_mm_lymph_hsc"

## We know that the colony_info file and the mutcounts_matrix samples are in the same order.
# Otherwise, must correctly order
mutcounts_matrix = read.table(file="/lustre/scratch116/casm/cgp/users/hm8/lymphocyteWGS/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/mutcounts_matrix_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
colonyinfo_all = read.table(file="/lustre/scratch116/casm/cgp/users/hm8/lymphocyteWGS/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/colonyinfo_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")

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
myinds = which(colonyinfo_all$Cell.type2 == "HSC" & colonyinfo_all$Donor %in% c("ARG","KX001","KX002","KX003") )  # 73 HSCs
mutcounts_matrix_both = mutcounts_matrix[,myinds]
samp_type = samp_typeA[myinds, ]
samp_names = samp_type$samp_names


####### testing
#mutcounts_matrix_both = mutcounts_matrix_both[,1:10]
#samp_type = samp_type[1:10,]
#samp_names = samp_names[1:10]


##### adding blood signature to cosmic signatures
#bloodsig = read.table("/lustre/scratch116/casm/cgp/users/hm8/lymphocyteWGS/sigprofiler_analysis/bloodsig_SBS96D.txt", stringsAsFactors = F, header=T)
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
