# ---
#   title: "run_hdp_lustre_pcawg_lymph_hsc"
# ---
# Feb 2019

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


args = commandArgs(trailingOnly = TRUE)
chrom = args[1]

###### load in ARG380 data
groupname="memoryB_1Mb_IgDenovo"

## We know that the colony_info file and the mutcounts_matrix samples are in the same order.
# Otherwise, must correctly order

# mutmax, mutmaxinfo # chr  start    end
## is the SNVs from all 74 normal memory B cells
load("../data/mutmats_byregion_bin1000000.Rdata")
mutcounts_matrix = mutmax
colonyinfo_all = mutmaxinfo

############# Make samp_type object
Nsamp = ncol(mutcounts_matrix)
samp_type = colonyinfo_all
samp_type$ppindex = 1
samp_type$cpindex = 1+1
samp_type$dpindex = 2:(Nsamp+1)
samp_type$samp_names = paste(samp_type$chr,  samp_type$start, sep="_")
samp_typeA = samp_type[,c("samp_names","chr", "start", "end")]



####### Combine two datasets

focalinds = which(mutmaxinfo$chr == chrom)
mutcounts_matrix_both = mutcounts_matrix[,focalinds]
samp_type = samp_typeA[focalinds,]

#samp_names = samp_type$samp_names


####### testing
#mutcounts_matrix_both = mutcounts_matrix_both[,1:10]
#samp_type = samp_type[1:10,]
#samp_names = samp_names[1:10]


##### adding blood signature to cosmic signatures (they are in the same order)
bloodsig = read.table("../data/sigfit_cosmic3_bloodsig_Aug2020.txt", stringsAsFactors = F, header=T)
cosmic_signatures_v3_blood = rbind(cosmic_signatures_v3, bloodsig$Signature.Blood)
rownames(cosmic_signatures_v3_blood)[68] = "Signature.blood"

##### selecting the signatures that are the union of hdp and sigprofiler results
hdp_sigprofiler = read.table(file="../data/signatures_union_hdp_sigprofiler_min10percent_memB.txt", header=F, stringsAsFactors = F)
selectsigsA = cosmic_signatures_v3_blood[rownames(cosmic_signatures_v3_blood) %in% hdp_sigprofiler[,1],]
Signature.Ig = read.table("../../mutsig_byregion_hdp_denovo/comp_trinuc_SIg.txt")
selectsigs = rbind(selectsigsA, Signature.Ig=Signature.Ig[,1])
#[1] "SBS1"            "SBS8"            "SBS9"            "SBS17b"
#[5] "SBS18"           "Signature.blood" "Signature.Ig"



###### fit cosmic + blood signatures
mcmc_samples_fit = fit_signatures(counts = t(mutcounts_matrix_both),
    signatures = selectsigs,
    iter = 10000,
    warmup = 1000,
    chains = 1,
    seed = 1113)


## refit to abundant signatures
exposures <- retrieve_pars(mcmc_samples_fit, par = "exposures", hpd_prob = 0.90)
save(samp_type, exposures, file=paste("exposures_",groupname, ".sigfit_union_hdp_sigprofiler_cosmic3.chr", chrom, ".Rdata", sep="") )
