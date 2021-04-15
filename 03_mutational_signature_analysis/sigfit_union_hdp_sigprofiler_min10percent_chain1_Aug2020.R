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


############## Read in pcawg focal7 June2020 matrix
mutpcawg = read.table("/lustre/scratch116/casm/cgp/users/hm8/pcawg/mutcounts_matrix_pcawgfocal7_mm_Aug2020.txt", stringsAsFactors = T, header=T)
samppcawg = read.table("/lustre/scratch116/casm/cgp/users/hm8/pcawg/samplesgroups_mutcounts_matrix_pcawgfocal7_mm_Aug2020.txt", stringsAsFactors = T, header=F)
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
#mutcounts_matrix_both = mutcounts_matrix_both[,1:10]
#samp_type = samp_type[1:10,]
#samp_names = samp_names[1:10]


##### adding blood signature to cosmic signatures (they are in the same order)
bloodsig = read.table("/lustre/scratch116/casm/cgp/users/hm8/lymphocyteWGS/metaAnalysis_pcawg_mm_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg_Aug2020/sigfit_cosmic3_bloodsig_Aug2020.txt", stringsAsFactors = F, header=T)
cosmic_signatures_v3_blood = rbind(cosmic_signatures_v3, bloodsig$Signature.Blood)
rownames(cosmic_signatures_v3_blood)[68] = "Signature.blood"


##### selecting the signatures that are the union of hdp and sigprofiler results
hdp_sigprofiler = read.table(file="signatures_union_hdp_sigprofiler_min10percent.txt", header=F, stringsAsFactors = F)
selectsigs = cosmic_signatures_v3_blood[rownames(cosmic_signatures_v3_blood) %in% hdp_sigprofiler[,1],]
#[1] "SBS1"            "SBS6"            "SBS7a"           "SBS8"
#[5] "SBS9"            "SBS17b"          "SBS18"           "SBS36"
#[9] "SBS56"           "Signature.blood"
## 10 signatures


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
# save(samp_type, mcmc_samples_fit, file=paste(groupname, ".sigfit_union_hdp_sigprofiler_min10percent.chain1.Rdata", sep="") )
