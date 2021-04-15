## Aug 2020


### preparing a matrix of mutation counts for sigprofilier signature extraction
mutcounts_matrix = read.table(file="../data/mutcounts_matrix_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
colonyinfo_all = read.table(file="../data/colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
rownames(colonyinfo_all) = colonyinfo_all$colony
colonyinfo_all2 = colonyinfo_all[colnames(mutcounts_matrix),]

# ############# Make samp_type object
Nsamp = ncol(mutcounts_matrix)
samp_type = colonyinfo_all2
samp_type$ppindex = 1
samp_type$cpindex = 1+1
samp_type$dpindex = 2:(Nsamp+1)
samp_type$samp_names = samp_type$colony
samp_typeA = samp_type[,c("samp_names","Cell.type2","Donor","Nmut_hsc_as")]
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
mutcounts_matrix = cbind(mutcounts_matrix, mutpcawg2)
samp_type = rbind(samp_typeA, samppcawg2)





## Load blood signature from the sigprofiler analysis
bloodsig = read.table("../data/sigfit_cosmic3_bloodsig_Aug2020.txt", stringsAsFactors = F, header=T)

colnames(bloodsig)[2] = "SBSblood"


## load cosmic V2 signatures (pcawg)
cosmicsigs = read.table("/Users/hm8/volumes/hm8_network/reference_files/mutational_signatures/signatures_probabilities.txt", header=T, stringsAsFactors = F, sep="\t") # 2583
# put cosmicsigs in sigprofiler order
rownames(cosmicsigs) = cosmicsigs$Somatic.Mutation.Type
#allsigs = cbind(bloodsig, cosmicsigs[bloodsig$MutationsType,4:33])
#write.table(allsigs, file="cosmicV2_SBS_signatures_bloodsig_SBS96D.txt", sep="\t", quote=F, col.names = T, row.names = F)
cosmicsigs2 = cosmicsigs[order(cosmicsigs$Substitution.Type, cosmicsigs$Trinucleotide),]

## my mut mat, new order
rownames(mutcounts_matrix) = cosmicsigs2$Somatic.Mutation.Type
mutcounts_matrix_new = data.frame(MutationsType=bloodsig[,1], mutcounts_matrix[bloodsig$MutationsType,])
colnames(mutcounts_matrix_new) = c("Mutation Types", colnames(mutcounts_matrix))

#/lustre/scratch116/casm/cgp/users/hm8/
write.table(mutcounts_matrix_new, file="../data/mutcounts_matrix_pcawg_mm_lymph_hsc_sigprofilerOrder.txt", sep="\t", quote=F, col.names = T, row.names = F)

write.table(mutcounts_matrix, file="../data/mutcounts_matrix_pcawg_mm_lymph_hsc.txt", sep="\t", quote=F, col.names = T, row.names = F)

write.table(samp_type, file="../data/samp_type_pcawg_mm_lymph_hsc.txt", sep="\t", quote=F, col.names = T, row.names = F)


