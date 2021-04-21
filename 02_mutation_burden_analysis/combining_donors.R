############################################################################################
## File: combining_donors.R
## Project: lymphocyte_somatic_mutation
## Description: Combine datasets
##
## Date: April 2021
## Author: Heather Machado
############################################################################################


## Combining files from each donor
library(dplyr)
library(plyr)

setwd("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_AX001_KX001_KX002_KX003_TX001_TX002_CB001_March2021")
projects = c("KX001hsc_lymph245", "KX002hsc_lymph", "KX003hsc19_lymph", "CB001", "TX001", "TX002","AX001")


## Combining mutation matrices
mutmatrix_all_list = list()
sample_names_all_list = list()
colonyinfo_all_list = list()
for (i in 1:length(projects)){
  load(paste("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/", projects[i],"_filtered/", projects[i],"_mutcounts_matrix.Rdata", sep="") )
  mutmatrix_all_list[[i]] = mutcounts_matrix
  sample_names_all_list[[i]] = sample_names
  colonyinfo_all_list[[i]] = read.csv(paste("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/", projects[i],"_filtered/", "meandp_cellnumber_colony_", projects[i],"_indel.csv", sep="") )
}
colonyinfo_all_list[[2]]$Donor = "KX002"
colonyinfo_all_list[[2]]$Age = 38
colonyinfo_all_list[[2]]$Tissue = colonyinfo_all_list[[2]]$tissue
colonyinfo_all_list[[2]]$CellType = colonyinfo_all_list[[2]]$Cell.type

## Age of KX001 at sampling was 29
colonyinfo_all_list[[1]]$Age = 29

focal_cols = c("colony","meandp","Donor","CellType","Cell.type2","Tissue","Age","Project","Nmut","Nmut_adj_as","Nmut_indel","Nmut_adj_as_indel")
colonyinfo_all_list2 = list()
for (i in 1:length(projects)){
  colonyinfo_all_list2[[i]] = colonyinfo_all_list[[i]][,focal_cols]
}

mutmatrix_all = do.call(cbind, mutmatrix_all_list)
sample_names_all = unlist(sample_names_all_list)
colnames(mutmatrix_all) = sample_names_all
colonyinfo_all_tmp = do.call(rbind, colonyinfo_all_list2)
rownames(colonyinfo_all_tmp) = colonyinfo_all_tmp$colony
colonyinfo_all = colonyinfo_all_tmp[sample_names_all,]

# rename cell types (column Cell.type2)
colonyinfo_all$Cell.type2 = revalue(colonyinfo_all$Cell.type2, c("B Memory"="Memory B", "T Memory"="Memory T", "B Naive" = "Naive B", "T Naive" = "Naive T"))

# write files
write.table(mutmatrix_all, file="mutcounts_matrix_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt", col.names = T, row.names = F, quote=F, sep="\t")
write.table(colonyinfo_all, file="colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt", col.names = T, row.names = F, quote=F, sep="\t")



############### Combining telomere lengths
projects = c("KX001hsc_lymph", "KX001lymph25", "KX002hsc_lymph", "KX003hsc19_lymph", "stemcellCB2hsc_lymph22", "tonsilMS", "tonsilPS","ARGhscPB_ARGlymph380Treg")
telomere_list = list()
for (i in 1:length(projects)){
    telomere_list[[i]] = read.csv(paste("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/", projects[i],"/telomerecat_length_", projects[i],".csv", sep=""), header=T, stringsAsFactors=F)
}
telomere_list[[9]] = read.csv("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/KX003hsc_lymph/telomerecat_length_PD43974ec.csv", header=T, stringsAsFactors=F)
## PD43974ec ran separately- check that it is not anomalous

telomereDF = do.call(rbind, telomere_list)
telomereDF$colony = unlist(lapply(telomereDF[,1], FUN=function(X) unlist(strsplit(X, split=".", fixed=TRUE))[1]   ))

colonyinfo_all = read.table(file="colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt",header=TRUE, stringsAsFactors=FALSE, sep="\t")

colonyinfo_allTelo = merge(colonyinfo_all, telomereDF[,c("colony","Length")], by="colony", all.x=TRUE)
colnames(colonyinfo_allTelo)[ncol(colonyinfo_allTelo)] = "Telomere"
write.table(colonyinfo_allTelo, file="colonyinfo_telo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt", col.names = T, row.names = F, quote=F, sep="\t")






######## writing table of ALL mutations (including indels) for dnds analysis
setwd("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_AX001_KX001_KX002_KX003_TX001_TX002_CB001")
projects = c("KX001hsc_lymph245", "KX002hsc_lymph", "KX003hsc19_lymph", "CB001", "TX001", "TX002","AX001")
genomeFile = "/Users/hm8/sanger/genomes/human/GRCh37d5/genome.fa"
donors = c("KX001", "KX002", "KX003", "CB001", "TX001", "TX002","AX001")


## Read in data
som2_all_list = list()
som2_indel_all_list = list()
sample_names_all_list = list()
sample_names_indel_all_list = list()

colonyinfo_all_list = list()
for (i in 1:length(projects)){
  load(paste("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/", projects[i],"_filtered/", projects[i],"_somsample_list2.Rdata", sep="") )
  som2_all_list[[i]] = somsample_list2
  sample_names_all_list[[i]] = sample_names
  
  load(paste("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/", projects[i],"_filtered/", projects[i],"_indel_somsample_list2.Rdata", sep="") )
  som2_indel_all_list[[i]] = somsample_list2
  sample_names_indel_all_list[[i]] = sample_names
}

som2_all_list2 = list()
som2_indel_all_list2 = list()
for (i in 1:length(projects)){
  for (j in 1:length(som2_all_list[[i]])){
    if (nrow(som2_all_list[[i]][[j]])>0){
      som2_all_list[[i]][[j]]$sample = sample_names_all_list[[i]][j]  
      som2_all_list[[i]][[j]]$donor = donors[i]
    }
    if (nrow(som2_indel_all_list[[i]][[j]])>0 ){
      som2_indel_all_list[[i]][[j]]$sample = sample_names_indel_all_list[[i]][j]
      som2_indel_all_list[[i]][[j]]$donor = donors[i]
    }
    
  }
  som2_all_list2[[i]] = do.call(rbind, som2_all_list[[i]])
  som2_indel_all_list2[[i]] = do.call(rbind, som2_indel_all_list[[i]])
}
som2_all = do.call(rbind, som2_all_list2)
som2_indel_all = do.call(rbind, som2_indel_all_list2)
write.table(rbind(som2_all, som2_indel_all), file="som2_indel_SNV_allDF.txt", col.names=T, row.names=F, quote=F, sep="\t")
