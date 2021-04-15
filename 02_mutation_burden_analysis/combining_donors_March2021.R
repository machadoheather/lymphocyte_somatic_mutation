## July 2019
## HEM

## Combining files from each donor

# knitr::opts_chunk$set(echo = TRUE)
# knitr::opts_chunk$set(cache=TRUE, fig.path="./graphics/plot",autodep=TRUE,fig.width=5, fig.height=4)
# knitr::opts_knit$set(root.dir = '/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_AX001_KX001_KX002_KX003_TX001_TX002_CB001_March2021/')
# library(ggplot2)
# library(cowplot)
# library(caTools)
# library(reshape2)
# library(dplyr)
# library(hdp)
# library(RColorBrewer)
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



############### Combining VDJ productive/non-productive strand (and mut#) info. 
## On lustre
for project in immunerepARGtreg immunerepKX001lymph immunerepKX001lymph25 immunerepKX002lymph immunerepKX003lymph immunerepstemcellCB2lymph22 immunereptonsilMS immunereptonsilPS; do
  for region in IGH IGL TRA TRB; do
    cat ../$project/results/*/*$region.txt | grep -v "File" - | cat /lustre/scratch116/casm/cgp/users/hm8/immune_rep/VDJsequence/header.$region.txt - > results_bydonor_July2019/$project.$region.VDJsequence.txt
  done
done


# splitting productive and non-productive strands
for project in immunerepARGtreg immunerepKX001lymph immunerepKX001lymph25 immunerepKX002lymph immunerepKX003lymph immunerepstemcellCB2lymph22 immunereptonsilMS immunereptonsilPS ARG380; do
  awk '{if ($7=="Yes") print}' results_bydonor_July2019/$project.TRA.VDJsequence.txt > results_bydonor_July2019/$project.TRA.VDJsequence.prod.txt
  awk '{if ($8=="Yes") print}' results_bydonor_July2019/$project.TRB.VDJsequence.txt > results_bydonor_July2019/$project.TRB.VDJsequence.prod.txt
  awk '{if ($8=="Yes") print}' results_bydonor_July2019/$project.IGH.VDJsequence.txt > results_bydonor_July2019/$project.IGH.VDJsequence.prod.txt
  awk '{if ($7=="Yes") print}' results_bydonor_July2019/$project.IGL.VDJsequence.txt > results_bydonor_July2019/$project.IGL.VDJsequence.prod.txt

  awk '{if ($7=="No") print}' results_bydonor_July2019/$project.TRA.VDJsequence.txt > results_bydonor_July2019/$project.TRA.VDJsequence.nonprod.txt
  awk '{if ($8=="No") print}' results_bydonor_July2019/$project.TRB.VDJsequence.txt > results_bydonor_July2019/$project.TRB.VDJsequence.nonprod.txt
  awk '{if ($8=="No") print}' results_bydonor_July2019/$project.IGH.VDJsequence.txt > results_bydonor_July2019/$project.IGH.VDJsequence.nonprod.txt
  awk '{if ($7=="No") print}' results_bydonor_July2019/$project.IGL.VDJsequence.txt > results_bydonor_July2019/$project.IGL.VDJsequence.nonprod.txt
done

# combining projects
cat results_bydonor_July2019/*.TRA.VDJsequence.prod.txt  | cat /lustre/scratch116/casm/cgp/users/hm8/immune_rep/VDJsequence/header.TRA.txt - > results_bydonor_July2019/all.TRA.VDJsequence.prod.txt
cat results_bydonor_July2019/*.TRB.VDJsequence.prod.txt | cat /lustre/scratch116/casm/cgp/users/hm8/immune_rep/VDJsequence/header.TRB.txt -  > results_bydonor_July2019/all.TRB.VDJsequence.prod.txt
cat results_bydonor_July2019/*.IGH.VDJsequence.prod.txt | cat /lustre/scratch116/casm/cgp/users/hm8/immune_rep/VDJsequence/header.IGH.txt -  > results_bydonor_July2019/all.IGH.VDJsequence.prod.txt
cat results_bydonor_July2019/*.IGL.VDJsequence.prod.txt | cat /lustre/scratch116/casm/cgp/users/hm8/immune_rep/VDJsequence/header.IGL.txt -  > results_bydonor_July2019/all.IGL.VDJsequence.prod.txt

cat results_bydonor_July2019/*.TRA.VDJsequence.nonprod.txt | cat /lustre/scratch116/casm/cgp/users/hm8/immune_rep/VDJsequence/header.TRA.txt -  > results_bydonor_July2019/all.TRA.VDJsequence.nonprod.txt
cat results_bydonor_July2019/*.TRB.VDJsequence.nonprod.txt | cat /lustre/scratch116/casm/cgp/users/hm8/immune_rep/VDJsequence/header.TRB.txt -  > results_bydonor_July2019/all.TRB.VDJsequence.nonprod.txt
cat results_bydonor_July2019/*.IGH.VDJsequence.nonprod.txt | cat /lustre/scratch116/casm/cgp/users/hm8/immune_rep/VDJsequence/header.IGH.txt -  > results_bydonor_July2019/all.IGH.VDJsequence.nonprod.txt
cat results_bydonor_July2019/*.IGL.VDJsequence.nonprod.txt | cat /lustre/scratch116/casm/cgp/users/hm8/immune_rep/VDJsequence/header.IGL.txt -  > results_bydonor_July2019/all.IGL.VDJsequence.nonprod.txt


# all.IGH.VDJsequence.nonprod.txt                                                              100%   12KB 662.7KB/s   00:00    
# all.IGH.VDJsequence.prod.txt                                                                 100%   67KB   1.5MB/s   00:00    
# all.IGL.VDJsequence.nonprod.txt                                                              100% 9714   709.5KB/s   00:00    
# all.IGL.VDJsequence.prod.txt                                                                 100%   34KB   1.4MB/s   00:00    
# all.TRA.VDJsequence.nonprod.txt                                                              100%  140KB   1.6MB/s   00:00    
# all.TRA.VDJsequence.prod.txt                                                                 100%  216KB   1.6MB/s   00:00    
# all.TRB.VDJsequence.nonprod.txt                                                              100%   17KB   1.0MB/s   00:00    
# all.TRB.VDJsequence.prod.txt                                                                 100%   65KB   1.6MB/s   00:00    
# 
# in dir: /Users/hm8/sanger/immune_rep/VDJsequence/results_bydonor_July2019














########################## Filtering out shared mutations (within a donor) and then creating mutation matrices and mutation lists
library("GenomicRanges")
library("Rsamtools")
library("MASS")
setwd("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_AX001_KX001_KX002_KX003_TX001_TX002_CB001_March2021")
projects = c("KX001hsc_lymph245", "KX002hsc_lymph", "KX003hsc19_lymph", "CB001", "TX001", "TX002","AX001")
genomeFile = "/Users/hm8/sanger/genomes/human/GRCh37d5/genome.fa"


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

## Filter out sites in more than one individual (per donor)
## Create mutation matrix
# som2_all_list_nodup = list()
# for (i in 1:length(projects)){
#   som2_all_list_focal = som2_all_list[[i]]
#   som2_all_list_nodup[[i]] = som2_all_list_focal[!(duplicated(som2_all_list_focal[,1])),]
# }
# 
# mutmatrix_all_list = list()
# for (i in 1:length(projects)){
#   focal = som2_all_list[[i]]
#   allmuts = do.call(rbind, focal)
#   mytable = table(allmuts$ID)
#   filterout = names(mytable)[mytable>1]
#   
#   mutcounts_matrix = matrix(ncol=length(focal), nrow=96)
#   for (focal.sample in 1:length(focal)){
#     focal2 = focal[[focal.sample]]
#     focal_filtered = focal2[!(focal2$ID %in% filterout),]
#     subs_only = focal_filtered[,2:6]
#     colnames(subs_only) = c("chr","pos","ref","mut","freq")
#     subs_only = subs_only[(subs_only$ref %in% c("A","C","G","T")) & (subs_only$mut %in% c("A","C","G","T")) & subs_only$chr %in% c(1:22,"X","Y"),]
#     subs_only$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(subs_only$chr, IRanges(subs_only$pos-1, subs_only$pos+1))))
#     
#     # 2. Annotating the mutation from the pyrimidine base
#     ntcomp = c(T="A",G="C",C="G",A="T")
#     subs_only$sub = paste(subs_only$ref,subs_only$mut,sep=">")
#     subs_only$trinuc_ref_py = subs_only$trinuc_ref
#     for (j in 1:nrow(subs_only)) {
#       if (subs_only$ref[j] %in% c("A","G")) { # Purine base
#         subs_only$sub[j] = paste(ntcomp[subs_only$ref[j]],ntcomp[subs_only$mut[j]],sep=">")
#         subs_only$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(subs_only$trinuc_ref[j],split="")[[1]])],collapse="")
#       }
#     }
#     
#     # 3. Counting subs
#     freqs = table(paste(subs_only$sub,paste(substr(subs_only$trinuc_ref_py,1,1),substr(subs_only$trinuc_ref_py,3,3),sep="-"),sep=","))
#     sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
#     ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
#     full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
#     freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
#     xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
#     y = freqs_full; maxy = max(y)
#     mutcounts_matrix[,focal.sample] = y
#   }
#   mutmatrix_all_list[[i]] = mutcounts_matrix
# }
# 
# 
# ## Combine mutation matrices
# mutmatrix_all = do.call(cbind, mutmatrix_all_list)
# sample_names_all = unlist(sample_names_all_list)
# colnames(mutmatrix_all) = sample_names_all
# 
# 
# ## Write output
# write.table(mutmatrix_all, file="mutcounts_matrix_privateOnly_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt", col.names = T, row.names = F, quote=F, sep="\t")
# 





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
