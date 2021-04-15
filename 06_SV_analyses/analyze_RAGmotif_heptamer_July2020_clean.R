### Identifying RAG mediated deletions
# testing if we can just use the heptamer

### July, 2020
#library(tidyverse)
library(stringr)
library(magrittr)
library(rtracklayer)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
#library(dplyr)
library(cowplot)
library(dplyr)

setwd("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/brass_metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/RAGmotifs")

# 
# ## Reading in filtered data
# myfiltered = read.table("../brassfiltered_all_July2019.txt", header=T, stringsAsFactors = F, sep="\t")
# # 1142 SVs (including some un-curated ones)
# # mydeletions = myfiltered[myfiltered$svclass=="deletion",]
# # 895 deletions
# 
# 
# ## Creating a control set of sites
# # Add a random offset of between 200-2000bp 
# # Going upstream for breakpoint 1 and downstream for breakpoint2 so can still do internal/external RSS motif analysis
# mycontrol = myfiltered
# offset1 = sample(c(-2000:-200), size=nrow(mycontrol), replace = T)
# offset2 = sample(c(200:2000), size=nrow(mycontrol), replace = T)
# mycontrol$start1 = myfiltered$start1 + offset1
# mycontrol$end1 = myfiltered$end1 + offset1
# mycontrol$start2 = myfiltered$start2 + offset2
# mycontrol$end2 = myfiltered$end2 + offset2
# myfiltered = mycontrol
# 
# 
# ## Identify SVs found in multiple colonies (exact breakpoints)
# # This is important for MEME analysis (for FIMO, include all)
# myfiltered$ID = paste(myfiltered$chr1, myfiltered$start1, myfiltered$chr2, myfiltered$start2, myfiltered$svclass, sep="_")
# myfiltered_nodup = myfiltered[!(duplicated(myfiltered$ID)), ]
# #myfiltered = myfiltered_nodup
# 
# 
# ## Fetching the flanking sequences
# mygranges1 = makeGRangesFromDataFrame(myfiltered, seqnames.field="chr1", start.field="start1", end.field="end1", strand.field="strand1", keep.extra.columns=TRUE)
# mygranges2 = makeGRangesFromDataFrame(myfiltered, seqnames.field="chr2", start.field="start2", end.field="end2", strand.field="strand2", keep.extra.columns=TRUE)
# 
# # _c. Get breakpoint flanking sequence ------------------------------------
# # __i. Get 50bp offsets ---------------------------------------------------
# # Make Grange around each breakpoint with flank 50bp
# brkpt.1.ext.gr <- promoters(mygranges1, upstream = 50, downstream = 0)
# brkpt.2.ext.gr <- promoters(mygranges2, upstream = 0, downstream = 50)
# brkpt.1.int.gr <- promoters(mygranges1, upstream = 0, downstream = 50)
# brkpt.2.int.gr <- promoters(mygranges2, upstream = 50, downstream = 0)
# 
# # __ii. Get Sequence ------------------------------------------------------
# # Make list of sequences and fasta style name
# brkpt.1.ext.gr.seq <- as.list(getSeq(hs37d5, brkpt.1.ext.gr, as.character = T))
# names(brkpt.1.ext.gr.seq) <- paste0(">bp1_ext_", myfiltered$id.name, "_", myfiltered$sample)
# 
# brkpt.2.ext.gr.seq <- as.list(getSeq(hs37d5, brkpt.2.ext.gr, as.character = T))
# names(brkpt.2.ext.gr.seq) <- paste0(">bp2_ext_", myfiltered$id.name, "_", myfiltered$sample)
# 
# brkpt.1.int.gr.seq <- as.list(getSeq(hs37d5, brkpt.1.int.gr, as.character = T))
# names(brkpt.1.int.gr.seq) <- paste0(">bp1_int_", myfiltered$id.name, "_", myfiltered$sample)
# 
# brkpt.2.int.gr.seq <- as.list(getSeq(hs37d5, brkpt.2.int.gr, as.character = T))
# names(brkpt.2.int.gr.seq) <- paste0(">bp2_int_", myfiltered$id.name, "_", myfiltered$sample)
# 
# # Add fasta name to original file bp1 & bp2
# myfiltered$bp1_ext_uniID <- names(brkpt.1.ext.gr.seq)
# myfiltered$bp2_ext_uniID <- names(brkpt.2.ext.gr.seq)
# myfiltered$bp1_int_uniID <- names(brkpt.1.int.gr.seq)
# myfiltered$bp2_int_uniID <- names(brkpt.2.int.gr.seq)
# 
# ## Annotate VDJ regions in original SV file
# ighst = 106304735 # 14
# ighend = 107283226 # 14
# iglst = 22385390 # 22
# iglend = 23263607 # 22
# tcrhst = 22090055 # 14 TRA
# tcrhend = 23014042 # 14
# tcrlst = 142000819 # 7  TRB
# tcrlend = 142510972 # 7
# 
# igkst = 89160078 # 2
# igkend = 90274237 # 2
# tcrgst = 38292979 #7
# tcrgend = 38407656 #7
# tcrdst = 22907537 #14   ## inside the TRA coordinates, so will just be labelled "TRA"
# tcrdend = 22938606 #14
# 
# # class switching genes
# # 14      106053274       106054731       IGHA2
# # 14      106066403       106068064       IGHE
# # 14      106090813       106092402       IGHG4
# # 14      106109540       106111126       IGHG2
# # 14      106173505       106175001       IGHA1
# # 14      106207810       106209407       IGHG1
# # 14      106232251       106237742       IGHG3
# # 14      106304737       106312010       IGHD
# # 14      106318298       106322322       IGHM
# cs_start = 106053274
# cs_end = 106322322
# 
# VDJlocus = vector()
# for (i in 1:nrow(myfiltered)){
#   X = myfiltered[i,]
#   length(X)
#   if ( (X$chr1 == 14 & X$start1 > ighst-1000 & X$start1 < ighend+1000) | (X$chr2 == 14 & X$end2 > ighst-1000 & X$end2 < ighend+1000 ) ) {VDJlocus[i] = "igh"} else
#     if ( (X$chr1 == 22 & X$start1 > iglst-1000 & X$start1 < iglend+1000) |  (X$chr2 == 22 & X$end2 > iglst-1000 & X$end2 < iglend+1000) ) {VDJlocus[i] = "igl"} else
#       if ( (X$chr1 == 14 & X$start1 > tcrhst-1000 & X$start1 < tcrhend+1000) | (X$chr2 == 14 & X$end2 > tcrhst-1000 & X$end2 < tcrhend+1000) )  {VDJlocus[i] = "tra"} else
#         if ( (X$chr1 == 7 & X$start1 > tcrlst-1000 & X$start1 < tcrlend+1000) | (X$chr2 == 7 & X$end2 > tcrlst-1000 & X$end2 < tcrlend+1000) ){VDJlocus[i] = "trb"} else 
#           if  ( (X$chr1 == 2 & X$start1 > igkst-1000 & X$start1 < igkend+1000) | (X$chr2 == 2 & X$end2 > igkst-1000 & X$end2 < igkend+1000) ) {VDJlocus[i] = "igk"} else 
#             if ( (X$chr1 == 7 & X$start1 > tcrgst-1000 & X$start1 < tcrgend+1000) | (X$chr2 == 7 & X$end2 > tcrgst-1000 & X$end2 < tcrgend+1000) ) {VDJlocus[i] = "trg"} else 
#               if ( (X$chr1 == 14 & X$start1 > tcrdst-1000 & X$start1 < tcrdend+1000) | (X$chr2 == 14 & X$end2 > tcrdst-1000 & X$end2 < tcrdend+1000) ) {VDJlocus[i] = "trd"} else 
#                 if ( (X$chr1 == 14 & X$start1 > cs_start-1000 & X$start1 < cs_end+1000) | (X$chr2 == 14 & X$end2 > cs_start-1000 & X$end2 < cs_end+1000) ) {VDJlocus[i] = "classS"} else 
#                 {VDJlocus[i] = FALSE}
# }
# myfiltered$VDJlocus = VDJlocus
# 
# write.table(myfiltered, file="motifnames_brassfiltered_all_July2019.txt", quote=FALSE, col.names = T, row.names = F, sep="\t")
# #write.table(myfiltered, file="motifnames_brassfiltered_all_controls_July2019.txt", quote=FALSE, col.names = T, row.names = F, sep="\t")
# #save(brass.del.assem.df, brass.del.assem.df, file="brass.del.assem.df.Rdata")
# 
# 
# 
# ## Write out fasta for MEME of FIMO
# # only write unique lines
# meme.fasta <- c(rbind(c(names(brkpt.1.ext.gr.seq), names(brkpt.2.ext.gr.seq),names(brkpt.1.int.gr.seq), names(brkpt.2.int.gr.seq) ), c(unlist(brkpt.1.ext.gr.seq), unlist(brkpt.2.ext.gr.seq), unlist(brkpt.1.int.gr.seq), unlist(brkpt.2.int.gr.seq) )))
# writeLines(meme.fasta,  "20190712_RAGfullmotif_MEME_filtered_updown_split.fasta")
# #writeLines(meme.fasta,  "20190712_RAGfullmotif_MEME_filtered_updown_split_controls.fasta")
# #writeLines(meme.fasta,  "20190712_RAGfullmotif_MEME_filtered_updown_split_nocolonydup.fasta")
# 
# 
# ## MEME running intructions from Dan:
# # # Submit o MEME with max motif set to 15 same as Elli paper
# # # Download HTML, txt and XML results
# # # /Users/dl8/Documents/002_ALL_Project/002_Analysis/000_MEME/000_MEME_5.0.11532383863623-1327004615
# # # I haven't built a full parser for the results files so I download all significant motif FASTAs results files    
# # # e. MEME results ---------------------------------------------------------
# # # _a. Count matches for top motif  ----------------------------------------
# # # Read in FASTA of motif 1 matches
# # MEME_motif01 <- readLines("000_MEME/000_MEME_5.0.11532383863623-1327004615/motif_1_fasta.txt")
# # 
# # # Match uniq ID back to input file and add header which includes match poistion and whether it is reverse complement
# # 
# # brass.del.assem.df$motif1.bp1 <-
# #   unlist(lapply(brass.del.assem.df$bp1_uniID, function(x) {
# #     ifelse(length(grep(x, MEME_motif01, value = T)) == 0, NA, grep(x, MEME_motif01, value = T) )
# #   }))
# # 
# # brass.del.assem.df$motif1.bp2 <-
# #   unlist(lapply(brass.del.assem.df$bp2_uniID, function(x) {
# #     ifelse(length(grep(x, MEME_motif01, value = T)) == 0, NA, grep(x, MEME_motif01, value = T))
# #   }))
# 
# 
# ## FIMO options (website)
# # OLD:  fimo --oc . --verbosity 1 --thresh 1.0E-4 RAG-motif_hepnon_combined.meme.txt 20190128_RAGfullmotif_MEME_filtered_updown_split.fasta
# # fimo --oc . --verbosity 1 --thresh 1.0E-4 RAG-motif_hepnon_combined.meme.txt 20190712_RAGfullmotif_MEME_filtered_updown_split.fasta
# 
# ## MEME options (website)
# # Options:
# # - Zero or one occurence per sequence
# # - 15 motifs
# # - Classic mode
# # - width: 6-50
# # - sites per motif: 2-600 (don't know what this means)
# # file: 20190712_RAGfullmotif_MEME_filtered_updown_split.fasta
# 


#########################################
#####   Check for RAG using FIMO
# Using the motifs file: RAG-motif_hepnon_combined.meme.txt
myfiltered = read.table(file="motifnames_brassfiltered_all_July2019.txt", header = T, stringsAsFactors = F, sep="\t")
fimoALL = read.table("fimo_20200710_RAGheptamer_fullmotif_MEME_filtered_updown_split.txt", header=TRUE, stringsAsFactors = F, sep="\t")
#pvalue = 1e-5
pvalue = 1e-4
fimo = fimoALL[fimoALL$p.value<pvalue,]
fimo_names = unlist(lapply(fimo$sequence_name, FUN=function(X) paste(">", X, sep="")))
myfilteredcont = read.table(file="motifnames_brassfiltered_all_controls_July2019.txt", header = T, stringsAsFactors = F, sep="\t")
fimoALLcont = read.table("fimo_20200710_control_RAGheptamer_fullmotif_MEME_filtered_updown_split.txt", header=TRUE, stringsAsFactors = F, sep="\t")
fimocont = fimoALLcont[fimoALLcont$p.value<pvalue,]
fimo_names_cont = unlist(lapply(fimocont$sequence_name, FUN=function(X) paste(">", X, sep="")))

g1=ggplot(fimoALL)+
  geom_histogram(aes(p.value))+
  geom_vline(xintercept=1e-4, col="red")+
  geom_vline(xintercept=1e-5, col="blue")+
  geom_vline(xintercept=1e-6, col="purple")+
  ggtitle("Observed")+
  theme_light()

g2=ggplot(fimoALLcont)+
  geom_histogram(aes(p.value))+
  geom_vline(xintercept=1e-4, col="red")+
  geom_vline(xintercept=1e-5, col="blue")+
  geom_vline(xintercept=1e-6, col="purple")+
  ggtitle("Control")+
  theme_light()
gArrange = grid.arrange(g1, g2, nrow = 1)
#ggsave(gArrange, file="filtered20190721_RSSmotif_pvalue_hist.pdf", width=6,height=3)


#####  Retrieving names
myfiltered$bp1_ext <-
  unlist(lapply(myfiltered$bp1_ext_uniID, function(x) {
    ifelse(length(grep(x, fimo_names, value = T)) == 0, NA, grep(x, fimo_names, value = T) )
  }))
myfiltered$bp2_ext <-
  unlist(lapply(myfiltered$bp2_ext_uniID, function(x) {
    ifelse(length(grep(x, fimo_names, value = T)) == 0, NA, grep(x, fimo_names, value = T) )
  }))
myfiltered$bp1_int <-
  unlist(lapply(myfiltered$bp1_int_uniID, function(x) {
    ifelse(length(grep(x, fimo_names, value = T)) == 0, NA, grep(x, fimo_names, value = T) )
  }))
myfiltered$bp2_int <-
  unlist(lapply(myfiltered$bp2_int_uniID, function(x) {
    ifelse(length(grep(x, fimo_names, value = T)) == 0, NA, grep(x, fimo_names, value = T) )
  }))

# for control sites (names are the same as original dataset, so they match up)
myfiltered$bp1_ext_cont <-
  unlist(lapply(myfiltered$bp1_ext_uniID, function(x) {
    ifelse(length(grep(x, fimo_names_cont, value = T)) == 0, NA, grep(x, fimo_names_cont, value = T) )
  }))
myfiltered$bp2_ext_cont <-
  unlist(lapply(myfiltered$bp2_ext_uniID, function(x) {
    ifelse(length(grep(x, fimo_names_cont, value = T)) == 0, NA, grep(x, fimo_names_cont, value = T) )
  }))
myfiltered$bp1_int_cont <-
  unlist(lapply(myfiltered$bp1_int_uniID, function(x) {
    ifelse(length(grep(x, fimo_names_cont, value = T)) == 0, NA, grep(x, fimo_names_cont, value = T) )
  }))
myfiltered$bp2_int_cont <-
  unlist(lapply(myfiltered$bp2_int_uniID, function(x) {
    ifelse(length(grep(x, fimo_names_cont, value = T)) == 0, NA, grep(x, fimo_names_cont, value = T) )
  }))


###### Plotting proportions of RAG hits for various categories
myfiltered$VDJlocusTF = !(myfiltered$VDJlocus == FALSE)
myfiltered$svclass2 = myfiltered$svclass
myfiltered$svclass2[myfiltered$svclass=="tandem-duplication"] = "tandem-dup."
myfiltered$non.template.T = myfiltered$non.template != "_" & myfiltered$non.template != "." 

# original sites
myfiltered$raghit = apply(myfiltered[,c("bp1_ext","bp2_ext","bp1_int","bp2_int")], MARGIN=1, FUN=function(X) sum(is.na(X)) <4)
myfiltered$raghit2 = apply(myfiltered[,c("bp1_ext","bp2_ext","bp1_int","bp2_int")], MARGIN=1, FUN=function(X) sum(is.na(X)) <=2)
myfiltered$raghit2_int = apply(myfiltered[,c("bp1_int","bp2_int")], MARGIN=1, FUN=function(X) sum(is.na(X)) ==0)
myfiltered$raghit2_ext = apply(myfiltered[,c("bp1_ext","bp2_ext")], MARGIN=1, FUN=function(X) sum(is.na(X)) ==0)
myfiltered$raghit_int = apply(myfiltered[,c("bp1_int","bp2_int")], MARGIN=1, FUN=function(X) sum(is.na(X)) <=1)
myfiltered$raghit_ext = apply(myfiltered[,c("bp1_ext","bp2_ext")], MARGIN=1, FUN=function(X) sum(is.na(X)) <=1)
myfiltered$raghit_intext = myfiltered$raghit_int & myfiltered$raghit_ext
myfiltered$raghit_intonly =  myfiltered$raghit_int & !myfiltered$raghit_ext
myfiltered$raghit_extonly =  !myfiltered$raghit_int & myfiltered$raghit_ext

# control sites
myfiltered$raghit_cont = apply(myfiltered[,c("bp1_ext_cont","bp2_ext_cont","bp1_int_cont","bp2_int_cont")], MARGIN=1, FUN=function(X) sum(is.na(X)) <4)
myfiltered$raghit2_cont = apply(myfiltered[,c("bp1_ext_cont","bp2_ext_cont","bp1_int_cont","bp2_int_cont")], MARGIN=1, FUN=function(X) sum(is.na(X)) <=2)
myfiltered$raghit2_int_cont = apply(myfiltered[,c("bp1_int_cont","bp2_int_cont")], MARGIN=1, FUN=function(X) sum(is.na(X)) ==0)
myfiltered$raghit2_ext_cont = apply(myfiltered[,c("bp1_ext_cont","bp2_ext_cont")], MARGIN=1, FUN=function(X) sum(is.na(X)) ==0)
myfiltered$raghit_int_cont = apply(myfiltered[,c("bp1_int_cont","bp2_int_cont")], MARGIN=1, FUN=function(X) sum(is.na(X)) <=1)
myfiltered$raghit_ext_cont = apply(myfiltered[,c("bp1_ext_cont","bp2_ext_cont")], MARGIN=1, FUN=function(X) sum(is.na(X)) <=1)
myfiltered$raghit_intext_cont = myfiltered$raghit_int_cont & myfiltered$raghit_ext_cont
myfiltered$raghit_intonly_cont =  myfiltered$raghit_int_cont & !myfiltered$raghit_ext_cont
myfiltered$raghit_extonly_cont =  !myfiltered$raghit_int_cont & myfiltered$raghit_ext_cont


## Writing results table
write.table(myfiltered, file="results_p10e4_obs_control_fimo_RAGheptamer_fullmotif_20200710_brassfiltered_all_July2020.txt", quote=F, col.names = T, row.names = F, sep="\t")


#########################################
##### RE-START HERE
#########################################
#### RAG for VDJ and non-VDJ, by SV class
myfiltered = read.table(file="/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/brass_metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/RAGmotifs/results_p10e4_obs_control_fimo_RAGheptamer_fullmotif_20200710_brassfiltered_all_July2020.txt", header = T, sep="\t") # 1142

### for internal VDJ hits (likely all TRUE positives), here is the distribution of hepatmers
# CACAGTG CACTGTG CACATTG CACAGCC CACAGTA CACAGTC CACAATG CACAGCG CACTCTG 
# 635     105      91      63      50      49      47      45      33 
# CACACTG CACAGAG CACAGCA CACGGTG AACTGTC CACAGAA CACTATA CACTGTA CACCCTG 
# 30      30      21      19      18      17      16      11      10 
# CACTATG CACAATC CACTCTA CACTGCA CACACCG CACAGGG CACGGCC CACAGAC AACTGTT 
# 10       9       9       8       7       7       7       5       4 
# CACCATG CACGGCT CACAAAA CACAGGA CACCATA CACTTTC CACACAG CACACTC CACCGTG 
# 4       4       3       3       3       3       2       2       1 
# CATAGCA CTCTCCA 
# 1       1 
## only half the time does it have the canonical heptamer


### For the MISSED VDJ, are there motifs? Use only T cells (b cells have CSR)
# Read in m3 table and use to include only colonies that made it into the main analysis
m3 = read.table("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/colonyinfo_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg.txt", sep="\t", header=T) # 734
unique(myfiltered$sample)[!(unique(myfiltered$sample) %in% m3$colony)]
# [1] B11_G2            BMH57             PD40521c          PD40521d          PD40521e         
# [6] PD40521g          PD40521te,PDv37is PD40667vr         PD43974ay         PD43974az        
# [11] PD43974bf         PD43974bn         PD43974cp         T1_A3 
# Remove the 14 colonies in SV/RAG analysis NOT included in main analysis (filtered out)
myfiltered2 = subset(myfiltered, colony %in% subset(m3, Cell.type2 %in% c("HSC","B Memory","T Memory","B Naive","T Naive"))$colony)   #  1118
myfiltered2$IDsample = paste(myfiltered2$chr1, myfiltered2$start1, myfiltered2$start2, myfiltered2$sample, sep="_")
#myfiltered_nodup = myfiltered[!(duplicated(myfiltered$ID)),]
#myfiltered2 = merge(myfiltered, m3, by.x="sample", by.y="colony") 
#write.table(myfiltered2, file="/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/brass_metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/RAGmotifs/results_p10e4_obs_control_fimo_RAGheptamer_fullmotif_20200710_brassfiltered_nodup_analyzed_info_July2020.txt", col.names = T, row.names = F, quote=F, sep="\t") # 1142
#myfiltered_type = merge(myfiltered_nodup, subset(m3, Cell.type2 %in% c("HSC","B Memory","T Memory","B Naive","T Naive","Treg"))[,c("colony","Cell.type2")], by="colony")   # 1074  116
#falseneg = subset(myfiltered_type, raghit==FALSE & VDJlocusTF==TRUE & svclass == "deletion" & Cell.type2 %in% c("T Memory","T Naive") )  ## 16
# myseqs = read.table(file="20190712_RAGfullmotif_MEME_filtered_updown_split.fasta", header = F, stringsAsFactors = F, sep="\t")
# myseqs2 = matrix(myseqs[,1], ncol=2, byrow=T)
# 
# falsenegID = c(as.character(falseneg$bp1_int_uniID), as.character(falseneg$bp2_int_uniID))
# falsenegseq = myseqs2[myseqs2[,1] %in% falsenegID,]
# grep(falsenegseq[,2], pattern="CACAGTG") # 1  2  3  4  7 11 13 15
# grep(falsenegseq[,2], pattern="CACTGTG") # 24 27
# grep(falsenegseq[,2], pattern="CACATTG") # 
# grep(falsenegseq[,2], pattern="CTCTCCA") # 


# remove False Positives discovered by manual curation (but does not remove IBD)
brass_ascat_curated = read.csv("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/brass_metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/brass_ascat_curated.csv")
brass_ascat_curated$IDsample = paste(brass_ascat_curated$chr1, brass_ascat_curated$start1, brass_ascat_curated$start2, brass_ascat_curated$sample, sep="_")
myfilteredTP = subset(myfiltered2, !(IDsample %in% brass_ascat_curated$IDsample[brass_ascat_curated$Correct==FALSE]) )  #  996
myfilteredTPinfo = merge(myfilteredTP, m3[,c("colony","Cell.type2")], by="colony", all.y=TRUE)

# What percent of lymphocytes have an off-target RAG?
rag_by_colony = subset(myfilteredTPinfo, VDJlocusTF == FALSE | is.na(VDJlocusTF)) %>% 
  group_by(colony) %>% 
    summarise(raghit.T = sum(raghit==TRUE, na.rm=TRUE),
              raghit.F = sum(raghit==FALSE, na.rm=TRUE),
              n = n())
mean(rag_by_colony$raghit.T > 0)
# 0.09705882


## including only the TP and non-dup (unique mutational events)
brassfilteredTPnodupNonASCAT = read.table(file="../brassfilteredTPnodupNonASCAT_July2020.txt", header=T, sep="\t", stringsAsFactors = F)
ragfinal = merge(myfiltered2, brassfilteredTPnodupNonASCAT[,c("IDsample", "Correct","Comment")], by="IDsample")  #  986 out of 1019 (33 missing - colonies without any SVs?- CHECK- only need to include those for SOME measurements)
## already checked- these are just lymphocytes

sum(brassfilteredTPnodupNonASCAT$sample %in% subset(m3, Cell.type2 %in% c("HSC","B Memory","T Memory","B Naive","T Naive"))$colony )  # 1019- no, these are only real SVs, so none should be missing.
tmpallSV = subset(brassfilteredTPnodupNonASCAT, sample %in% subset(m3, Cell.type2 %in% c("HSC","B Memory","T Memory","B Naive","T Naive"))$colony)
myfiltered2$ID2 = paste(myfiltered2$chr1, myfiltered2$start1, myfiltered2$start2, sep="_")
sum(tmpallSV$ID %in% myfiltered2$ID2)
# perhaps some were missing in the brass file, but that existed in the curated one (TRUE! Peter added some variants...)
# Consider re-running analyses using the manually curated SV list of 1019


## read in the background rate
mybackground = read.table("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/genomic_controls/RAG/results_RAGmotif_genomiccontrol.txt", header=T, stringsAsFactors = F)


## perform bootstrap resampling (for VDJ, non-VDJ, non-VDJ B and non-VDJ T)
bootstrap_summary = function(x, n){
  # x: sv results table
  # n: number of bootstrap re-samplings
  #x=ragfinal
  outlist = list()
  for (i in 1:n){
    newsubset = x[sample(1:nrow(x), replace = TRUE, size=nrow(x)),]
    m.raghit = mean(newsubset$raghit==TRUE)
    m.raghit_int = mean(newsubset$raghit_int==TRUE)
    m.raghit_intonly = mean(newsubset$raghit_intonly==TRUE)
    outlist[[i]] = c(m.raghit, m.raghit_int, m.raghit_intonly)
  }
  outDF = data.frame(do.call(rbind, outlist))
  colnames(outDF) = c("m.raghit", "m.raghit_int", "m.raghit_intonly")
  out95CI = apply(outDF, MARGIN=2, FUN=quantile, probs=c(0.025, 0.5, 0.975) )
  list(out95CI, outDF)
}
bootstrap_VDJ = bootstrap_summary(x=subset(ragfinal, VDJlocusTF==TRUE), n=1000)
bootstrap_nonVDJ = bootstrap_summary(x=subset(ragfinal, VDJlocusTF==FALSE), n=1000)
bootstrap_nonVDJ_B = bootstrap_summary(x=subset(ragfinal, VDJlocusTF==FALSE & sample %in% subset(m3, Cell.type2 %in% c("B Memory","B Naive"))$colony ), n=1000)  # 37 SVs
bootstrap_nonVDJ_T = bootstrap_summary(x=subset(ragfinal, VDJlocusTF==FALSE & sample %in% subset(m3, Cell.type2 %in% c("T Memory","T Naive"))$colony ), n=1000) #113 SVs


# Summary figures
ragfinal %>% 
   group_by(raghit, VDJlocusTF) %>% 
   summarise(nontemplate.T = sum(non.template.T==TRUE),
             nontemplate.F = sum(non.template.T==FALSE),
       m.nontemplate = mean(non.template.T==TRUE),
       n = n())
# raghit VDJlocusTF nontemplate.T nontemplate.F m.nontemplate     n
# <lgl>  <lgl>              <int>         <int>         <dbl> <int>
#   1 FALSE  FALSE                 24            90         0.211   114
# 2 FALSE  TRUE                  19            15         0.559    34
# 3 TRUE   FALSE                 16            20         0.444    36
# 4 TRUE   TRUE                 704            98         0.878   802

# for non-Ig/TCR , proportion of nontemplate in raghit vs nonraghit
fisher.test(matrix(c(24,90,16,20), nrow=2))
# p-value = 0.008943
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.1401806 0.8055242
# sample estimates:
#   odds ratio 
# 0.336165 

# for Ig/TCR, proportion of nontemplate in raghit vs nonraghit
fisher.test(matrix(c(704,98,19,15), nrow=2))
# p-value = 7.579e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   2.581727 12.178294
# sample estimates:
#   odds ratio 
# 5.652196 


gd <- ragfinal %>% 
  group_by(svclass2, VDJlocusTF) %>% 
  summarise(m.raghit = mean(raghit),
            m.raghit_int = mean(raghit_int),
            m.raghit_intonly = mean(raghit_intonly),
            n.raghit = sum(raghit),
            sd.raghit = sd(raghit),
            nontemplate = sum(non.template.T==TRUE),
            m.nontemplate = mean(non.template.T==TRUE),
            n = n())
# svclass2      VDJlocusTF m.raghit n.raghit sd.raghit     n
# <fct>         <lgl>         <dbl>    <int>     <dbl> <int>
# 1 deletion      FALSE         0.276       21     0.450    76
# 2 deletion      TRUE          0.966      714     0.181   739
# 3 inversion     FALSE         0.286        4     0.469    14
# 4 inversion     TRUE          0.948       73     0.223    77
# 5 tandem-dup.   FALSE         0.242        8     0.435    33
# 6 tandem-dup.   TRUE          0.818        9     0.405    11
# 7 translocation FALSE         0.111        3     0.320    27
# 8 translocation TRUE          0.667        6     0.5       9
gd_cont <- ragfinal %>% 
  group_by(svclass2, VDJlocusTF) %>% 
  summarise(m.raghit_cont = mean(raghit_cont),
            m.raghit_int_cont = mean(raghit_int_cont),
            m.raghit_intonly_cont = mean(raghit_intonly_cont),
            n.raghit = sum(raghit_cont),
            sd.raghit_cont = sd(raghit_cont),
            n = n())
# svclass2      VDJlocusTF m.raghit_cont n.raghit sd.raghit_cont     n
# <fct>         <lgl>              <dbl>    <int>          <dbl> <int>
# 1 deletion      FALSE             0.145        11          0.354    76
# 2 deletion      TRUE              0.120        89          0.326   739
# 3 inversion     FALSE             0.0714        1          0.267    14
# 4 inversion     TRUE              0.169        13          0.377    77
# 5 tandem-dup.   FALSE             0.0606        2          0.242    33
# 6 tandem-dup.   TRUE              0.0909        1          0.302    11
# 7 translocation FALSE             0.148         4          0.362    27
# 8 translocation TRUE              0.111         1          0.333     9
gd$control = "obs"
gd_cont$control = "control"
colnames(gd_cont) = colnames(gd)
gd_all = rbind(gd, gd_cont)

gd_total <- ragfinal %>% 
  group_by(VDJlocusTF) %>% 
  summarise(m.raghit = mean(raghit),
            m.raghit_int = mean(raghit_int),
            m.raghit_intonly = mean(raghit_intonly),
            n.raghit = sum(raghit),
            sd.raghit = sd(raghit),
            n = n())
# VDJlocusTF m.raghit n.raghit sd.raghit     n
# <lgl>         <dbl>    <int>     <dbl> <int>
#   1 FALSE         0.24        36     0.429   150
# 2 TRUE          0.959      802     0.198   836
gd_cont_total <- ragfinal %>% 
  group_by(VDJlocusTF) %>% 
  summarise(m.raghit_cont = mean(raghit_cont),
            m.raghit_int_cont = mean(raghit_int_cont),
            m.raghit_intonly_cont = mean(raghit_intonly_cont),
            n.raghit = sum(raghit_cont),
            sd.raghit_cont = sd(raghit_cont),
            n = n())
# VDJlocusTF m.raghit_cont n.raghit sd.raghit_cont     n
# <lgl>              <dbl>    <int>          <dbl> <int>
#   1 FALSE              0.12        18          0.326   150
# 2 TRUE               0.124      104          0.330   836
gd_total$control = "obs"
gd_cont_total$control = "control"
colnames(gd_cont_total) = colnames(gd_total)
gd_all_total = rbind(gd_total, gd_cont_total)
gd_all_total$svclass2 = "total"
gd_all_total = gd_all_total[,colnames(gd)]

gd_all_both = rbind(data.frame(gd_all), data.frame(gd_all_total))
gd_all_both$n.nonraghit = gd_all_both$n - gd_all_both$n.raghit

myclasses = unique(gd_all_both$svclass2)
mychi2 = mychi2vdj = vector()
for (i in 1:length(myclasses)){
  focalclass = subset(gd_all_both, svclass2==myclasses[i])
  mychi2[i] = chisq.test(matrix(c(subset(focalclass,VDJlocusTF==FALSE & control=="obs")$n.raghit,
               subset(focalclass,VDJlocusTF==FALSE & control=="obs")$n.nonraghit,
               subset(focalclass,VDJlocusTF==FALSE & control=="control")$n.raghit,
               subset(focalclass,VDJlocusTF==FALSE & control=="control")$n.nonraghit),ncol=2) )$p.value
  mychi2vdj[i] = chisq.test(matrix(c(subset(focalclass,VDJlocusTF==TRUE & control=="obs")$n.raghit,
                                  subset(focalclass,VDJlocusTF==TRUE & control=="obs")$n.nonraghit,
                                  subset(focalclass,VDJlocusTF==TRUE & control=="control")$n.raghit,
                                  subset(focalclass,VDJlocusTF==TRUE & control=="control")$n.nonraghit),ncol=2) )$p.value
}
chi2pvalues = data.frame(svclass = myclasses, mychi2, mychi2vdj)
#write.table(chi2pvalues, file="chi2pvalues_p10e4.txt", quote=F)
# svclass     mychi2     mychi2vdj
# 1      deletion 0.04453074 1.239396e-232
# 2     inversion 1.00000000  1.025133e-21
# 3   tandem-dup. 0.05142226  3.992470e-04
# 4 translocation 1.00000000  5.311619e-02
# 5         total 0.02151843 4.230487e-257






## By cell type (nonVDJ)
gd <- subset(ragfinal, VDJlocusTF==FALSE)  %>% 
  group_by(Cell.type2) %>% 
  summarise(m.raghit = mean(raghit),
            n.raghit = sum(raghit),
            n = n(),
            n.nonraghit = n-n.raghit)
# Cell.type2 m.raghit n.raghit     n
#   1 B Memory      0.258        8    31
# 2 B Naive       0.167        1     6
# 3 T Memory      0.159        7    44
# 4 T Naive       0.290       20    69
chisq.test(matrix(c(8+7, 23+37,1+20,5+49), ncol=2)) # mem rag, mem non-rag, naive rag, naive non-rag
# p=0.3391


## Just the RAG hit for VDJ or non-VDJ
gd_all_total$mygroup = c("non-Ig/TCR","Ig/TCR","non-Ig/TCR","Ig/TCR")
gd_all_total$sd = sqrt(gd_all_total$m.raghit*(1-gd_all_total$m.raghit)/gd_all_total$n)

# ggplot(subset(gd_all_total, control=="obs"), aes(mygroup, m.raghit, fill=mygroup))+
#   geom_col(position = position_dodge())+
#   geom_errorbar(aes(ymin=m.raghit-2*sd, ymax=m.raghit+2*sd),col="grey", width=0.5 , size=0.3,   position = position_dodge2())+
#   geom_col(data=subset(gd_all_total, control=="control"), fill="grey", aes(mygroup, m.raghit), position = position_dodge())+
#   scale_fill_manual("",values=c("#1b9e77","#e7298a"), guide=F)+
#   theme_light()+
#   xlab("")+
#   ylab("RAG motif associated (prop.)")+
#   coord_flip()+
#   theme(panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       strip.background = element_blank(),
#       panel.border = element_rect(colour = "black"))
# ggsave("fimoRAGheptamer_fullmotif_20200710_brassfiltered_barplot_horizontal_Jul2020.pdf", width=8,height=2)      


# with genomic background and bootstrap CI
bootstrap_VDJ = bootstrap_summary(x=subset(ragfinal, VDJlocusTF==TRUE), n=1000)
bootstrap_nonVDJ = bootstrap_summary(x=subset(ragfinal, VDJlocusTF==FALSE), n=1000)
bootstrap_vdj_nonvdj = rbind( c(bootstrap_nonVDJ[[1]][1,],bootstrap_nonVDJ[[1]][3,]),c(bootstrap_VDJ[[1]][1,],bootstrap_VDJ[[1]][3,])) 
colnames(bootstrap_vdj_nonvdj) = c("m.raghit.CIupper", "m.raghit_int.CIupper", "m.raghit_intonly.CIupper",  "m.raghit.CIlower", "m.raghit_int.CIlower", "m.raghit_intonly.CIlower")
obs_bootstrap = cbind(subset(gd_all_total, control=="obs"), bootstrap_vdj_nonvdj)

# ggplot(obs_bootstrap, aes(mygroup, m.raghit, fill=mygroup))+
#   geom_col(position = position_dodge())+
#   geom_errorbar(aes(ymin=m.raghit.CIlower, ymax=m.raghit.CIupper),col="black", width=0.5 , size=0.3,   position = position_dodge2())+
#   geom_hline(yintercept = median(mybackground$m.raghit), lty=2)+
#   scale_fill_manual("",values=c("#1b9e77","#e7298a"), guide=F)+
#   theme_light()+
#   xlab("")+
#   ylab("RAG motif associated (prop.)")+
#   coord_flip()+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_blank(),
#         panel.border = element_rect(colour = "black"))
# ggsave("fimoRAGheptamer_fullmotif_20200710_TPnodupNonASCAT_genomicbackground_bootstrap_barplot_horizontal_Jul2020.pdf", width=5,height=1.5)      

ggplot(obs_bootstrap, aes(mygroup, m.raghit, fill=mygroup))+
  geom_col(width=0.7, position = position_dodge())+
  geom_errorbar(aes(ymin=m.raghit.CIlower, ymax=m.raghit.CIupper),col="black", width=0.3 , size=0.3,   position = position_dodge2())+
  geom_hline(yintercept = median(mybackground$m.raghit), lty=2)+
  scale_fill_manual("",values=c("#df65b0","#005a32"), guide=F)+
  theme_bw()+
  xlab("")+
  ylab("RAG motif (prop.)")+
  scale_y_continuous(expand=c(0,0), limits=c(0,1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("fimoRAGheptamer_fullmotif_20200710_TPnodupNonASCAT_genomicbackground_bootstrap_barplot_Jan2021.pdf", width=1.3,height=2.6)  


### VDJ/nonVDJ for at least one internal RAG motif
ggplot(obs_bootstrap, aes(mygroup, m.raghit_int, fill=mygroup))+
  geom_col(position = position_dodge())+
  geom_errorbar(aes(ymin=m.raghit_int.CIlower, ymax=m.raghit_int.CIupper),col="black", width=0.5 , size=0.3,   position = position_dodge2())+
  geom_hline(yintercept = median(mybackground$m.raghit_int), lty=2)+
  scale_fill_manual("",values=c("#1b9e77","#e7298a"), guide=F)+
  theme_light()+
  xlab("")+
  ylab("RAG motif internal to bp (prop.)")+
  coord_flip()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
ggsave("fimoRAGheptamer_fullmotif_internal_20200710_TPnodupNonASCAT_genomicbackground_bootstrap_barplot_horizontal_Jul2020.pdf", width=5,height=1.5)  
  

## Proportion of SVs with RSS motif
#g1=
  ggplot(gd, aes(x = svclass2, y=m.raghit, group=VDJlocusTF, fill=VDJlocusTF) )+
  geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
  geom_text(aes(y=0.03, label = n) ,  position=position_dodge(width=1), size=2) + 
  scale_fill_manual(name="",values=c("TRUE"="lightgrey", "FALSE"="darkgrey"), labels=c("TRUE"="VDJ region", "FALSE"="non-VDJ region"))+
  ylab("RSS motif (prop.)")+
  #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
  ylim(c(0,1))+
  xlab("SV class")+
  theme_light()+
  ggtitle("Observed")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("fimoRAGheptamer_fullmotif_20200710_brassfiltered_barplot_Jul2020.pdf", width=7,height=5)


ggplot(gd_cont, aes(x = svclass2, y=m.raghit, group=VDJlocusTF, fill=VDJlocusTF) )+
  geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
  geom_text(aes(y=0.03, label = n) ,  position=position_dodge(width=1), size=2) + 
  scale_fill_manual(name="",values=c("TRUE"="lightgrey", "FALSE"="darkgrey"), labels=c("TRUE"="VDJ region", "FALSE"="non-VDJ region"))+
  ylab("RSS motif (prop.)")+
  #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
  ylim(c(0,1))+
  xlab("SV class")+
  theme_light()+
  ggtitle("Control")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("fimo_control_RAGheptamer_fullmotif_20200710_brassfiltered_barplot_Jul2020.pdf", width=7,height=5)


#### Proportion of SVs with RSS motif- obs and control
ggplot(subset(gd_all, VDJlocusTF==FALSE), aes(x = svclass2, y=m.raghit, group=control, fill=control) )+
  geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
  geom_text(aes(y=0.01, label = n) ,  position=position_dodge(width=1), size=2) + 
  scale_fill_manual(name="",values=c("control"="lightgrey", "obs"="darkgrey"), labels=c("control"="Control", "obs"="Observed"))+
  ylab("RSS motif (prop.)")+
  #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
  #ylim(c(0,1))+
  xlab("SV class")+
  theme_light()+
  ggtitle("Non-VDJ")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("fimoRAGheptamer_fullmotif_control_obs_brassfiltered_nonVDJ_barplot_Jul2020.pdf", width=5,height=4)

ggplot(subset(gd_all, VDJlocusTF==TRUE), aes(x = svclass2, y=m.raghit, group=control, fill=control) )+
  geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
  geom_text(aes(y=0.04, label = n) ,  position=position_dodge(width=1), size=2) + 
  scale_fill_manual(name="",values=c("control"="lightgrey", "obs"="darkgrey"), labels=c("control"="Control", "obs"="Observed"))+
  ylab("RSS motif (prop.)")+
  #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
  #ylim(c(0,1))+
  xlab("SV class")+
  theme_light()+
  ggtitle("VDJ")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("fimoRAGheptamer_fullmotif_control_obs_brassfiltered_VDJ_barplot_Jul2020.pdf", width=5,height=4)






#### Proportion of SVs with RSS motif by cell subset- obs and control

bootstrap_nonVDJ_B = bootstrap_summary(x=subset(ragfinal, VDJlocusTF==FALSE & sample %in% subset(m3, Cell.type2 %in% c("B Memory","B Naive"))$colony ), n=1000)  # 37 SVs
bootstrap_nonVDJ_T = bootstrap_summary(x=subset(ragfinal, VDJlocusTF==FALSE & sample %in% subset(m3, Cell.type2 %in% c("T Memory","T Naive"))$colony ), n=1000) #113 SVs


gd_celltype <- myfiltered_type %>% 
  group_by(svclass2, VDJlocusTF,Cell.type2) %>% 
  summarise(m.raghit = mean(raghit),
            n.raghit = sum(raghit),
            sd.raghit = sd(raghit),
            n = n())
gd_celltype_cont <- myfiltered_type %>% 
  group_by(svclass2, VDJlocusTF,Cell.type2) %>% 
  summarise(m.raghit_cont = mean(raghit_cont),
            n.raghit = sum(raghit_cont),
            sd.raghit_cont = sd(raghit_cont),
            n = n())
gd_celltype$control = "obs"
gd_celltype_cont$control = "control"
colnames(gd_celltype_cont) = colnames(gd_celltype)
gd_celltype_all = rbind(gd_celltype, gd_celltype_cont)

ggplot(subset(gd_celltype_all, VDJlocusTF==FALSE), aes(x = svclass2, y=m.raghit, group=control, fill=control) )+
  geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
  facet_wrap(.~Cell.type2)+
  geom_text(aes(y=0.01, label = n) ,  position=position_dodge(width=1), size=2) + 
  scale_fill_manual(name="",values=c("control"="lightgrey", "obs"="darkgrey"), labels=c("control"="Control", "obs"="Observed"))+
  ylab("RSS motif (prop.)")+
  #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
  #ylim(c(0,1))+
  xlab("SV class")+
  theme_light()+
  ggtitle("Non-VDJ")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("fimoRAGheptamer_fullmotif_control_obs_allSV_celltype_nonVDJ_barplot_Jul2020.pdf", width=5,height=4)

ggplot(subset(gd_celltype_all, VDJlocusTF==TRUE), aes(x = svclass2, y=m.raghit, group=control, fill=control) )+
  geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
  geom_text(aes(y=0.04, label = n) ,  position=position_dodge(width=1), size=2) + 
  scale_fill_manual(name="",values=c("control"="lightgrey", "obs"="darkgrey"), labels=c("control"="Control", "obs"="Observed"))+
  facet_wrap(.~Cell.type2)+
  ylab("RSS motif (prop.)")+
  #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
  #ylim(c(0,1))+
  xlab("SV class")+
  theme_light()+
  ggtitle("VDJ")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("fimoRAGheptamer_fullmotif_control_obs_allSV_celltype_VDJ_barplot_Jul2020.pdf", width=5,height=4)



# #### Proportion of RAG deletions with RSS internal to deletion 
# myfiltered = ragfinal  ## using the 
# 
# gd_rag_intext <- myfiltered %>% 
#   group_by(svclass2, VDJlocusTF, raghit) %>% 
#   summarise(m.raghit_int = mean(raghit_intonly),
#             m.raghit_ext = mean(raghit_extonly),
#             m.raghit_intext = mean(raghit_intext),
#             sd.raghit_int = sd(raghit_intonly),
#             sd.raghit_ext = sd(raghit_extonly),
#             sd.raghit_intext = sd(raghit_intext),
#             n = n())
# 
# gd_rag_intext_total <- myfiltered %>% 
#   group_by(VDJlocusTF, raghit) %>% 
#   summarise(m.raghit_int = mean(raghit_intonly),
#             m.raghit_ext = mean(raghit_extonly),
#             m.raghit_intext = mean(raghit_intext),
#             sd.raghit_int = sd(raghit_intonly),
#             sd.raghit_ext = sd(raghit_extonly),
#             sd.raghit_intext = sd(raghit_intext),
#             n = n())
# 
# gd_rag_intext_cont <- myfiltered %>% 
#   group_by(svclass2, VDJlocusTF, raghit_cont) %>% 
#   summarise(m.raghit_int_cont = mean(raghit_intonly_cont),
#             m.raghit_ext_cont = mean(raghit_extonly_cont),
#             m.raghit_intext_cont = mean(raghit_intext_cont),
#             sd.raghit_int_cont = sd(raghit_intonly_cont),
#             sd.raghit_ext_cont = sd(raghit_extonly_cont),
#             sd.raghit_intext_cont = sd(raghit_intext_cont),
#             n = n())
# 
# gd_rag_intext_total_cont <- myfiltered %>% 
#   group_by(VDJlocusTF, raghit_cont) %>% 
#   summarise(m.raghit_int_cont = mean(raghit_intonly_cont),
#             m.raghit_ext_cont = mean(raghit_extonly_cont),
#             m.raghit_intext_cont = mean(raghit_intext_cont),
#             sd.raghit_int_cont = sd(raghit_intonly_cont),
#             sd.raghit_ext_cont = sd(raghit_extonly_cont),
#             sd.raghit_intext_cont = sd(raghit_intext_cont),
#             n = n())
# 

# ## proportion of RAG deletions with RSS motif internal only
# g1= 
#   ggplot(subset(gd_rag_intext, raghit==TRUE), aes(x = svclass2, y=m.raghit_int, group=VDJlocusTF, fill=VDJlocusTF) )+
#   geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
#   geom_text(aes(y=0.03, label = n) ,  position=position_dodge(width=1), size=2) + 
#   scale_fill_manual(name="",values=c("TRUE"="lightgrey", "FALSE"="darkgrey"), labels=c("TRUE"="VDJ region", "FALSE"="non-VDJ region"))+
#   ylab("Only internal RSS motif (prop.)")+
#   #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
#   ylim(c(0,1))+
#   xlab("SV class")+
#   theme_light()+
#   ggtitle("Observed")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# #ggsave("fimoRAG_20190712_allSV_intRSSonly_barplot_ggplotJul2019.pdf", width=7,height=5)
# g2= 
# ggplot(subset(gd_rag_intext_cont, raghit_cont==TRUE), aes(x = svclass2, y=m.raghit_int_cont, group=VDJlocusTF, fill=VDJlocusTF) )+
#   geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
#   geom_text(aes(y=0.03, label = n) ,  position=position_dodge(width=1), size=2) + 
#   scale_fill_manual(name="",values=c("TRUE"="lightgrey", "FALSE"="darkgrey"), labels=c("TRUE"="VDJ region", "FALSE"="non-VDJ region"))+
#   ylab("Only internal RSS motif (prop.)")+
#   #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
#   ylim(c(0,1))+
#   xlab("SV class")+
#   theme_light()+
#   ggtitle("Control")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# gArrange = grid.arrange(g1, g2, nrow = 1)
# ggsave(gArrange, file="fimoRAG_obs_control_20190721_allSV_intRSSonly_barplot_ggplotJul2019.pdf", width=9,height=3.5)
# #ggsave(gArrange, file="fimoRAG_obs_control_p10e05_20190721_allSV_intRSSonly_barplot_ggplotJul2019.pdf", width=9,height=3.5)
# 
# 
# 
# ## proportion of RAG deletions with RSS motif external only 
# g1=
#   ggplot(subset(gd_rag_intext, raghit==TRUE), aes(x = svclass2, y=m.raghit_ext, group=VDJlocusTF, fill=VDJlocusTF) )+
#   geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
#   geom_text(aes(y=0.03, label = n) ,  position=position_dodge(width=1), size=2) + 
#   scale_fill_manual(name="",values=c("TRUE"="lightgrey", "FALSE"="darkgrey"), labels=c("TRUE"="VDJ region", "FALSE"="non-VDJ region"))+
#   ylab("Only external RSS motif (prop.)")+
#   ylim(c(0,1))+
#   xlab("SV class")+
#   theme_light()+
#   ggtitle("Observed")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# #ggsave("fimoRAG_20190712_allSV_extRSSonly_barplot_ggplotJul2019.pdf", width=7,height=5)
# g2=
#   ggplot(subset(gd_rag_intext_cont, raghit_cont==TRUE), aes(x = svclass2, y=m.raghit_ext_cont, group=VDJlocusTF, fill=VDJlocusTF) )+
#   geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
#   geom_text(aes(y=0.03, label = n) ,  position=position_dodge(width=1), size=2) + 
#   scale_fill_manual(name="",values=c("TRUE"="lightgrey", "FALSE"="darkgrey"), labels=c("TRUE"="VDJ region", "FALSE"="non-VDJ region"))+
#   ylab("Only external RSS motif (prop.)")+
#   ylim(c(0,1))+
#   xlab("SV class")+
#   theme_light()+
#   ggtitle("Control")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# gArrange = grid.arrange(g1, g2, nrow = 1)
# ggsave(gArrange, file="fimoRAG_obs_control_20190721_allSV_extRSSonly_barplot_ggplotJul2019.pdf", width=9,height=3.5)
# #ggsave(gArrange, file="fimoRAG_obs_control_p10e05_20190721_allSV_extRSSonly_barplot_ggplotJul2019.pdf", width=9,height=3.5)
# 
# 
# ## proportion of RAG deletions with RSS motif both internal and external
# g1=
#   ggplot(subset(gd_rag_intext, raghit==TRUE), aes(x = svclass2, y=m.raghit_intext, group=VDJlocusTF, fill=VDJlocusTF) )+
#   geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
#   geom_text(aes(y=0.03, label = n) ,  position=position_dodge(width=1), size=2) + 
#   scale_fill_manual(name="",values=c("TRUE"="lightgrey", "FALSE"="darkgrey"), labels=c("TRUE"="VDJ region", "FALSE"="non-VDJ region"))+
#   ylab("Both int. and ext. RSS motif (prop.)")+
#   ylim(c(0,1))+
#   xlab("SV class")+
#   theme_light()+
#   ggtitle("Observed")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# #ggsave("fimoRAG_20190712_allSV_intextRSS_barplot_ggplotJul2019.pdf", width=7,height=5)
# g2=
#   ggplot(subset(gd_rag_intext_cont, raghit_cont==TRUE), aes(x = svclass2, y=m.raghit_intext_cont, group=VDJlocusTF, fill=VDJlocusTF) )+
#   geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
#   geom_text(aes(y=0.03, label = n) ,  position=position_dodge(width=1), size=2) + 
#   scale_fill_manual(name="",values=c("TRUE"="lightgrey", "FALSE"="darkgrey"), labels=c("TRUE"="VDJ region", "FALSE"="non-VDJ region"))+
#   ylab("Both int. and ext. RSS motif (prop.)")+
#   ylim(c(0,1))+
#   xlab("SV class")+
#   theme_light()+
#   ggtitle("Control")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# gArrange = grid.arrange(g1, g2, nrow = 1)
# ggsave(gArrange, file="fimoRAG_obs_control_20190721_allSV_intextRSS_barplot_ggplotJul2019.pdf", width=9,height=3.5)
# #ggsave(gArrange, file="fimoRAG_obs_control_p10e05_20190721_allSV_intextRSS_barplot_ggplotJul2019.pdf", width=9,height=3.5)
# 
# 
# 
# ### Non-template insertions (control not relevant here- is prop of non-temp ins for the SV)
# nt <- myfiltered %>% 
#   group_by(svclass2, VDJlocusTF, raghit) %>% 
#   summarise(m.non.template.T = mean(non.template.T),
#             sd.non.template.T = sd(non.template.T),
#             n = n())
# nt$raghit[nt$raghit==FALSE] = "no RSS"
# nt$raghit[nt$raghit==TRUE] = "RSS motif"
# nt$raghit = factor(nt$raghit, levels=c("RSS motif","no RSS") )
# 
# nt_cont <- myfiltered %>% 
#   group_by(svclass2, VDJlocusTF, raghit_cont) %>% 
#   summarise(m.non.template.T_cont = mean(non.template.T),
#             sd.non.template.T_cont = sd(non.template.T),
#             n = n())
# nt_cont$raghit_cont[nt_cont$raghit_cont==FALSE] = "no RSS"
# nt_cont$raghit_cont[nt_cont$raghit_cont==TRUE] = "RSS motif"
# nt_cont$raghit_cont = factor(nt_cont$raghit_cont, levels=c("RSS motif","no RSS") )
# 
# # Proportion of SVs with RSS motif
# ggplot(nt, aes(x = raghit, y=m.non.template.T, group=VDJlocusTF, fill=VDJlocusTF) )+
#   geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
#   geom_text(aes(y=0.03, label = n) ,  position=position_dodge(width=1), size=2) + 
#   scale_fill_manual(name="",values=c("TRUE"="lightgrey", "FALSE"="darkgrey"), labels=c("TRUE"="VDJ region", "FALSE"="non-VDJ region"))+
#   facet_wrap(~svclass2)+
#   ylab("Non-template insertion (prop.)")+
#   #geom_errorbar(aes(ymin=m.raghit-sd.raghit, ymax=m.raghit+sd.raghit), width=.2, position=position_dodge(1)) +
#   ylim(c(0,1))+
#   theme_light()+
#   #ggtitle("Observed")+
#   xlab("")
# ggsave(file="fimoRAG_obs_20190721_nontemplateI_barplot_ggplotJul2019.pdf", width=6,height=6)
# # interestingly, very high rate of non-template insertions in VDJ deletion w/o RSS motif
# # wonder if I'm missing some RSS motifs there, and they really are RAG. ~5% false neg rate? 
# 
# nt2 = nt
# nt2$VDJlocus = "VDJ"
# nt2$VDJlocus[nt2$VDJlocusTF==FALSE] = "non-VDJ"
# ggplot(nt2, aes(x = VDJlocus, y=m.non.template.T, group=raghit, fill=raghit) )+
#   geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
#   geom_text(aes(y=0.03, label = n) ,  position=position_dodge(width=1), size=2) + 
#   scale_fill_manual(name="",values=c("no RSS"="lightgrey", "RSS motif"="darkgrey"), labels=c("no RSS"="No motif", "RSS motif"="RSS motif"))+
#   facet_wrap(~svclass2)+
#   ylab("Non-template insertion (prop.)")+
#   #geom_errorbar(aes(ymin=m.raghit-sd.raghit, ymax=m.raghit+sd.raghit), width=.2, position=position_dodge(1)) +
#   ylim(c(0,1))+
#   theme_light()+
#   #ggtitle("Observed")+
#   xlab("")
# ggsave(file="fimoRAG_obs_20190721_nontemplateI_barplot_ggplotJul2019.pdf", width=6,height=6)
# 
# 
# 



# ###### What are the off-target RAG deletions
# Read in meta data to get cell type
#m3 = read.csv(file="/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/ARGhsc_ARGlymph_all380/meandp_cellnumber_colony_all380.csv",header=T, stringsAsFactors = F)
#myfiltered_m3 = merge(myfiltered, m3, by="colony")

offtarget = subset(myfiltered, raghit ==T & VDJlocus==FALSE)
# one deletion per colony (none with 2+)
# 10 T cells and one B cell (could be due to # colonies + depth)
# write.table(myfiltered_m3, file="e04_filtered20190206_myfiltered_RAG_nonRAG_SVs_ARG.txt", col.names = T, row.names = F, quote=F)
# write.table(myfiltered_m3[myfiltered_m3$VDJlocus==FALSE & myfiltered_m3$raghit==TRUE,] , file="e04_filtered20190206_nonVDJdeletions_RAG_SVs_ARG.txt", col.names = T, row.names = F, quote=F, sep="\t")


offtargetS = offtarget[order(offtarget$chr1, offtarget$start1), ]
plot()

manhattan.plot(offtargetS$chr1, offtargetS$start1, rep(1,times=nrow(offtarget)) , chromsizes)
ggsave("filtered20190206_offtargetRAGdeletions.pdf", width=7,height=4)

offtargetS2 = data.frame(offtargetS[,c(1:6,11)])

####### CODE TO RUN
chromsizes = read.table("/Users/hm8/sanger/genomes/human/GRCh37d5/hg19.chrom_sizes.txt")
manhattan.plot = function(chr, pos, height, chromsizes){
  if (length(chr) != length(pos) | length(pos) != length(height)) {
    warning("vectors not of same length")
  }
  #chr=offtargetS$X..chr1
  #pos=offtargetS$start1
  #height=rep(1,times=nrow(offtarget)) 
  # Calculate total genomic positions for each chrom
  base1_genomic = vector()
  for (i in 1:(nrow(chromsizes)-1) ){
    base1_genomic[i] = sum(as.numeric(chromsizes[1:i,3]) )
  }
  chromsizes$base1_genomic = c(0,base1_genomic)
  posgenomic = vector()
  for (i in 1:length(chr)){
    chromsizes[which(chromsizes[,1] == chr[i]), "base1_genomic"]
    posgenomic[i] = chromsizes[which(chromsizes[,1] == chr[i]),]$base1_genomic + pos[i]
  }
  df1 = data.frame(chr, pos, posgenomic, height)
  ggplot(df1, aes(posgenomic, height))+
    geom_jitter(height=0.05, alpha=0.5) +
    #geom_dotplot()+
    scale_x_continuous(name ="Genomic location", breaks = chromsizes$base1_genomic,
                       labels=chromsizes$V1) +
    ylim(0.5,1.5)+
    theme_light()+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}



# ############################## Plotting the RAG motif:
# #source("http://bioconductor.org/biocLite.R")
# #biocLite("seqLogo")
# library(seqLogo)
# # heptamer
# A = scan()
# C = scan()
# G = scan()
# T = scan()
# motif1 <- data.frame(A, C, G, T)
# pdf("RSSheptamer_logo.pdf", width=5,height=2.5)
# seqLogo( t(motif1))
# dev.off()
# 
# # nonamer
# A = scan()
# C = scan()
# G = scan()
# T = scan()
# motif2 <- data.frame(A, C, G, T)
# pdf("RSSnonamer_logo.pdf", width=5.2,height=2.5)
# seqLogo( t(motif2))
# dev.off()
# 
# # both
# A = scan()
# C = scan()
# G = scan()
# T = scan()
# motif3 <- data.frame(A, C, G, T)
# pdf("RSSheptamer_nonamer_logo.pdf", width=8,height=2.5)
# seqLogo( t(motif3))
# dev.off()







