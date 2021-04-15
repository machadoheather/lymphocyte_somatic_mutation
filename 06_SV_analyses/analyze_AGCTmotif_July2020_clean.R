### Identifying RAG mediated deletions

### July, 2019
library(tidyverse)
library(stringr)
library(magrittr)
library(rtracklayer)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(dplyr)
library(cowplot)
library(reshape2)

setwd("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/brass_metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/AGCTmotifs")


## Reading in filtered data
myfiltered = read.table("../brassfiltered_all_July2019.txt", header=T, stringsAsFactors = F, sep="\t")
# 1142 SVs (including some un-curated ones)
# mydeletions = myfiltered[myfiltered$svclass=="deletion",]
# 895 deletions


## Creating a control set of sites
# Add a random offset of between 200-2000bp 
# Going upstream for breakpoint 1 and downstream for breakpoint2 so can still do internal/external RSS motif analysis
mycontrol = myfiltered
offset1 = sample(c(-30000:-20000), size=nrow(mycontrol), replace = T)
offset2 = sample(c(20000:30000), size=nrow(mycontrol), replace = T)
mycontrol$start1 = myfiltered$start1 + offset1
mycontrol$end1 = myfiltered$end1 + offset1
mycontrol$start2 = myfiltered$start2 + offset2
mycontrol$end2 = myfiltered$end2 + offset2
myfiltered = mycontrol


## Identify SVs found in multiple colonies (exact breakpoints)
# This is important for MEME analysis (for FIMO, include all)
myfiltered$ID = paste(myfiltered$chr1, myfiltered$start1, myfiltered$chr2, myfiltered$start2, myfiltered$svclass, sep="_")
myfiltered_nodup = myfiltered[!(duplicated(myfiltered$ID)), ]
#myfiltered = myfiltered_nodup


## Fetching the flanking sequences
mygranges1 = makeGRangesFromDataFrame(myfiltered, seqnames.field="chr1", start.field="start1", end.field="end1", strand.field="strand1", keep.extra.columns=TRUE)
mygranges2 = makeGRangesFromDataFrame(myfiltered, seqnames.field="chr2", start.field="start2", end.field="end2", strand.field="strand2", keep.extra.columns=TRUE)

# _c. Get breakpoint flanking sequence ------------------------------------
# __i. Get 50bp offsets ---------------------------------------------------
# Make Grange around each breakpoint with flank 50bp
brkpt.1.ext.gr <- promoters(mygranges1, upstream = 1000, downstream = 0)
brkpt.2.ext.gr <- promoters(mygranges2, upstream = 0, downstream = 1000)
brkpt.1.int.gr <- promoters(mygranges1, upstream = 0, downstream = 1000)
brkpt.2.int.gr <- promoters(mygranges2, upstream = 1000, downstream = 0)

# __ii. Get Sequence ------------------------------------------------------
# Make list of sequences and fasta style name
brkpt.1.ext.gr.seq <- as.list(getSeq(hs37d5, brkpt.1.ext.gr, as.character = T))
names(brkpt.1.ext.gr.seq) <- paste0(">bp1_ext_", myfiltered$id.name, "_", myfiltered$sample)

brkpt.2.ext.gr.seq <- as.list(getSeq(hs37d5, brkpt.2.ext.gr, as.character = T))
names(brkpt.2.ext.gr.seq) <- paste0(">bp2_ext_", myfiltered$id.name, "_", myfiltered$sample)

brkpt.1.int.gr.seq <- as.list(getSeq(hs37d5, brkpt.1.int.gr, as.character = T))
names(brkpt.1.int.gr.seq) <- paste0(">bp1_int_", myfiltered$id.name, "_", myfiltered$sample)

brkpt.2.int.gr.seq <- as.list(getSeq(hs37d5, brkpt.2.int.gr, as.character = T))
names(brkpt.2.int.gr.seq) <- paste0(">bp2_int_", myfiltered$id.name, "_", myfiltered$sample)

# Add fasta name to original file bp1 & bp2
myfiltered$bp1_ext_uniID <- names(brkpt.1.ext.gr.seq)
myfiltered$bp2_ext_uniID <- names(brkpt.2.ext.gr.seq)
myfiltered$bp1_int_uniID <- names(brkpt.1.int.gr.seq)
myfiltered$bp2_int_uniID <- names(brkpt.2.int.gr.seq)

## Annotate VDJ regions in original SV file
ighst = 106304735 # 14
ighend = 107283226 # 14
iglst = 22385390 # 22
iglend = 23263607 # 22
tcrhst = 22090055 # 14 TRA
tcrhend = 23014042 # 14
tcrlst = 142000819 # 7  TRB
tcrlend = 142510972 # 7

igkst = 89160078 # 2
igkend = 90274237 # 2
tcrgst = 38292979 #7
tcrgend = 38407656 #7
tcrdst = 22907537 #14   ## inside the TRA coordinates, so will just be labelled "TRA"
tcrdend = 22938606 #14

# class switching genes
# 14      106053274       106054731       IGHA2
# 14      106066403       106068064       IGHE
# 14      106090813       106092402       IGHG4
# 14      106109540       106111126       IGHG2
# 14      106173505       106175001       IGHA1
# 14      106207810       106209407       IGHG1
# 14      106232251       106237742       IGHG3
# 14      106304737       106312010       IGHD
# 14      106318298       106322322       IGHM
cs_start = 106053274
cs_end = 106322322

VDJlocus = vector()
for (i in 1:nrow(myfiltered)){
  X = myfiltered[i,]
  length(X)
  if ( (X$chr1 == 14 & X$start1 > ighst-1000 & X$start1 < ighend+1000) | (X$chr2 == 14 & X$end2 > ighst-1000 & X$end2 < ighend+1000 ) ) {VDJlocus[i] = "igh"} else
    if ( (X$chr1 == 22 & X$start1 > iglst-1000 & X$start1 < iglend+1000) |  (X$chr2 == 22 & X$end2 > iglst-1000 & X$end2 < iglend+1000) ) {VDJlocus[i] = "igl"} else
      if ( (X$chr1 == 14 & X$start1 > tcrhst-1000 & X$start1 < tcrhend+1000) | (X$chr2 == 14 & X$end2 > tcrhst-1000 & X$end2 < tcrhend+1000) )  {VDJlocus[i] = "tra"} else
        if ( (X$chr1 == 7 & X$start1 > tcrlst-1000 & X$start1 < tcrlend+1000) | (X$chr2 == 7 & X$end2 > tcrlst-1000 & X$end2 < tcrlend+1000) ){VDJlocus[i] = "trb"} else 
          if  ( (X$chr1 == 2 & X$start1 > igkst-1000 & X$start1 < igkend+1000) | (X$chr2 == 2 & X$end2 > igkst-1000 & X$end2 < igkend+1000) ) {VDJlocus[i] = "igk"} else 
            if ( (X$chr1 == 7 & X$start1 > tcrgst-1000 & X$start1 < tcrgend+1000) | (X$chr2 == 7 & X$end2 > tcrgst-1000 & X$end2 < tcrgend+1000) ) {VDJlocus[i] = "trg"} else 
              if ( (X$chr1 == 14 & X$start1 > tcrdst-1000 & X$start1 < tcrdend+1000) | (X$chr2 == 14 & X$end2 > tcrdst-1000 & X$end2 < tcrdend+1000) ) {VDJlocus[i] = "trd"} else 
                if ( (X$chr1 == 14 & X$start1 > cs_start-1000 & X$start1 < cs_end+1000) | (X$chr2 == 14 & X$end2 > cs_start-1000 & X$end2 < cs_end+1000) ) {VDJlocus[i] = "classS"} else 
                {VDJlocus[i] = FALSE}
}
myfiltered$VDJlocus = VDJlocus

#write.table(myfiltered, file="motifnames_brassfiltered_all_1000bp_July2020.txt", quote=FALSE, col.names = T, row.names = F, sep="\t")
#write.table(myfiltered, file="motifnames_brassfiltered_all_controls_1000bp_July2020.txt", quote=FALSE, col.names = T, row.names = F, sep="\t")


## Write out fasta for MEME of FIMO
# only write unique lines
meme.fasta <- c(rbind(c(names(brkpt.1.ext.gr.seq), names(brkpt.2.ext.gr.seq),names(brkpt.1.int.gr.seq), names(brkpt.2.int.gr.seq) ), c(unlist(brkpt.1.ext.gr.seq), unlist(brkpt.2.ext.gr.seq), unlist(brkpt.1.int.gr.seq), unlist(brkpt.2.int.gr.seq) )))
#writeLines(meme.fasta,  "MEME_filtered_updown_split_1000bp_July2020.fasta")
#writeLines(meme.fasta,  "MEME_filtered_updown_split_1000bp_controls_July2020.fasta")


## MEME running intructions from Dan:
# # Submit o MEME with max motif set to 15 same as Elli paper
# # Download HTML, txt and XML results
# # /Users/dl8/Documents/002_ALL_Project/002_Analysis/000_MEME/000_MEME_5.0.11532383863623-1327004615
# # I haven't built a full parser for the results files so I download all significant motif FASTAs results files    
# # e. MEME results ---------------------------------------------------------
# # _a. Count matches for top motif  ----------------------------------------
# # Read in FASTA of motif 1 matches
# MEME_motif01 <- readLines("000_MEME/000_MEME_5.0.11532383863623-1327004615/motif_1_fasta.txt")
# 
# # Match uniq ID back to input file and add header which includes match poistion and whether it is reverse complement
# 
# brass.del.assem.df$motif1.bp1 <-
#   unlist(lapply(brass.del.assem.df$bp1_uniID, function(x) {
#     ifelse(length(grep(x, MEME_motif01, value = T)) == 0, NA, grep(x, MEME_motif01, value = T) )
#   }))
# 
# brass.del.assem.df$motif1.bp2 <-
#   unlist(lapply(brass.del.assem.df$bp2_uniID, function(x) {
#     ifelse(length(grep(x, MEME_motif01, value = T)) == 0, NA, grep(x, MEME_motif01, value = T))
#   }))


## FIMO options (website)
# OLD:  fimo --oc . --verbosity 1 --thresh 1.0E-4 RAG-motif_hepnon_combined.meme.txt 20190128_RAGfullmotif_MEME_filtered_updown_split.fasta
# fimo --oc . --verbosity 1 --thresh 1.0E-4 RAG-motif_hepnon_combined.meme.txt 20190712_RAGfullmotif_MEME_filtered_updown_split.fasta

## MEME options (website)
# Options:
# - Zero or one occurence per sequence
# - 15 motifs
# - Classic mode
# - width: 6-50
# - sites per motif: 2-600 (don't know what this means)
# file: 20190712_RAGfullmotif_MEME_filtered_updown_split.fasta



#########################################
#####   Check for RAG using FIMO
# Using the motifs file: RAG-motif_hepnon_combined.meme.txt
myfiltered = read.table(file="motifnames_brassfiltered_all_1000bp_July2020.txt", header = T, stringsAsFactors = F, sep="\t")
#fimoALL = read.table("mcast_obs.tsv", header=TRUE, stringsAsFactors = F, sep="\t")
## new pvalue consistent with pcwag analysis
fimoALL = read.table("/Users/hm8/volumes/hm8_network/lymphocyteWGS/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/brass_metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/AGCTmotifs/lymph/mcast.tsv", header=TRUE, stringsAsFactors = F, sep="\t")
#pvalue = 1e-5
pvalue = 1
fimo = fimoALL[fimoALL$p.value<pvalue,]
fimo_names = unlist(lapply(fimo$sequence_name, FUN=function(X) paste(">", X, sep="")))
myfilteredcont = read.table(file="motifnames_brassfiltered_all_controls_1000bp_July2020.txt", header = T, stringsAsFactors = F, sep="\t")
fimoALLcont = read.table("/Users/hm8/volumes/hm8_network/lymphocyteWGS/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/brass_metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/AGCTmotifs/lymph_controls/mcast.tsv", header=TRUE, stringsAsFactors = F, sep="\t")
fimocont = fimoALLcont[fimoALLcont$p.value<pvalue,]
fimo_names_cont = unlist(lapply(fimocont$sequence_name, FUN=function(X) paste(">", X, sep="")))

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

myfiltered = myfiltered[!(duplicated(myfiltered$ID)),]  # remove 16 lines that are duplicated

## Writing results table
#write.table(myfiltered, file="results_mcastWGCW_p01e1000space100_brassfiltered_all_1000bp_July2020.txt", quote=F, col.names = T, row.names = F, sep="\t")
#write.table(myfiltered, file="results_mcastWGCW_p01e100space50_brassfiltered_all_1000bp_July2020.txt", quote=F, col.names = T, row.names = F, sep="\t")
#write.table(myfiltered, file="results_mcastWGCW_p1_brassfiltered_all_1000bp_July2020.txt", quote=F, col.names = T, row.names = F, sep="\t")


#########################################
##### RE-START HERE
#########################################
#### RAG for VDJ and non-VDJ, by SV class
myfiltered = read.table(file="results_mcastWGCW_p1_brassfiltered_all_1000bp_July2020.txt", header = T, sep="\t") # 1126
brassfilteredTPnodupNonASCAT = read.table(file="/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/brass_metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/brassfilteredTPnodupNonASCAT_July2020.txt", header=T, sep="\t", stringsAsFactors = F)
myfiltered$IDsample = paste(myfiltered$chr1, myfiltered$start1, myfiltered$start2, myfiltered$sample, sep="_")
agctfinal = myfiltered[myfiltered$IDsample %in% brassfilteredTPnodupNonASCAT$IDsample, ]  #  986 116 (not sure 

# Read in m3 table and use to include only colonies that made it into the main analysis
# m3 = read.table("../../colonyinfo_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg.txt", sep="\t", header=T) # 734
# unique(myfiltered$sample)[!(unique(myfiltered$sample) %in% m3$colony)]
# # [1] B11_G2            BMH57             PD40521c          PD40521d          PD40521e         
# # [6] PD40521g          PD40521te,PDv37is PD40667vr         PD43974ay         PD43974az        
# # [11] PD43974bf         PD43974bn         PD43974cp         T1_A3 
# # Remove the 14 colonies in SV/RAG analysis NOT included in main analysis (filtered out)
# #myfiltered = subset(myfiltered, colony %in% subset(m3, Cell.type2 %in% c("HSC","B Memory","T Memory","B Naive","T reg"))$colony)   #  1118
# myfiltered_type = merge(myfiltered, subset(m3, Cell.type2 %in% c("HSC","B Memory","T Memory","B Naive","T Naive"))[,c("colony","Cell.type2")], by="colony")   #  1118
myfiltered_type = agctfinal


# Summary figures
gd <- myfiltered_type %>% 
  group_by(svclass2, VDJlocusTF) %>% 
  summarise(m.raghit = mean(raghit),
            n.raghit = sum(raghit),
            sd.raghit = sd(raghit),
            n = n())
gd_cont <- myfiltered_type %>% 
  group_by(svclass2, VDJlocusTF) %>% 
  summarise(m.raghit_cont = mean(raghit_cont),
            n.raghit = sum(raghit_cont),
            sd.raghit_cont = sd(raghit_cont),
            n = n())
gd$control = "obs"
gd_cont$control = "control"
colnames(gd_cont) = colnames(gd)
gd_all = rbind(gd, gd_cont)

gd_total <- myfiltered_type %>% 
  group_by(VDJlocusTF) %>% 
  summarise(m.raghit = mean(raghit),
            n.raghit = sum(raghit),
            sd.raghit = sd(raghit),
            n = n())
gd_cont_total <- myfiltered_type %>% 
  group_by(VDJlocusTF) %>% 
  summarise(m.raghit_cont = mean(raghit_cont),
            n.raghit = sum(raghit_cont),
            sd.raghit_cont = sd(raghit_cont),
            n = n())
gd_total$control = "obs"
gd_cont_total$control = "control"
colnames(gd_cont_total) = colnames(gd_total)
gd_all_total = rbind(gd_total, gd_cont_total)
gd_all_total$svclass2 = "total"
gd_all_total = gd_all_total[,colnames(gd)]
gd_all_both = rbind(data.frame(gd_all), data.frame(gd_all_total))
gd_all_both$n.nonraghit = gd_all_both$n - gd_all_both$n.raghit
#svclass2 VDJlocusTF   m.raghit n.raghit sd.raghit   n control n.nonraghit
# 1       deletion      FALSE 0.00000000        0 0.0000000  76     obs          76
# 2       deletion       TRUE 0.01768707       13 0.1319011 735     obs         722
# 3      inversion      FALSE 0.00000000        0 0.0000000  14     obs          14
# 4      inversion       TRUE 0.04054054        3 0.1985695  74     obs          71
# 5    tandem-dup.      FALSE 0.00000000        0 0.0000000  33     obs          33
# 6    tandem-dup.       TRUE 0.00000000        0 0.0000000  11     obs          11
# 7  translocation      FALSE 0.00000000        0 0.0000000  27     obs          27
# 8  translocation       TRUE 0.22222222        2 0.4409586   9     obs           7
# 9       deletion      FALSE 0.00000000        0 0.0000000  76 control          76
# 10      deletion       TRUE 0.00000000        0 0.0000000 735 control         735
# 11     inversion      FALSE 0.00000000        0 0.0000000  14 control          14
# 12     inversion       TRUE 0.00000000        0 0.0000000  74 control          74
# 13   tandem-dup.      FALSE 0.00000000        0 0.0000000  33 control          33
# 14   tandem-dup.       TRUE 0.00000000        0 0.0000000  11 control          11
# 15 translocation      FALSE 0.00000000        0 0.0000000  27 control          27
# 16 translocation       TRUE 0.00000000        0 0.0000000   9 control           9
# 17         total      FALSE 0.00000000        0 0.0000000 150     obs         150
# 18         total       TRUE 0.02171291       18 0.1458325 829     obs         811
# 19         total      FALSE 0.00000000        0 0.0000000 150 control         150
# 20         total       TRUE 0.00000000        0 0.0000000 829 control         829

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
# 1      deletion 0.02624899 2.592287e-255
# 2     inversion 0.74278686  2.233821e-21
# 3   tandem-dup. 0.67181432  3.992470e-04
# 4 translocation 1.00000000  2.669060e-02
# 5         total 0.01535235 2.052748e-279

## Proportion of SVs with RSS motif
#g1=
  ggplot(gd, aes(x = svclass2, y=m.raghit, group=VDJlocusTF, fill=VDJlocusTF) )+
  geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
  geom_text(aes(y=0.03, label = n) ,  position=position_dodge(width=1), size=2) + 
  scale_fill_manual(name="",values=c("TRUE"="lightgrey", "FALSE"="darkgrey"), labels=c("TRUE"="VDJ region", "FALSE"="non-VDJ region"))+
  ylab("WGCW clusters (prop.)")+
  #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
#  ylim(c(0,1))+
  xlab("SV class")+
  theme_light()+
  ggtitle("Observed")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("allSV_barplot_AGCT_p1_July2020.pdf", width=7,height=5)


# bytype
gd_celltype <- myfiltered_type %>% 
group_by(svclass2, VDJlocusTF,Cell.type2) %>% 
  summarise(m.raghit = mean(raghit),
            n.raghit = sum(raghit),
            sd.raghit = sd(raghit),
            n = n())
ggplot(gd_celltype, aes(x = svclass2, y=m.raghit, group=VDJlocusTF, fill=VDJlocusTF) )+
  geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
  geom_text(aes(y=0.03, label = n) ,  position=position_dodge(width=1), size=2) + 
  scale_fill_manual(name="",values=c("TRUE"="lightgrey", "FALSE"="darkgrey"), labels=c("TRUE"="VDJ region", "FALSE"="non-VDJ region"))+
  facet_wrap(.~Cell.type2)+
  ylab("WGCW clusters (prop.)")+
  #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
  #  ylim(c(0,1))+
  xlab("SV class")+
  theme_light()+
  ggtitle("Observed")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("celltype_allSV_barplot_AGCT_p1_July2020.pdf", width=7,height=5)


  

#### Proportion of SVs with RSS motif- obs and control
ggplot(subset(gd_all, VDJlocusTF==FALSE), aes(x = svclass2, y=m.raghit, group=control, fill=control) )+
  geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
  geom_text(aes(y=0.01, label = n) ,  position=position_dodge(width=1), size=2) + 
  scale_fill_manual(name="",values=c("control"="lightgrey", "obs"="darkgrey"), labels=c("control"="Control", "obs"="Observed"))+
  ylab("WGCW clusters (prop.)")+
  #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
  #ylim(c(0,1))+
  xlab("SV class")+
  theme_light()+
  ggtitle("Non-VDJ")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("control_obs_nonVDJ_barplot_AGCT_p1_July2020.pdf", width=5,height=4)

ggplot(subset(gd_all, VDJlocusTF==TRUE), aes(x = svclass2, y=m.raghit, group=control, fill=control) )+
  geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
  geom_text(aes(y=0.04, label = n) ,  position=position_dodge(width=1), size=2) + 
  scale_fill_manual(name="",values=c("control"="lightgrey", "obs"="darkgrey"), labels=c("control"="Control", "obs"="Observed"))+
  ylab("WGCW clusters (prop.)")+
  #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
  #ylim(c(0,1))+
  xlab("SV class")+
  theme_light()+
  ggtitle("VDJ")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("control_obs_VDJ_barplot_AGCT_p1_July2020.pdf", width=5,height=4)




#### Plotting with RAG results
myRAG = read.table(file="/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/brass_metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/RAGmotifs/results_p10e4_obs_control_fimo_RAGheptamer_fullmotif_20200710_brassfiltered_nodup_analyzed_info_July2020.txt", header=T, stringsAsFactors = F, sep="\t")


myATGC = myfiltered_type[,c("ID","raghit")]
colnames(myATGC) = c("ID","AGCT")
myboth = merge(myATGC, myRAG, by="ID")  # 979 sites
sum(myboth$AGCT==TRUE & myboth$raghit==TRUE) ## none that have both- awesome
gdboth <- myboth %>% 
  group_by(Cell.type2, VDJlocusTF) %>% 
  summarise(m.raghit = mean(raghit),
            m.AGCT = mean(AGCT),
            n = n(),
            n.raghit = sum(raghit),
            n.AGCT = sum(AGCT))
gdboth = data.frame(gdboth)
gdboth$Cell.type2 = c("Memory B","Memory B","Naive B","Naive B","Memory T","Memory T","Naive T","Naive T")
gdboth$Cell.type2 = factor(gdboth$Cell.type2, levels=c("Naive B","Memory B","Naive T","Memory T"))
mytab = melt(subset(gdboth, VDJlocusTF==TRUE)[,c(1,3,4)])
mytab$variable = factor(mytab$variable, levels=c("m.AGCT","m.raghit"))

# prop of Ig/TCR SVs accounted for by RAG OR AID
( sum(subset(gdboth, VDJlocusTF==TRUE)$n.raghit) + sum(subset(gdboth, VDJlocusTF==TRUE)$n.AGCT) ) / sum(subset(gdboth, VDJlocusTF==TRUE)$n)
# 98.2%

ggplot(mytab, aes(x = Cell.type2, y=value, fill=variable) )+
  geom_bar(stat = "identity", col="black", width=0.6, size=0.2)+
  scale_fill_manual(name="", values=c("m.raghit"="#F0BADB", "m.AGCT"="#DF65B0"), labels=c("m.AGCT"="CSR","m.raghit"="RAG"))+
  scale_y_continuous(expand=c(0,0), limits=c(0,1))+
  ylab("RAG/CSR motif (prop.)")+
  #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
  xlab("")+
  theme_bw()+
  #geom_hline(yintercept = 1, lty=2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text=element_text(size=8),
        legend.key.size = unit(0.35, 'cm'),
        legend.position="top")
ggsave("RAG_AGCT_VDJ_barplot_July2020.pdf", width=2,height=2.7)
