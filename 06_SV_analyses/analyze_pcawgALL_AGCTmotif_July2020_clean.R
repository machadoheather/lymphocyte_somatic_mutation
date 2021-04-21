############################################################################################
## File: analyze_pcawgALL_AGCTmotif_July2020_clean.R
## Project: lymphocyte_somatic_mutation
## Description: CSR motif analysis- cancer datasets
##
## Date: April 2021
## Author: Heather Machado
############################################################################################

### July, 2019
library(tidyverse)
library(stringr)
library(magrittr)
library(rtracklayer)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(dplyr)
library(cowplot)



## Reading in filtered data: only blood cell relevant types
myfiltered = read.table("../data/pcawg_ALL_SVs_postanalysis_July2020.txt", header=T, stringsAsFactors = F, sep="\t")
# 6519 SVs   2227
myfiltered2 = myfiltered[c(1:2226,2228:nrow(myfiltered)),]
#myfiltered2 = subset(myfiltered, sample=="82b8cda8-fbff-455e-b0db-7ff6528bd6c8" & start1==154770320 & start2==133854722)
myfiltered = myfiltered2[!(duplicated(myfiltered2)),]

## removing the one SV that is too close to the edge to  
# myfiltered[2213,]
# sample                     celltype group chr1    start1
# 2227 82b8cda8-fbff-455e-b0db-7ff6528bd6c8 Diffse large B-cell lymphoma     8    4 154770320
# end1 chr2    start2      end2 svtype VDJlocus VDJlocusTF bp1_ext_uniID bp2_ext_uniID
# 2227 154770321   12 133854722 133854723    ITX    FALSE      FALSE >bp1_ext_2227 >bp2_ext_2227
# bp1_int_uniID bp2_int_uniID
# 2227 >bp1_int_2227 >bp2_int_2227

## Fetching the flanking sequences
mygranges1 = makeGRangesFromDataFrame(myfiltered, seqnames.field="chr1", start.field="start1", end.field="end1", keep.extra.columns=TRUE)
mygranges2 = makeGRangesFromDataFrame(myfiltered, seqnames.field="chr2", start.field="start2", end.field="end2", keep.extra.columns=TRUE)

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
names(brkpt.1.ext.gr.seq) <- paste0(">bp1_ext_", rownames(myfiltered))

brkpt.2.ext.gr.seq <- as.list(getSeq(hs37d5, brkpt.2.ext.gr, as.character = T))
names(brkpt.2.ext.gr.seq) <- paste0(">bp2_ext_", rownames(myfiltered))

brkpt.1.int.gr.seq <- as.list(getSeq(hs37d5, brkpt.1.int.gr, as.character = T))
names(brkpt.1.int.gr.seq) <- paste0(">bp1_int_",rownames(myfiltered))

brkpt.2.int.gr.seq <- as.list(getSeq(hs37d5, brkpt.2.int.gr, as.character = T))
names(brkpt.2.int.gr.seq) <- paste0(">bp2_int_", rownames(myfiltered))

# Add fasta name to original file bp1 & bp2
myfiltered$bp1_ext_uniID <- names(brkpt.1.ext.gr.seq)
myfiltered$bp2_ext_uniID <- names(brkpt.2.ext.gr.seq)
myfiltered$bp1_int_uniID <- names(brkpt.1.int.gr.seq)
myfiltered$bp2_int_uniID <- names(brkpt.2.int.gr.seq)
write.table(myfiltered, file="motifnames_pcawgALL_controls_1000bp_July2020.txt", quote=FALSE, col.names = T, row.names = F, sep="\t")

## Write out fasta for MEME of FIMO
# only write unique lines
meme.fasta <- c(rbind(c(names(brkpt.1.ext.gr.seq), names(brkpt.2.ext.gr.seq),names(brkpt.1.int.gr.seq), names(brkpt.2.int.gr.seq) ), c(unlist(brkpt.1.ext.gr.seq), unlist(brkpt.2.ext.gr.seq), unlist(brkpt.1.int.gr.seq), unlist(brkpt.2.int.gr.seq) )))
writeLines(meme.fasta,  "MEME_pcawgALL_updown_split_1000bp_July2020.fasta")




## Creating a control set of sites
# Add a random offset of between 200-2000bp 
# Going upstream for breakpoint 1 and downstream for breakpoint2 so can still do internal/external RSS motif analysis
## Reading in filtered data: only blood cell relevant types
myfiltered = read.table("../data/pcawg_ALL_SVs_postanalysis_July2020.txt", header=T, stringsAsFactors = F, sep="\t")
# 6519 SVs   2227
myfiltered2 = myfiltered[c(1:2226,2228:nrow(myfiltered)),]
#myfiltered2 = subset(myfiltered, sample=="82b8cda8-fbff-455e-b0db-7ff6528bd6c8" & start1==154770320 & start2==133854722)
myfiltered = myfiltered2[!(duplicated(myfiltered2)),]
mycontrol = myfiltered
offset1 = sample(c(-30000:-20000), size=nrow(mycontrol), replace = T)
offset2 = sample(c(20000:30000), size=nrow(mycontrol), replace = T)
mycontrol$start1 = myfiltered$start1 + offset1
mycontrol$end1 = myfiltered$end1 + offset1
mycontrol$start2 = myfiltered$start2 + offset2
mycontrol$end2 = myfiltered$end2 + offset2
myfiltered = mycontrol

## Fetching the flanking sequences
mygranges1 = makeGRangesFromDataFrame(myfiltered, seqnames.field="chr1", start.field="start1", end.field="end1", keep.extra.columns=TRUE)
mygranges2 = makeGRangesFromDataFrame(myfiltered, seqnames.field="chr2", start.field="start2", end.field="end2", keep.extra.columns=TRUE)

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
names(brkpt.1.ext.gr.seq) <- paste0(">bp1_ext_", rownames(myfiltered))

brkpt.2.ext.gr.seq <- as.list(getSeq(hs37d5, brkpt.2.ext.gr, as.character = T))
names(brkpt.2.ext.gr.seq) <- paste0(">bp2_ext_", rownames(myfiltered))

brkpt.1.int.gr.seq <- as.list(getSeq(hs37d5, brkpt.1.int.gr, as.character = T))
names(brkpt.1.int.gr.seq) <- paste0(">bp1_int_",rownames(myfiltered))

brkpt.2.int.gr.seq <- as.list(getSeq(hs37d5, brkpt.2.int.gr, as.character = T))
names(brkpt.2.int.gr.seq) <- paste0(">bp2_int_", rownames(myfiltered))

# Add fasta name to original file bp1 & bp2
myfiltered$bp1_ext_uniID <- names(brkpt.1.ext.gr.seq)
myfiltered$bp2_ext_uniID <- names(brkpt.2.ext.gr.seq)
myfiltered$bp1_int_uniID <- names(brkpt.1.int.gr.seq)
myfiltered$bp2_int_uniID <- names(brkpt.2.int.gr.seq)
write.table(myfiltered, file="motifnames_pcawgALL_controls_1000bp_July2020.txt", quote=FALSE, col.names = T, row.names = F, sep="\t")

## Write out fasta for MEME of FIMO
# only write unique lines
meme.fasta <- c(rbind(c(names(brkpt.1.ext.gr.seq), names(brkpt.2.ext.gr.seq),names(brkpt.1.int.gr.seq), names(brkpt.2.int.gr.seq) ), c(unlist(brkpt.1.ext.gr.seq), unlist(brkpt.2.ext.gr.seq), unlist(brkpt.1.int.gr.seq), unlist(brkpt.2.int.gr.seq) )))
writeLines(meme.fasta,  "MEME_pcawgALL_updown_split_1000bp_controls_July2020.fasta")



###########################################################################################################################
######################################### RUN MCAST
## MCAST options (website)
# motif p-value threshold 0.01
# output threshold <= 100
# spacing max 50bp
#
# type in 2 motifs:
# AGCT
# TGCA
###########################################################################################################################



#########################################
#####   Check for switch motif using MCAST
myfiltered = read.table(file="../data/motifnames_pcawgALL_1000bp_July2020.txt", header = T, stringsAsFactors = F, sep="\t")
fimoALL = read.table("../data/pcawgALL_mcast.tsv", header=TRUE, stringsAsFactors = F, sep="\t")
#pvalue = 1e-5
pvalue = 1
fimo = fimoALL[fimoALL$p.value<pvalue,]
fimo_names = unlist(lapply(fimo$sequence_name, FUN=function(X) paste(">", X, sep="")))
myfilteredcont = read.table(file="../data/motifnames_pcawgALL_controls_1000bp_July2020.txt", header = T, stringsAsFactors = F, sep="\t")
fimoALLcont = read.table("../data/pcawgALL_controls_mcast.tsv", header=TRUE, stringsAsFactors = F, sep="\t")
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
myfiltered$svtype2 = myfiltered$svtype
myfiltered$svtype2[myfiltered$svtype=="tandem-duplication"] = "tandem-dup."
#myfiltered$non.template.T = myfiltered$non.template != "_" & myfiltered$non.template != "." 

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
#myfiltered = myfiltered[!(duplicated(myfiltered$ID)),]  # remove 16 lines that are duplicated

## Writing results table
#write.table(myfiltered, file="results_mcastWGCW_p1e10000space100_pcawgALL_1000bp_July2020.txt", quote=F, col.names = T, row.names = F, sep="\t")


#########################################
##### RE-START HERE
#########################################
#### RAG for VDJ and non-VDJ, by SV class
myfiltered = read.table(file="../data/results_mcastWGCW_p1e10000space100_pcawgALL_1000bp_July2020.txt", header = T, sep="\t") # 1126

# Summary figures
gd <- myfiltered %>% 
  group_by(celltype, VDJlocusTF) %>% 
  summarise(m.raghit = mean(raghit),
            n.raghit = sum(raghit),
            sd.raghit = sd(raghit),
            n = n())
gd_cont <- myfiltered %>% 
  group_by(celltype, VDJlocusTF) %>% 
  summarise(m.raghit_cont = mean(raghit_cont),
            n.raghit = sum(raghit_cont),
            sd.raghit_cont = sd(raghit_cont),
            n = n())
gd$control = "obs"
gd_cont$control = "control"
colnames(gd_cont) = colnames(gd)
gd_all = rbind(gd, gd_cont)


## Proportion of SVs with CSR motif
#g1=
  ggplot(gd, aes(x = celltype, y=m.raghit, group=VDJlocusTF, fill=VDJlocusTF) )+
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
ggsave("pcawgall_bycelltype_barplot_AGCT_July2020.pdf", width=7,height=5)


# bytype
gd_celltype <- myfiltered %>% 
group_by(svtype2, VDJlocusTF,celltype) %>% 
  summarise(m.raghit = mean(raghit),
            n.raghit = sum(raghit),
            sd.raghit = sd(raghit),
            n = n())
ggplot(gd_celltype, aes(x = svtype2, y=m.raghit, group=VDJlocusTF, fill=VDJlocusTF) )+
  geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
  geom_text(aes(y=0.03, label = n) ,  position=position_dodge(width=1), size=2) + 
  scale_fill_manual(name="",values=c("TRUE"="lightgrey", "FALSE"="darkgrey"), labels=c("TRUE"="VDJ region", "FALSE"="non-VDJ region"))+
  facet_wrap(.~celltype)+
  ylab("WGCW clusters (prop.)")+
  #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
  #  ylim(c(0,1))+
  xlab("SV class")+
  theme_light()+
  ggtitle("Observed")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("pcawgall_celltype_allSV_barplot_AGCT_July2020.pdf", width=7,height=10)


#### Proportion of SVs with RSS motif- obs and control
ggplot(subset(gd_all, VDJlocusTF==FALSE), aes(x = celltype, y=m.raghit, group=control, fill=control) )+
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
ggsave("pcawgall_control_obs_nonVDJ_barplot_AGCT_July2020.pdf", width=5,height=4)



ggplot(subset(gd_all, VDJlocusTF==TRUE), aes(x = celltype, y=m.raghit, group=control, fill=control) )+
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
ggsave("pcawgall_control_obs_VDJ_barplot_AGCT_July2020.pdf", width=5,height=4)




#### Proportion of SVs with RSS motif by cell subset- obs and control
gd_celltype <- myfiltered_type %>% 
  group_by(svtype2, VDJlocusTF,Cell.type2) %>% 
  summarise(m.raghit = mean(raghit),
            n.raghit = sum(raghit),
            sd.raghit = sd(raghit),
            n = n())
gd_celltype_cont <- myfiltered_type %>% 
  group_by(svtype2, VDJlocusTF,Cell.type2) %>% 
  summarise(m.raghit_cont = mean(raghit_cont),
            n.raghit = sum(raghit_cont),
            sd.raghit_cont = sd(raghit_cont),
            n = n())
gd_celltype$control = "obs"
gd_celltype_cont$control = "control"
colnames(gd_celltype_cont) = colnames(gd_celltype)
gd_celltype_all = rbind(gd_celltype, gd_celltype_cont)

ggplot(subset(gd_celltype_all, VDJlocusTF==FALSE), aes(x = svtype2, y=m.raghit, group=control, fill=control) )+
  geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
  facet_wrap(.~Cell.type2)+
  geom_text(aes(y=0.01, label = n) ,  position=position_dodge(width=1), size=2) + 
  scale_fill_manual(name="",values=c("control"="lightgrey", "obs"="darkgrey"), labels=c("control"="Control", "obs"="Observed"))+
  ylab("WGCW clusters  (prop.)")+
  #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
  #ylim(c(0,1))+
  xlab("SV class")+
  theme_light()+
  ggtitle("Non-VDJ")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("celltype_control_obs_nonVDJ_barplot_AGCT_July2020.pdf", width=5,height=4)

ggplot(subset(gd_celltype_all, VDJlocusTF==TRUE), aes(x = svtype2, y=m.raghit, group=control, fill=control) )+
  geom_bar(stat = "identity", position=position_dodge(), col="black", size=0.2)+
  geom_text(aes(y=0.04, label = n) ,  position=position_dodge(width=1), size=2) + 
  scale_fill_manual(name="",values=c("control"="lightgrey", "obs"="darkgrey"), labels=c("control"="Control", "obs"="Observed"))+
  facet_wrap(.~Cell.type2)+
  ylab("WGCW clusters  (prop.)")+
  #geom_errorbar(aes(ymin=m.raghit_int-sd.raghit_int, ymax=m.raghit_int+sd.raghit_int), width=.2, position=position_dodge(1)) #+
  #ylim(c(0,1))+
  xlab("SV class")+
  theme_light()+
  ggtitle("VDJ")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("celltype_control_obs_VDJ_barplot_AGCT_July2020.pdf", width=5,height=4)






