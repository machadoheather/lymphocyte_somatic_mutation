############################################################################################
## File: analyze_RAGmotif_heptamer_July2020_clean.R
## Project: lymphocyte_somatic_mutation
## Description: RAG motif analysis, signal decay with distance from bp
##
## Date: April 2021
## Author: Heather Machado
############################################################################################


### Identifying RAG mediated deletions
library(tidyverse)
library(stringr)
library(magrittr)
library(rtracklayer)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(dplyr)
library(cowplot)
library(reshape2)



mydistances = c(0,5,10,15,20,25,50,75,100,125,150,175,200,300,400,500,600,700,800,900,1000,5000,10000)
for (d in 1:length(mydistances)){
  ## Reading in filtered data
  brassfiltered = read.csv(file="../data/brass_ascat_curated.csv", stringsAsFactors = F, header=T)  # 1160   91
  myfiltered = subset(brassfiltered, Correct != FALSE & !(is.na(start1)) )
  myfiltered$ID = paste(myfiltered$chr1, myfiltered$start1, myfiltered$chr2, myfiltered$start2, myfiltered$svclass, sep="_")
  myfiltered$strand1 = "+"
  myfiltered$strand2 = "+"

  ## Fetching the flanking sequences
  mygranges1 = makeGRangesFromDataFrame(myfiltered, seqnames.field="chr1", start.field="start1", end.field="end1", strand.field="strand1", keep.extra.columns=TRUE)
  mygranges2 = makeGRangesFromDataFrame(myfiltered, seqnames.field="chr2", start.field="start2", end.field="end2", strand.field="strand2", keep.extra.columns=TRUE)

  focaldist = mydistances[d]

  # _c. Get breakpoint flanking sequence ------------------------------------
  # __i. Get 50bp offsets ---------------------------------------------------
  # Make Grange around each breakpoint with flank 50bp
  brkpt.1.ext.gr <- promoters( shift(mygranges1, -focaldist), upstream = 50, downstream = 0)
  brkpt.2.ext.gr <- promoters(shift(mygranges2, focaldist), upstream = 0, downstream = 50)
  brkpt.1.int.gr <- promoters( shift(mygranges1, focaldist), upstream = 0, downstream = 50)
  brkpt.2.int.gr <- promoters(shift(mygranges2, -focaldist), upstream = 50, downstream = 0)

  # # __ii. Get Sequence ------------------------------------------------------
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

  write.table(myfiltered, file=paste("motifnames_brassfiltered_all_", focaldist,"bp_July2019.txt", sep=""), quote=FALSE, col.names = T, row.names = F, sep="\t")

  ## Write out fasta for MEME of FIMO
  # only write unique lines
  meme.fasta <- c(rbind(c(names(brkpt.1.ext.gr.seq), names(brkpt.2.ext.gr.seq),names(brkpt.1.int.gr.seq), names(brkpt.2.int.gr.seq) ), c(unlist(brkpt.1.ext.gr.seq), unlist(brkpt.2.ext.gr.seq), unlist(brkpt.1.int.gr.seq), unlist(brkpt.2.int.gr.seq) )))
  writeLines(meme.fasta, paste("MEME_brassfiltered_all_", focaldist, "bp_July2019.fasta", sep="")  )
}



###########################################################################################################################
############################# Run FIMO to detect RSS motifs in the sequences from the fasta's just written
## FIMO options (website)
# fimo --oc . --verbosity 1 --thresh 1.0E-4 RAG-motif_hepnon_combined.meme.txt 20190712_RAGfullmotif_MEME_filtered_updown_split.fasta
###########################################################################################################################




########## Read in data
mydistances = c(0,5,10,15,20,25,50,75,100,125,150,175,200,300,400,500,600,700,800,900,1000,5000,10000)
mean.vdj = mean.vdj.int = mean.vdj.ext = mean.nonvdj = mean.nonvdj.int = mean.nonvdj.ext = mean.total = mean.total.int = mean.total.ext = n.vdj = n.nonvdj = n.total = vector()

for (d in 1:length(mydistances)){
  
  focaldist = mydistances[d]
  
  myfiltered = read.table(file=paste("motifnames_brassfiltered_all_", focaldist,"bp_July2019.txt", sep=""), header = T, stringsAsFactors = F, sep="\t")
  fimoALL = read.table(paste("decay_",focaldist,"bp/fimo.tsv", sep=""), header=TRUE, stringsAsFactors = F, sep="\t")
  pvalue = 1e-4
  fimo = fimoALL[fimoALL$p.value<pvalue,]
  fimo_names = unlist(lapply(fimo$sequence_name, FUN=function(X) paste(">", X, sep="")))
  
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
  
  
  ### For the MISSED VDJ, are there motifs? Use only T cells (b cells have CSR)
  # Read in m3 table and use to include only colonies that made it into the main analysis
  m3 = read.table("../../colonyinfo_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg.txt", sep="\t", header=T) # 734
  myfiltered = subset(myfiltered, colony %in% subset(m3, Cell.type2 %in% c("HSC","B Memory","T Memory","B Naive","T Naive"))$colony)   #  1118
  myfiltered_nodup = myfiltered[!(duplicated(myfiltered$ID)),]
  myfiltered = myfiltered_nodup
  
  # Summary figures
  
  gd_byvdj <- myfiltered %>% 
    group_by(VDJlocusTF) %>% 
    summarise(m.raghit = mean(raghit),
              m.raghit_int = mean(raghit_int),
              m.raghit_ext = mean(raghit_ext),
              n.raghit = sum(raghit),
              n = n())

  gd_all <- myfiltered %>% 
    summarise(m.raghit = mean(raghit),
              n.raghit = sum(raghit),
              m.raghit_int = mean(raghit_int),
              m.raghit_ext = mean(raghit_ext),
              n = n())
  
  mean.vdj[d] = gd_byvdj$m.raghit[gd_byvdj$VDJlocusTF==TRUE]
  mean.nonvdj[d] = gd_byvdj$m.raghit[gd_byvdj$VDJlocusTF==FALSE]
  mean.total[d] = gd_all$m.raghit
  mean.vdj.int[d] = gd_byvdj$m.raghit_int[gd_byvdj$VDJlocusTF==TRUE]
  mean.vdj.ext[d] = gd_byvdj$m.raghit_ext[gd_byvdj$VDJlocusTF==TRUE]
  mean.nonvdj.int[d] = gd_byvdj$m.raghit_int[gd_byvdj$VDJlocusTF==FALSE]
  mean.nonvdj.ext[d] = gd_byvdj$m.raghit_ext[gd_byvdj$VDJlocusTF==FALSE]
  mean.total.int[d] = gd_all$m.raghit_int
  mean.total.ext[d] = gd_all$m.raghit_ext
    
  n.vdj[d] = gd_byvdj$n[gd_byvdj$VDJlocusTF==TRUE]
  n.nonvdj[d] = gd_byvdj$n[gd_byvdj$VDJlocusTF==FALSE]
  n.total[d] = gd_all$n
}
decay_summary = data.frame(dist=mydistances, mean.vdj, mean.vdj.int, mean.vdj.ext, mean.nonvdj, mean.nonvdj.int, mean.nonvdj.ext,  mean.total,  mean.total.int, mean.total.ext, n.vdj, n.nonvdj, n.total)

decay_summary_sym = data.frame(dist=c(mydistances, -mydistances), mean.vdj=c(mean.vdj.int, mean.vdj.ext), mean.nonvdj=c(mean.nonvdj.int, mean.nonvdj.ext),  mean.total= c(mean.total.int, mean.total.ext))
decay_summary_sym2 = melt(decay_summary_sym[,c(1,3,4)], id.vars="dist")


#### include genomic background rate
mybackground = read.table("../data/results_RAGmotif_genomiccontrol.txt", header=T, stringsAsFactors = F)
median(mybackground$m.raghit_int) ## ext is equivalent: median(mybackground$m.raghit_ext)
# 0.05103042

### starting 25bp off (midpoint)
decay_summary_sym_25bp = data.frame(dist=c(mydistances+25, -mydistances-25), mean.vdj=c(mean.vdj.int, mean.vdj.ext), mean.nonvdj=c(mean.nonvdj.int, mean.nonvdj.ext),  mean.total= c(mean.total.int, mean.total.ext))

ggplot(subset(decay_summary_sym_25bp, dist>0), aes(dist, mean.nonvdj))+
  geom_vline(xintercept = 0, lty=1 ,col="grey")+
  geom_hline(yintercept = median(mybackground$m.raghit_int), lty=2 ,col="black")+
  geom_line(color="#005a32")+
  geom_point(color="#005a32", pch=16, size=0.8)+
  geom_line(data=subset(decay_summary_sym_25bp, dist<0), aes(dist, mean.nonvdj), color="#005a32")+
  geom_point(data=subset(decay_summary_sym_25bp, dist<0), aes(dist, mean.nonvdj), color="#005a32", pch=16, size=0.8)+
  scale_x_continuous(limits = c(-325,325))+
  scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.2)) )+
  theme_light()+
  xlab("Distance from breakpoint (center of 50bp bin)")+
  ylab("RAG motif (prop.)")+
  ggtitle("non-Ig/TCR SVs")+
  theme(#strip.text = element_text(hjust = 0, size=10, face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 10))
ggsave("RAG_decay_nonVDJ_25bp_offset.pdf", width=4,height=1.5)

ggplot(subset(decay_summary_sym_25bp, dist>0), aes(dist, mean.vdj))+
  geom_vline(xintercept = 0, lty=1 ,col="grey")+
  geom_hline(yintercept = median(mybackground$m.raghit_int), lty=2 ,col="black")+
  geom_line(color="#df65b0")+  
  geom_point(color="#df65b0", pch=16, size=0.8)+
  geom_line(data=subset(decay_summary_sym_25bp, dist<0), aes(dist, mean.vdj), color="#df65b0")+
  geom_point(data=subset(decay_summary_sym_25bp, dist<0), aes(dist, mean.vdj), color="#df65b0", pch=16, size=0.8)+
  scale_x_continuous(limits = c(-325,325))+
  scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.2)) )+
  theme_light()+
  xlab("Distance from breakpoint (center of 50bp bin)")+
  ggtitle("Ig/TCR SVs")+
  ylab("RAG motif (prop.)")+
  theme(#strip.text = element_text(hjust = 0, size=10, face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 10) )
ggsave("RAG_decay_VDJ_25bp_offset.pdf",  width=4,height=1.5)


save(decay_summary_sym_25bp, decay_summary_sym, file="RAG_decay_summary_sym.Rdata")
