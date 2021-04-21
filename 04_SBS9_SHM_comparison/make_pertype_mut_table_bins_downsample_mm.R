############################################################################################
## File: make_pertype_mut_table_bins_downsample_mm.R
## Project: lymphocyte_somatic_mutation
## Description: Downsample datasets and create 1Mb mutation matrices (malignancy)
##
## Date: April 2021
## Author: Heather Machado
############################################################################################

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

library("GenomicRanges")
library("Rsamtools")
library("MASS")
library(stringr)




# ######## per MB
# ## Calculate mutmat per 1MB window, 1 chrom at a time (ALL samples combined- almost zero data per sample)
# # 1 row per window
# # columns: chrom, st, end, Nmut, Nsamples, meanS9, sumS9, 1se=sd(S9)/sqrt(Nmut)


#########
# 2
#########
genomeFile = paste(root, "/reference_files/GRCh37d5/genome.fa", sep="")
chromsize = read.table(paste(root, "/reference_files/GRCh37d5/hg19.chrom_sizes.txt", sep=""), stringsAsFactors = FALSE)
mychroms = chromsize[,1] #c(1:22, "X","Y")
my_groups = "Multiple_myeloma"

samplesize=45844
  focalsample = read.table(file=paste(root, "/pcawg/", "chr_pos_ref_alt.",my_groups, ".txt", sep=""), sep="\t", stringsAsFactors = F, header=T)
    if (nrow(focalsample > samplesize)){
        focalsample = focalsample[sample(1:nrow(focalsample), size=samplesize),]
    }
  Nregions = sum(chromsize[,3]/1000000) + length(mychroms)
  mydf = matrix(nrow=96, ncol=Nregions)
  regionstart = regionend = regionchrom = vector()
  count=0
  for (chrom in 1:length(mychroms)){
    focalchrom = focalsample[focalsample$chr==mychroms[chrom],]
    bins = seq(from=1, to=chromsize[chrom,3], by=1000000) # add last bin
    bins[length(bins)+1] = bins[length(bins)] + 1000000 # add last bin
    for (j in 1:(length(bins)-1)){
      count = count+1
      regionstart[count] = bins[j]
      regionend[count] = bins[j+1]
      regionchrom[count] = mychroms[chrom]
      focalbin = focalchrom[focalchrom$pos >= bins[j] & focalchrom$pos < bins[j+1], ]
      if (nrow(focalbin)==0){
        mydf[ ,count] = rep(0,times=96) # skip
      } else {
        # calc mut mat
        subs_only = focalbin
        colnames(subs_only) = c("chr","pos","ref","mut")
        subs_only = subs_only[(subs_only$ref %in% c("A","C","G","T")) & (subs_only$mut %in% c("A","C","G","T")) & subs_only$chr %in% c(1:22,"X","Y"),]
        subs_only$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(subs_only$chr, IRanges(subs_only$pos-1, subs_only$pos+1))))
        ntcomp = c(T="A",G="C",C="G",A="T")
        subs_only$sub = paste(subs_only$ref,subs_only$mut,sep=">")
        subs_only$trinuc_ref_py = subs_only$trinuc_ref
        for (k in 1:nrow(subs_only)) {
          if (subs_only$ref[k] %in% c("A","G")) { # Purine base
            subs_only$sub[k] = paste(ntcomp[subs_only$ref[k]],ntcomp[subs_only$mut[k]],sep=">")
            subs_only$trinuc_ref_py[k] = paste(ntcomp[rev(strsplit(subs_only$trinuc_ref[k],split="")[[1]])],collapse="")
          }
        }
        # Counting subs
        freqs = table(paste(subs_only$sub,paste(substr(subs_only$trinuc_ref_py,1,1),substr(subs_only$trinuc_ref_py,3,3),sep="-"),sep=","))
        sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
        ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
        full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
        freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
        mydf[ ,count] = freqs_full
        }
    }
  }
  
  mutmax = mydf[,1:count]  # 3113
  mutsum = apply(mutmax, MARGIN=2, FUN=sum)
  mutmaxinfo = data.frame(chr=regionchrom, start=regionstart, end=regionend, Nmut=mutsum)
  
  save(mutmax, mutmaxinfo, file=paste(root, "/lymphocyteWGS/mutsig_byregion/pcawg/mutmats_byregion_", my_groups, "_1Mb_N45844_Sep2020.Rdata", sep="") )
  



