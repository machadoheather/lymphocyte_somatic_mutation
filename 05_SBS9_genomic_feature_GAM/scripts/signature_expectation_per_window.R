## Feb 2020
# test_signature_expectation_per_window.R

library("GenomicRanges")
library("Rsamtools")
library("MASS")

## Designed to be run from one of the two following directories
#setwd('/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/mutsig_byregion/regression_analysis')


## Feb 2020

# ## Calculating per window expected S9 contribution
# # 1) fetch S9 profile
# # 2) fetch window trinucleotide list
# # 3) match window trinucleotide and S9 profile
# # 4) calculate statistic over window

## Example files
# chrom 14 windows
# infile="/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/trinuc_sigprob_attribution.ARG380_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/S9attributed_1KB_test.Rdata"
# outfile="int_files/windows1kb_S9attribute_S9exp_chr14.txt"

# args = c("int_files/S1attributed_window1000bp_chr1.Rdata", "int_files/S1exp_S1attributed_window1000bp_chr1.txt", "S1")
args = commandArgs(trailingOnly = TRUE)
infile=args[1]
outfile=args[2]
signature=args[3]

# reads in the functions: fetch_sig_trinuc_profile and per_window_sig_exp1
source("scripts/calculate_signature_expectation_per_window.R")

# load in persig_tri and mut_count objects (from main analysis)
load(file="int_files/persig_tri.metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg.Rdata")
# signatures to choose: ("X0","S1","Sblood","S9","S8","S7","N1","N2","S17")



# Load in windows info
load(infile) # df1
# start  end Nmut Nsamples meanS9 sumS9 chr
#1     1 1000    0        0      0     0   1
#2  1001 2000    0        0      0     0   1



# Annotating Ig regions
#### annotate ig regions +/- 1kb
ighst = 106304735 # 14
ighend = 107283226 # 14
iglst = 22385390 # 22
iglend = 23263607 # 22
igkst = 89160078 # 2
igkend = 90274237 # 2
# encompassing all class switching genes
cs_start = 106053274
cs_end = 106322322

# call as Ig if any of the window overlaps an Igh Igl or Igk gene, or class switch region
VDJlocus = vector()
for (i in 1:nrow(df1)){
  X = df1[i,]
  if ( (X$chr == 14 & X$start > ighst-1000 & X$end < ighend+1000) | (X$chr == 14 & X$end > ighst-1000 & X$end < ighend+1000 ) ) {VDJlocus[i] = "igh"} else
    if ( (X$chr == 22 & X$start > iglst-1000 & X$start < iglend+1000) |  (X$chr == 22 & X$end > iglst-1000 & X$end < iglend+1000) ) {VDJlocus[i] = "igl"} else
          if  ( (X$chr == 2 & X$start > igkst-1000 & X$start < igkend+1000) | (X$chr == 2 & X$end > igkst-1000 & X$end < igkend+1000) ) {VDJlocus[i] = "igk"} else 
                if ( (X$chr == 14 & X$start > cs_start-1000 & X$start < cs_end+1000) | (X$chr == 14 & X$end > cs_start-1000 & X$end < cs_end+1000) ) {VDJlocus[i] = "classS"} else 
                {VDJlocus[i] = FALSE}
}
df1$Ig = VDJlocus


# subsetting the datasets
df1_nonIg = subset(df1, Ig==FALSE)
nonzero_mut = subset(df1_nonIg, Nmut>0)
zero_mut = subset(df1_nonIg, Nmut==0)


# calculate expected number of S9 mutations per trinucleotide
sigprop_trinuc = fetch_sig_trinuc_profile(signature=signature)


# calculate expected number of S9 mutations per window (compared to 1 mutation in a window of trinucs of equal proportions)
# include windows with no mutations as well: 12x as many as with mutations for a final 10x
windows_S9exp = per_window_sig_exp1(nonzero_mut, sigprop_trinuc)
if ( nrow(zero_mut) <= nrow(nonzero_mut)*14 ){
    windows_S9exp_zero = per_window_sig_exp1(zero_mut, sigprop_trinuc)
} else {
    windows_S9exp_zero = per_window_sig_exp1(zero_mut[sample(1:nrow(zero_mut), size=nrow(nonzero_mut)*14),], sigprop_trinuc)
}


# there are NA's- some windows have N's: drop from further analysis
windows_S9expNA = windows_S9exp[!(is.na(windows_S9exp$sig_exp1)) ,]
windows_S9exp_zeroNA = windows_S9exp_zero[!(is.na(windows_S9exp_zero$sig_exp1)) ,]

# downsample zero mut windows to be 10x that of the windows with mutations
if (nrow(windows_S9exp_zeroNA) > 10*nrow(windows_S9expNA) ) {
  windows_S9exp_zeroNA_10x =  windows_S9exp_zeroNA[sample(1:nrow(windows_S9exp_zeroNA), size=10*nrow(windows_S9expNA)), ]
} else {
  windows_S9exp_zeroNA_10x =  windows_S9exp_zeroNA
}


## write table of both
windows_S9exp_bothNA = rbind(windows_S9expNA, windows_S9exp_zeroNA_10x)
write.table(windows_S9exp_bothNA, file=outfile, col.names = T, quote=F, sep="\t", row.names = F)
