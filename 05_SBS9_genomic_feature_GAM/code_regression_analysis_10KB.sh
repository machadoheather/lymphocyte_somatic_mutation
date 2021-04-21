#!/usr/bin/env bash

############################################################################################
## File: code_regression_analysis_10KB.sh
## Project: lymphocyte_somatic_mutation
## Description: Helper script to run/combine genomic feature extraction and attribution
##
## Date: April 2021
## Author: Heather Machado
############################################################################################

module load R/3.6.1

## Real fetching of genomic attributes

# 1a) Get S9 profile (done locally)
#   hdp analysis

# 1) Get attributes per window (general to all chromosomes)
#   script: trinuc_sigprob_attribution_by_window.r
#   input: attribution file: "/lustre/scratch116/casm/cgp/users/hm8/lymphocyteWGS/mutsig_byregion/mutsig_byregion_hdp_denovo/attribution_perbp_hdp_denovo_Memory_B_1Mb.txt"
#           window size: 10000
#           signature: S9
#   output: per_window_attributes.Rdata
signature=S9
signature=SIg
signature=X1
bsub -J trinuc_win -q normal -R 'select[mem>=8000] rusage[mem=8000]' -M8000 -e bsub.%J.err -o bsub.%J.out Rscript scripts/trinuc_sigprob_attribution_by_window.R /lustre/scratch116/casm/cgp/users/hm8/lymphocyteWGS/mutsig_byregion/mutsig_byregion_hdp_denovo/attribution_perbp_hdp_denovo_Memory_B_1Mb_colnames.txt 10000 $signature


# 2) Get per window signature expectation (run for each chromosome and each signature)
#   script: scripts/signature_expectation_per_window.R
#   input:  per_window_attributes.Rdata (int_files/S9attributed_window1000bp_chr$chrom.Rdata)
#           signature: S9
#   output: per_window_attributes_S9exp.txt (int_files/windows1kb_S9attribute_S9exp_chr$chrom.txt)
#signature=S9
#signature=S1
#while read chrom; do
#    bsub -J S9exp_win -n1 -R "select[mem>16000]" -R "rusage[mem=16000]" -M 16000 -R "span[hosts=1]" -q long -e error.%J -o output.%J Rscript scripts/signature_expectation_per_window.R int_files/$signature\attributed_window10000bp_chr$chrom.Rdata int_files/$signature\exp_$signature\attributed_window10000bp_chr$chrom.txt $signature
#done < ../regression_analysis/chrom_list.txt


3) Get per window genomic features (run for each chromosome, maybe for each signature)
#   script: scripts/get_hg19_properties_HEM.R
#   input: int_files/windows1kb_S9attribute_S9exp_chr$chrom.txt
#           location: "farm"
#   output: results/windows1kb_S9attribute_S9exp_chr$chrom_hg19properties.txt
while read chrom; do
    bsub -J attr_win -n1 -R "select[mem>4000]" -R "rusage[mem=4000]" -M 4000 -R "span[hosts=1]" -q normal -e error.%J -o output.%J "/software/R-3.6.1/bin/Rscript scripts/get_hg19_properties_HEM.R int_files/S9exp_S9attributed_window10000bp_chr$chrom.txt  results/hg19properties_S9exp_S9attributed_window10000bp_chr$chrom.txt"
done < ../regression_analysis/chrom_list.txt

# without the expected
signature=S9
signature=SIg
signature=X1
while read chrom; do
    bsub -J attr_win -n1 -R "select[mem>4000]" -R "rusage[mem=4000]" -M 4000 -R "span[hosts=1]" -q normal -e bsub.%J.err -o bsub.%J.out Rscript scripts/get_hg19_properties_HEM.R int_files/$signature\attributed_window10000bp_chr$chrom.Rdata results/hg19properties_$signature\attributed_window10000bp_chr$chrom.txt
done < ../regression_analysis/chrom_list.txt


## combine chromosomes and save as R object
signature="S9"
signature="SIg"
signature="X1"
library(BSgenome.Hsapiens.UCSC.hg19)
library(selectiveInference)
allfiles = list.files("results")
S9files = allfiles[grep(allfiles, pattern="S9")]
s9_results_list = lapply(S9files, FUN=function(X) read.table(file=paste("results/", X, sep=""), header=TRUE, stringsAsFactors = F, sep="\t") )
s9_results = do.call(rbind, s9_results_list)
s9_results$start = s9_results$start-5000
s9_results$end = s9_results$end+4999

sigfiles = allfiles[grep(allfiles, pattern="SIg")]
sig_results_list = lapply(sigfiles, FUN=function(X) read.table(file=paste("results/", X, sep=""), header=TRUE, stringsAsFactors = F, sep="\t") )
sig_results = do.call(rbind, sig_results_list)

X1files = allfiles[grep(allfiles, pattern="X1")]
X1_results_list = lapply(X1files, FUN=function(X) read.table(file=paste("results/", X, sep=""), header=TRUE, stringsAsFactors = F, sep="\t") )
X1_results = do.call(rbind, X1_results_list)

s9_results$meanX1 = X1_results$meanX1
s9_results$sumX1 = X1_results$sumX1
s9_results$meanSIg = sig_results$meanSIg
s9_results$sumSIg = sig_results$sumSIg
s9_resultsGR = makeGRangesFromDataFrame(s9_results, keep.extra.columns=TRUE)

save(s9_resultsGR, file=paste("results/resultsGR10Kb_hdp_denovo_1Mb.Rdata", sep="") )
#load(file="results/s9_resultsGR10Kb_hdp_denovo_1Mb.Rdata")


### annotation the results
 bsub -J S9exp_win -n8 -R "select[mem>32000]" -R "rusage[mem=32000]" -M 32000 -R "span[hosts=1]" -q normal -e error.%J -o output.%J Rscript scripts/annotate_10Kb_bins.R
