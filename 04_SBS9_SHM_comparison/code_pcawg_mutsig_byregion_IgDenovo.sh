############################################################################################
## File: sigfit_pcawg_1Mb_Sep2020_min10percent_IgDenovo_clean.R
## Project: lymphocyte_somatic_mutation
## Description: Wrapper script for per-Mb and per-bp attribution (malignancy)
##
## Date: April 2021
## Author: Heather Machado
############################################################################################

# Sept 2020
# calling signatures in the pcawg data by 1Mb bins
# call both S9 and cAID (using sigfit)

# 1) combine all mutations from a given cancer type into one table. Write table.
# 2) for each cancer type, split mutations into 1Mb bins. Save R object.
# 3) for each bin, calculate mutational matrix. Save table per cancer type with info table.
# 4) run sigfit per cancer type using tables from (3)
# 5) Run per-nucleotide attribution using sigfit results (from 4) and mutations per bin (from 2)


# 1) combine all mutations from a given cancer type into one table. Write table.
# 2) for each cancer type, split mutations into 1Mb bins. Save R object.
# 3) for each bin, calculate mutational matrix. Save table per cancer type with info table.
module load R/3.6.1

## done in folder ../pcawg
#bsub -J mutbins -q normal -R 'select[mem>=20000] rusage[mem=20000]' -M20000 -e bsub.%J.err -o bsub.%J.out Rscript make_pertype_mut_table_bins_downsample_mm.R


# 4) run sigfit per cancer type using tables from (3)
#export LD_LIBRARY_PATH=/lustre/scratch116/casm/cgp/users/hm8/software/libv8-3.14/usr/lib:$LD_LIBRARY_PATH
module load R/3.6.1

while read group; do
    mkdir -p $group

    while read line; do
        chrom=$(echo $line | cut -d " " -f1)
        bsub -J sigfit -n1 -R "select[mem>16000]" -R "rusage[mem=16000]" -M 16000 -R "span[hosts=1]" -q normal -e bsub.error.%J -o bsub.output.%J Rscript sigfit_pcawg_1Mb_Sep2020_min10percent_IgDenovo.R $chrom $group
    done < ../data/hg19.chrom_sizes.txt
done < ../data/my_groups_pcawg_mm.txt

while read group; do
    bsub -J mutbins -q normal -R 'select[mem>=32000] rusage[mem=32000]' -M32000 -e bsub.%J.err -o bsub.%J.out Rscript signature_prob_per_trinuc_attribution_per1MB_pcawg_IgDenovo.R $group
done < ../data/my_groups_pcawg_mm.txt
