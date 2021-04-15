# Aug 2020

# running sigfit on 7 select pcawg cancer types, mm plus all lymph and hsc

#export LD_LIBRARY_PATH=/lustre/scratch116/casm/cgp/users/hm8/software/libv8-3.14/usr/lib:$LD_LIBRARY_PATH
module load R/3.6.1

while read line; do
    chrom=$(echo $line | cut -d " " -f1)
    bsub -J sigfit -n1 -R "select[mem>32000]" -R "rusage[mem=32000]" -M 32000 -R "span[hosts=1]" -q normal -e bsub.error.%J -o bsub.output.%J Rscript sigfit_memoryB_1Mb_Feb2021_min10percent_IgDenovo.R $chrom
done < /lustre/scratch116/casm/cgp/users/hm8/reference_files/GRCh37d5/hg19.chrom_sizes.txt

bsub -J sigfit -n1 -R "select[mem>32000]" -R "rusage[mem=32000]" -M 32000 -R "span[hosts=1]" -q normal -e bsub.error.%J -o bsub.output.%J Rscript signature_prob_per_trinuc_attribution_per1MB_memoryB.R

#bsub -J mutbins -q normal -R 'select[mem>=64000] rusage[mem=64000]' -M64000 -e bsub.%J.err -o bsub.%J.out Rscript genomic_distribution_10Kb_memoryB_sigfit_1Mb_min10percent.R
#bsub -J mutbins -q normal -R 'select[mem>=32000] rusage[mem=32000]' -M32000 -e bsub.%J.err -o bsub.%J.out Rscript genomic_distribution_10Kb_annotation.R
