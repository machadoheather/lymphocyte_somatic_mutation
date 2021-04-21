############################################################################################
## File: code_sigprofiler_pcwag_mm_lymph_hsc_clean.sh
## Project: lymphocyte_somatic_mutation
## Description: Sigprofiler signature extraction
##
## Date: April 2019
## Author: Heather Machado
############################################################################################

# make vcfs for each sample
# run sigprofiler

### First make mutational matrix
R prep_mutmat.R


## run sigprofiler
bsub -o log/bsub.%J.out -e log/bsub.%J.err -q basement -n16 -R 'select[mem>=32000] rusage[mem=32000]' -M32000 -J sigprof20 "/software/singularity-v3.5.1/bin/singularity exec \
/software/CASM/singularity/sigprofiler/sigprofiler-docker_0.17-GRCh37.sif \
docker-sigprofiler \
-i sigprofiler/mutcounts_matrix_pcawg_mm_lymph_hsc_sigprofilerOrder.txt \
-it table \
-n 1000 \
-o sigprofiler/results20 \
-r GRCh37 \
-ms 3 \
-ts 14 \
-m SBS96"
