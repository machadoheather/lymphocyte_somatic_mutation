## July 2020\

# run sigprofiler

# first make vcfs for each sample

### First make mutational matrix
R prep_mutmat.R

## set working directory
dir=/lustre/scratch116/casm/cgp/users/hm8/lymphocyteWGS/metaAnalysis_pcawg_mm_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg_Aug2020

## run sigprofiler
bsub -o log/bsub.%J.out -e log/bsub.%J.err -q basement -n16 -R 'select[mem>=32000] rusage[mem=32000]' -M32000 -J sigprof20 "/software/singularity-v3.5.1/bin/singularity exec \
/software/CASM/singularity/sigprofiler/sigprofiler-docker_0.17-GRCh37.sif \
docker-sigprofiler \
-i $dir/sigprofiler/mutcounts_matrix_pcawg_mm_lymph_hsc_sigprofilerOrder.txt \
-it table \
-n 1000 \
-o $dir/sigprofiler/results20 \
-r GRCh37 \
-ms 3 \
-ts 14 \
-m SBS96"
