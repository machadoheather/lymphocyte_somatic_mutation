# lymphocyte_somatic_mutation
 
 ## General info
 The page provides code for analyses conducted in Machado et al. 2021.
 "Genome-wide mutagenesis during immunological diversification in normal lymphocytes"
 
 Author: Heather E. Machado \
 Contact: heather.machado@sanger.ac.uk
 
## Table of contents
* [Variant filtering](01_variant_filtering)
* [Mutation burden](02_mutation_burden_analysis)
* [Mutational signature extraction](03_mutational_signature_analysis)
* [Mutational signature plotting and statistical tests](03_mutational_signature_analysis)
* [SBS9 and SHM mutations](04_SBS9_SHM_comparison)
* [GAM regression of genomic features and mutations](05_SBS9_genomic_feature_GAM)
* [SV analyses](06_SV_analyses)

 ## [Variant filtering](01_variant_filtering)
 * [Filtering SNVs](01_variant_filtering/snv_results)
 * [Filtering indels](01_variant_filtering/indel_results)
 
 ## [Mutation burden](02_mutation_burden_analysis)
 * [Depth-based burden correction](02_mutation_burden_analysis/analysis_HSC_mutburden_correction.html)
 * [Combine datasets](02_mutation_burden_analysis/combining_donors.R)
 * [SNV burden mixed model](02_mutation_burden_analysis/analyses_SNVmixedmodel_March2021.html)
 * [Indel burden mixed model](02_mutation_burden_analysis/analyses_indelmixedmodel_March2021.html)

## [Mutational signature extraction](03_mutational_signature_analysis)
* [hdp signature extraction](03_mutational_signature_analysis/run_hdp_lustre_pcawg_lymph_hsc_clean.R)
* [Sigprofiler signature extraction](03_mutational_signature_analysis/code_sigprofiler_pcwag_mm_lymph_hsc_clean.sh)
* [Sigprofiler signature analysis](03_mutational_signature_analysis/sigprofiler_match_signatures_cosmic.html)
* [SBSblood signature extraction using sigfit](03_mutational_signature_analysis/sigfit_hsc_sig1_Aug2020_clean.R)
* [Comparison of SBSblood](03_mutational_signature_analysis/comparesigs_BloodSig_analysis_Aug2020_clean.html)
* [Analyze hdp and sigprofiler signature extraction](03_mutational_signature_analysis/analyses_hdp_pcawg_mm_AX001_KX001_KX002_KX003_TX001_TX002_CB001_Aug2020_clean.html)
* [sigfit per-genome signature attribution: all signatures](03_mutational_signature_analysis/sigfit_union_hdp_sigprofiler_cosmic3_chain1_Aug2020_clean.R)
* [sigfit per-genome signature attribution: signatures min10%](03_mutational_signature_analysis/sigfit_union_hdp_sigprofiler_min10percent_chain1_Aug2020_clean.R)

## [Mutational signature plotting and statistical tests](03_mutational_signature_analysis)
* [Signature proportion per genome plotting and statistical tests](03_mutational_signature_analysis/analyses_pcawg_mm_AX001_KX001_KX002_KX003_TX001_TX002_CB001_sigfit_union_hdp_sigprofiler_min10percent_clean.html)
* [Correlations with SHM](03_mutational_signature_analysis/analyses_immuno_pcawg_mm_AX001_KX001_KX002_KX003_TX001_TX002_CB001_sigfit_union_hdp_sigprofiler_min10percent_clean.html)

## [SBS9 and SHM mutations](04_SBS9_SHM_comparison)
* [hdp denovo extraction of the SHM signature](04_SBS9_SHM_comparison/mutsig_byregion_hdp_denovo_Oct2020_clean.R)
* [Wrapper script for per-Mb and per-bp attribution (normal memory B cells)](04_SBS9_SHM_comparison/code_memoryB_100kb_Sep2020_1Mb_min10percent_IgDenovo.sh)
* [sigfit per-Mb signature attribution (normal memory B cells)](04_SBS9_SHM_comparison/sigfit_memoryB_1Mb_Feb2021_min10percent_IgDenovo_clean.R)
* [Per-bp signature attribution (normal memory B cells)](04_SBS9_SHM_comparison/signature_prob_per_trinuc_attribution_per1MB_memoryB.R)
* [Downsample datasets and create 1Mb mutation matrices (malignancy)](04_SBS9_SHM_comparison/make_pertype_mut_table_bins_downsample_mm.R)
* [Wrapper script for per-Mb and per-bp attribution (malignancy)](04_SBS9_SHM_comparison/code_pcawg_mutsig_byregion_IgDenovo.sh)
* [sigfit per-Mb signature attribution (malignancy)](04_SBS9_SHM_comparison/sigfit_pcawg_1Mb_Sep2020_min10percent_IgDenovo_clean.R)
* [Per-bp signature attribution (malignancy)](04_SBS9_SHM_comparison/signature_prob_per_trinuc_attribution_per1MB_pcawg_IgDenovo_clean.R)
* [Enrichment per gene](04_SBS9_SHM_comparison/analyze_pergene_sigfit_pcawg_1Mb_min10percent_IgDenovo_ttest_clean.html)

## [GAM regression of genomic features and mutations](05_SBS9_genomic_feature_GAM)
* [Extract genomic feature value per 10Kb](05_SBS9_genomic_feature_GAM/scripts/get_hg19_properties_HEM.R)
* [Helper script to run/combine genomic feature extraction and attribution](05_SBS9_genomic_feature_GAM/code_regression_analysis_10KB.sh)
* [Full model and individual GAM regressions, plotting](05_SBS9_genomic_feature_GAM/S9_majority_gam_reptimingGm_Jan2021_clean.R)

## [SV analyses](06_SV_analyses)
* [SV burden per cell subset](06_SV_analyses/analyze_brass_metaAnalysis_AX001_KX001_KX002_KX003_TX001_TX002_CB001_clean.html)
* [RAG motif analysis](06_SV_analyses/analyze_RAGmotif_heptamer_July2020_clean.R)
* [RAG motif analysis, signal decay with distance from bp](06_SV_analyses/analyze_RAGmotif_heptamer_July2020_decay_clean.R)
* [CSR motif analysis and RAG/CSR plotting](06_SV_analyses/analyze_AGCTmotif_July2020_clean.R)
* [CSR motif analysis (malignancy)](06_SV_analyses/analyze_pcawgALL_AGCTmotif_July2020_clean.R)


## Notes
* Data used in analyses but not supplied here (e.g. file paths external to the github repo) are available on request.
