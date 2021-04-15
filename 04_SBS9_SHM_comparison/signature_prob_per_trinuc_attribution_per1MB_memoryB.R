## Aug 2020 


## For sigfit results, for 1Mb bins
# 1) load in signatures
# 2) extraction attribution function
# 3) create mutation list, with each element being the mutations in one 1Mb bin
# 4) apply attribution function per 1Mb bin

## Making an attribution function to assign per-mutation+trinucelotide type signature probabilities
# then merge theses attributions with per-genome variants.

#Probability of a particular mutation type being from a particular signature for genome1, eg. probability of A(A>T)G being from Signature1 for genome1:
  
#  (Proportion of Sig1 in genome1) * (Proportion of A(A>T)G in Sig1) /
#  Sum over all signatures(    (Proportion of SigX in genome1) * (Proportion of A(A>T)G in SigX)  )


## if working on the farm, might need these:
#export LD_LIBRARY_PATH=/lustre/scratch116/casm/cgp/users/hm8/software/libv8-3.14/usr/lib:$LD_LIBRARY_PATH
#module load R/3.6.1

#args = commandArgs(trailingOnly = TRUE)


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


# 1) load in signatures
mysigs = read.table(paste(root, "/lymphocyteWGS/mutsig_byregion/mysigs_min10percent_Sep2020.txt", sep=""), header=T, stringsAsFactors = F)


# 2a) load in sigfit results, combine chromosomes
files <- list.files(paste(root, "/lymphocyteWGS/mutsig_byregion/sigfit/sigfit_memoryB_1Mb_min10percent_IgDenovo/", sep=""), pattern = "exposures_memoryB_1Mb_IgDenovo.sigfit_union_hdp_sigprofiler_cosmic")
inlist = list()
infolist = list()
for (i in 1:length(files)){
  load(paste(root, "/lymphocyteWGS/mutsig_byregion/sigfit/sigfit_memoryB_1Mb_min10percent_IgDenovo/", files[i], sep=""))
  inlist[[i]] = exposures$mean
  infolist[[i]] = samp_type
}
sigbybin = data.frame(do.call(rbind, inlist))
sampbybin = data.frame(do.call(rbind, infolist))


# 2) extraction attribution function
extract_attribution = function(mysigs, sigfit_means, sampleinfo){
  persig_tri = t(mysigs)
  # labels for the trinucs
  # trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
  # group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
  #                               each=16))
  # sub_trinuc_context = paste(group_factor, trinuc_context, sep="_")
  # write.table(sub_trinuc_context, paste(root, "/lymphocyteWGS/mutsig_byregion/sub_trinuc_context.txt", sep=""), col.names = F, row.names = F, quote=F)
  sub_trinuc_context = read.table(paste(root, "/lymphocyteWGS/mutsig_byregion/sub_trinuc_context.txt", sep=""), stringsAsFactors = F)[,1]
  
  ## Calculating per-mutation+trinucelotide type signature probabilities
  pergenome_attribution_prob = list()
  for (i in 1:nrow(sigfit_means)){ ## for each genome
    focal_genome=sigfit_means[i,]
    prob_tri_mat = matrix(nrow=96, ncol=length(focal_genome))
    for (j in 1:96){  ## for each trinucleotide prop of each trinucleotide per sig
      focal_tri = persig_tri[j,]
      prob_tri_mat[j,] = unlist((focal_genome*focal_tri)/sum(focal_genome*focal_tri))
    }
    rownames(prob_tri_mat) = sub_trinuc_context
    colnames(prob_tri_mat) = names(focal_genome)
    pergenome_attribution_prob[[i]] = prob_tri_mat
  }
  names(pergenome_attribution_prob) = sampleinfo$samp_names
  pergenome_attribution_prob
}
my_attribution = extract_attribution(mysigs, sigfit_means=sigbybin, sampleinfo=sampbybin)



# 3) create mutation list, with each element being the mutations in one 1Mb bin
# each element should be named the SAME as in the attribution
focalsample = read.table(paste(root, "/lymphocyteWGS/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/chr_pos_ref_alt.memoryB_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2.txt", sep=""), stringsAsFactors = F, header=T)

binsize=1000000
chromsize = read.table(paste(root, "/reference_files/GRCh37d5/hg19.chrom_sizes.txt", sep=""), stringsAsFactors = FALSE)
mychroms = chromsize[,1] # alternatively, c(1:22, "X","Y")
Nregions = sum(chromsize[,3]/binsize) + length(mychroms)
#mydf = matrix(nrow=96, ncol=Nregions)
regionstart = regionend = regionchrom = vector()
mymut_list = list()
mynames = c()

count=0
for (chrom in 1:length(mychroms)){
  focalchrom = focalsample[focalsample$chr==mychroms[chrom],]
  bins = seq(from=1, to=chromsize[chrom,3], by=binsize) # add last bin
  bins[length(bins)+1] = bins[length(bins)] + binsize # add last bin
  for (j in 1:(length(bins)-1)){
    count = count+1
    regionstart[count] = bins[j]
    regionend[count] = bins[j+1]
    regionchrom[count] = mychroms[chrom]
    focalbin = focalchrom[focalchrom$pos >= bins[j] & focalchrom$pos < bins[j+1], ]
    mymut_list[[count]] = focalbin
    mynames[count] = paste(mychroms[chrom], bins[j], sep="_")
  }
}
names(mymut_list) = mynames



# 4) apply attribution function per 1Mb bin
apply_attribution = function(attribution, mutlist){
  genomeFile = paste(root, "/reference_files/GRCh37d5/genome.fa", sep="")
  require("GenomicRanges")
  require("Rsamtools")
  require("MASS")
  #attribution = my.attribution
  #mutlist = somlist2
  
  mutlistnames = names(mutlist)
  attributionnames = names(attribution)
  
  muts_prob_per_sig_list = list()
  for (i in 1:length(mutlist)){ ## for each genome
    if ( (mutlistnames[i] %in% attributionnames) & nrow(mutlist[[i]])>0 ){ ## if this genome is in the attribution list (by name)
      focal_mut = mutlist[[i]]
      focal_att = attribution[[mutlistnames[i]]]
      
      # get trinucleotide context
      subs_only = focal_mut
      colnames(subs_only) = c("chr","pos","ref","mut")
      subs_only = subs_only[(subs_only$ref %in% c("A","C","G","T")) & (subs_only$mut %in% c("A","C","G","T")) & subs_only$chr %in% c(1:22,"X","Y"),]
      subs_only$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(subs_only$chr, IRanges(subs_only$pos-1, subs_only$pos+1))))
      
      # 2. Annotating the mutation from the pyrimidine base
      ntcomp = c(T="A",G="C",C="G",A="T")
      subs_only$sub = paste(subs_only$ref,subs_only$mut,sep=">")
      subs_only$trinuc_ref_py = subs_only$trinuc_ref
      for (j in 1:nrow(subs_only)) {
        if (subs_only$ref[j] %in% c("A","G")) { # Purine base
          subs_only$sub[j] = paste(ntcomp[subs_only$ref[j]],ntcomp[subs_only$mut[j]],sep=">")
          subs_only$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(subs_only$trinuc_ref[j],split="")[[1]])],collapse="")
        }
      }
      subs_only$sub_trinuc_ref_py = paste(subs_only$sub, subs_only$trinuc_ref_py, sep="_")
      muts_prob_per_sig = merge(subs_only, data.frame(sub_trinuc_ref_py=rownames(focal_att), focal_att), by="sub_trinuc_ref_py")
      muts_prob_per_sig_list[[i]] = muts_prob_per_sig[order(muts_prob_per_sig$chr, muts_prob_per_sig$pos),c(2:8,1,9:ncol(muts_prob_per_sig))]
    }
  }
  muts_prob_per_sigDF = do.call(rbind, muts_prob_per_sig_list)
  muts_prob_per_sigDF
}
my_attribution_perbp = apply_attribution(attribution=my_attribution, mutlist=mymut_list)


## write final per bp attribution 
write.table(my_attribution_perbp, paste(root, "/lymphocyteWGS/mutsig_byregion/sigfit/sigfit_memoryB_1Mb_min10percent_IgDenovo/attribution_perbp_sigfit_memoryB_1Mb_min10percent_IgDenovo.txt", sep=""), col.names = T, row.names = F, quote=F, sep="\t")
