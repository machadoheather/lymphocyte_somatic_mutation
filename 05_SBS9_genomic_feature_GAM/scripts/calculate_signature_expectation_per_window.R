## Feb 2020

## Calculating per window expected S9 contribution
# 1) fetch S9 profile
# 2) fetch window trinucleotide list
# 3) match window trinucleotide and S9 profile
# 4) calculate statistic over window

# library("GenomicRanges")
# library("Rsamtools")
# library("MASS")
# 
# # load in luad_multi object (from main analysis)
# load(file="/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/luad_multi_first.noprior.ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg.Rdata")
# 
# # get focal window: chrom, start, end
# windows = data.frame(c("1","1"), c(100001, 101001), c(101000, 102000) )
# colnames(windows) = c("chr","start","end")
# windows$size = windows$end - windows$start + 1
# 
# sigprop_trinuc = fetch_sig_trinuc_profile(luad_multi)
# windows_S9exp = per_window_sig_exp1(windows, sigprop_trinuc)
# 

####################################################################################
############################ FUNCTIONS
####################################################################################
fetch_sig_trinuc_profile = function(signature){
  # Uses object persig_tri and mut_count (pre-load)
  # 1) fetch S9 profile
  ###############################################
  ## Extract components first
  #luad_multi2 <- hdp_extract_components(luad_multi)
  
  ## proportion of each trinuc per signature 
  #comp_trinuc = comp_categ_distn(luad_multi2)
  #persig_tri = t(comp_trinuc$mean)
  #colnames(persig_tri) = c("X0","S1","Sblood","S9","S8","S7","N1","N2","S17")
  focal_sig = which(colnames(persig_tri) == signature)
  
  # labels for the trinucs
  trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
  group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                                each=16))
  sub_trinuc_context = paste(group_factor, trinuc_context, sep="_")
  rownames(persig_tri) = sub_trinuc_context
  sigprop = data.frame(trinuc_context = trinuc_context, S9prop = persig_tri[,focal_sig])
  
  # each trinucleotide appears three times (three possible mutations away)
  # average over those three
  sigprop_trinuc_res = tapply(sigprop[,2], INDEX=sigprop[,1], FUN=sum)
  sigprop_trinuc = data.frame(trinuc = names(sigprop_trinuc_res), prop = sigprop_trinuc_res)
  
  # normalize this to an average of 1 per trinucleotide (rather sum of 1 across trinucleotides)
  n_trinuc = length(unique(sigprop_trinuc$trinuc))
  if (n_trinuc != 32) {warning(paste(n_trinuc, " unique trinucleotides (should be 32)", sep=""))}
  sigprop_trinuc$prop_norm = sigprop_trinuc$prop*32
  sigprop_trinuc
}


# windows=nonzero_mut_ig
# sigprop_trinuc

per_window_sig_exp1 = function(windows, sigprop_trinuc){
# 2) fetch window trinucleotide list
  ###############################################
  # windows object (df) must have columns (named):
  # chr
  # start
  # end
  require("GenomicRanges")
  require("Rsamtools")
  require("MASS")
  genomeFile = "/lustre/scratch116/casm/cgp/users/hm8/reference_files/GRCh37d5/genome.fa"
  
  windows$size = windows$end - windows$start + 1
  # get trinucleotide for each position
  trinuc_list = list()
  for(w in 1:nrow(windows)){
    focal_window = windows[w,]
    focal_window_df = suppressWarnings(data.frame(chr=focal_window$chr, pos=as.numeric(focal_window$start):as.numeric(focal_window$end), stringsAsFactors=FALSE))
    trinuc = as.vector(scanFa(genomeFile, GRanges(focal_window_df$chr, IRanges(focal_window_df$pos-1, focal_window_df$pos+1))))
    nuc = as.vector(scanFa(genomeFile, GRanges(focal_window_df$chr, IRanges(focal_window_df$pos, focal_window_df$pos))))
    
    # annotating from the pyrimidine base
    ntcomp = c(T="A",G="C",C="G",A="T")
    trinuc_py = trinuc
    for (j in 1:length(trinuc)) {
      if (nuc[j] %in% c("A","G")) { # Purine base
        trinuc_py[j] = paste(ntcomp[rev(strsplit(trinuc[j],split="")[[1]])],collapse="")
      }
    }
    trinuc_list[[w]] = trinuc_py
  }
  
  
  ###############################################
  # 3) match window trinucleotide and S9 profile
  ###############################################
  sig_exp1 = vector()
  for (i in 1:length(trinuc_list)){
    # merges trinuc of the window with sig probabilities
    # will add a line for each position of the window that does not match a trinuc and list the prop as "NA"
    m1 = merge(data.frame(trinuc_list[[i]]), sigprop_trinuc, by=1, all.x = TRUE)
    
    # divide by the window size to get the expected number of mutations (based on 1 mutation in a window of equal trinuc props)
    sig_exp1[i] = sum(m1$prop_norm)/windows$size[i]
  }
  windows$sig_exp1 = sig_exp1
  windows
}

