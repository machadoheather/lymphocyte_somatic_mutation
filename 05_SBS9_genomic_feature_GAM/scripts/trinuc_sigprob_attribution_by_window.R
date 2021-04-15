## Feb 2020
## Calculate S9 prob per window (e.g. 1kb)

## Can run on farm or locally
#setwd("/lustre/scratch116/casm/cgp/users/hm8/lymphocyteWGS/mutsig_byregion/regression_analysis")
#args = c("/lustre/scratch116/casm/cgp/users/hm8/lymphocyte_somatic_mutations_final/results/trinuc_sigprob_attribution.ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg.nopriors_July2019.ALL.txt", "1000", "S9")

args = commandArgs(trailingOnly = TRUE)
trinuc_file = args[1]
window_size = as.numeric(args[2])
signature = args[3]

# Read in data
#trinuc = read.table("/lustre/scratch116/casm/cgp/users/hm8/lymphocyte_somatic_mutations_final/results/trinuc_sigprob_attribution.ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg.nopriors_July2019.ALL.txt", header=TRUE, stringsAsFactors = F, sep="\t")

# Only use memory B cells in the S9 calculation (other cells just introduce noise)
trinuc = read.table(trinuc_file, header=TRUE, stringsAsFactors = F)
bmem = trinuc


#################### 1KB Jan2020
## Calculate S9 prob per 1KB window, 1 chrom at a time
# 1 row per window
# columns: chrom, st, end, Nmut, Nsamples, meanS9, sumS9

mychroms = c(1:22, "X","Y")
#allMean = mean(bmem[,colnames(bmem) ])
focalsig_column = which(colnames(bmem)==signature)

for (i in 1:length(mychroms)){
    #i= 14
  focalchrom = bmem[bmem$chr==mychroms[i],]
  bins = seq(from=1, to=max(focalchrom$pos), by=window_size) # ignore last bin
  chromlist = list()
  for (j in 1:(length(bins)-1)){
    mystart = bins[j]
    myend = bins[j+1]-1
    focalbin = focalchrom[focalchrom$pos >=mystart & focalchrom$pos < myend,  ]
    Nmut = nrow(focalbin)
    Nsamples = length(unique(focalbin$sample))
    if (Nmut==0){
      meanS9 = 0
      sumS9 = 0
    } else{
      meanS9 = mean(focalbin[,focalsig_column])
      sumS9 = sum(focalbin[,focalsig_column])
    }
    # if (Nmut < 2){
    #   se = NA
    #   pS9 = NA
    # } else{
    #   se = sd(focalbin$S9)/sqrt(Nmut)
    #   pS9 = try(t.test(focalbin$S9, alternative = c("greater"), mu=allMean)$p.value, silent =TRUE)
    # }
    # if (typeof(pS9)=="character"){
    #   pS9 = NA
    # }
    # #chromlist[[j]] = c(mystart, myend, Nmut, Nsamples, meanS9, sumS9, se, pS9)
    chromlist[[j]] = c(mystart, myend, Nmut, Nsamples, meanS9, sumS9)
  }
  #df1 = data.frame(na.omit(do.call(rbind, chromlist)))
  df1 = data.frame(do.call(rbind, chromlist))
  df1$chr = mychroms[i]
  colnames(df1) = c("start", "end", "Nmut",  "Nsamples",paste("mean",signature, sep=""), paste("sum",signature, sep=""),"chr")
  if (nrow(df1)==0){ warning("Zero rows resulting") }
  save(df1, file=paste("int_files/", signature, "attributed_window", window_size, "bp_chr", mychroms[i], ".Rdata", sep="") )
}
