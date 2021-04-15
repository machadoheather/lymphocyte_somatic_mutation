## Aug 2020

# comparing blood signatures
# 1) one extracted from sigprofiler,  (date: )
# 2) one extracted from sigfit, Aug 28, 2020

setwd("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/bloodsignature_analysis/BloodSig_analysis_Aug2020")
library("lsa")

######## Loading PCAWG signatures (to get S1)
cosmic.sigs <- read.table('/Users/hm8/sanger/example_code/nicola_hdp/signatures_probabilities.txt', header=TRUE, sep='\t')
cosmic.sigs <- cosmic.sigs[order(cosmic.sigs$Substitution.Type, cosmic.sigs$Trinucleotide),]
pcawg.sigs <- read.csv('/Users/hm8/sanger/example_code/nicola_hdp/sigProfiler_SBS_signatures_2018_03_28.csv', header=TRUE)
#  sort by Substitution Type and Trinucleotide
pcawg.sigs <- pcawg.sigs[order(pcawg.sigs$Type, pcawg.sigs$SubType),]
prior_sigsP1 <- as.matrix(pcawg.sigs[,grep('SB', colnames(pcawg.sigs))])
# number of prior signatures to condition on (65)
# each col (signatures) must sum to 1 (over all mutation types)
colsum1 = apply(prior_sigsP1, MARGIN=2, FUN=sum)
prior_sigsP = prior_sigsP1
for (i in 1:length(colsum1)){
  prior_sigsP[,i] = prior_sigsP1[,i]/colsum1[i]
}
S1 = prior_sigsP[,1]
barplot(S1)
# pdf("S1_cosmic3.pdf", width=6,height=4)
# plotspectrum(prior_sigsP[,1])
# dev.off()
# pdf("S1_cosmic2.pdf", width=6,height=4)
# plotspectrum(cosmic.sigs[,"Signature.1"])
# dev.off()
pdf("S5_cosmic2.pdf", width=6,height=4)
plotspectrum(cosmic.sigs[,"Signature.5"])
dev.off()

# 1) sigprofiler
bloodsig = read.table("/Users/hm8/volumes/hm8_network/lymphocyteWGS/metaAnalysis_pcawg_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg_June2020/sigprofiler/sig5_12/bloodsig_SBS96D.txt", header=T, stringsAsFactors = F, sep="\t")
colnames(bloodsig)[2] = "SBSblood"
# put in order of cosmic sigs
rownames(bloodsig) = bloodsig$MutationsType
bloodsigO = bloodsig[cosmic.sigs$Somatic.Mutation.Type,]


# 2a) sigfit
load("/Users/hm8/volumes/hm8_network/lymphocyteWGS/bloodsignature_analysis/BloodSig_analysis_Aug2020/OLD/exposures_hsc_S1_extra1.sigfit.Rdata")
head(myexposures$mean)
head(mysigs$mean)
mysigs2 = t(mysigs$mean)
mysigsO = mysigs2[cosmic.sigs$Somatic.Mutation.Type, ]
mysigsO = mysigs2
# save new blood sig from sigfit
newbloodsig = data.frame(Somatic.Mutation.Type = bloodsigO[,1], Signature.Blood = mysigsO[,2])
write.table(newbloodsig, file="sigfit_cosmic3_bloodsig_Aug2020.txt", col.names=T, row.names = F, quote=F, sep="\t")
pdf("Sblood_sigfitextra1_cosmic3S1.pdf", width=6,height=4)
plotspectrum(newbloodsig$Signature.Blood)
dev.off()


# 2b) sigfit
load("/Users/hm8/volumes/hm8_network/lymphocyteWGS/bloodsignature_analysis/BloodSig_analysis_Aug2020/exposures_hsc_cosmicv2_S1_extra1.sigfit.Rdata")
head(myexposures$mean)
head(mysigs$mean)
mysigs2 = t(mysigs$mean)
#mysigsO = mysigs2[cosmic.sigs$Somatic.Mutation.Type, ]
mysigsO = mysigs2
# save new blood sig from sigfit
newbloodsigv2 = data.frame(Somatic.Mutation.Type = bloodsigO[,1], Signature.Blood = mysigsO[,2])
write.table(newbloodsigv2, file="sigfit_cosmic2_bloodsig_Aug2020.txt", col.names=T, row.names = F, quote=F, sep="\t")
pdf("Sblood_sigfitextra1_cosmic2S1.pdf", width=6,height=4)
plotspectrum(newbloodsigv2$Signature.Blood)
dev.off()

## blood sigs are VERY similar regardless of using S1 from cosmic2 or cosmic3
cosine(newbloodsigv2$Signature.Blood, newbloodsig$Signature.Blood) #0.9937693
par(mfrow=c(2,1))
barplot(newbloodsig[,2])
barplot(newbloodsigv2[,2])


# 3b) sigfit S1 + S5
load("/Users/hm8/volumes/hm8_network/lymphocyteWGS/bloodsignature_analysis/BloodSig_analysis_Aug2020/exposures_hsc_S1_S5_extra1.sigfit.Rdata")
head(myexposures$mean)
head(mysigs$mean)
mysigs2 = t(mysigs$mean)
#mysigsO = mysigs2[cosmic.sigs$Somatic.Mutation.Type, ]
mysigsO = mysigs2
# save new blood sig from sigfit
newbloodsigv2S5v3 = data.frame(Somatic.Mutation.Type = bloodsigO[,1], Signature.Blood = mysigsO[,3])
write.table(newbloodsigv2S5v3, file="sigfit_cosmic3S1S5_bloodsig_Aug2020.txt", col.names=T, row.names = F, quote=F, sep="\t")
pdf("Sblood_sigfitextra1_cosmic3S1S5.pdf", width=6,height=4)
plotspectrum(newbloodsigv2S5v3$Signature.Blood)
dev.off()


# 3a) sigfit S1 + S5
load("/Users/hm8/volumes/hm8_network/lymphocyteWGS/bloodsignature_analysis/BloodSig_analysis_Aug2020/exposures_hsc_cosmicv2_S1_S5_extra1.sigfit.Rdata")
head(myexposures$mean)
head(mysigs$mean)
mysigs2 = t(mysigs$mean)
#mysigsO = mysigs2[cosmic.sigs$Somatic.Mutation.Type, ]
mysigsO = mysigs2
# save new blood sig from sigfit
newbloodsigv2S5 = data.frame(Somatic.Mutation.Type = bloodsigO[,1], Signature.Blood = mysigsO[,3])
write.table(newbloodsigv2S5, file="sigfit_cosmic2S1S5_bloodsig_Aug2020.txt", col.names=T, row.names = F, quote=F, sep="\t")
pdf("Sblood_sigfitextra1_cosmic2S1S5.pdf", width=6,height=4)
plotspectrum(newbloodsigv2S5$Signature.Blood)
dev.off()

par(mfrow=c(3,1))
plotspectrum(mysigsO[,1])
plotspectrum(mysigsO[,2])
plotspectrum(mysigsO[,3])

# comparison
cosine(bloodsigO$SBSblood, mysigsO[,2]) # 0.9366
cosine(newbloodsigv2S5[,2], newbloodsigv2[,2]) # 0.9727061

# plotting the signatures S1, S5, and Sblood for HSCs
myexposuresS1S5 = myexposures$mean
colnames(myexposuresS1S5) = c("S1","S5","Sblood")
myexposuresS1S5$colony = rownames(myexposuresS1S5)
m3 = read.table(file="/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis//metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/colonyinfo_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
sig_prop_type = merge(m3, myexposuresS1S5, by="colony")
sig_prop_type_stack = sig_prop_type[,c(1,3,7,13:(ncol(sig_prop_type)))]
sig_prop_type_stack$Age[sig_prop_type_stack$Age==27] = 29


sig_prop_type_stackO = sig_prop_type_stack[order(sig_prop_type_stack[,4], decreasing = TRUE),]
sig_prop_type_stackO$colony = factor(sig_prop_type_stackO$colony, levels=unique(sig_prop_type_stackO$colony))
sig_prop_type2_stack = melt(sig_prop_type_stackO, id.vars=c("colony", "Donor", "Age")) # 
sig_prop_type2_stack$variable = factor(sig_prop_type2_stack$variable, levels=unique(sig_prop_type2_stack$variable))
#sig_prop_type2_stack$variable = factor(sig_prop_type2_stack$variable, levels = levels(sig_prop_type2_stack$variable)[c(length(levels(sig_prop_type2_stack$variable)),1:(length(levels(sig_prop_type2_stack$variable))-1) )])

my10 = c(brewer.pal(10, "Set1"),"black") # max for this set is 9 
sigcol_values = c("S1" = my10[2],  "Sblood"=my10[3], "S5" = my10[1], "SBS7a"=my10[6], "SBS8"=my10[5], "SBS9"=my10[4],"SBS17b"=my10[8],"SBS18"=my10[9],  "C0"="white")

ggplot(subset(sig_prop_type2_stack ), aes(x=colony, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  facet_wrap(.~Age, scales = "free_x", ncol=5)+
  xlab("") +
  ylab("Proportion per genome") +
  scale_fill_manual("",values=sigcol_values)+
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        #axis.text = element_tex9t(size = 10),
        axis.line.x=element_blank(),
        strip.text = element_text(size = 10),
        strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")) 
  )
#ggsave(file=paste("figures/signatures_", name_sigs, "_", groupname, "_cellType_stackedBar_cancercompare_Prop_names_labels_allsigs.pdf", sep=""), width=12,height=4)


###### plotting
plotspectrum = function(freqs, samplename=NULL){
  if (is.null(samplename)){
    samplename = ""
  }
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  #freqs_full = freqs[full_vec]
  freqs_full = freqs
  freqs_full[is.na(freqs_full)] = 0
  names(freqs_full) = full_vec
  y = freqs_full; maxy = max(y)
  h=barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="# SNVs")
  mtext(side=4, text=samplename)
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j-1)*16+1,j*16)]
    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colvec[j*16])
    text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
  } 
}


pdf("S1_Sblood_compare.pdf", width=12,height=12)
par(mfrow=c(4,2))
plotspectrum(mysigsO[,1], "sigfit S1")
plotspectrum(mysigsO[,2], "sigfit denovo 1")
plotspectrum(S1, "pcawg S1")
plotspectrum(bloodsigO$SBSblood, "sigprofiler blood")
barplot(mysigsO[,1], main="sigfit S1 + denovo 1", border=NA, col=rgb(0,0,1,0.5))
barplot(mysigsO[,2], border=NA, col=rgb(1,0,0,0.5), add=TRUE)
barplot(S1, main="pcawg S1 + sigprofiler blood", border=NA, col=rgb(0,0,1,0.5))
barplot(bloodsigO$SBSblood, border=NA, col=rgb(1,0,0,0.5), add=TRUE)
barplot(bloodsigO$SBSblood, border=NA, col=rgb(1,0,0,0.5))
barplot(mysigsO[,2], main="pcawg S1 + sigprofiler blood", border=NA, col=rgb(0,0,1,0.5), add=TRUE)
frame()
legend("center",legend=c("sigfit denovo 1", "sigprofiler blood"), fill=c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)) )
dev.off()


######## RANDOM OTHER PLOTS
pdf("S9_cosmic2.pdf", width=6,height=4)
plotspectrum(cosmic.sigs[,"Signature.9"])
dev.off()

pdf("S7_cosmic2.pdf", width=6,height=4)
plotspectrum(cosmic.sigs[,"Signature.7"])
dev.off()

pdf("S9_cosmic3.pdf", width=6,height=4)
plotspectrum(pcawg.sigs[,"SBS9"])
dev.off()

load(file="/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/mutsig_byregion/mutmats_byregion.Rdata")
Signature.myIG=mutmaxIG[,1]/sum(mutmaxIG[,1])
pdf("myIg.pdf", width=6,height=4)
plotspectrum(Signature.myIG)
dev.off()
## IG defined by: 
# focalbin = focalsample[  (focalsample$chr == "14" & focalsample$pos >106032614 &  focalsample$pos< 107288051) | 
#                            (focalsample$chr == "22" & focalsample$pos >22380474 &  focalsample$pos< 23265085) | 
#                            (focalsample$chr == "2" & focalsample$pos >89890568 &  focalsample$pos< 90274235) |
#                            (focalsample$chr == "2" & focalsample$pos >89156874 &  focalsample$pos< 89630436), ]

library(dplyr)
table(subset(m3, Donor=="KX003")$Cell.type2)
