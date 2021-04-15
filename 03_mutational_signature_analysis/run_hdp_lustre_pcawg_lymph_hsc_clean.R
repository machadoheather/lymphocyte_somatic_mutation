# ---
#   title: "run_hdp_lustre_pcawg_lymph_hsc"
# ---
# Feb 2019

#setwd('/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_pcawg_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg_June2020')

library(hdp)
library(ggplot2)
library(foreach)
library(doParallel)
library(reshape2)


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




###### load in ARG380 data
#setwd("/Users/hm8/sanger/lymphocyteExpansionSequenceAnalysis/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg")
groupname="pcawg_lymph_hsc"

## We know that the colony_info file and the mutcounts_matrix samples are in the same order.
# Otherwise, must correctly order
mutcounts_matrix = read.table(file="/lustre/scratch116/casm/cgp/users/hm8/lymphocyteWGS/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/mutcounts_matrix_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
colonyinfo_all = read.table(file="/lustre/scratch116/casm/cgp/users/hm8/lymphocyteWGS/metaAnalysis_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg/colonyinfo_ARGhscPBonly_KX001_KX002_KX003_tonsilMS_tonsilPS_stemcellCB2_ARGtreg.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")

############# Make samp_type object
Nsamp = ncol(mutcounts_matrix)
samp_type = colonyinfo_all
samp_type$ppindex = 1
samp_type$cpindex = 1+1
samp_type$dpindex = 2:(Nsamp+1)
samp_type$samp_names = samp_type$colony
samp_typeA = samp_type[,c("samp_names","Cell.type2","Donor","Nmut_adj_as")]
colnames(samp_typeA) = c("samp_names","celltype","group","Nmut")
  
############## Read in pcawg focal7 June2020 matrix
mutpcawg = read.table(paste(root, "/pcawg/mutcounts_matrix_pcawgfocal7_mm_Aug2020.txt",sep=""), stringsAsFactors = T, header=T)
samppcawg = read.table(paste(root, "/pcawg/samplesgroups_mutcounts_matrix_pcawgfocal7_mm_Aug2020.txt", sep=""), stringsAsFactors = T, header=F)
colnames(samppcawg) = c("samp_names","celltype")
samppcawg$group = samppcawg$celltype
samppcawg$Nmut = apply(mutpcawg, MARGIN=2, FUN=sum)

# exclude those with more than 100K mutations:
keepsamples = which(samppcawg$Nmut < 60000 & samppcawg$celltype %in% c("Lymph-BNHL","Lymph-CLL","mm","Myeloid-AML") )
#samp_names         celltype            group
#102 f9dc999f-6dde-448d-9cf1-2897ddcf7b0b       Lymph-BNHL       Lymph-BNHL
#216 04aa6b77-8074-480c-872e-a1a47afa5314    Skin-Melanoma    Skin-Melanoma
#219 0ab4d782-9a50-48b9-96e4-6ce42b2ea034    Skin-Melanoma    Skin-Melanoma
#220 0dd0718d-5ddf-4c59-8c47-0f51303daeb5    Skin-Melanoma    Skin-Melanoma
#221 108749d2-5c62-4ef1-92df-aec6941ba53b    Skin-Melanoma    Skin-Melanoma
#236 00aa769d-622c-433e-8a8a-63fb5c41ea42 ColoRect-AdenoCA ColoRect-AdenoCA
#241 0980e7fd-051d-45e9-9ca6-2baf073da4e8 ColoRect-AdenoCA ColoRect-AdenoCA
#244 14c5b81d-da49-4db1-9834-77711c2b1d38 ColoRect-AdenoCA ColoRect-AdenoCA
#245 154f80bd-984c-4792-bb89-20c4da0c08e0 ColoRect-AdenoCA ColoRect-AdenoCA
#Nmut
#102  112529
#216  744553
#219  214810
#220  159349
#221  105876
#236  234231
#241  849282
#244 2426470
#245  259867

mutpcawg2 = mutpcawg[,keepsamples]
samppcawg2 = samppcawg[keepsamples,]


####### Combine two datasets
mutcounts_matrix_both = cbind(mutcounts_matrix, mutpcawg2)
samp_type = rbind(samp_typeA, samppcawg2)
samp_names = samp_type$samp_names


####### testing
# mutcounts_matrix_both = mutcounts_matrix_both[,1:10]
# samp_type = samp_type[1:10,]
# samp_names = samp_names[1:10]



######### Running hdp without priors
####### Set the indicies for the samples
groupname="pcawg_lymph_hsc"

samp_type$ppindex = 1
samp_type$cpindex = 1+1

####### Run the program
### All samples
hdp_mut <- hdp_init(ppindex = c(0, samp_type$ppindex),  # index of parental nodes
                    cpindex = c(1, samp_type$cpindex), # index of concentration param
                    hh=rep(1, 96), # prior is uniform over the 96 mutation categories
                    alphaa=rep(1, max(samp_type$cpindex)), # shape hyperparams for different CPs
                    alphab=rep(1, max(samp_type$cpindex))) # rate hyperparams for different CPs

hdp_mut <- hdp_setdata(hdp_mut,
                       dpindex = 2:numdp(hdp_mut),   # index of nodes to add data to
                       t(mutcounts_matrix_both[,1:ncol(mutcounts_matrix_both)]) )   ## flip so samples are in rows

samp_type$dpindex = 2:numdp(hdp_mut)



# Run multiple posterior sampling chains
cl <- makeCluster(4) 
registerDoParallel(cl)


#chlist <- vector("list", 4)

chlist = foreach(i=1:4) %dopar% {
  require(hdp)
      ###### run in parallel, remember to request enough cores when running (plus one extra?)
  hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10, seed=i*200)
  chlist_part <- hdp_posterior(hdp_activated,
                               burnin=20000,
                               n=50,
                               space=200,
                               cpiter=3,
                               seed=i*1e3)
  chlist_part
}


luad_multi <- hdp_multi_chain(chlist)
name_sigs = "noprior"
my_prior=""
save(luad_multi, samp_type, file=paste("luad_multi_first.", name_sigs, ".", groupname, ".Rdata", sep="") )



############# More extraction
luad_multi <- hdp_extract_components(luad_multi)
dp_distn = comp_dp_distn(luad_multi)
comp_trinuc = comp_categ_distn(luad_multi)
Nsig = ncol(dp_distn$mean)
trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))
inds = samp_type$dpindex
sig_prop = data.frame(dp_distn$mean[samp_type$dpindex, ])
samp_names = samp_type$samp_names
sig_prop$samples = samp_names
Nsamp = length(samp_names)
dp_distn_lowerCI = data.frame(matrix(ncol=ncol(dp_distn$cred.int[[1]]), nrow=Nsamp))
dp_distn_upperCI = data.frame(matrix(ncol=ncol(dp_distn$cred.int[[1]]), nrow=Nsamp))
for (i in inds){
  focal = dp_distn$cred.int[[i]]
  dp_distn_lowerCI[i-(min(inds)-1), ] = focal[1,] 
  dp_distn_upperCI[i-(min(inds)-1), ] = focal[2,] 
}
dp_distn_lowerCI$samples = samp_names
dp_distn_upperCI$samples = samp_names
colnames(dp_distn_lowerCI) = colnames(sig_prop)
colnames(dp_distn_upperCI) = colnames(sig_prop)

s1 = melt(sig_prop[,2:ncol(sig_prop)], id=c("samples") )
s2 = melt(dp_distn_lowerCI[,2:ncol(sig_prop)], id=c("samples") )
s3 = melt(dp_distn_upperCI[,2:ncol(sig_prop)], id=c("samples") )
sig_propCI = cbind(s1, s2$value, s3$value)
colnames(sig_propCI) = c("samples","variable","freq","lower","upper")
sig_propCI$samples = factor(sig_propCI$samples, levels(as.factor(sig_propCI$samples)))

#### Looking at the contributions of each signature to the different cell types
#sig_prop_type = merge(m3, sig_prop, by.x="colony", by.y="samples")
save(luad_multi, samp_type, sig_propCI, dp_distn, comp_trinuc, Nsig, file=paste("luad_multi.", name_sigs, ".", groupname, ".Rdata", sep="") )
save(sig_propCI,samp_type,sig_prop, file=paste("misc_samp_prop_objects.", groupname, ".Rdata", sep="") )

# ### Checking burn-in, etc
# par(mfrow=c(2,2))
# p1 <- lapply(chains(luad_multi), plot_lik, bty='L', start=10) # check for convergence (if not, increase burnin)- needs 20K
# p2 <- lapply(chains(luad_multi), plot_numcluster, bty='L')
# p3 <- lapply(chains(luad_multi), plot_data_assigned, bty='L')
# # Only needs 10K, for future reference
# 
# # Extract components (mutational signatures)
# mut_colours <- c(RColorBrewer::brewer.pal(10, 'Paired')[seq(1,10,2)], 'grey70')
# par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
# plot_comp_size(luad_multi, bty="L")
# trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
# group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
#                               each=16))
# Ncolumns = ceiling(sqrt(Nsig))
# Nrows = ceiling(Nsig/Ncolumns)
# pdf(file=paste("signatures",name_sigs,"_prior", my_prior, "_", groupname, "_componentsSpectra.pdf", sep=""), width=Ncolumns*3, height=Nrows*2)
# par(mfrow=c(Nrows,Ncolumns))
# par(mar=c(3,3,3,1))
# plot_comp_distn(luad_multi, cat_names=trinuc_context,
#                   grouping=group_factor, col=mut_colours,
#                   col_nonsig="grey80", show_group_labels=TRUE)
# dev.off()
# 
# # Save mutational spectrum for components
# comp_distn <- comp_categ_distn(luad_multi)
# write.table(t(comp_distn$mean), file=paste("mutspectrum_components.comp_distn_mean.", name_sigs, ".", groupname, ".Rdata", sep=""), col.names = T, row.names = F, quote=F)
#     
# 




# ############################################################################
# ############ RE-START HERE ##############
# groupname="ARGmemory_pcawgfocal7_mm_ptcl"
# name_sigs="noprior"
# #save(luad_multi, samp_type, sig_propCI, file=paste("luad_multi.", name_sigs, ".", groupname, ".Rdata", sep="") )
# load(file=paste("luad_multi.", name_sigs, ".", groupname, ".Rdata", sep="") )
# #save(Nrows,Ncolumns,samp_type,sig_prop, file=paste("misc_samp_prop_objects.", groupname, ".Rdata", sep="") )
# load(file=paste("misc_samp_prop_objects.", groupname, ".Rdata", sep="") )
# ############################################################################
# ############################################################################
# 
# 
# #### Looking at the contributions of each signature to the different cell types
# sig_prop_type = merge(samp_type, sig_prop, by.x="samp_names", by.y="samples")
# mygroups = unique(sig_prop_type$celltype)
# Ngroups = length(mygroups)
# 
# # Re-ordering the signatures and re-naming based on cosmic similarity
# colnames(sig_prop_type)[9:22] = c("S7","S5","X3","S1","S17","X6","Sblood","S9","X9","S18","X11","S3","X13","S2")
# 
# 
# 
# ## Add subtype info for BNHL and CLL
# IDconversion = read.table("/Users/hm8/sanger/pcawg/release_may2016.v1.4.tsv", sep="\t", header=TRUE)
# bnhl_types = read.table("/Users/hm8/sanger/pcawg/pcawg_specimen_histology_August2016_v9_BCL.txt", sep="\t", header=TRUE)
# # max1 = max(apply(IDconversion, 2, FUN=function(x) mean(sig_prop_type$samp_names %in% x) ) )
# # which(apply(IDconversion, 2, FUN=function(x) mean(sig_prop_type$samp_names %in% x) ) ==max1)
# # # IDconversion$tumor_wgs_aliquot_id
# # which(apply(IDconversion, 2, FUN=function(x) mean(bnhl_types$icgc_specimen_id %in% x) ) ==1)
# # # IDconversion$tumor_wgs_icgc_specimen_id is the column that matches the 
# merge1 = merge(bnhl_types[,c("icgc_specimen_id","histology_tier4")],  IDconversion[,c("tumor_wgs_icgc_specimen_id", "tumor_wgs_aliquot_id")], by=1)
# merge1B = merge1[merge1$histology_tier4 %in% c("Burkitt lymphoma","Diffse large B-cell lymphoma","Follicular lymphoma") ,c("tumor_wgs_aliquot_id","histology_tier4")]
# # omitting Marginal zone B-cell lymphoma (2) and Post-transplant lymphoproliferative disorder (1)
# colnames(merge1B) = c("samp_names","histology_tier4")
# 
# # CLL
# cll_types = read.table("/Users/hm8/sanger/francesco_MM_PTCL/CLL/sample_with_VDJ.txt", stringsAsFactors = F, sep="\t", header=T)
# cll_typesB = cll_types[,c("pancan_ID","IGHV_status_BCN")]
# colnames(cll_typesB) = c("samp_names","histology_tier4")
# 
# # Both
# mergeboth = rbind(cll_typesB, merge1B)
# merge2 = merge(sig_prop_type, mergeboth, by="samp_names", all.x=T)
# merge2$celltype = as.character(merge2$celltype)
# merge2$histology_tier4 = as.character(merge2$histology_tier4)
# 
# # Remove CLL without mutated status (7 genomes)
# merge3 = merge2[!(merge2$celltype=="Lymph-CLL" & is.na(merge2$histology_tier4)), ]
# table(merge3$histology_tier4)
# # Burkitt lymphoma Diffse large B-cell lymphoma          Follicular lymphoma 
# # 17                           51                           36 
# # MUT                        UNMUT 
# # 38                           45 
# merge3$celltype[merge3$celltype=="Lymph-BNHL" | merge3$celltype=="Lymph-CLL"] = merge3$histology_tier4[merge3$celltype=="Lymph-BNHL" | merge3$celltype=="Lymph-CLL"]
# sig_prop_type = merge3[ ,1:37]
# 
# 
# ############# Plotting contributions
# sig_prop_type2 = melt(sig_prop_type[,c(2,(ncol(samp_type)+2):ncol(sig_prop_type))]) # not unattributed
# sig_prop_type2$variable = factor(sig_prop_type2$variable, levels=unique(sig_prop_type2$variable)[c(7,4,2,8,5,1,14,12,3,10, 6,9,11,13,15:29)] )
# #sig_prop_type2$celltype = factor(sig_prop_type2$celltype, levels=c("B Memory","Lymph-BNHL","Lymph-CLL","mm","Myeloid-AML","Myeloid-MPN","Skin-Melanoma","Eso-AdenoCa","T Memory","ptcl","Lymph-NOS") )
# sig_prop_type2$celltype = factor(sig_prop_type2$celltype, levels=c("B Memory","Burkitt lymphoma","Diffse large B-cell lymphoma", "Follicular lymphoma", "MUT","UNMUT","mm","Myeloid-AML","Myeloid-MPN","Skin-Melanoma","Eso-AdenoCa","T Memory","ptcl","Lymph-NOS") )
# 
# sig_prop_type2_known = sig_prop_type2[sig_prop_type2$variable %in% c("S7","S5","X3","S1","S17","Sblood","S9","S2","S3","S18"), ]
# sig_prop_type2_known$variable = factor(sig_prop_type2_known$variable, levels=c("Sblood","S1","S5","S9","S17","S7","S2","S3","S18","X3") )
# #sig_prop_type2_known$celltype = factor(sig_prop_type2_known$celltype, levels=c("B Memory","Lymph-BNHL","Lymph-CLL","mm","Myeloid-AML","Myeloid-MPN","T Memory","ptcl","Lymph-NOS","Skin-Melanoma","Eso-AdenoCa") )
# 
# my_median = data.frame(sig_prop_type2_known) %>% group_by(variable, celltype) %>% summarise(median = median(value))
# 
# # Edited scatterplot or proportions with names and medians
# library(RColorBrewer)
# my10 = c(brewer.pal(9, "Set1"),"black") # max for this set is 9 
# my_colors = c(my10[2:9],my10[1],"black")
# # ggplot(sig_prop_type2_known, aes(x=variable, y=value, fill=variable)) +
# #   geom_point(pch=21, alpha=0.5) +
# #   facet_wrap(.~celltype, nrow=2) +
# #   #geom_crossbar(data=my_median,aes(x=variable,ymin=median, ymax=median,y=median,group=celltype), width = 0.5, fatten = 1) +
# #   xlab("Mutational signature") +
# #   ylab("Proportion per genome") +
# #   scale_fill_manual(values=my_colors) +
# #   #scale_fill_brewer(palette = "Set1") +
# #   guides(fill=FALSE)+
# #   #scale_x_discrete(labels=c("Unattributed","Novel1 (HSC)", "S1 (aging)", "S9 (AID)", "S7 (UV)", "Novel2 (S8)", "Novel3","Novel4","S17")) +
# #   theme_light() +
# #   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
# #         axis.ticks.x=element_blank(),
# #         axis.text = element_text(size = 10),
# #         axis.line.x=element_blank(),
# #         legend.justification=c(1,0), legend.position=c(1,0.5),
# #         panel.spacing = unit(1, "lines"),
# #         strip.text = element_text(size = 14),
# #         strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")) )
# # ggsave(file=paste("figures/signaturesno_prior_", groupname, "_cellType_scatterProp_names_knownsigs.pdf", sep=""), width=10,height=4)
# 
# # ggplot(sig_prop_type2, aes(x=variable, y=value, fill=variable)) +
# #   geom_point(pch=21, alpha=0.5) +
# #   facet_rep_wrap(.~celltype, nrow=3, repeat.tick.labels = TRUE) +
# #   #geom_crossbar(data=my_median,aes(x=variable,ymin=median, ymax=median,y=median,group=celltype), width = 0.5, fatten = 1) +
# #   xlab("Mutational signature") +
# #   ylab("Proportion per genome") +
# #   #scale_fill_brewer(palette = "Set1") +
# #   guides(fill=FALSE)+
# #   #scale_x_discrete(labels=c("Unattributed","Novel1 (HSC)", "S1 (aging)", "S9 (AID)", "S7 (UV)", "Novel2 (S8)", "Novel3","Novel4","S17")) +
# #   theme_light() +
# #   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
# #         axis.ticks.x=element_blank(),
# #         axis.text = element_text(size = 8),
# #         axis.line.x=element_blank(),
# #         legend.justification=c(1,0), legend.position=c(1,0.5),
# #         panel.spacing = unit(1, "lines"),
# #         strip.text = element_text(size = 14),
# #         strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")) )
# # ggsave(file=paste("figures/signaturesno_prior_", groupname, "_cellType_scatterProp_names_allsigs.pdf", sep=""), width=14,height=7)
# 
# 
# 
# #### Plotting contributions combining signatures that were split, and just focal ones (for poster)
# # combine C8 and C3 = sig9
# # C6 and C14 = sig2
# sig_prop_type$S9all = sig_prop_type[,16] + sig_prop_type[,11]
# 
# posterfocal = c("Sblood","S1","S9all","S7","S17","S2","S3","S5","S18")
# posterother = c(8,14,17,19,21,23:37)
# sumfocal = apply(sig_prop_type[,posterfocal], MARGIN=1, FUN=sum)
# sumother = apply(sig_prop_type[,posterother], MARGIN=1, FUN=sum)
# #sumfocal+sumother # all sum to 1
# sig_prop_type$other = sumother
# posterfocal2 = c(posterfocal,"other")
# sig_prop_type_poster = sig_prop_type[,c("samp_names","celltype",posterfocal2)]
# sig_prop_type_posterO = sig_prop_type_poster[order(sig_prop_type_poster$Sblood, decreasing = TRUE),]
# sig_prop_type_posterO$samp_names = factor(sig_prop_type_posterO$samp_names, levels=unique(sig_prop_type_posterO$samp_names) )
# sig_prop_type2_poster = melt(sig_prop_type_posterO) # not unattributed
# #sig_prop_type2_poster_focal = subset(sig_prop_type2_poster,celltype %in%  c("B Memory","Burkitt lymphoma","Diffse large B-cell lymphoma", "Follicular lymphoma","Lymph-CLL","T Memory","ptcl","Skin-Melanoma"))
# sig_prop_type2_poster_focal = subset(sig_prop_type2_poster,celltype %in%  c("B Memory","Burkitt lymphoma","Diffse large B-cell lymphoma", "Follicular lymphoma","MUT","UNMUT","T Memory","Skin-Melanoma"))
# sig_prop_type2_poster_focal$celltype = factor(sig_prop_type2_poster_focal$celltype, levels = c("B Memory","Burkitt lymphoma","Diffse large B-cell lymphoma", "Follicular lymphoma","MUT","UNMUT","T Memory","Skin-Melanoma"))
# sig_prop_type2_poster_focal$celltype = revalue(sig_prop_type2_poster_focal$celltype, c("Skin-Melanoma"="Melanoma","MUT"="CLL mutated","UNMUT"="CLL unmutated","Diffse large B-cell lymphoma"="DLBC lymphoma"))
# 
# 
# my10 = brewer.pal(9, "Set1") 
# my_colors = c(my10[2:9],my10[1],"white")
# ggplot(sig_prop_type2_poster_focal, aes(x=samp_names, y=value, fill=variable)) +
#   geom_bar(stat="identity") +
#   facet_wrap(.~celltype, nrow=2, scales = "free_x")+
#   xlab("") +
#   ylab("Proportion per genome") +
#   #scale_fill_manual(values=my_colors) +
#   #scale_fill_brewer(palette = "Set1") +
#   #guides(fill=FALSE)+
#   scale_fill_manual("",values=my_colors, labels=c("HSC", "S1 (aging)", "S9 (ncAID)", "S7 (UV)", "S17", "S2 (APOBEC)","S3 (ds break repair)","S5 (aging)","S18","other") ) +
#   theme_light() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x=element_blank(),
#         #axis.text = element_text(size = 10),
#         axis.line.x=element_blank(),
#         #legend.justification=c(1,0), legend.position=c(1,0.5),
#         #panel.spacing = unit(1, "lines"),
#         strip.text = element_text(size = 14),
#         strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")) )
# #ggsave(file=paste("figures/signaturesno_prior_", groupname, "_cellType_stackedBar_Prop_names_knownsigs_poster.pdf", sep=""), width=10.5,height=6)
# ggsave(file=paste("figures/signaturesno_prior_", groupname, "_cellType_stackedBar_Prop_names_knownsigs_poster_cellsubset.pdf", sep=""), width=10.5,height=6)
# 
# 
# #### Plotting the number of Sig9 mutations B cells vs BNHL, CLL, Myeloma?
# sig_prop_type$S9allN = sig_prop_type$S9all*sig_prop_type$Nmut
# sig_prop_type_posterN = sig_prop_type[,c("samp_names","celltype","S9allN")]
# sig_prop_type2_posterN_focal = subset(sig_prop_type_posterN,celltype %in%  c("B Memory","Burkitt lymphoma","Diffse large B-cell lymphoma", "Follicular lymphoma","MUT","UNMUT","Myeloid-MPN"))
# sig_prop_type2_posterN_focal$celltype = factor(sig_prop_type2_posterN_focal$celltype, levels = c("B Memory","Burkitt lymphoma","Diffse large B-cell lymphoma", "Follicular lymphoma","MUT","UNMUT","Myeloid-MPN"))
# sig_prop_type2_posterN_focal$celltype = revalue(sig_prop_type2_posterN_focal$celltype, c("Diffse large B-cell lymphoma"= "DLBC lymphoma", "MUT"="CLL mutated","UNMUT"="CLL unmutated","Myeloid-MPN"="MPN"))
# 
# #my10 = brewer.pal(9, "Set1") # max for this set is 9 
# #my_colors = c(my10[2:9],my10[1],"white")
# ggplot(sig_prop_type2_posterN_focal, aes(x=celltype, y=S9allN)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(alpha=0.3, width=0.1, size=0.6) +
#   #facet_wrap(.~celltype, nrow=2, scales = "free_x")+
#   xlab("") +
#   #ylab("# ncAID mutations / genome") +
#   ylab("# S9 mutations / genome") +
#   coord_cartesian(ylim=c(0, 6050)) + # remove 1 outlier of DLBC lymphoma at 8000
# #  theme_light() +
#   theme(text = element_text(size=14),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")) )
# #ggsave(file=paste("figures/signaturesno_prior_", groupname, "_cellType_boxplot_Nsig9_poster.pdf", sep=""), width=4,height=4)
# ggsave(file=paste("figures/signaturesno_prior_", groupname, "_cellType_boxplot_Nsig9_cellsubset_poster.pdf", sep=""), width=6,height=7)
# 
# 
# s9aov = aov(S9allN~celltype, data=sig_prop_type2_posterN_focal)
# TukeyHSD(s9aov)
# ##                  diff        lwr       upr     p adj
# # BNHL-B Memory   811.6638   217.2431 1406.0846 0.0027644
# # CLL-B Memory   -279.4100  -882.3506  323.5306 0.6280518
# # MPN-B Memory   -451.2300 -1207.2427  304.7827 0.4127530
# # CLL-BNHL      -1091.0738 -1455.2399 -726.9078 0.0000000
# # MPN-BNHL      -1262.8939 -1846.5331 -679.2546 0.0000004
# # MPN-CLL        -171.8200  -764.1342  420.4942 0.8763579
# 
# #### Plotting the proportion of Sig9 mutations B cells vs BNHL, CLL, Myeloma?
# sig_prop_type_posterP = sig_prop_type[,c("samp_names","celltype","S9all")]
# sig_prop_type2_posterP_focal = subset(sig_prop_type_posterP,celltype %in%  c("B Memory","Lymph-BNHL","Lymph-CLL","Myeloid-MPN"))
# sig_prop_type2_posterP_focal$celltype = factor(sig_prop_type2_posterP_focal$celltype, levels = c("B Memory","Lymph-BNHL","Lymph-CLL","Myeloid-MPN"))
# sig_prop_type2_posterP_focal$celltype = revalue(sig_prop_type2_posterP_focal$celltype, c("Lymph-BNHL"= "BNHL", "Lymph-CLL"="CLL","Myeloid-MPN"="MPN"))
# 
# #my10 = brewer.pal(9, "Set1") # max for this set is 9 
# #my_colors = c(my10[2:9],my10[1],"white")
# ggplot(sig_prop_type2_posterP_focal, aes(x=celltype, y=S9all)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(alpha=0.3, width=0.1, size=0.6) +
#   #facet_wrap(.~celltype, nrow=2, scales = "free_x")+
#   xlab("") +
#   ylab("Prop. AID mutations / genome") +
#   #  theme_light() +
#   theme(text = element_text(size=14),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")) )
# ggsave(file=paste("figures/signaturesno_prior_", groupname, "_cellType_boxplot_propsig9_poster.pdf", sep=""), width=4,height=4)
# s9aov = aov(S9allN~celltype, data=sig_prop_type2_posterN_focal)
# TukeyHSD(s9aov)
# 
# # scatter of # mut by prop Sig9
# ggplot(subset(sig_prop_type, celltype %in%  c("B Memory","Lymph-BNHL","Lymph-CLL","Myeloid-MPN")), aes(x=Nmut, y=S9all)) +
#   geom_point(alpha=0.5) +
#   facet_wrap(.~celltype, nrow=2, scales = "free_x")+
#   xlab("# SNVs") +
#   ylab("Prop. AID mutations / genome") +
#   #theme_light() +
#   theme(text = element_text(size=14),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")) )
# ggsave(file=paste("figures/signaturesno_prior_", groupname, "_cellType_scatter_propsig9_Nmut.pdf", sep=""), width=4,height=5)
# 
# 
# 
# # Scatterplot of number of mutations per genome
# Nmuts_sigs = sig_prop_type[,c((ncol(samp_type)+2):ncol(sig_prop_type))]*sig_prop_type$Nmut
# sig_prop_type2N = melt(data.frame(celltype=sig_prop_type[,2],Nmuts_sigs)) # not unattributed
# sig_prop_type2N$variable = factor(sig_prop_type2N$variable, levels=unique(sig_prop_type2N$variable)[c(7,4,2,8,5,1,14,12,3,10, 6,9,11,13,15:29)] )
# sig_prop_type2N$Ncelltype = factor(sig_prop_type2N$celltype, levels=c("B Memory","Lymph-BNHL","Lymph-CLL","mm","Myeloid-AML","Myeloid-MPN","Skin-Melanoma","Eso-AdenoCa","T Memory","ptcl","Lymph-NOS") )
# 
# sig_prop_type2N_known = sig_prop_type2N[sig_prop_type2N$variable %in% c("S7","S5","X3","S1","S17","Sblood","S9","S2","S3","S18"), ]
# sig_prop_type2N_known$variable = factor(sig_prop_type2N_known$variable, levels=c("Sblood","S1","S5","S9","S17","S7","S2","S3","S18","X3") )
# sig_prop_type2N_known$celltype = factor(sig_prop_type2N_known$celltype, levels=c("B Memory","Lymph-BNHL","Lymph-CLL","mm","Myeloid-AML","Myeloid-MPN","T Memory","ptcl","Lymph-NOS","Skin-Melanoma","Eso-AdenoCa") )
# 
# ggplot(sig_prop_type2N_known, aes(x=variable, y=value, fill=variable)) +
#   geom_point(pch=21, alpha=0.5) +
#   facet_wrap(.~celltype, nrow=2) +
#   geom_crossbar(data=my_median,aes(x=variable,ymin=median, ymax=median,y=median,group=celltype), width = 0.5, fatten = 1) +
#   xlab("Mutational signature") +
#   ylab("Proportion per genome") +
#   scale_fill_manual(values=my_colors) +
#   #scale_fill_brewer(palette = "Set1") +
#   guides(fill=FALSE)+
#   #scale_x_discrete(labels=c("Unattributed","Novel1 (HSC)", "S1 (aging)", "S9 (AID)", "S7 (UV)", "Novel2 (S8)", "Novel3","Novel4","S17")) +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
#         axis.ticks.x=element_blank(),
#         axis.text = element_text(size = 10),
#         axis.line.x=element_blank(),
#         legend.justification=c(1,0), legend.position=c(1,0.5),
#         panel.spacing = unit(1, "lines"),
#         strip.text = element_text(size = 14),
#         strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")) )
# ggsave(file=paste("figures/signaturesno_prior_", groupname, "_cellType_scatterNmut_names_knownsigs.pdf", sep=""), width=10,height=4)
# 
# # ylim of 10K?
# ggplot(sig_prop_type2N_known, aes(x=variable, y=value, fill=variable)) +
#   geom_point(pch=21, alpha=0.5) +
#   facet_wrap(.~celltype, nrow=2) +
#   geom_crossbar(data=my_median,aes(x=variable,ymin=median, ymax=median,y=median,group=celltype), width = 0.5, fatten = 1) +
#   xlab("Mutational signature") +
#   ylab("Proportion per genome") +
#   scale_fill_manual(values=my_colors) +
#   #scale_fill_brewer(palette = "Set1") +
#   guides(fill=FALSE)+
#   ylim(0,10000)+
#   #scale_x_discrete(labels=c("Unattributed","Novel1 (HSC)", "S1 (aging)", "S9 (AID)", "S7 (UV)", "Novel2 (S8)", "Novel3","Novel4","S17")) +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
#         axis.ticks.x=element_blank(),
#         axis.text = element_text(size = 10),
#         axis.line.x=element_blank(),
#         legend.justification=c(1,0), legend.position=c(1,0.5),
#         panel.spacing = unit(1, "lines"),
#         strip.text = element_text(size = 14),
#         strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")) )
# ggsave(file=paste("figures/signaturesno_prior_", groupname, "_cellType_scatterNmut_10Kmax_names_knownsigs.pdf", sep=""), width=10,height=4)
# 
# 
# ggplot(sig_prop_type2N, aes(x=variable, y=value, fill=variable)) +
#   geom_point(pch=21, alpha=0.5) +
#   facet_rep_wrap(.~celltype, nrow=3, repeat.tick.labels = TRUE) +
#   geom_crossbar(data=my_median,aes(x=variable,ymin=median, ymax=median,y=median,group=celltype), width = 0.5, fatten = 1) +
#   xlab("Mutational signature") +
#   ylab("Proportion per genome") +
#   #scale_fill_brewer(palette = "Set1") +
#   guides(fill=FALSE)+
#   ylim(0,10000)+
#   #scale_x_discrete(labels=c("Unattributed","Novel1 (HSC)", "S1 (aging)", "S9 (AID)", "S7 (UV)", "Novel2 (S8)", "Novel3","Novel4","S17")) +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
#         axis.ticks.x=element_blank(),
#         axis.text = element_text(size = 8),
#         axis.line.x=element_blank(),
#         legend.justification=c(1,0), legend.position=c(1,0.5),
#         panel.spacing = unit(1, "lines"),
#         strip.text = element_text(size = 14),
#         strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")) )
# ggsave(file=paste("figures/signaturesno_prior_", groupname, "_cellType_scatterNmut_10Kmax_names_allsigs.pdf", sep=""), width=14,height=7)
# 
# 
# ##########Heatmap of contributions 
# sig_prop_type_mat = data.matrix(sig_prop_type[,9:ncol(sig_prop_type)])
# rownames(sig_prop_type_mat) = paste(sig_prop_type$colony, sig_prop_type$celltype, sep=":")
# #library(gplots)
# my.rowcols = as.character(sig_prop_type$celltype)
# my.rowcols[sig_prop_type$celltype=="HSC"] = "#006837"
# my.rowcols[sig_prop_type$celltype=="B Naive"] = "#980043"
# my.rowcols[sig_prop_type$celltype=="B Memory"] = "#dd1c77"
# my.rowcols[sig_prop_type$celltype=="T Naive"] = "#253494"
# my.rowcols[sig_prop_type$celltype=="T Memory"] = "#2c7fb8"
# 
# {pdf(file=paste("figures/signatures",name_sigs,"_prior", my_prior, "_", groupname, "_heatmap.pdf", sep=""), width=5*Nsig/3,height=9)
#   heatmap.2(sig_prop_type_mat, Colv=NA, scale="none", col= colorRampPalette(brewer.pal(8, "Blues"))(25), colRow=my.rowcols, dendrogram="row", trace="none", key.title=NA, key.xlab="Proportion", xlab=paste(name_sigs, "signature", sep=" "), keysize=0.5, cexRow=0.2, key=F, lwid=c(1,0.05,1), lhei=c(0.03,1), lmat=rbind(c(5,0,4),c(3,1,2)),RowSideColors=my.rowcols )
#   legend("topleft",legend=c("HSC","B Naive","B Memory","T Naive","T Memory"), fill=c("#006837","#980043","#dd1c77","#253494","#2c7fb8"))
# }
# 
# 
# ##########  
# #### Plotting contributions combining signatures that were split, and all cancer types
# 
# sig_prop_type2_poster_all = subset(sig_prop_type2_poster, celltype != "Lymph-NOS")
# sig_prop_type2_poster_all$celltype = factor(sig_prop_type2_poster_all$celltype, levels = c("B Memory","Lymph-BNHL","Lymph-CLL","mm","Eso-AdenoCa","T Memory","ptcl","Myeloid-AML","Myeloid-MPN","Skin-Melanoma") )
# sig_prop_type2_poster_all$celltype = revalue(sig_prop_type2_poster_all$celltype, c("Lymph-BNHL"= "BNHL", "Lymph-CLL"="CLL","ptcl"="PTCL","Skin-Melanoma"="Melanoma", "mm"="MM","Myeloid-AML"="AML","Myeloid-MPN"="MPN" ) )
# 
# my10 = brewer.pal(9, "Set1") 
# my_colors = c(my10[2:9],my10[1],"white")
# ggplot(sig_prop_type2_poster_all, aes(x=samp_names, y=value, fill=variable)) +
#   geom_bar(stat="identity") +
#   facet_wrap(.~celltype, nrow=2, scales = "free_x")+
#   xlab("") +
#   ylab("Proportion per genome") +
#   #scale_fill_manual(values=my_colors) +
#   #scale_fill_brewer(palette = "Set1") +
#   #guides(fill=FALSE)+
#   scale_fill_manual("",values=my_colors, labels=c("HSC", "S1 (aging)", "S9 (AID)", "S7 (UV)", "S17", "S2 (APOBEC)","S3 (ds break repair)","S5 (aging)","S18","other") ) +
#   theme_light() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x=element_blank(),
#         #axis.text = element_text(size = 10),
#         axis.line.x=element_blank(),
#         #legend.justification=c(1,0), legend.position=c(1,0.5),
#         #panel.spacing = unit(1, "lines"),
#         strip.text = element_text(size = 14),
#         strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")) )
# ggsave(file=paste("figures/signaturesno_prior_", groupname, "_cellType_stackedBar_Prop_names_knownsigs_allcancer.pdf", sep=""), width=12,height=5)
