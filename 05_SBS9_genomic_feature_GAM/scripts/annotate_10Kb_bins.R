## Oct 2020

## annotate genomic regions

# read in data (10kb bins of signature attributions from hdp denovo extraction on memory B cells broken into 1MB bins)

# object: s9_resultsGR
load(file=paste("results/resultsGR10Kb_hdp_denovo_1Mb.Rdata", sep="") )


require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(bumphunter)
require(org.Hs.eg.db)

genes = annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene, annotation="org.Hs.eg.db")
mymatch = matchGenes(s9_resultsGR, genes)
myatt_geneinfo = cbind(s9_resultsGR, mymatch[,c("name" , "description", "region", "distance", "insideDistance", "exonnumber", "nexons", "UTR", "strand", "geneL")])
write.table(myatt_geneinfo, file=paste("results/resultsGR10Kb_hdp_denovo_1Mb_geneannotation.txt", sep=""), quote=T, col.names = T, row.names = F, sep="\t")



