############################################################################################
## File: get_hg19_properties_HEM.R
## Project: lymphocyte_somatic_mutation
## Description: Extract genomic feature value per 10Kb
##
## Date: April 2021
## Author: Heather Machado
##
## Adapted from scripts by Nicola Roberts, for analyses in:
## Y. Li, N. D. Roberts, J. A. Wala, O. Shapira, S. E. Schumacher, K. Kumar, E. Khurana, S. Waszak, J. O. Korbel, J. E. Haber, M. Imielinski, J. Weischenfeldt, R. Beroukhim, P. J. Campbell, Patterns of somatic structural variation in human cancer genomes. Nature. 578, 112–121 (2020).
############################################################################################


## Fetching hg19 genomic properties per window

### Referencing lustre mounted locally: /Users/hm8/volumes/hm8_network
### Lustre dir: /lustre/scratch116/casm/cgp/users/hm8/lymphocyteWGS/mutsig_byregion/regression_analysis
#df1 = read.table("int_files/windows1kb_S9attribute_S9exp_chr14.txt", header=TRUE, stringsAsFactors = F)

#location="local"
library(GenomicFeatures)
library(stringr)
library("BSgenome.Hsapiens.UCSC.hg19")
# library(nrmisc) ## not loading
library(GenomicRanges)

args = commandArgs(trailingOnly = TRUE)
load(args[1]) # df1
outfile=args[2]

# requires a position that is a single nucleotide. Use the middle (e.g. 501st bp, assuming the 10KB is centered around it)
df2 = df1
df2$chr = paste("chr",df1$chr, sep="")  # chrom format has "chr" prefix in the genomic data files
df2$start = df1$start+5000
df2$end = df1$start+5000
gpos_testset = makeGRangesFromDataFrame(df2, keep.extra.columns=TRUE)
gpos_testset$hist = "Lymphoid"  

get_hg19_properties <- function(gpos, location='farm'){
    require(GenomicFeatures)
    require(stringr)
    require("BSgenome.Hsapiens.UCSC.hg19")
    # require(nrmisc)
  
    # all gpos elements must have width 1.
    if(any(width(gpos)!=1)){stop('All gpos elements must have width 1.')}
    # gpos must have mcol called "hist"
    if(!"hist" %in% colnames(mcols(gpos))) {stop("gpos must have column 'hist'")}

    gpos <- sort(gpos, ignore.strand=TRUE)

    # check location value, and set root accordingly
    location <- tolower(location)
    if (!location %in% c('farm', 'local')){stop("location must be one of 'farm' or 'local'")}
    if (location=='farm'){
        root <- '/lustre/scratch116/casm/cgp/users/hm8/lymphocyteWGS/mutsig_byregion/regression_analysis'
    } else {
        root <- '/Users/hm8/volumes/hm8_network/lymphocyteWGS/mutsig_byregion/regression_analysis'
    }

    # genome property files - list by property type
    in_dir <- file.path(root, 'get_hg19_properties/PanCan_genome_properties_lustre/results/GRanges/')
    f_in <- list.files(in_dir)
    f_in <- split(f_in, sapply(strsplit(f_in, "_"), `[`, 1))

    # don't include plain g_quadruplex_GR.RData, just use g4 subset with loops <=4
    f_in <- f_in[-which(names(f_in)=='g')]

    # load callable genome for ECDF calc
    load(file.path(root, 'get_hg19_properties/PanCan_callable_genome_lustre/results/callable_genome.RData'))


    for (i in 1:length(f_in)){

        if (length(f_in[[i]])==1) {
            load(file.path(in_dir, f_in[[i]]))
            pname <- gsub(".RData", "", f_in[[i]])
            prop <- get(pname)
            rm(list=pname)

            #suppressWarnings(seqinfo(prop) <- seqinfo(gpos)) ## HEM: unnecessary??
            prop <- trim(prop)

            ovl <- findOverlaps(gpos, prop)

            # only use one metric
            if (ncol(mcols(prop))==1) j <- 1 else {
              if (pname %in% c("LAD_GR", "gene_GR")) {
                j <- which(names(mcols(prop))=="dens_1e6")
              } else if (pname %in% c("cruciform_inverted_rep_GR", "short_tandem_rep_GR")){
                j <- which(names(mcols(prop))=="dens_3e3")
              } else  {
                j <- which(names(mcols(prop))=="dist_log10")
              }
            }

            cname <- gsub(" ", "_", paste(gsub("_GR", "", pname), names(mcols(prop))[j]))
            mcols(gpos)[,cname] <- mcols(prop)[subjectHits(ovl),j]

            # get quantile value - callable overlap (force 10kb tiles)
            kbs <- unlist(tile(prop, width=1e4))
            ovl <- findOverlaps(kbs, prop)
            kbs$value <- mcols(prop)[subjectHits(ovl),j]
            prop <- subsetByOverlaps(kbs, callable_genome)

            p_ecdf <- ecdf(jitter(prop$value))
            mcols(gpos)[,paste0('q_', cname)] <- p_ecdf(jitter(mcols(gpos)[,cname]))

            rm(prop, kbs, j, p_ecdf, ovl, cname, pname)

        } else {
            for(k in seq_along(f_in[[i]])){
                fname <- f_in[[i]][k]
                pname <- gsub(".RData", "", fname)
                tname <- strsplit(pname, "_")[[1]][2]

                if (!tname %in% gpos$hist) next

                load(file.path(in_dir, fname))
                prop <- get(pname)
                #suppressWarnings(seqinfo(prop) <- seqinfo(gpos)) ## HEM: unnecessary??
                prop <- trim(prop)
                rm(list=pname)

                ovl <- findOverlaps(gpos[gpos$hist==tname], prop)
                cname <- names(f_in)[i]

                mcols(gpos)[gpos$hist==tname,cname] <- prop$value[subjectHits(ovl)]

                # get quantile value - callable overlap (force 1kb tiles)
                kbs <- unlist(tile(prop, width=1e4))
                ovl <- findOverlaps(kbs, prop)
                kbs$value <- prop$value[subjectHits(ovl)]
                prop <- subsetByOverlaps(kbs, callable_genome)

                p_ecdf <- ecdf(jitter(prop$value))
                mcols(gpos)[gpos$hist==tname,paste0('q_', cname)] <- p_ecdf(jitter(mcols(gpos)[gpos$hist==tname,cname]))


                rm(prop, kbs, p_ecdf, ovl, cname, pname, tname, fname)
            }
        }

    }


    # # nucleosome occupancy - use RAW data (at 1e5 random callable pos): HEM: this analysis is SLOW! omitting for now
    # rpos <- nrmisc::sample_GRanges(callable_genome, 1e5)
    # 
    # mnase <- import.bw(file.path(root, "PanCan_genome_properties_lustre/raw/ENCFF000VNN.bigWig"),
    #                    which=c(rpos, granges(gpos)))
    # mnase <- keepSeqlevels(mnase, seqlevels(gpos))
    # # seqinfo(mnase) <- seqinfo(gpos)  ## HEM: unnecessary??
    # fo <- findOverlaps(gpos, mnase)
    # gpos$nucleosome_occupancy_raw <- as.numeric(NA)
    # gpos$nucleosome_occupancy_raw[queryHits(fo)] <- mnase$score[subjectHits(fo)]
    # 
    # # get quantile value - the rpos from callable
    # p_ecdf <- ecdf(jitter(subsetByOverlaps(mnase, rpos)$score))
    # gpos$q_nucleosome_occupancy_raw <- p_ecdf(jitter(gpos$nucleosome_occupancy_raw))

    return(gpos)

}

out1 = get_hg19_properties(gpos_testset)
write.table(out1, file=outfile, row.names = F, col.names = T, quote=F, sep="\t")
