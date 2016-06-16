
# format variants from mutation summary as ABSOLUTE input

#---------------
# LOAD LIBRARIES
#---------------

pacman::p_load( dplyr,lazyeval,readr,tidyr,magrittr,purrr,stringr,rlist,openxlsx,  # base
                crayon,colorspace,RColorBrewer,                                    # coloring
                ggplot2,grid,gridExtra,gplots,scales,                              # plot layout
                parallel,
                ABSOLUTE )


#-----------
# PARAMETERS
#-----------

build.maf       = TRUE
build.cncf      = FALSE

absolute.dir    = 'absolute'
platform        = 'Illumina_WES'  # possible values: SNP_250K_STY | SNP_6.0 | Illumina_WES
cncf.dir        = 'facets/cncf'


#----------
# FUNCTIONS
#----------

CncfHeader <- function(cncf.file, cncf.dir){

   # read cncf files
   cncf.file %>%
   str_c(cncf.dir, ., sep='/') %>%
   read.delim(stringsAsFactors=FALSE, sep='\t') %>%
   tbl_df
}

CncfWrite <- function(cncf, absolute.dir) {

   out.file <- str_c(absolute.dir, 'segment', str_c(cncf.file %>% str_split('.cncf.') %>% list.map(.[[1]]) %>% unlist, '.cncf.txt'), sep='/')

   # write copies using correct column headers
   cncf %>%
   rename(Chromosome=chrom, Start=loc.start, End=loc.end, Num_Probes=num.mark, cncfment_Mean=cnlr.median) %>%
   write_tsv(out.file)

}


#-------
# STEP 1
#-------

# load mutations, etc
source('modules/summary/variantMaps.R')

# create directories
MakeDirs( str_c(absolute.dir, '/cncfment'))


# process cncfs
cncfs <-
    list.files('facets', pattern='*.cncf.txt') %>%
    list.map(CncfHeader(., cncf.dir)) %>%
    list.map(CncfWrite(., absolute.dir)) %>%
    bind_rows %>%
    mutate(sample=ID) %>%
    separate(sample, into=c('sample', 'normal'), sep='_')


# process mafs
muts.maf <-
    muts.su %>%
    mutate(dbSNP_Val_Status='validated') %>%
    mutate(ref.depth=cov-val_al_count) %>%
    arrange(sample, chrom, pos) %>%
    select(Tumor_Sample_Barcode   = sample,
           Hugo_Symbol            = gene,
           t_ref_count            = ref.depth,
           t_alt_count            = most_common_al_count,
           dbSNP_Val_Status,
           Chromosome             = chrom,
           Start_position         = pos) %>%
    split(.$Tumor_Sample_Barcode)


map2(.x=muts.maf, .y=names(muts.maf), ~ {
    write_tsv(.x, str_c(absolute.dir, '/maf/', .y, '.maf.txt'))
}) %>% invisible


#------------
# EXAMPLE RUN
#------------

# > make -f modules/clonality/absoluteSeq.mk

sigma.p          = 0
max.sigma.h      = 0.07
min.ploidy       = 0.95
max.ploidy       = 7
primary.disease  = "disease_name"
sample.name      = "TUMOR_NORMAL"
platform         = "Illumina_WES"
max.as.cncf.count = 3500
copynum.type     = "total"
max.neg.genome   = 0
max.non.clonal   = 0
min.mut.af       = 0
cncf.dat.fn       = "absolute/cncfment/TUMOR_NORMAL.cncf.txt"
results.dir      = "absolute/results"
output.fn.base   = "TUMOR_NORMAL"
maf.fn           = "absolute/maf/TUMOR_NORMAL.maf.txt"

#> RunAbsolute(cncf.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy, primary.disease, platform,sample.name, results.dir, max.as.cncf.count, max.non.clonal, max.neg.genome, copynum.type, maf.fn = maf.fn, min.mut.af = min.mut.af, output.fn.base = output.fn.base, verbose = T)

#----------------------
# STEP 2 (after review)
#----------------------

# > ExtractReviewedResults("absolute/review/all.PP-calls_tab.reviewed.txt", "shb", "absolute/review/all.PP-modes.data.RData", "absolute", "all", verbose = T, copy_num_type = "total")

#------------------
# facets histograms
#------------------

library()

MakeDirs <- function(dir.list, debug=TRUE) {
    dir.list %>% list.map({
        dir <- .
        if(debug){ message(blue(str_c('- making directory: ', dir))) }
        dir.create(dir, recursive=TRUE, showWarnings=FALSE)
    }) %>% invisible
}

MakeDirs('summary/cn_hist')

PlotHist <- function(cncf.h, sample.name) {

    message(blue(str_c('- plotting: ', sample.name)))

    hp <-
        ggplot(cncf.h, aes(cn, ..density.., fill=cn.type)) +
        geom_histogram(aes(weight=weight), binwidth=0.2) +
        stat_density(bw=0.1, alpha=0.2, color='black') +
        xlab('copy number') +
        ylab('density') +
        scale_x_continuous(breaks=pretty_breaks()) +
        facet_wrap(~ cn.type, nrow=2)

    pdf(str_c('summary/cn_hist/', sample.name, '.pdf'), 5, 5)
        plot(hp)
    dev.off()
}


cncfs.h <-
    cncfs %>%
    mutate(weight=loc.end - loc.start) %>%
    select(ID, weight, tcn.em, lcn.em) %>%
    gather(cn.type, cn, tcn.em:lcn.em) %>%
    filter(!is.na(cn)) %>%
    split(.$ID)

map2(cncfs.h, names(cncfs.h), ~ { PlotHist(.x, .y) }) %>% invisible


#---------------------------------
# sufam all variants all bam files
#---------------------------------

library(parallel)

MakeDirs(c('sumamry/log', 'summary/sufam'))

muts <-
    read.delim('summary/tsv/mutation_summary.tsv', stringsAsFactors=FALSE, sep='\t') %>%
    tbl_df %>%
    FormatEvents(allosome='keep') %>%
    mutate(ID=str_c(chrom, pos, ref, alt, gene, sep=',')) %>%
    #mutate(chrom=str_c('chr', chrom)) %>%
    unique %>%
    select(`#CHROM`=chrom, POS=pos, ID, REF=ref, ALT=alt)

write_tsv(muts, 'summary/all_sufam_in.vcf')

bams <- list.files('bam', pattern='.bam$')

system('rm summary/sufam/*')
system('rm summary/log/*')

mclapply(bams, function(bam) {
    cmd <- str_c('~/share/usr/bin/sufam /ifs/e63data/reis-filho/reference/ucsc_gatk_bundle_2.8/ucsc.hg19.fasta summary/all_sufam_in.vcf bam/', bam, ' 2> summary/log/', bam, '.sufam.log > summary/sufam/', bam, '.sufam.tsv')
    message(cmd)
    system(cmd)
}, mc.cores=20) %>% invisible



muts.suf <-
   list.files('summary/sufam', full.names=TRUE) %>%
   map(~ { read.delim(.x, sep='\t', stringsAsFactors=FALSE) %>% tbl_df }) %>%
   bind_rows %>%
   filter(val_maf != 0) %>%
   mutate(pos=as.numeric(pos)) %>%
   mutate(cov=as.numeric(cov)) %>%
   mutate(val_maf=as.numeric(val_maf)) %>%
   mutate(val_al_count=as.numeric(val_al_count)) %>%
   mutate(most_common_al_count=as.numeric(most_common_al_count))




###############


muts.gene <-
    muts %>%
    select(chrom, pos, gene) %>%
    mutate(chrom=as.character(ifelse(chrom==23, 'X', chrom)))

muts.sufam <-
   read.delim('recurrent_mutations/sufam/all_sufam.txt', stringsAsFactors=FALSE, sep='\t') %>%
   tbl_df %>%
   mutate(chrom=as.character(ifelse(chrom==23, 'X', chrom))) %>%
   mutate(pos=as.numeric(pos)) %>%
   mutate(cov=as.numeric(cov)) %>%
   mutate(val_maf=as.numeric(val_maf)) %>%
   mutate(val_al_count=as.numeric(val_al_count)) %>%
   mutate(most_common_al_count=as.numeric(most_common_al_count)) %>%
   filter(!is.na(pos) & !is.nan(pos)) %>%
   left_join(muts.gene, copy=TRUE) %>%
   select(sample, chrom, pos, gene, everything()) %>%
   filter(!is.na(gene)) %>%
   filter(val_maf != 0)

muts.id <- muts %>% select(sample, normal) %>% mutate(id=str_c(sample, normal, sep='_')) %>% unique


# process mafs
muts.maf <-
    muts.sufam %>%
    left_join(muts.id, by='sample') %>%
    filter(!is.na(id)) %>%
    mutate(chrom=ifelse(chrom=='X', 23, chrom)) %>%
    unique %>%
    mutate(dbSNP_Val_Status='validated') %>%
    select(id, everything()) %>%
    mutate(ref.depth=cov-val_al_count) %>%
    arrange(sample, chrom, pos) %>%
    select( Tumor_Sample_Barcode   = id,
            Hugo_Symbol            = gene,
            t_ref_count            = ref.depth,
            t_alt_count            = most_common_al_count,
            dbSNP_Val_Status,
            Chromosome             = chrom,
            Start_position         = pos ) %>%
    split(.$Tumor_Sample_Barcode)


map2(.x=muts.maf, .y=cncfs$ID %>% unique %>% sort, ~ {
    write_tsv(.x, str_c(absolute.dir, '/maf/', .y, '.maf.txt'))
}) %>% invisible





sigma.p <- 0
max.sigma.h <- 0.07
min.ploidy <- 0.95
max.ploidy <- 7
primary.disease <- "breast"
sample.name <- "R06P_R06N"
platform <- "Illumina_WES"
max.as.cncf.count <- 3500
copynum.type <- "total"
max.neg.genome <- 0
max.non.clonal <- 0
min.mut.af <- 0
cncf.dat.fn <- "absolute/cncfment/R06P_R06N.cncf.txt"
results.dir <- "absolute/results"
output.fn.base = "R06P_R06N"
maf.fn = "absolute/maf/R06P_R06N.maf.txt"

RunAbsolute(cncf.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy, primary.disease, platform,sample.name, results.dir, max.as.cncf.count, max.non.clonal, max.neg.genome, copynum.type, maf.fn = maf.fn, min.mut.af = min.mut.af, output.fn.base = output.fn.base, verbose = T)
