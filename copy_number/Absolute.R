
# format variants from mutation summary as ABSOLUTE input

#---------------
# LOAD LIBRARIES
#---------------

pacman::p_load( dplyr,lazyeval,readr,tidyr,magrittr,purrr,stringr,rlist,openxlsx,  # base
                crayon,colorspace,RColorBrewer,                                    # coloring
                ggplot2,grid,gridExtra,gplots,                                     # plot layout
                parallel,
                ABSOLUTE )

#-----------
# PARAMETERS
#-----------

project.name    = 'cancer_type'
platform        = 'Illumina_WES'  # possible values: SNP_250K_STY | SNP_6.0 | Illumina_WES
cncf.dir        = 'facets/cncf'
absolute.dir    = 'absolute'


#----------
# FUNCTIONS
#----------

NewHeader <- function(cncf.file, seg.dir, absolute.dir){

    # read seg files
    cncf <-
        cncf.file %>% str_c(seg.dir, ., sep='/') %>%
        read.delim(stringsAsFactors=FALSE, sep='\t') %>%
        tbl_df

    # write copies using correct column headers
    cncf %>%
    rename(Chromosome=chrom, Start=loc.start, End=loc.end, Num_Probes=num.mark, Segment_Mean=cnlr.median) %>%
    write_tsv(str_c(absolute.dir, 'segment', str_c(cncf.file %>% str_split('.cncf.') %>% list.map(.[[1]]) %>% unlist, '.seg.txt'), sep='/'))

    return(cncf)
}


#-------
# STEP 1
#-------

# create step 1 directories
dir.create(str_c(absolute.dir, 'segment', sep='/'), recursive=TRUE, showWarnings=FALSE)
dir.create(str_c(absolute.dir, 'maf', sep='/'), recursive=TRUE, showWarnings=FALSE)
dir.create(str_c(absolute.dir, 'results', sep='/'), recursive=TRUE, showWarnings=FALSE)

# load mutations, etc
source('modules/summary/variantMaps.R')


# process cncfs
cncfs <-
    list.files('facets/cncf', pattern='*.cncf.txt') %>%
    list.map(NewHeader(., cncf.dir, absolute.dir)) %>%
    bind_rows %>%
    mutate(sample=ID) %>%
    separate(sample, into=c('sample', 'normal'), sep='_')


# process mafs
muts.maf <-
    muts %>%
    mutate(alt.depth=round(maf.t * depth.t, 0)) %>%
    mutate(ref.depth=depth.t - alt.depth) %>%
    mutate(dbSNP_Val_Status='validated') %>%
    arrange(sample, chrom, pos) %>%
    select( Tumor_Sample_Barcode   = sample,
            Hugo_Symbol            = gene,
            t_ref_count            = ref.depth,
            t_alt_count            = alt.depth,
            dbSNP_Val_Status,
            Chromosome             = chrom,
            Start_position         = pos ) %>%
    split(.$Tumor_Sample_Barcode)

map2(.x=muts.maf, .y=cncfs$ID %>% unique %>% sort, ~ {
    write_tsv(.x, str_c(absolute.dir, '/maf/', .y, '.maf.txt'))
}) %>% invisible
 

#------------
# EXAMPLE RUN
#------------

# > make -f modules/clonality/absoluteSeq.mk

# sigma.p          = 0
# max.sigma.h      = 0.07
# min.ploidy       = 0.95
# max.ploidy       = 7
# primary.disease  = "disease_name"
# sample.name      = "TUMOR_NORMAL"
# platform         = "Illumina_WES"
# max.as.seg.count = 3500
# copynum.type     = "total"
# max.neg.genome   = 0
# max.non.clonal   = 0
# min.mut.af       = 0
# seg.dat.fn       = "absolute/segment/TUMOR_NORMAL.seg.txt"
# results.dir      = "absolute/results"
# output.fn.base   = "TUMOR_NORMAL"
# maf.fn           = "absolute/maf/TUMOR_NORMAL.maf.txt"

# RunAbsolute(seg.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy, primary.disease, platform,sample.name, results.dir, max.as.seg.count, max.non.clonal, max.neg.genome, copynum.type, maf.fn = maf.fn, min.mut.af = min.mut.af, output.fn.base = output.fn.base, verbose = T)


