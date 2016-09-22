#!/usr/bin/env Rscript

# variant processing + cascade, heatmap & tree plots

#----------
#
# LIBRARIES
#
#----------


# // 

suppressPackageStartupMessages(library('plyr'))
suppressPackageStartupMessages(library('data.table'))
suppressPackageStartupMessages(requireNamespace('phytools'))

# source('http://bioconductor.org/biocLite.R'); biocLite(c('RBGL', 'graph'))

pacman::p_load( Vennerable,                                                               # venn diagram
                ggdendro, dendextend, dendextendRcpp, dynamicTreeCut, gclus, phangorn,    # dentrograms
                gtools, nnet,                                                             # permutation
                crayon, colorspace, RColorBrewer,                                         # coloring
                ggplot2, grid, gridExtra, gplots,                                         # plotting
                dplyr, lazyeval, tidyr, magrittr, purrr, stringr, rlist,                  # dplyr
                readr, openxlsx,                                                          # I/O
                parallel, optparse )                                                      # utilities


source('~/pascal/R/variant_maps/variantFishers.R')
source('~/pascal/R/lib/fun.R')

# // -- libraries


#--------
#
# OPTIONS
#
#--------


# // 

#----------------
# program control
#----------------

# commandArgs <- function() 1:3

run.input.parameters = 1
run.data.processing  = 1
run.sub.sets         = 1
run.cn.heatmap.plots = 1
run.cascade.plots    = 1
run.trees            = 1
run.venn             = 0
run.metrics          = 0

run.fishers.plots    = 0
run.experimental     = 0


#-------------------
# subgroup selection
#-------------------

sub.sub              = 1  # run samples individually
sub.run              = 'syn_nonsyn'  # syn_nonsyn, nonsyn
sub.nums             = 0


#-------------------
# annotation options
#-------------------

call.loh             = 1
loh.closest          = 1  # should copy number / loh assignments be made according to closest segment if variant does not fall within segment: bool
call.patho           = 1
call.abs             = 1


#----------------------
# global run parameters
#----------------------

verbose              = FALSE
excel.summary.sheet  = 'MUTATION_SUMMARY'
tab.summary.sheet    = NULL #'summary/mutation_summary.tsv'  # NULL
allosome             = 'merge'
sub.set.combinations = FALSE
use.sufam            = TRUE


#-------------------
# clustering & trees
#-------------------

dist.method          = 'hamming'   # 'hamming', euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'     | see ?dist for details [hamming method implemented manually]
clust.method         = 'complete'  # 'complete', ward', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid' | see ?hclust for details
sort.method          = 'distance'  # for tree sorting: ladder 'ladderize'
tree.labels          = TRUE
min.cluster          = 3


#-----------------
# coloring options
#-----------------

random.pheno.color   = TRUE
color.seed           = 0
pheno.palette        = c('#2d4be0', '#20e6ae', '#ccb625', '#969696')
cn.cols              = 'threshold'  # 'threshold': reserved string for selecting columns with threshold sufix, 'all' = reserved string for using all columns, else a string specifing columns to use


#---------------
# Fisher's plots
#---------------

gene.names           = FALSE
targets.file         = NULL
variant.plots        = FALSE

# // -- options


#--------------
#
# I/O VARIABLES
#
#--------------


# // 

if(interactive()) {  # set detaults for interactive use

    opts = list( # input
                 project_config            = 'project_config.yaml',
                 samples_config            = 'samples.yaml',
                 subsets_config            = 'subsets.yaml',
                 muts_table_in             = 'summary/mutation_summary.xlsx',
                 gene_cn_in                = 'facets/geneCN.txt',
                 seg_maf_path              = 'absolute/reviewed/SEG_MAF',
                 cncf_path                 = 'facets/cncf',
                 # output
                 muts_table_silent_out     = 'summary/map_tables/mutations_silent.tsv',
                 muts_table_out            = 'summary/map_tables/mutations.tsv',
                 cnas_table_out            = 'summary/map_tables/copy_number_alterations.tsv',
                 variants_table_silent_out = 'summary/map_tables/variants_silent.tsv',
                 variants_table_out        = 'summary/map_tables/variants.tsv' )

} else {  # define input options

    opts.list <- list( # input
                       make_option('--project_config',             default='', help='project configuration file'),
                       make_option('--samples_config',             default='', help='sample configuration file'),
                       make_option('--subsets_config',             default='', help='subsets configuration file'),
                       make_option('--muts_table_in',              default='', help='mutation summary file'),
                       make_option('--gene_cn_in',                 default='', help='per-gene copy number file'),
                       make_option('--cncf_path',                  default='', help='location of facets cncf files'),
                       make_option('--seg_maf_path',               default='', help='location of absolute segmented mafs'),
                       # output
                       make_option('--muts_table_silent_out',      default='', help='mutation output summary'),
                       make_option('--muts_table_out',             default='', help='mutation output summary'),
                       make_option('--cnas_table_out',             default='', help='copy number aberration output summary'),
                       make_option('--variants_table_silent_out',  default='', help='variant output with silent mutations summary'),
                       make_option('--variants_table_out',         default='', help='all-variant output summary') )

    # parse input
    parser <- OptionParser(usage="%prog [options] [project_config] [samples_config] [subsets_config] [muts_table_in] [gene_cn_in] [cncf_path] [seg_maf_path] [muts_table_silent_out] [muts_table_out] [cnas_table_out] [variants_table_silent_out] [variants_table_out]", option_list=opts.list)

    # build options table
    opts <-
        parse_args(parser, positional_arguments=TRUE)$options %>%
        head(-1) %>%
        stack %>%
        mutate(ind=as.character(ind)) %>%
        rowwise %>%
        mutate(pass=(if(values == ''){  # throw error if argument is missing
                print_help(parser)
                stop(str_c('missing ', ind, ' file'))
            } else {TRUE}))

    opts <- opts$values %>% set_names(opts$ind) %>% as.list

}

# // -- I/O variables


#-----------------------
#
H1('INPUT & PARAMETERS')
#
#-----------------------


# // 

#-------------
# run settings
#-------------

debug <- interactive() & verbose==TRUE


#---------------------
# assign config params
#---------------------

# load config yaml file
config <- list.load(opts$project_config)


# load sample list
samples <-
   list.load(opts$samples_config) %>%
   list.map(., as.data.frame(., stringsAsFactors=FALSE)) %>%
   bind_rows %>%
   rowwise %>%
   mutate(normal=gsub("^\\s+|\\s+$", '', normal)) %>%
   mutate(tumor=gsub("^\\s+|\\s+$", '', tumor)) %>%
   ungroup


if('name' %in% names(samples)) {
   samples <-
        samples %>%
        rowwise %>%
        mutate(id=sub(name, '', tumor)) %>%
        mutate(name=gsub("^[[:punct:]]+|[[:punct:]]+$", '', name)) %>%
        ungroup
} else {
   samples <-
      samples %>%
      rowwise %>%
      mutate(name=strsplit(tumor, 'T') %>% .[[1]] %>% head(1) %>% unlist) %>%
      ungroup
}

samples %<>% mutate(id=name %>% substring(., 1, nchar(.)-1) %>% gsub('-', '', .))


# sample key values
keys <- config$keys %>% unlist

if(is.null(keys)) {
    Warn('no keys found in config file, using defaults')

    if(!all(is.na(samples$name))) {
        #keys <- samples$name %>% set_names(samples$tumor)
        keys <- samples$tumor %>% set_names(samples$tumor)
    } else {
        keys <- samples$tumor %>% unique %>% set_names(samples$tumor %>% unique)
    }
}


# load cancer gene list
cancer.gene.list <- read_tsv('~/share/reference/cancer+kandoth+lawrence_gene_lists.tsv') %>% unlist %>% list.filter(!is.na(.)) %>% unique %>% sort


# define sub.sets
sub.sets <- list(all=names(keys) %>% set_names(keys))
sub.sets <-
    c(sub.sets, samples %>% arrange(name) %>%
        group_by(name) %>%
        select(name, tumor) %>%
        unique %>%
        do(name=data.frame(select(.,tumor))) %>%
        lapply(FUN=function(x) { unlist(x, recursive=FALSE) }) %>%
        .[[1]] %>%
        set_names(samples$name %>% unique %>% sort) )


if(file.exists(opts$subsets_config)) {
    sub.sets <- c( sub.sets, list.load(opts$subsets_config) %>% map(~ { .x %>% KeyMod(keys, debug)}) )
}


# sub.set groups for pairwise comparisons
if(!is.null(config$subset_groups)) {
    sub.groups <-
        config$subset_groups %>%
        list.map(data_frame(.) %>% t %>% as_data_frame) %>%
        plyr::rbind.fill(.) %>%
        set_names(letters[1:(config$subset_groups %>% list.mapv(length(.)) %>% max)]) %>%
        tbl_df %>%
        gather(col,sub.set) %>%
        group_by(col) %>%
        mutate(group.id=row_number()) %>%
        ungroup %>%
        group_by(group.id) %>%
        mutate(subv=sub.sets[sub.set]) %>%
        mutate(overlap=anyDuplicated(unlist(subv))) %>%
        ungroup %>%
        select(-subv) %>%
        spread(col,sub.set) %>%
        mutate(comparison=names(config$subset_groups))

        if(any(sub.groups$overlap>0)){ message(red('warning: sub.set groups contain overlapping samples')) }

} else if(length(sub.sets) > 1 & sub.set.combinations == TRUE) {
    sub.groups <-
        combinations(n=length(sub.sets), r=2, v=names(sub.sets)) %>%
        as_data_frame %>%
        set_names(c('a', 'b')) %>%
        mutate(group.id=row_number()) %>%
        select(group.id, a, b) %>%
        rowwise %>%
        mutate(overlap=length(intersect(unlist(sub.sets[a]), unlist(sub.sets[b])))) %>%
        ungroup %>%
        mutate(comparison=str_c(a, ' x ', b))
} else if(length(sub.sets) == 1) {
    sub.groups <- data_frame(overlap=NA, group.id=1, extension=names(sub.sets[[1]]), a=names(sub.sets[[1]]), b=NA)
} else {
    sub.groups <- data_frame(overlap=NA, group.id=1, extension=NA, comparison=NA, a=NA, b=NA)
}


# stop if sub.set specifications absent
if(sub.groups %>% select(a, b) %>% unlist %in% c(names(sub.sets), NA) %>% all == FALSE & length(sub.sets) > 1) {
    print((sub.groups %>% select(a, b) %>% unlist)[!sub.groups %>% select(a, b) %>% unlist %in% names(sub.sets)] %>% unname %>% unique)
    stop('missing sub.sets specified in sub.set groups')
}


# define sub groups
sub.groups <-
    sub.groups %>%
    filter(overlap!=0 | !is.na(overlap)) %>%
    select(-overlap) %>%
    bind_rows(data_frame(extension=names(sub.sets) %>% FixName, comparison=NA, a=names(sub.sets) %>% FixName), .) %>%
    mutate(extension=ifelse(is.na(extension), str_c(comparison, sep='_') %>% FixName, extension)) %>%
    mutate(group.id=row_number()-1) %>%
    rowwise %>%
    do({ n.samples <- .[!names(.) %in% c('extension', 'comparison', 'group.id')] %>%
        unlist %>%
        sub.sets[.] %>%
        unlist %>%
        length
        c(., n.samples=n.samples) %>% as_data_frame }) %>%
    select(group.id, n.samples, extension, comparison, everything())

# add all samples as sub.set
#sub.sets <- c( sub.sets, samples$tumor %>% KeyMod(keys, debug) %>% set_names(KeyMod(., keys, debug)) %>% as.list )


if(run.data.processing != TRUE) { QuietStop('- exiting after: input & parameters //') }

# // -- input & parameters


#--------------------
#
H1('DATA PROCESSING')
#
#--------------------


# // 

#------------------------
H2('pheno palette colors')
#------------------------

# pheno color selection
if(random.pheno.color==TRUE){
    pheno.palette <-
        grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] %>%
        sample(.,length(sub.sets)) %>%
        setNames(names(sub.sets))
} else {
    pheno.palette <-
        colorRampPalette(pheno.palette)(length(sub.sets)) %>%
        setNames(names(sub.sets))  # static palette
}

# build pheno bar lookup tables
sub.sets.pheno <-
    map2(sub.sets, pheno.palette, ~{ data_frame(sample=.x, .y=.y) %>% set_names(c('sample', .y)) }) %>%
    plyr::join_all(by='sample', type='full', match='all') %>%
    plyr::rename(replace=names(pheno.palette) %>% set_names(pheno.palette)) %>%
    tbl_df %>%
    KeyMod(keys, debug)


#-------------------------
H2('mutations processing')
#-------------------------

if(!is.null(tab.summary.sheet)) {
   muts.file <- tab.summary.sheet
} else {
   muts.file <- opts$muts_table_in
}

# muts <-
#     ReadMuts(muts.file) %>%
#     FormatEvents(keys=keys)

muts <-
   ReadMuts(muts.file) %>%
   FormatEvents(keys=keys, col.names='all') %>%
   FormatEvents(keys=keys) %>%
   select(-band)

muts %<>% RmCols(c('purity'))


if(intersect.bed == TRUE) {

  bed1 <- read.delim('', sep='\t', stringsAsFactors=FALSE, skip=1, header=FALSE) %>% .[,1:3] %>% tbl_df %>% set_names(c('chrom', 'start.target', 'end.target')) %>% rowwise %>% mutate(chrom=substr(chrom, 4, nchar(chrom))) %>% ungroup %>% ChromMod(allosome='merge')
  bed2 <- read.delim('', sep='\t', stringsAsFactors=FALSE, header=FALSE) %>% .[,1:3] %>% tbl_df %>% set_names(c('chrom', 'start.target', 'end.target')) %>% ChromMod(allosome='merge')

  muts %<>%
  left_join(bed1, by='chrom') %>%
  filter(start.target <= pos & end.target >= pos) %>%
  select(-start.target, -end.target) %>%
  left_join(bed1, by='chrom') %>%
  filter(start.target <= pos & end.target >= pos) %>%
  select(-start.target, -end.target) %>%
  unique

}


#-------------------------------
H2('optinoally use sufam calls')
#-------------------------------

if(use.sufam ==TRUE) {

   H3('using sufam calls')

   # build sufam table
   muts.suf <-
      list.files('recurrent_mutations/sufam/', pattern='*_validated_sufam.txt', full.names=TRUE) %>%
      map(~
        read.delim(.x, sep='\t', stringsAsFactors=FALSE) %>%
        tbl_df %>%
        TypeCol(c('sample', 'normal', 'chrom', 'pos', 'ref',  'alt',  'cov', 'val_al_count', 'val_maf', 'X..2', 'X..1', 'A',    'C',    'G',    'T',    'X.'),
                c('char',   'char',   'num',   'num', 'char', 'char', 'num', 'num',          'num',     'char', 'char', 'char', 'char', 'char', 'char', 'int'))
        ) %>%
      bind_rows %>%
      ChromMod(allosome) %>%
      mutate(ref=substr(ref, 1, 1)) %>%
      mutate(alt=substr(val_alt, 1, 1))


   # break out normal calls
   muts.suf.n <-
      muts.suf %>%
      filter(sample %in% samples$normal) %>%
      select(normal=sample, chrom, pos, ref, alt, depth.n=cov, val.al.count.n=val_al_count, maf.n=val_maf) %>%
      unique


   # break out tumor calls
   muts.suf.t <-
      muts.suf %>%
      filter(sample %in% samples$tumor) %>%
      select(tumor=sample, chrom, pos, ref, alt, depth.t=cov, val.al.count.t=val_al_count, maf.t=val_maf) %>%
      unique %>%
      left_join(samples, by='tumor')


   # pair tumor & normal calls
   muts.suf.tn <-
      muts.suf.t %>%
      left_join(muts.suf.n, by=c('normal', 'chrom', 'pos', 'ref', 'alt')) %>%
      select(tumor, normal, chrom, pos, ref, alt, depth.n, depth.t, val.al.count.n, val.al.count.t, maf.n, maf.t) %>%
      filter(maf.t != 0 & !is.na(chrom) & !is.na(pos))


   muts.pipe <-
      muts %>%
      select(sample, chrom, pos, ref, alt) %>%
      mutate(ref=substr(ref, 1, 1)) %>%
      mutate(alt=substr(alt, 1, 1)) %>%
      mutate(pipeline=TRUE)


   # select muts columns to carry over with sufam join
   muts.all <-
      muts %>%
      mutate(ref=substr(ref, 1, 1)) %>%
      mutate(alt=substr(alt, 1, 1)) %>%
      select(chrom, pos, gene, effect, normal, ref, alt, hgvs.p, hgvs.c, variant, pathogenic, chasm, mut.taster, fathmm, provean, cancer.gene.census, kandoth, lawrence, impact, haplo, loh) %>%  # absent: sample, maf.n, depth.n, maf.t, depth.t, ex.af, mut.id,
      unique %>%
      full_join(muts.suf.tn, c('chrom', 'pos', 'normal', 'ref', 'alt'), copy=TRUE) %>%
      rename(sample=tumor) %>%
      filter(!is.na(gene)) %>%
      arrange(normal, sample, chrom, pos, ref, alt) %>%
      mutate(mut.id=str_c(sample, ':', row_number())) %>%
      full_join(muts.pipe, by=c('sample', 'chrom', 'pos', 'ref', 'alt'), copy=TRUE) %>%
      mutate(pipeline=ifelse(is.na(pipeline), FALSE, TRUE))

} else { H3('skipping sufam calls') }


#--------------
H2('LOH calls')
#--------------

if(call.loh == TRUE) {

    H3('calling loh')
    # read cnf & make calls
    segs.loh <-
        list.files(opts$cncf_path, pattern='.cncf.txt', full.names=TRUE) %>%
        map(~ { CallLOH(.x) }) %>%
        bind_rows %>%
        mutate(sample=keys[sample])

    # assign loh calls to each mutation
    muts.loh <-
        muts %>%
        select(sample, chrom, pos) %>%
        group_by(sample, chrom, pos) %>%
        full_join(segs.loh, by=c('sample', 'chrom')) %>%
        slice(which.min(pos-loc.mid)) %>%
        ungroup %>%
        mutate(in.seg=(loc.start<=pos & pos<=loc.end))

    # filter loh calls if specified
    if(loh.closest != TRUE) {
        muts.loh %<>% filter(in.seg==TRUE)
    }

    if('loh' %in% colnames(muts)) {
        Warn('removing existing loh column')
        muts %<>% select(-loh)
    } else {
        muts.loh %<>% select(-loh)
    }

    # add loh & seg columns to mutation tables
    muts %<>% left_join(muts.loh, by=c('sample', 'chrom', 'pos'))

} else {
    H3('skipping loh calling')
    if('loh' %in% names(muts)) {
        Warn('warning: no prexisting loh column')
    }
}


#------------------------
H2('pathogenicity calls')
#------------------------

if(call.patho == TRUE) {

    muts <-
        muts %>%
            rowwise %>%
            mutate(k.l.c=ifelse(any(strsplit(gene,'\\|') %in% cancer.gene.list),TRUE,FALSE)) %>%
            ungroup %>%
            mutate(pathogenic =
                            ifelse(effect == 'Missense SNV' & (mut.taster %in% c('D', 'A') | chasm <= 0.3) & (fathmm == 'CANCER' | chasm <= 0.3) & k.l.c == TRUE, 'Pathogenic',
                            ifelse(effect == 'Missense SNV' & (mut.taster %in% c('D', 'A') | chasm <= 0.3) & (fathmm == 'CANCER' | chasm <= 0.3), 'Potentially Pathogenic',

                            ifelse(effect == 'Inframe In-Del' & mut.taster %in% c('D', 'A') & provean == TRUE & (haplo == TRUE | loh == 'LOH' | k.l.c == TRUE), 'Pathogenic',
                            ifelse(effect == 'Inframe In-Del' & mut.taster %in% c('D', 'A') & provean == TRUE, 'Potentially Pathogenic',

                            ifelse(effect %in% c('Frameshift In-Del', 'Splice site variant', 'Truncating SNV', 'Upstream, start/stop, or de novo modification') & (haplo == TRUE | loh == 'LOH') & k.l.c == TRUE, 'Pathogenic',
                            ifelse(effect %in% c('Frameshift In-Del', 'Splice site variant', 'Truncating SNV', 'Upstream, start/stop, or de novo modification') & (haplo == TRUE | loh == 'LOH' | k.l.c == TRUE), 'Potentially Pathogenic',

                            'Passenger')))))))
}


#-----------------------
H2('ABSOLUTE CCF calls')
#-----------------------

if(call.abs == TRUE) {

    # retreive absolute ccf calls
    abs.ccf <-
        list.files(opts$seg_maf_path, pattern='_ABS_MAF.txt', full.names=TRUE) %>%
        map(~ { read.delim(.x, sep='\t',stringsAsFactors=FALSE) %>% tbl_df }) %>%
        bind_rows %>%
        separate(sample, into=c('sample','normal'), sep='_', fill='right') %>%
        rename(pos=Start_position, alt.depth=alt) %>%
        FormatEvents %>%
        tbl_df %>%
        mutate(clonality=ifelse(pr.somatic.clonal>=0.5 | ci95.low >=0.9, 'Clonal', 'Subclonal')) %>%
        select(sample, gene, chrom, pos, ccf, clonality, purity, pr.somatic.clonal, ci95.low) %>%
        mutate(sample=KeyMod(sample, keys=keys, debug=debug))

    if('ccf' %in% names(muts)) {
        muts %<>% select(-ccf)
    }
    if('pr.somatic.clonal' %in% names(muts)) {
        muts %<>% select(-pr.somatic.clonal)
    }
    if('clonality' %in% names(muts)) {
        muts %<>% select(-clonality)
    }
    if('ci95.low' %in% names(muts)) {
        muts %<>% select(-ci95.low)
    }

    # add absolute calls to mutation table
    muts <-
        muts %>%
        left_join(abs.ccf, by=c('sample', 'chrom', 'pos', 'gene')) %>%
        DummyCols('cf', debug) %>%
        mutate(ccf.replace=ifelse(is.na(ccf), TRUE, FALSE)) %>%
        mutate(ccf=ifelse(ccf.replace == TRUE, cf, ccf)) %>%
        unique

}


#-------------------------
H2('per-gene copy number')
#-------------------------

# select threshold columns
gene.cn <- read.delim(opts$gene_cn_in, sep='\t',stringsAsFactors=FALSE) %>% tbl_df %>% arrange(chrom, start)
gene.cn <- set_names(gene.cn, gsub('\\.', '-', names(gene.cn)))

if(cn.cols=='threshold') {
    gene.cn %<>% select(gene=hgnc,chrom,start,end,band,matches('threshold'))
    cn.sample.keys <- grep('threshold',colnames(gene.cn), value=TRUE) %>% list.map(strsplit(.,'_') %>% unlist %>% head(1) %>% KeyMod(keys, debug) %>% unname)
    names(gene.cn)[which(names(gene.cn) %in% names(cn.sample.keys))] <- cn.sample.keys[which(names(cn.sample.keys) %in% names(gene.cn))]
} else if(cn.cols!='all') {
    gene.cn %<>% select(gene=hgnc,chrom,start,end,band,one_of(cn.cols))
    cn.sample.keys <- colnames(gene.cn) %>% list.filter(!. %in% c('gene','hgnc','chrom','start','end','band'))
} else {
    cn.sample.keys <- colnames(gene.cn) %>% list.filter(!. %in% c('gene','hgnc','chrom','start','end','band'))
}


# geneCN column selection vector
gene.cn.cols <- sub.sets %>% unlist %>% unique %>% unname %>% KeyMod(keys, debug)


# sub.set amp / del rows
cnas <- 
    gene.cn %>%
    mutate(amp=rowSums(.[gene.cn.cols]==2)) %>%
    mutate(del=rowSums(.[gene.cn.cols]==-2)) %>%
    filter(amp==TRUE|del==TRUE) %>%
    mutate_each(funs(ifelse(.==1,0,.)), one_of(gene.cn.cols)) %>%
    mutate_each(funs(ifelse(.==-1,0,.)), one_of(gene.cn.cols)) %>%
    unique %>%
    arrange(chrom, start) %>%
    group_by(chrom) %>%
    mutate(max.end=max(end)) %>%
    mutate(lag.end=lag(end)) %>%
    (function(df=., col.v=c("chrom", "band", unname(gene.cn.cols))) {
        DistinctCalls <- sapply(col.v, function(col) {
            lazyeval::interp(~ col.name, col.name=as.name(col))
        })
        df %>% distinct_(.dots=DistinctCalls, .keep_all=TRUE)
    }) %>%
    mutate(lead.lag.end=lead(lag.end)) %>%
    mutate(end=ifelse(!is.na(lead.lag.end), lead.lag.end, max.end)) %>%
    select(chrom, start, end, band, one_of(gene.cn.cols)) %>%
    gather(sample, effect, -band, -chrom, -start, -end) %>%
    filter(effect==2|effect==-2)

if(nrow(cnas) != 0) {
    cnas <-
        cnas %>%
        rowwise %>%
        mutate(sample=strsplit(sample,'_') %>% unlist %>% head(1)) %>%
        ungroup %>%
        unique %>%
        FormatEvents %>%
        select(sample, band, chrom, start, end, effect)

    # collapse consecutive bands
    cnas %<>%
        arrange(sample, chrom, start) %>%
        separate(band, into=c('arm', 'band'), sep='\\.', fill='right') %>%
        mutate(band=as.numeric(band)) %>%
        group_by(sample, chrom, arm, effect) %>%
        mutate(lead.band=lead(band)) %>%
        mutate(consec=band==lead.band-1 | band==lead.band+1) %>%
        mutate(no.consec=ifelse(is.na(consec),TRUE,FALSE)) %>%
        mutate(no.consec=ifelse(row_number()==n(),TRUE,no.consec)) %>%
        mutate(consec=ifelse(row_number()==n(),TRUE,consec)) %>%
        ungroup %>%
        mutate(cu.consec=cumsum(.$no.consec)) %>%
        mutate(cu.consec=ifelse(lag(no.consec)==FALSE & no.consec==TRUE, lag(cu.consec), cu.consec)) %>%
        ungroup %>%
        group_by(sample, chrom, arm, effect, cu.consec) %>%
        slice(c(1,n())) %>%
        unique %>%
        mutate(lead.jump.band=lead(band)) %>%
        mutate(lead.jump.end=lead(end)) %>%
        mutate(band.test=ifelse(!is.na(lead.jump.band), lead.jump.band, band)) %>%
        mutate(end=ifelse(!is.na(lead.jump.end), lead.jump.end, end)) %>%
        mutate(band=as.character(band)) %>%
        mutate(band.rep=str_c(lead.jump.band, '-', band)) %>%
        mutate(band=ifelse(!is.na(band.rep),band.rep,band)) %>%
        mutate(n.group=n()) %>%
        filter(!is.na(lead.jump.end) | is.na(band) | n.group==1) %>%
        top_n(1) %>%
        ungroup %>%
        mutate(arm=as.character(arm)) %>%
        mutate(band=str_c('.',band)) %>%
        mutate(band=ifelse(is.na(band),'',band)) %>%
        mutate(band=str_c(arm,band)) %>%
        select(sample, band, chrom, start, end, effect)
}


#------------------------
H2('combine muts & cnas')
#------------------------

message(green('  combining variants'))

# all mutations + cnas, add pheno columns
variants <-
    bind_rows( muts %>%  # abbreviated mutations
               mutate(class='muts') %>%
               DummyCols(c('clonality', 'loh', 'ccf', 'cancer.gene.census'), debug) %>%
               select(class, sample, chrom, pos, span=gene, effect, ccf, cancer.gene.census, loh, haplo, clonality, everything()),
               cnas %>%  # abbreviated cnas
               mutate(class='cnas') %>%
               select(class, sample, chrom, pos=start, span=band, effect, everything()) ) %>%
    left_join(sub.sets.pheno, by='sample') %>%
    arrange(sample, effect) %>%
    rowwise %>%
    mutate(k.l.c=ifelse(any(strsplit(span,'\\|') %in% cancer.gene.list),TRUE,FALSE)) %>%
    ungroup


#-------------------------
H2('write variants files')
#-------------------------

MakeDirs('summary/map_tables')

variants %>%
filter(effect!='Silent') %>%
write_tsv(opts$variants_table_silent_out)

variants %>%
write_tsv(opts$variants_table_out)


#---------------------
H2('write muts files')
#---------------------

muts %>%
DummyCols('ccf', debug) %>%
arrange(sample, desc(ccf)) %>%
filter(effect!='Silent') %>%
write_tsv(opts$muts_table_silent_out)

muts %>%
DummyCols('ccf', debug) %>%
arrange(sample, desc(ccf)) %>%
write_tsv(opts$muts_table_out)


#--------------------
H2('write CNAS file')
#--------------------

cnas %>%
spread(sample, effect) %>%
arrange(chrom, start) %>%
write_tsv(opts$cnas_table_out)


if(run.sub.sets != TRUE) { QuietStop('- exiting after: data processing //') }

# // -- data processing


#-------------
#
H1('SUB SETS')
#
#-------------

if(sub.run == 'nonsyn') {

   H2('NONSYNONYMOUS')
   run.variants <-
      variants %>%
      filter(effect != 'Silent') %>%
      filter(ccf > 0)

 } else {

   H2('NONSYNONYMOUS + SYNONYMOUS')
   run.variants <-
      variants %>% filter(ccf > 0)

 }


if(sub.sub != TRUE) {

   for(sub.run in c('nonsyn', 'syn_nonsyn')) {

      if(sub.run == 'nonsyn') {

         H2('NONSYNONYMOUS')
         run.variants <-
            variants %>%
            filter(effect != 'Silent') %>%
            filter(ccf > 0)

       } else {

         H2('NONSYNONYMOUS + SYNONYMOUS')
         run.variants <-
            variants %>% filter(ccf > 0)

       }

   mclapply(0:(nrow(sub.groups)-1), function(sub.num) { SubGroup(sub.num) }, mc.silent=FALSE, mc.cores=20)

   }

} else if(length(sub.nums) != 1) {

   lapply(sub.nums, function(sub.num) { SubGroup(sub.num) })
} else {
   SubGroup(sub.nums)
}
