

#---------------------------
# add dummy column if absent
#---------------------------

DummyCols <- function(events, col.names, debug) {

    for(col.name in col.names) {
        if(!col.name %in% colnames(events)) {
            if(debug) { message(green(str_c('adding dummy column: ', col.name))) }
            events.names <- colnames(events) %>% replace(which(is.na(.)),'empty') %>% list.filter(.!='add.col')
            events$add.col <- NA
            events %<>% setNames(c(events.names,col.name))
        }
    }
    return(events)

}


#-----------------------------
# remove columns if they exist
#-----------------------------

RmCols <- function(df, cols) {

  rm.cols <- which(names(df) %in% cols)

  if(length(rm.cols > 0)) {
    df[, -rm.cols]
  } else {
    df
  }
}


#--------------------------
# convert chromosome format
#--------------------------

ChromMod <- function(events, allosome){

    if('X' %in% events$chrom & allosome != 'keep') { Warn('X chromosome labelling converted to numeric') }

    if(allosome == 'distinct') {
        if('Y' %in% events$chrom) { Warn('Y chromosome labelling converted to 24-based numeric') }
        events %>%
        mutate(.,chrom=as.numeric(ifelse(chrom=='X', 23,chrom))) %>%
        mutate(chrom=as.numeric(ifelse(chrom=='Y', 24,chrom)))
    } else if(allosome == 'merge') {
        if('Y' %in% events$chrom) { Warn('Y chromosome labelling converted to 23-based numeric') }
        events %>%
        mutate(chrom=as.numeric(ifelse(chrom %in% c('X','Y'), 23, chrom)))
    } else if(allosome == 'keep') {
        events
    } else {
        events %>%
        filter(!chrom %in% c('X','Y'))
    }
}


#-----------------------------------
# soft sample key:value modification
#-----------------------------------

KeyMod <- function(sample.object, keys, debug=FALSE, force.keys=FALSE) {

    if(is.vector(sample.object)) {  # convert strings
        keys.index <- which(sample.object %in% names(keys))
        if(length(keys.index) == 0 & force.keys==FALSE) {

            if(debug) {
                message(green('no keys to convert'))
            }

        } else if(force.keys==FALSE) {

            conversions <- data.frame(pre=as.character(sample.object[keys.index]), post=keys[as.character(sample.object[keys.index])], row.names=NULL) %>% unique
            if(debug) {
                message(green(str_c('key conversions: ',length(conversions$pre),'/',length(unique(sample.object)))))
                print(conversions)
            }
            sample.object <- keys[sample.object]

        } else {

            conversions <- data.frame(pre=as.character(unique(sample.object)), post=keys[unique(as.character(sample.object))], in.keys=unique(as.character(events)) %in% names(keys), row.names=NULL)
            if(debug) {
                Warn('forcing key conversions')
                print(conversions)
            }
            sample.object[keys.index] <- keys[sample.object[keys.index]]
        }
    } else {  # convert data.frame

        keys.index <- which(sample.object$sample %in% names(keys))

        if(length(keys.index) == 0 & force.keys==FALSE) {
            if(debug) { message(green('no keys to convert')) }
        } else if(force.keys==FALSE) {
            conversions <- data.frame(pre=as.character(sample.object$sample[keys.index]), post=keys[as.character(sample.object$sample[keys.index])], row.names=NULL) %>% unique
            if(debug) {
                message(green(str_c('key conversions: ',length(conversions$pre),'/',length(unique(sample.object$sample)))))
                print(conversions)
            }
            sample.object %<>% mutate(sample=ifelse(row_number() %in% keys.index, keys[sample], sample))
        } else {
            conversions <- data.frame(pre=as.character(unique(sample.object$sample)), post=keys[unique(as.character(sample.object$sample))], in.keys=unique(as.character(events$sample)) %in% names(keys), row.names=NULL)
            if(debug) {
                Warn('forcing key conversions')
                print(conversions)
            }
            sample.object %<>% mutate(sample=keys[sample])
        }
    }
    sample.object
}


#-----------------------
# coerce columns to type
#-----------------------

TypeCol <- function(events, columns, types) {
    for (col in 1:length(columns)) {
        if(columns[col] %in% colnames(events)) {
            if(types[col] == 'char') {
                events[,columns[col]] <- events[,columns[col]] %>% unlist %>% as.character
            } else if(types[col] == 'num') {
                events[,columns[col]] <- events[,columns[col]] %>% unlist %>%  ifelse(.=='.',NA,.) %>% as.numeric
            } else if(types[col] == 'logic') {
                events[,columns[col]] <- events[,columns[col]] == 'TRUE' %>% as.vector
            }
        }
    }
    return(events)
}


#---------------------------------
# multiple pheno columns to single
#---------------------------------

MergePheno <- function(events, merge.cols, col.name='pheno') {
    events %>%
    do({ df <- .
        pheno <- select(df, one_of(merge.cols)) %>%
        apply(1, function(pheno){
            pheno %<>% na.omit
            if(length(pheno)==0){
                pheno=NA
            } else {
                pheno
            }
        })
        df <- mutate(df, col.name=pheno)
        colnames(df)[colnames(df) == 'col.name'] <- col.name
        return(df)
    })
}


#-----------------------------
# read mutations file & format
#-----------------------------

ReadMuts <- function(muts.file, keys=NULL, allosome='merge') {

    # read mutations file
    if(length(grep('.xlsx$', muts.file))){
        muts <- read.xlsx(muts.file)
    } else {
        muts <- read.delim(muts.file, sep='\t', stringsAsFactors=FALSE)
    }

    if(anyDuplicated(colnames(muts))) {
        message(yellow(str_c('duplicated columns:', colnames(foo)[duplicated(names(foo))], sep=' ')))
    } else {
        muts %<>% tbl_df
    }

    return(muts)
}


#------------------------------------
# rename columns & clean nomenclature
#------------------------------------

FormatEvents <- function(events, col.names=NULL, drop=FALSE, allosome='merge', keys=NULL, force.keys=FALSE) {

    # col.names:
    #
    #   NULL = use only the columns available in input table
    #   character string = dumy column names to be added if absent
    #   'all' = reserved character string to specify addition of all maintained columns
    #
    #   drop = TRUE / FALSE, should additional columns be dropped
    #
    # allosome:
    #
    #   'none'      = exclude X & Y chromosome
    #   'merge'     = treat X & Y coordinates as homologous pair
    #   'distinct'  = treat X & Y as seperate, sequential chromosomes

    # lookup table
    col.keys <- c( 'alt'='alt', 'ALT'='alt', 'Alternate.allele'='alt',
                   'alt.count.t'='alt.count.t', 't_alt_count'='alt.count.t', 
                   'band'='band', 'Band'='band',
                   'cancer.gene.census'='cancer.gene.census', 'Cancer.Gene.Census'='cancer.gene.census', 'Cancer Gene Census'='cancer.gene.census','Cancer5000.S.genes..Lawrence.et.al.'='cancer.gene.census', 'cancer_gene_census'='cancer.gene.census', 'cancer.gene.census'='cancer.gene.census',
                   'cancer.gene'='cancer.gene', 'Cancer.Gene'='cancer.gene',
                   'ccf'='ccf', 'Cancer_Cell_Fraction'='ccf', 'CCF'='ccf', 'cancer_cell_frac'='ccf', 'Cancer cell fraction'='ccf', 'Cancer.cell.fraction'='ccf',
                   'chasm'='chasm', 'CHASM'='chasm', 'Breast_chasm_score'='chasm','Chasm'='chasm',
                   'chrom'='chrom','Chrom'='chrom','Chromosome'='chrom','CHROM'='chrom',
                   'Clonal'='clonal', 'clonal'='clonal', 'Clonal_Status'='clonal', 'Clonality'='clonal',
                   'cn'='cn', 'CN'='cn',
                   'ci95.low'='ci95.low', 'ccf_CI95_low'='ci95.low', 'CI95_low'='ci95.low',
                   'NORMAL.DP'='depth.n', 'depth.n'='depth.n',
                   'TUMOR.DP'='depth.t', 'depth.t'='depth.t',
                   'Effect'='variant', 'Variant_Classification'='variant', 'ANN....EFFECT'='variant', 'ANN[*].EFFECT'='variant',
                   'effect'='effect',
                   'ExAC_AF'='ex.af', 'ex.af'='ex.af', 'dbNSFP_ExAC_Adj_AF'='ex.af',
                   'end'='end', 'stop'='end', 'End_position'='end',
                   'fathmm'='fathmm', 'FATHMM'='fathmm', 'fathmm_pred'='fathmm',
                   'gene'='gene', 'Gene'='gene', 'Hugo_Symbol'='gene', 'GENE'='gene', 'hgnc'='gene', 'ANN....GENE'='gene', 'ANN[*].GENE'='gene', 'Gene.symbol'='gene', 'Hugo_Symbol'=='gene',
                   'haploinsufficient'='haplo', 'hap_insuf'='haplo', 'haplo'='haplo', 'Haploinsufficient'='haplo',
                   'hotspot'='hotspot', 'HOTSPOT'='hotspot','Hotspot'='hotspot',
                   'ANN[*].HGVS_C'='hgvs.c', 'ANN....HGVS_P'='hgvs.c',
                   'ANN[*].HGVS_P'='hgvs.p', 'ANN....HGVS_C'='hgvs.p',
                   'ANN[*].IMPACT'='impact', 'impact'='impact', 'ANN....IMPACT'='impact',
                   'kandoth'='kandoth', '127 significantly mutated genes (Kandoth et al)'='kandoth', '127 significantly.mutated genes.(Kandoth et al)'='kandoth','X127.significantly.mutated.genes..Kandoth.et.al.'='kandoth', 'Kandoth'='kandoth',
                   'lawrence'='lawrence', 'Cancer5000-S genes (Lawrence et al)'='lawrence', 'Cancer5000-S genes.(Lawrence et al)'='lawrence',
                   'loh'='loh', 'LOH'='loh', 'Loss.of.heterozygocity.(LOH)'='loh', 'Loss of heterozygocity (LOH)'='loh','Loss.of.heterozygocity..LOH.'='loh',
                   'NORMAL_SAMPLE'='normal', 'normal'='normal',
                   'NORMAL_MAF'='maf.n', 'maf.n'='maf.n',
                   'maf.t'='maf.t', 'Mutant allele fraction', 'Mutant.allele.fraction'='maf.t', 'TUMOR_MAF'='maf.t', 'maf'='maf.t',
                   'maf.n'='maf.n',
                   'mut.taster'='mut.taster', 'dbNSFP_MutationTaster_pred'='mut.taster', 'Mutation.Taster'='mut.taster',
                   'pathogenic'='pathogenic', 'Pathogenic'='pathogenic', 'pathogenicity'='pathogenic', 'Pathogenicity'='pathogenic',
                   'pheno'='pheno', 'pheno.bar'='pheno',
                   'pos'='pos','POS'='pos', 'Position'='pos', 'position'='position',
                   'pr.somatic.clonal'='pr.somatic.clonal', 'Pr_somatic_clonal'='pr.somatic.clonal',
                   'provean'='provean', 'Provean'='provean', 'dbNSFP_PROVEAN_pred'='provean',
                   'purity'='purity',
                   'ref'='ref', 'REF'='ref', 'Reference.allele'='ref',
                   'ref.count.t'='ref.count.t', 't_ref_count'='ref.count.t', 
                   'sample'='sample', 'Sample'='sample', 'Tumor_Sample_Barcode'='sample', 'TUMOR_SAMPLE'='sample', 'Sample.ID'='sample',
                   'start'='start', 'Start'='start', 'Start_position'='start' )

    # save column names to fill unhandled
    fallback.cols <- colnames(events)

    # rename columns
    names(events) <- col.keys[names(events)]

    # specify output columns
    if(is.null(col.names)) {
        col.names <- colnames(events)
        cols.match <- col.names %>% list.filter(!is.na(.))
        col.names[which(is.na(col.names))] <- fallback.cols[which(is.na(col.names))]
        names(events) <- col.names
    } else if('all' %in% col.names) {
        col.names <- c( 'alt','band','cancer.gene.census','ccf','chasm','chrom','clonal','cn','ci95.low','depth.t', 'depth.n', 'variant',
                        'end','fathmm','gene','haplo','hgvs.p','hotspot', 'kandoth','lawrence','loh','maf.n', 'maf.t','mut.taster',
                        'pathogenic','pheno','pos','pr.somatic.clonal','provean','purity','ref','sample','start' )
    }

    # df typing to avoid row name conflicts
    events %<>% as.data.frame

    # add columns if absent
    events %<>% DummyCols(col.names, debug)

    # select columns
    events <- events[colnames(events) %>% list.filter(.!='empty')]

    # chromosome type conversion
    events %<>% ChromMod(allosome)

    # column typing
    events %<>% TypeCol( c('sample', 'gene', 'effect', 'variant', 'ref',  'alt',  'cancer.gene.census', 'depth.t', 'depth.n', 'hgvs.p', 'kandoth', 'lawrence', 'maf.n', 'maf.t', 'mut.taster', 'haplo', 'hotspot', 'fathmm', 'chasm', 'provean'),
                         c('char',   'char', 'char',   'char',    'char', 'char', 'logic',              'num',     'num',     'logic',  'logic',   'logic',    'num',   'num',   'char',       'logic', 'logic',   'char',   'num',   'logic') )

    if('variant' %in% colnames(events) & !'effect' %in% colnames(events)) {
        events %<>% mutate(effect=variant) %>% tbl_df
    }

    if('effect' %in% colnames(events)) {

        # clip effect names to first if pipe or '&' seperated
        events %<>% rowwise %>% mutate(effect=effect %>% str_split('\\|') %>% .[[1]] %>% head(1)) %>% ungroup
        events %<>% rowwise %>% mutate(effect=effect %>% str_split('&') %>% .[[1]] %>% head(1)) %>% ungroup

        # remove NA events
        if(!all(is.na(events$effect))) {
            events %<>% filter(!is.na(effect))
        }

        # rename variant classifications
        events <-
            events %>%
            unique %>%
            mutate(effect=
                ifelse(effect%in%c('STOP_GAINED','Nonsense_Mutation','stop_gained&splice_region_variant','stop_gained','Nonsense_Mutation','Stop_Codon_Ins','nonsense','truncating snv','Truncating snv','Truncating snv','Truncating SNV'),'Truncating SNV',
                ifelse(effect%in%c('FRAME_SHIFT','FRAME_SHIFT','Frame_Shift_Del','Frame_Shift_Ins','frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant','frameshift_variant','frameshift_variant&stop_gained','frameshift_variant&splice_region_variant','frameshift_variant&splice_acceptor_variant&splice_region_variant&splice_region_variant&intron_variant','Frame_Shift_Del','Frame_Shift_Ins','frame_shift_del','frame_shift_ins','frameshift indel','Frameshift indel','Frameshift In-Del','frameshift_variant'),'Frameshift In-Del',
                ifelse(effect%in%c('NON_SYNONYMOUS_CODING','STOP_LOST','Missense_Mutation','missense_variant','missense_variant&splice_region_variant','missense_variant|missense_variant','Missense_Mutation','missense','missense snv','Missense snv','Missense SNV'),'Missense SNV',
                ifelse(effect%in%c('CODON_CHANGE_PLUS_CODON_DELETION','CODON_DELETION','CODON_INSERTION','In_Frame_Ins','In_Frame_Del','disruptive_inframe_deletion','disruptive_inframe_insertion','inframe_deletion','inframe_insertion','disruptive_inframe_deletion&splice_region_variant','inframe_deletion&splice_region_variant','In_Frame_Del','In_Frame_Ins','in_frame_del','in_frame_ins','inframe indel','Inframe indel','Inframe In-Del'),'Inframe In-Del',
                ifelse(effect%in%c('splice_donor_variant','splice_region_variant','splice_acceptor_variant','SPLICE_SITE_DONOR','SPLICE_SITE_ACCEPTOR','SPLICE_SITE_REGION','Splice_Site','splice_donor_variant&intron_variant','splice_acceptor_variant&intron_variant','splicing','splice_donor_variant&splice_region_variant&intron_variant','splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&splice_region_variant&intron_variant','Splice_Site','splice','splice site variant','Splice site variant','missense_variant & splice_region_variant'),'Splice site variant',
                ifelse(effect%in%c('STOP_LOST','START_LOST','START_GAINED','UTR_5_PRIME','start_lost','stop_lost',"5'UTR","5'Flank",'De_novo_Start_InFrame','De_novo_Start_OutOfFrame','Stop_Codon_Del','Start_Codon_SNP','Start_Codon_Ins','Start_Codon_Del','Nonstop_Mutation','nonstop','upstream, start/stop, or de novo modification','Upstream, start/stop, or de novo modification'),'Upstream, start/stop, or de novo modification',
                ifelse(effect%in%c('synonymous_variant','splice_region_variant&synonymous_variant','splice_region_variant&synonymous_variant','non_coding_exon_variant','upstream_gene_variant','downstream_gene_variant','intron_variant','frameshift_variant&splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant','non_coding_exon_variant|synonymous_variant','SYNONYMOUS_CODING','synonymous_variant|synonymous_variant','splice_region_variant&synonymous_variant|splice_region_variant&non_coding_exon_variant','splice_acceptor_variant & intron_variant','intragenic_variant',"3'UTR",'IGR','lincRNA','RNA','Intron','silent','intron_exon','silent','Silent','intron_variant & missense_variant'),'Silent',
                ifelse(effect%in%c('Amplification','amplification','amp','2'),'Amplification',
                ifelse(effect%in%c('Deletion','deletion','del','-2'),'Deletion',
                ifelse(is.na(effect),NA,str_c('UNACCOUNTED: ', effect))))))))))))

        # warn on unknown effects
        if(all(is.na(events$effect))){
            message(yellow('warning: empty effect column, left as NA'))
        } else if(events %>% filter(is.na(effect)) %>% nrow > 0) {
            message(yellow('warning: some variants not accounted for'))
        }
    }

    # change chasm NA column to Inf
    if('chasm' %in% colnames(events)) {
        events %<>% mutate(chasm=ifelse(is.na(chasm), Inf, chasm))
    }

    # remove unmutated LOH if present & format
    if('loh' %in% colnames(events)) {
        events %<>% mutate(loh=ifelse(loh%in%c('LOH','loh','Loss of heterozygosity') & !is.na(effect),'Loss of heterozygosity','.'))
    }

    if(drop!=FALSE) {
        events %<>% select(one_of(cols.match))
    }

    if(!is.null(keys)) {
        events %<>% KeyMod(keys, debug)
    }

    events %<>% tbl_df

    return(events)
}



















#---------------------
# header message setup
#---------------------

# H1 message
H1 <- function(txt) {
    div <- str_c('  ', str_c(rep('-', nchar(txt)), collapse=''))
    message(cyan(str_c('\n', div, '\n\n  ', txt, '\n\n', div, '\n')))
}

# H2 message
H2 <- function(txt) {
    div <- str_c('  ', str_c(rep('-', nchar(txt)), collapse=''))
    message(green(str_c(div, '\n  ', txt, '\n', div)))
}

# H3 message
H3 <- function(txt) {
    message(blue(str_c('- ', txt)))
}

# warning
Warn <- function(txt) {
    message(yellow(str_c('- ', txt)))
}


#--------------
# graceful exit
#--------------

QuietStop <- function(txt) {
    message(blue(txt))
    stop(simpleError('\r      \r'))
}


#--------------------------
# create system directories
#--------------------------

MakeDirs <- function(dir.list, debug=TRUE) {
    dir.list %>% list.map({
        dir <- .
        if(debug){ message(blue(str_c('- making directory: ', dir))) }
        dir.create(dir, recursive=TRUE, showWarnings=FALSE)
    }) %>% invisible
}


#--------------------------
# Hamming distance function
#--------------------------

Hamming <- function(event.matrix) {
    D <- (1 - event.matrix) %*% t(event.matrix)
    D + t(D)
}


#---------------------------
# add dummy column if absent
#---------------------------

DummyCols <- function(events, col.names, debug) {

    for(col.name in col.names) {
        if(!col.name %in% colnames(events)) {
            if(debug) { message(green(str_c('adding dummy column: ', col.name))) }
            events.names <- colnames(events) %>% replace(which(is.na(.)),'empty') %>% list.filter(.!='add.col')
            events$add.col <- NA
            events %<>% setNames(c(events.names,col.name))
        }
    }
    return(events)

}


#--------------------------
# convert chromosome format
#--------------------------

ChromMod <- function(events, allosome){

    if('X' %in% events$chrom & allosome != 'keep') { Warn('X chromosome labelling converted to numeric') }

    if(allosome == 'distinct') {
        if('Y' %in% events$chrom) { Warn('Y chromosome labelling converted to 24-based numeric') }
        events %>%
        mutate(.,chrom=as.numeric(ifelse(chrom=='X', 23,chrom))) %>%
        mutate(chrom=as.numeric(ifelse(chrom=='Y', 24,chrom)))
    } else if(allosome == 'merge') {
        if('Y' %in% events$chrom) { Warn('Y chromosome labelling converted to 23-based numeric') }
        events %>%
        mutate(chrom=as.numeric(ifelse(chrom %in% c('X','Y'), 23, chrom)))
    } else if(allosome == 'keep') {
        events
    } else {
        events %>%
        filter(!chrom %in% c('X','Y'))
    }
}


#-----------------------------------
# soft sample key:value modification
#-----------------------------------

KeyMod <- function(sample.object, keys, debug, force.keys=FALSE) {

    if(is.vector(sample.object)) {  # convert strings
        keys.index <- which(sample.object %in% names(keys))
        if(length(keys.index) == 0 & force.keys==FALSE) {

            if(debug) {
                message(green('no keys to convert'))
            }

        } else if(force.keys==FALSE) {

            conversions <- data.frame(pre=as.character(sample.object[keys.index]), post=keys[as.character(sample.object[keys.index])], row.names=NULL) %>% unique
            if(debug) {
                message(green(str_c('key conversions: ',length(conversions$pre),'/',length(unique(sample.object)))))
                print(conversions)
            }
            sample.object <- keys[sample.object]

        } else {

            conversions <- data.frame(pre=as.character(unique(sample.object)), post=keys[unique(as.character(sample.object))], in.keys=unique(as.character(events)) %in% names(keys), row.names=NULL)
            if(debug) {
                Warn('forcing key conversions')
                print(conversions)
            }
            sample.object[keys.index] <- keys[sample.object[keys.index]]
        }
    } else {  # convert data.frame

        keys.index <- which(sample.object$sample %in% names(keys))

        if(length(keys.index) == 0 & force.keys==FALSE) {
            if(debug) { message(green('no keys to convert')) }
        } else if(force.keys==FALSE) {
            conversions <- data.frame(pre=as.character(sample.object$sample[keys.index]), post=keys[as.character(sample.object$sample[keys.index])], row.names=NULL) %>% unique
            if(debug) {
                message(green(str_c('key conversions: ',length(conversions$pre),'/',length(unique(sample.object$sample)))))
                print(conversions)
            }
            sample.object %<>% mutate(sample=ifelse(row_number() %in% keys.index, keys[sample], sample))
        } else {
            conversions <- data.frame(pre=as.character(unique(sample.object$sample)), post=keys[unique(as.character(sample.object$sample))], in.keys=unique(as.character(events$sample)) %in% names(keys), row.names=NULL)
            if(debug) {
                Warn('forcing key conversions')
                print(conversions)
            }
            sample.object %<>% mutate(sample=keys[sample])
        }
    }
    sample.object
}


#-----------------------
# coerce columns to type
#-----------------------

TypeCol <- function(events, columns, types) {
    for (col in 1:length(columns)) {
        if(columns[col] %in% colnames(events)) {
            if(types[col] == 'char') {
                events[,columns[col]] <- events[,columns[col]] %>% unlist %>% as.character
            } else if(types[col] == 'num') {
                events[,columns[col]] <- events[,columns[col]] %>% unlist %>%  ifelse(.=='.',NA,.) %>% as.numeric
            } else if(types[col] == 'logic') {
                events[,columns[col]] <- events[,columns[col]] == 'TRUE' %>% as.vector
            }
        }
    }
    return(events)
}


#-----------------------
# melted table to matrix
#-----------------------

# function to convert event table into sample-gene occurance matrix
MeltToMatrix <- function(event.melt) {
    # reshape melted events to wide matrix for distance calculation
    event.melt %>%
    select(sample, span) %>%
    unique %>%
    mutate(exists=1) %>%
    spread(sample,exists, fill=0) %>%
    data.frame(., row.names=1, stringsAsFactors=FALSE, check.names=FALSE) %>%
    as.matrix %>%
    t
}


#------------------------
# distance method wrapper
#------------------------

# wrap Hamming function for easy specification
DistExtra <- function(event.matrix, dist.method){
    if(dist.method=='hamming'){
        dist.matrix <-
            event.matrix %>%
            Hamming %>%
            as.dist
    }else{
        dist.matrix <-
            event.matrix %>%
            dist(method=dist.method)
    }
    return(dist.matrix)
}


#---------------
# fix file names
#---------------

FixName <- function(file.name) {
    file.name %>% gsub('\\s', '_', .) %>% gsub('/', '-', .)
}


#-------------------------------------
# melted event table into ordered tree
#-------------------------------------

MeltToTree <- function(events, dist.method='hamming', clust.method='complete', sort.method='distance', span='gene', tree.samples) {

    events %<>% rename_(span=span)

    # build master event matrix
    event.matrix <- events %>% MeltToMatrix

    missing.samples <- tree.samples[!tree.samples %in% row.names(event.matrix)]

    if(length(missing.samples) > 0) {
        message(yellow('found samples with zero events, adding common root'))
        missing.matrix <- matrix(0, ncol=length(missing.samples), nrow=ncol(event.matrix)) %>% as_data_frame %>% set_names(missing.samples) %>% as.matrix %>% t
        event.matrix %<>% rbind(missing.matrix)
        event.matrix %<>% cbind( data_frame(project=rep(1,length(tree.samples))) %>% as.matrix)
    }

    # compute distance matrix using method specified
    dist <- DistExtra(event.matrix, dist.method)

    # cluster using method specified
    hc <- dist %>% hclust(method=clust.method)

    # construct dendrogram
    dend <-
        hc %>%
        as.dendrogram %>%
        set('branches_lwd', 10)

    # rotate tree, order using method specified
    if(sort.method == 'distance') {
        dend %<>% rotate_DendSer(ser_weight=dist(x))
    } else {
        dend %<>% ladderize
    }

    # reorder rows using distance matrix min
    dist <- reorder(dend, dist)

    tree <- list( event.matrix=event.matrix,
                  dist=dist,
                  hc=hc,
                  dend=dend,
                  dist=dist )

    return(tree)
}


#---------------------------------
# multiple pheno columns to single
#---------------------------------

MergePheno <- function(events, merge.cols, col.name='pheno') {
    events %>%
    do({ df <- .
        pheno <- select(df, one_of(merge.cols)) %>%
        apply(1, function(pheno){
            pheno %<>% na.omit
            if(length(pheno)==0){
                pheno=NA
            } else {
                pheno
            }
        })
        df <- mutate(df, col.name=pheno)
        colnames(df)[colnames(df) == 'col.name'] <- col.name
        return(df)
    })
}


#---------
# call LOH
#---------

CallLOH <- function(cncf.file) {
    cncf <-
        read.delim(cncf.file, stringsAsFactors=FALSE, sep="\t") %>%
        tbl_df %>%
        separate(ID, into=c('sample','normal'), sep='_') %>%
        DummyCols(c('lcn', 'tcn', 'cf'), debug=debug) %>%
        rowwise %>%
        mutate(lcn.em=replace(lcn.em, is.na(lcn.em), -Inf)) %>%
        mutate(lcn=replace(lcn, is.na(lcn), -Inf)) %>%
        ungroup %>%
        # loh calclulation
        dplyr::mutate(loh=
            ifelse(lcn.em==0,'LOH',
            ifelse(lcn.em!=-Inf,'.',
            ifelse(lcn==0,'LOH',
            ifelse(lcn!=-Inf,'.',
            ifelse(lcn.em==-Inf & lcn==-Inf, '.',
            NA)))))) %>%
        mutate(lcn=replace(lcn, lcn==-Inf, NA)) %>%
        mutate(lcn.em=replace(lcn.em, lcn.em==-Inf, NA)) %>%
        mutate(loc.mid=(loc.start+loc.end)/2) %>%
        mutate(id=row_number())

    cncf <-
        cncf %>%
        group_by(id) %>%
        full_join(
            cncf %>% filter(!is.na(loh)) %>%
            select(chrom, loc.mid.adj=loc.mid, cnlr.median.clust.adj=cnlr.median.clust, adj.lcn.em=lcn.em, adj.lcn=lcn, adj.loh=loh) %>%
            bind_rows(data_frame(chrom=1:max(cncf$chrom), loc.mid.adj=-Inf, cnlr.median.clust.adj=1, adj.lcn.em=1, adj.lcn=1, adj.loh='.')),
        by='chrom') %>%
        slice(which.min(abs(loc.mid-loc.mid.adj))) %>%
        ungroup %>%
        mutate(loh=ifelse(is.na(loh), adj.loh, loh)) %>%
        #mutate(loh=ifelse(loh=='LOH' & cnlr.median.clust < 0 & mafR > 0.2, 'LOH', '.')) %>%
        #mutate(loh=ifelse(loh=='.',NA,loh)) %>%
        mutate(seg.id=str_c(sample, ':', row_number())) %>%
        select(sample,chrom,loc.start,loc.mid,loc.end,loh,lcn.em,lcn,tcn.em,tcn,mafR,cnlr.median,cnlr.median.clust,cf,seg,num.mark,seg.id) %>%
        invisible

    return(cncf)
}






#----------------------------------
# prepare melted table for plotting
#----------------------------------

OrgEvents <- function(events, sample.order=NULL, pheno.order=NULL, sub.sets.pheno=NULL, recurrence=1, allosome='merge', event.type='gene', run.type=NULL, debug) {

    # sample.order:
    #   NULL = arrange alphabetically
    #   character string = specify specific order
    #
    # recurrence:
    #   integer cutoff for recurrent mutations
    #
    # allosome:
    #
    #   'none'      = exclude X & Y chromosome
    #   'merge'     = treat X & Y coordinates as homologous pair
    #   'distinct'  = treat X & Y as seperate, sequential chromosomes

    # specify output columns
    col.names <- c('sample','gene','chrom','effect','pheno','band','span','pos','maf','ccf','loh','hotspot', 'cancer.gene.census','clonality','pathogenic', 'k.l.c')

    events %<>% mutate(ccf=ifelse(is.na(ccf),1,ccf))

    events %<>% filter(!is.na(pheno)) %>% filter(!is.na(effect) | !is.na(ccf))

    if(is.null(sub.sets.pheno)) {
        sub.sets.pheno <- events %>% select(sample, pheno) %>% unique %>% spread(pheno, pheno)
    }

    if(is.null(pheno.order)) {
        pheno.order <- events$pheno %>% unique %>% sort
    }

    # add dummy column names
    events %<>% DummyCols(col.names, debug)

    # clip gene names to first if pipe seperated
    events %<>% rowwise %>% mutate(gene=gene %>% str_split('\\|') %>% .[[1]] %>% head(1)) %>% ungroup

    # effect prescedence
    events <-
        events %>%
        select(one_of(col.names)) %>%
        unique %>%
        filter(!is.na(effect)) %>%
        # rank variant importance for plot overlay
        mutate(precedence=
            ifelse(effect=='Deletion',1,
            ifelse(effect=='Amplification',2,
            ifelse(effect=='Truncating SNV',3,
            ifelse(effect=='Frameshift In-Del',4,
            ifelse(effect=='Missense SNV',5,
            ifelse(effect=='Inframe In-Del',6,
            ifelse(effect=='Splice site variant',7,
            ifelse(effect=='Upstream, start/stop, or de novo modification',8,
            ifelse(effect=='Silent',9,
            NA))))))))))

    if(event.type=='gene') {

        events <-
            events %>%
            # remove genes with lower prescedence
            group_by(sample, gene) %>%
            arrange(gene, precedence) %>%
            top_n(1) %>%
            # count number of variants per gene
            unique %>%
            ungroup %>%
            group_by(gene) %>%
            mutate(n.gene=n()) %>%
            ungroup %>%
            { events <- .
                if(recurrence > 0) { events %<>% filter(n.gene >= recurrence) }
                return(events)
            } %>%
            unique

        if(is.null(sample.order)) {

            sample.pool <- events$sample %>% unique

            while(length(sample.order) < length(sample.pool)) {

                sub.samples <- sample.pool %>% list.filter(!. %in% sample.order)
                take  = 1
                pick.sample = 0
                n.rounds = 1

                while(pick.sample!=1) {

                    pick.df <- lapply(sub.samples, function(sub.sample) {

                        sub.events <-
                            events %>%
                            select(sample, gene, n.gene) %>%
                            arrange(desc(n.gene)) %>%
                            filter(sample==sub.sample) %>%
                            bind_rows(data_frame(n.gene=c(0,0,0)))

                        num.v = sub.events$n.gene
                        sum.n <- combn(num.v, take) %>% apply(2, sum) %>% max

                        data_frame(sample=sub.sample, sum.n=sum.n)
                    }) %>%
                    bind_rows %>%
                    filter(sum.n==max(sum.n))

                    pick.sample = nrow(pick.df)

                    if(pick.sample == 1) {

                        sample.order <- c(sample.order, pick.df$sample)

                    } else if(n.rounds == 2) {

                        pick.sample = 1

                        sample.order <- c(sample.order,
                            events %>%
                            filter(sample %in% sub.samples) %>%
                            group_by(sample) %>%
                            summarise(burden=n()) %>%
                            arrange(desc(burden)) %>%
                            top_n(1, burden) %>%
                            .$sample )

                    } else {

                        take = take + 1
                        sub.samples = pick.df$sample
                        n.rounds = n.rounds + 1

                    }
                }
            }
        }

        CheckDepth <- function(gene) {
            found.gene = FALSE
            while(found.gene == FALSE) {
                for (sample.n in 1:length(sample.order)) {
                    if(gene %in% (events %>% filter(sample %in% sample.order[sample.n]) %>% .$ gene)) {
                        found.gene = TRUE
                        return(sample.n)
                    }
                }
            }
        }

        events %<>%
        rowwise %>%
        mutate(sample.order.depth=CheckDepth(gene)) %>%
        ungroup %>%
        arrange(desc(n.gene), sample.order.depth, sample, desc(ccf), precedence, gene)

        if(!is.null(run.type)) {
           events %<>%
            rowwise %>%
            mutate(gene=ifelse(k.l.c==TRUE, str_c('* ',gene), gene)) %>%
            ungroup
        }
    }

    if(event.type=='band') {
        events %<>%
            group_by(band) %>%
            mutate(n.band=n()) %>%
            ungroup %>%
            rowwise %>%
            mutate(band=str_c('chr',chrom,': ',band))
    }

    if(event.type=='span') {
        events %<>%
            # remove spans with lower prescedence
            group_by(sample,span) %>%
            arrange(span,precedence) %>%
            top_n(1) %>%
            # count number of variants per gene
            unique %>%
            ungroup %>%
            group_by(span) %>%
            mutate(n.span=n()) %>%
            ungroup %>%
            { events <- .
                if(recurrence>0) {events %<>% filter(n.span>=recurrence)}
                return(events)
            } %>%
            unique %>%
            arrange(desc(n.span),sample,desc(ccf),precedence,span)
    }

    # define plot gene order
    events %<>%
        mutate(order=row_number()) %>%
        ungroup %>%
        select(-precedence) %>%
        unique

    missing.samples <- sample.order[!sample.order %in% events$sample]

    missing.pheno <-
        sub.sets.pheno %>%
        MergePheno(merge.cols=pheno.order, col.name='pheno') %>%
        select(sample, pheno) %>%
        filter(sample %in% missing.samples)


    if(length(missing.samples) > 0){

        message(yellow(str_c('missing specified samples: ',missing.samples,'\n')))

        missing.fill <-
            expand.grid(missing.samples,unique(unlist(events[,event.type])), stringsAsFactors=FALSE) %>%
            set_names(c('sample', event.type)) %>%
            tbl_df %>%
            left_join(missing.pheno) %>%
            invisible
    }

    if(event.type=='gene') {
        events.fill <-
            events %>%
            select(sample,effect,pheno,gene) %>%
            unique %>%
            group_by(pheno) %>%
            spread(gene,effect) %>%
            gather(gene,effect,-pheno,-sample) %>%
            filter(is.na(effect)) %>%
            ungroup
    } else if(event.type=='band') {
        events.fill <-
            events %>%
            select(sample,chrom,effect,pheno,band) %>%
            unique %>%
            group_by(pheno) %>%
            spread(band,effect) %>%
            gather(band,effect,-chrom,-pheno,-sample) %>%
            filter(is.na(effect)) %>%
            ungroup
    } else {
        events.fill <-
            events %>%
            select(sample,chrom,effect,pheno,span) %>%
            unique %>%
            group_by(pheno) %>%
            spread(span,effect) %>%
            gather(span,effect,-chrom,-pheno,-sample) %>%
            filter(is.na(effect)) %>%
            ungroup
    }

    if(length(missing.samples) > 0){
        events %<>% bind_rows(missing.fill)
    }

    if(nrow(events.fill) > 0) {
        events %<>% full_join(events.fill) %>% invisible
    }

    # plot aesthetics
    events %<>%
        filter(!is.na(sample)) %>%
        filter_(paste('!is.na(', event.type, ')' )) %>%
        filter(sample %in% sample.order) %>%
        filter(!is.na(pheno)) %>%
        # push NAs to bottom of stack
        mutate(order=ifelse(is.na(effect),-Inf,order)) %>%
        unique %>%
        # fix plot ordering & assign gene factor levels
        arrange(!is.na(effect),desc(order)) %>%
        mutate(gene=factor(gene,levels=filter(.,!is.na(effect)) %>% .$gene %>% unique)) %>%
        mutate(span=factor(span,levels=filter(.,!is.na(effect)) %>% .$span %>% unique)) %>%
        # ccf binning
        mutate(ccf=
            ifelse(ccf==0,               'CCF = 0%',
            ifelse(ccf>0.00 & ccf<=0.05, '0% < CCF <= 5%',
            ifelse(ccf>0.05 & ccf<=0.20, '5% < CCF <= 20%',
            ifelse(ccf>0.20 & ccf<=0.40, '20% < CCF <= 40%',
            ifelse(ccf>0.40 & ccf<=0.60, '40% < CCF <= 60%',
            ifelse(ccf>0.60 & ccf<=0.80, '60% < CCF <= 80%',
            ifelse(ccf>0.80,             '80% < CCF <= 100%',
            NA)))))))) %>%
        mutate(ccf=ifelse(is.na(effect),NA,ccf)) %>%
        mutate(ccf=factor(ccf, levels=c('CCF = 0%','0% < CCF <= 5%','5% < CCF <= 20%','20% < CCF <= 40%','40% < CCF <= 60%','60% < CCF <= 80%','80% < CCF <= 100%'))) %>%
        mutate(clonal=ifelse(clonality=='Clonal', clonality, NA)) %>%
        select(-order) %>%
        mutate(sample=factor(sample, levels=sample.order))

    # change working columns to NA
    events %<>% mutate(loh=ifelse(loh == '.', NA, loh))

    if(event.type=='gene') {
        events %<>%
            filter(!is.na(gene)) %>%
            mutate(effect=ifelse(hotspot==TRUE, 'Hotspot', effect))

    } else if(event.type=='band') {
        events %<>%
            arrange(!is.na(effect), desc(n.band), desc(chrom), desc(band)) %>%
            mutate(band=factor(band, levels=filter(.,!is.na(effect)) %>% .$band %>% unique)) %>%
            filter(!is.na(band))
    } else if( event.type=='span') {
        events %<>% filter(!is.na(span))
    }

    return(events)
}


#-----------------------
# main plotting function
#-----------------------

PlotVariants <- function(events, output.file, clonal=FALSE, cancer.gene.census=FALSE, pathogenic=FALSE, ccf=FALSE, loh=TRUE, width=15, height=7, text.size=6, event.type='gene'){

    # set graphics device
    options(device=pdf)

    # rename for plot output
    events %<>% plyr::rename(replace=c(sample='Sample', gene='Gene', band='Band', span='Span', variant='Variant', effect='Effect', cancer.gene.census='Carcinogenic', clonal='clonal', ccf='CCF', cn='CN'), warn_duplicated=FALSE)

    # plot aesthetic definitions
    palette <- c( 'Truncating SNV'='#C84DDD',
                  'Frameshift In-Del'='#C17200',
                  'Missense SNV'='#00A5A8',
                  'Inframe In-Del'='#E44988',
                  'Splice site variant'='#008AE9',
                  'Upstream, start/stop, or de novo modification'='#749000',
                  'Silent'='#666666',
                  'Hotspot'='#d64242',
                  'Amplification'='#333399',
                  'Deletion'='#e60000',
                  'CCF = 0%'='#e5e5e5',
                  `0% < CCF <= 5%`='#c7dbee',
                  `5% < CCF <= 20%`='#a0cae0',
                  `20% < CCF <= 40%`='#6eaed4',
                  `40% < CCF <= 60%`='#2772b3',
                  `60% < CCF <= 80%`='#10539a',
                  `80% < CCF <= 100%`='#0b3269' )

    # main plot params
    if(event.type=='gene') {
        hp <- ggplot(events, aes(Sample, Gene, drop=FALSE))
    } else if(event.type=='band') {
        hp <- ggplot(events, aes(Sample, Band, drop=FALSE))
    } else {
        hp <- ggplot(events, aes(Sample, Span, drop=FALSE))
    }

    # CCF coloring
    if(ccf == TRUE) {
        hp <- hp +
        geom_tile(data=events, aes(fill=CCF, drop=FALSE), colour='white') +
        scale_fill_manual(breaks=names(palette), values=palette, na.value='gray90', drop=FALSE)
    } else {
        # draw tiles and color
        hp <- hp +
        geom_tile(data=events, aes(fill=Effect, drop=FALSE), colour='white') +
        scale_fill_manual(breaks=names(palette), values=palette, na.value='gray90', drop=FALSE)
    }

    if(loh==TRUE & !all(is.na(events$loh))) {
        hp <- hp +
            # LOH slashes
            geom_segment(data=ggplot_build(hp)$data[[1]][which(events$loh == 'LOH'), ],
                         aes(x=xmin, xend=xmax, y=ymin, yend=ymax),
                         color='white',
                         size=1)
    }

    if(clonal == TRUE & !all(is.na(events$clonal))) {
        hp <- hp +
            geom_tile(data=events %>% filter(!is.na(Effect) & !is.na(clonal)), aes(colour=clonal), width=0.95, height=0.95, size=1, fill=NA) +
            scale_color_manual(values='#DCA43E')
    }

    hp <- hp +
        # specify legend
        guides(colour='white') +
        guides(colour=guide_legend(override.aes=list(alpha=1, fill=NA)))

    hp <- hp +
        # tile groups
        facet_wrap(~pheno, nrow=1, scales='free_x') +
        scale_x_discrete(expand=c(0, 0.5)) +
        scale_y_discrete(expand=c(0, 0.5))

    hp <- hp +
        # theme params
        theme( legend.title        = element_blank(),
               panel.grid.major    = element_blank(),
               panel.grid.minor    = element_blank(),
               text                = element_text(size=18),
               axis.title.x        = element_blank(),
               axis.title.y        = element_blank(),
               axis.text.x         = element_text(angle=90, vjust=0.5, hjust=1, margin=margin(0,12,0,0)),
               axis.text.y         = element_text(face='italic', hjust=1),
               axis.text           = element_text(size=text.size),
               axis.ticks.x        = element_blank(),
               axis.ticks.y        = element_blank(),
               legend.key          = element_rect(colour='black', fill=NULL, size=1),
               legend.key.size     = unit(1.8, 'lines'),
               legend.text         = element_text(size=text.size),
               strip.background    = element_rect(fill='grey'),
               strip.text.x        = element_text(colour='white', size=text.size*1.2),
               panel.background    = element_rect(fill=NA, color=NA),
               panel.border        = element_rect(fill=NA, colour='black', size=1),
               plot.margin         = unit(c(1,1,1,1), 'pt') )

    # build grob object
    hpg <- suppressWarnings(ggplotGrob(hp))

    # draw plot
    pdf(output.file, width, height, bg='white')
        grid.draw(hpg)
    dev.off()
}


#------------
# CNA heatmap
#------------

PlotCNHeatmap <- function(gene.cn, file.name, sample.names=NULL, threshold=FALSE) {

    # set graphics device
    options(device=pdf)

    if(is.null(sample.names) & threshold==TRUE) {
        sample.names <- gene.cn %>% select(matches('threshold')) %>% names %>% sort
        sample.names <- sample.names[sample.names %>% gsub('[^0-9]', '',.) %>% as.numeric %>% order]
    } else if(is.null(sample.names)) {
        sample.names <- gene.cn %>% names %>% list.filter(! . %in% c('hgnc', 'gene', 'chrom', 'start', 'mid', 'end', 'band'))
        sample.names <- sample.names[sample.names %>% gsub('[^0-9]', '',.) %>% as.numeric %>% order]
    }

    chr.rle <- gene.cn$chrom %>% rle
    chr.sep <- chr.rle$lengths %>% cumsum
    chr.mid <- c(0, chr.sep[-length(chr.sep)]) + chr.rle$lengths/2

    pdf(file.name, width=24, height=2+length(sample.names)/2)

        g.cn <- gene.cn %>% select(one_of(rev(sample.names)))

        layout(matrix(c(0,1),2,1,byrow=TRUE), c(length(sample.names),1), TRUE)

        par(mfrow=c(1,1), mar=c(8,5,1,1))
        image(as.matrix(g.cn), col=c('#CF3A3D', '#DC9493', '#FFFFFF', '#7996BA', '#2A4B94'), xaxt='n', yaxt='n', zlim=c(-2, 2))

        for (i in seq(-1, max(((2*(ncol(g.cn)-1))+1),1), 2)) {
            abline(h=i/(2*(ncol(g.cn)-1)), col='white', lwd=2)
        }

        for (i in (chr.sep*2)-1) {
            abline(v=i/((max(chr.sep)-1)*2), lwd=1.5, col='grey', lty='dashed')
        }

        box()

         axis( 1,
               at       = chr.mid/(max(chr.sep)-1),
               label    = chr.rle$values,
               cex.axis = 1.3,
               tick     = FALSE )

         axis( 2,
               at       = if(ncol(g.cn)==1){ 0.5 }else{seq(0, 1, 1/max((ncol(g.cn)-1),1))},
               label    = sub('T_.*', '', colnames(g.cn)),
               las      = 2,
               cex.axis = 1.1,
               tick     = FALSE )

         legend( 'bottom',
                 inset  = c(0, -0.5),
                 legend = c('Homozygous deletion', 'Loss', 'Gain', 'Amplification'),
                 fill   = c('#CF3A3D', '#DC9493', '#7996BA', '#2A4B94'),
                 xpd    = TRUE,
                 ncol   = 2,
                 cex    = 1.1 )
     dev.off()

}
























SubGroup <- function(sub.num) {

   #------------------------
   H2('setup sub group run')
   #------------------------

   # // 

   # specify sub group info
   sub.group <- sub.groups %>% filter(group.id == sub.num) %>% mutate(extension=str_c(sub.run, '_', extension))
   sub.group.sets <- sub.group[names(sub.group) %>% list.filter(!. %in% c('group.id', 'comparison', 'extension', 'n.samples'))] %>% list.filter(!is.na(.)) %>% unlist
   sub.group.samples <- sub.group.sets %>% list.map(sub.sets[.]) %>% unlist %>% unique %>% KeyMod(keys, debug) %>% sort
   sub.group.prefix <- str_c('summary/sub_group/', sub.group$extension, '/')


   # print sub group
   H1(str_c('SUB NUM: ', sub.num))
   print(sub.group)


   # make output directories
   c('mutation_heatmap', 'cna_heatmap', 'variant_heatmap', 'tree', 'venn') %>%
   str_c(sub.group.prefix, .) %>%
   MakeDirs(debug)


   #------------------------------
   H2('define subsetted variants')
   #------------------------------

   # create pheno column & filter samples to sub group set
   sub.variants <-
      run.variants %>%
      select(-pheno) %>%
      filter(sample %in% sub.group.samples) %>%
      MergePheno(merge.cols=sub.group.sets, col.name='pheno') %>%
      #MergePheno(merge.cols=sub.group.samples, col.name='sample.colors') %>%
      filter(!is.na(pheno))


   # subset mutations
   sub.muts <- sub.variants %>% filter(class=='muts') %>% rename(gene=span)


   # subset cnas
   sub.cnas <- sub.variants %>% filter(class=='cnas') %>% rename(band=span)


   #------------------------
   H2('write subset tables')
   #------------------------

   # write cnas file
   sub.variants %>% write_tsv(str_c(sub.group.prefix, 'variant_summary_', sub.group$extension, '.tsv'))


   # write mutations raw file
   sub.muts %>% write_tsv(str_c(sub.group.prefix, 'mutation_summary_', sub.group$extension, '.tsv'))


   # write mutations summary files
   sub.muts.summary <-
        sub.muts %>%
        arrange(sample, effect) %>%
        rowwise %>%
        mutate_each(funs(ifelse(is.character(.),ifelse(is.na(.),'.',.),.))) %>%
        mutate(chasm=ifelse(chasm <= 0.3,'Driver', ifelse(chasm != Inf,'Passenger','.'))) %>%
        mutate(fathmm=ifelse(fathmm == 'CANCER', 'Driver', ifelse(fathmm == 'PASSENGER/OTHER', 'Passenger', ifelse(fathmm == 'PASSENGER', 'Passenger', '.')))) %>%
        mutate(kandoth=ifelse(kandoth !=TRUE, '.', 'TRUE')) %>%
        mutate(lawrence=ifelse(lawrence !=TRUE, '.', 'TRUE')) %>%
        mutate(hotspot=ifelse(hotspot !=TRUE, '.', 'TRUE')) %>%
        mutate(haplo=ifelse(haplo !=TRUE, '.', 'TRUE')) %>%
        mutate(cancer.gene.census=ifelse(cancer.gene.census == TRUE, 'TRUE', '.')) %>%
        select(`Sample ID`=sample, Gene=gene, Chromosome=chrom, Position=pos, AA=hgvs.p, Effect=variant, LOH=loh, CCF=ccf, `Clonal Probability`=pr.somatic.clonal,
            `95% Confidence Interval Low`=ci95.low, Clonality=clonality, `Reference Allele`=ref, `Alternate Allele`=alt, `Tumor MAF`=maf.t, `Normal MAF`=maf.n,
            `Tumor Depth`=depth.t, `Normal Depth`=depth.n, `Mutation Taster`=mut.taster, FATHMM=fathmm, Chasm=chasm, `Cancer Gene Census`=cancer.gene.census,
            Kandoth=kandoth, Lawrence=lawrence, Haploinsufficient=haplo, Pathogenicity=pathogenic, Hotspot=hotspot, `Cancer Gene`=k.l.c )

   sub.muts.summary %>% write_tsv(str_c(sub.group.prefix, 'mutation_summary_', sub.group$extension, '.tsv'))
   sub.muts.summary %>% write.xlsx(str_c(sub.group.prefix, 'mutation_summary_', sub.group$extension, '.xlsx'))

   if(sub.muts$sample %>% unique %>% length == 1) {
        sub.muts %<>% mutate(gene=str_c(gene, hgvs.p, sep='-'))
   }


   # write master summary sheet sheet
   if(sub.num==0) {
      gene.cn.curated <-
         read.delim(opts$gene_cn_in, sep='\t', stringsAsFactors=FALSE) %>%
         mutate(chrom=as.character(chrom)) %>%
         mutate(chrom=ifelse(chrom==23, 'X', chrom)) %>%
         rowwise %>%
         mutate(mid=(start+end)/2) %>%
         ungroup %>%
         tbl_df

    if(run.metrics) {

          metrics <-
             read.delim('metrics/hs_metrics.tsv', sep='\t', stringsAsFactors=FALSE) %>%
             tbl_df

        if('SAMPLE' %in% names(metrics)) {
            metrics %<>%
                DummyCols('PCT_TARGET_BASES_40X', debug) %>%
                DummyCols('PCT_TARGET_BASES_100X', debug) %>%
                select(`Sample ID`=SAMPLE, `Target Territory`=TARGET_TERRITORY, `Percent Selected Bases`=PCT_SELECTED_BASES, `Mean Target Coverage`=MEAN_TARGET_COVERAGE,
                    `Percent Selected Bases 2X`=PCT_TARGET_BASES_2X, `Percent Selected Bases 10X`=PCT_TARGET_BASES_10X, `Percent Selected Bases 20X`=PCT_TARGET_BASES_20X,
                    `Percent Selected Bases 30X`=PCT_TARGET_BASES_30X, `Percent Selected Bases 40X`=PCT_TARGET_BASES_40X, `Percent Selected Bases 50X`=PCT_TARGET_BASES_50X,
                    `Percent Selected Bases 100X`=PCT_TARGET_BASES_100X )
        } else {
            metrics %<>%
                DummyCols('Percent.Target.Bases.40X', debug) %>%
                DummyCols('Percent.Target.Bases.100X', debug) %>%
                select(`Sample ID`=Sample.ID, `Target Territory`=Target.Territory, `Percent Selected Bases`=Percent.Selected.Bases, `Mean Target Coverage`= Mean.Target.Coverage,
                    `Percent Selected Bases 2X`=Percent.Target.Bases.2X, `Percent Selected Bases 10X`=Percent.Target.Bases.10X, `Percent Selected Bases 20X`=Percent.Target.Bases.20X,
                    `Percent Selected Bases 30X`=Percent.Target.Bases.30X, `Percent Selected Bases 40X`=Percent.Target.Bases.40X, `Percent Selected Bases 50X`=Percent.Target.Bases.50X,
                    `Percent Selected Bases 100X`=Percent.Target.Bases.100X )
        }

        metrics %<>%
            left_join(samples %>% select(normal, tumor) %>% gather(`Tissue Type`, `Sample ID`), by='Sample ID') %>%
            mutate(`Tissue Type`=ifelse(`Tissue Type` == 'tumor', 'Tumor', ifelse(`Tissue Type` == 'normal', 'Normal', `Tissue Type`))) %>%
            mutate(`Percent Selected Bases`=str_c(round(100 * `Percent Selected Bases`, 2), '%')) %>%
            mutate(`Mean Target Coverage`=round(`Mean Target Coverage`, 2)) %>%
            mutate(`Percent Selected Bases 2X`=str_c(round(100 * `Percent Selected Bases 2X`, 2), '%')) %>%
            mutate(`Percent Selected Bases 10X`=str_c(round(100 * `Percent Selected Bases 10X`, 2), '%')) %>%
            mutate(`Percent Selected Bases 20X`=str_c(round(100 * `Percent Selected Bases 20X`, 2), '%')) %>%
            mutate(`Percent Selected Bases 30X`=str_c(round(100 * `Percent Selected Bases 30X`, 2), '%')) %>%
            mutate(`Percent Selected Bases 40X`=str_c(round(100 * `Percent Selected Bases 40X`, 2), '%')) %>%
            mutate(`Percent Selected Bases 50X`=str_c(round(100 * `Percent Selected Bases 50X`, 2), '%')) %>%
            mutate(`Percent Selected Bases 100X`=str_c(round(100 * `Percent Selected Bases 100X`, 2), '%')) %>%
            select(`Sample ID`, `Tissue Type`, everything()) %>%
            arrange(`Sample ID`)

      list(`Sequencing_Statistics`=metrics, `Mutations`=sub.muts.summary, `geneCN`=gene.cn.curated) %>%
      write.xlsx(str_c(sub.group.prefix, 'master_summary_', sub.group$extension, '.xlsx'))

    } else {

      list(`Mutations`=sub.muts.summary, `geneCN`=gene.cn.curated) %>%
      write.xlsx(str_c(sub.group.prefix, 'master_summary_', sub.group$extension, '.xlsx'))

    }

   }


   # write cnas file
   sub.cnas %>%
   arrange(sample, effect) %>%
   write_tsv(str_c(sub.group.prefix, 'cna_summary_', sub.group$extension, '.tsv'))

   # write new geneCN file
   gene.cn %>%
   select(gene, chrom, start, end, band, one_of(sub.group.samples)) %>%
   arrange(chrom, start) %>%
   write.xlsx(str_c(sub.group.prefix, 'geneCN_', sub.group$extension, '.xlsx'))

   gene.cn %>%
   select(gene, chrom, start, end, band, one_of(sub.group.samples)) %>%
   arrange(chrom, start) %>%
   write_tsv(str_c(sub.group.prefix, 'geneCN_', sub.group$extension, '.tsv'))


   #-------------------------------
   H2('select event table columns')
   #-------------------------------

   sub.variants %<>% select(class, sample, chrom, pos, span, effect, hotspot, ccf, loh, clonality, cancer.gene.census, k.l.c, haplo, pheno) %>% mutate(sub.id=str_c(chrom, pos, span, sep='-'))
   sub.muts %<>% select(class, sample, chrom, pos, gene, effect, hotspot, ccf, loh, clonality, cancer.gene.census, k.l.c, haplo, pheno) %>% mutate(sub.id=str_c(chrom, pos, gene, sep='-'))
   sub.cnas %<>% select(class, sample, chrom, pos, band, effect, hotspot, ccf, loh, clonality, cancer.gene.census, k.l.c, haplo, pheno) %>% mutate(sub.id=str_c(chrom, pos, band, sep='-'))


   # // -- sub group initialization

   #---------------------
   H2('CN HEATMAP PLOTS')
   #---------------------

   # // 

   if(run.cn.heatmap.plots != FALSE) {

      for(sub.set.num in 1:length(sub.group.sets)) {

         sub.set.name <- sub.sets[sub.group.sets][sub.set.num] %>% names %>% FixName
         sub.set <- sub.sets[sub.group.sets][[sub.set.num]] %>% KeyMod(keys, debug)

         # cut down gene.cn table
         sub.gene.cn <- gene.cn %>% select(gene, chrom, start, end, band, one_of(sub.set))

         # call CN heatmap plotting function
         PlotCNHeatmap(sub.gene.cn, file.name=str_c(sub.group.prefix, 'cna_heatmap/genecn_heatmap_', sub.group$extension, '_',sub.set.name,'.pdf'), threshold=FALSE)
      }

   }

   # // -- cn heatmap plots

  #------------------
  H2('CASCADE PLOTS')
  #------------------

   # // 

   if(run.cascade.plots != FALSE) {

        if(nrow(sub.muts) > 0) {

          H3('muts [type] cascade plot')

          sub.muts %>%
          OrgEvents(sample.order=NULL, pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=1, allosome='merge', event.type='gene', run.type=NULL, debug) %>%
          PlotVariants( output.file = str_c(sub.group.prefix, 'mutation_heatmap/mutation_heatmap_type_', sub.group$extension, '.pdf'),
                        clonal      = FALSE,
                        pathogenic  = FALSE,
                        ccf         = FALSE,
                        loh         = FALSE,
                        event.type  = 'gene',
                        width       = (length(unique(.$sample))/2) + 4.5,
                        height      = (length(unique(.$gene))/6) + 3,
                        text.size   = 12 )


          H3('muts [type] [recurrent] cascade plot')

          if( sub.muts %>% select(sample, gene) %>% group_by(sample) %>% unique %>% ungroup %>% .$gene %>% duplicated %>% any ) {
             sub.muts %>%
             OrgEvents(sample.order=NULL, pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=2, allosome='merge', event.type='gene', run.type=NULL, debug) %>%
             PlotVariants( output.file = str_c(sub.group.prefix, 'mutation_heatmap/mutation_heatmap_type_recurrent_', sub.group$extension, '.pdf'),
                           clonal      = FALSE,
                           pathogenic  = FALSE,
                           ccf         = FALSE,
                           loh         = FALSE,
                           event.type  = 'gene',
                           width       = (length(unique(.$sample))/2) + 4.5,
                           height      = (length(unique(.$gene))/6) + 3,
                           text.size   = 14 )
          }

          H3('muts [CCF] cascade plot')

          sub.muts %>%
          OrgEvents(sample.order=NULL, pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=1, allosome='merge', event.type='gene', run.type='ccf', debug) %>%
          PlotVariants( output.file = str_c(sub.group.prefix, 'mutation_heatmap/mutation_heatmap_ccf_', sub.group$extension, '.pdf'),
                        clonal      = TRUE,
                        pathogenic  = FALSE,
                        ccf         = TRUE,
                        loh         = TRUE,
                        event.type  = 'gene',
                        width       = (length(unique(.$sample))/2) + 4.8,
                        height      = (length(unique(.$gene))/6) + 3,
                        text.size   = 12 )


          H3('muts [CCF] [recurrent] cascade plot')

          if( sub.muts %>% select(sample, gene) %>% group_by(sample) %>% unique %>% ungroup %>% .$gene %>% duplicated %>% any ) {
             sub.muts %>%
             OrgEvents(sample.order=NULL, pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=2, allosome='merge', event.type='gene', run.type='ccf', debug) %>%
             PlotVariants( output.file = str_c(sub.group.prefix, 'mutation_heatmap/mutation_heatmap_ccf_recurrent_', sub.group$extension, '.pdf'),
                           clonal      = TRUE,
                           pathogenic  = FALSE,
                           ccf         = TRUE,
                           loh         = TRUE,
                           event.type  = 'gene',
                           width       = (length(unique(.$sample))/2) + 4.8,
                           height      = (length(unique(.$gene))/6) + 3,
                           text.size   = 14 )
          }

    }

      H3('cnas cascade plot')

      if(nrow(sub.cnas) > 0) {

         # fix sample order
         sample.order.1 <-
             sub.muts %>%
             OrgEvents(sample.order=NULL, pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=1, allosome='merge', event.type='gene', run.type=NULL, debug) %>%
             .$sample %>%
             as.character %>%
             rev %>%
             unique

          sub.cnas %>%
          OrgEvents(sample.order=sample.order.1, pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=1, allosome='merge', event.type='band', run.type=NULL, debug) %>%
          PlotVariants( output.file = str_c(sub.group.prefix, 'cna_heatmap/cna_heatmap_', sub.group$extension, '.pdf'),
                        clonal      = FALSE,
                        pathogenic  = FALSE,
                        ccf         = FALSE,
                        loh         = FALSE,
                        event.type  = 'band',
                        width       = (length(unique(.$sample))/2) + 4.5,
                        height      = (length(unique(.$gene))/6) + 3,
                        text.size   = 12 )

         if( sub.muts %>% select(sample, gene) %>% group_by(sample) %>% unique %>% ungroup %>% .$gene %>% duplicated %>% any ) {
            # fix sample order [recurrent]
            sample.order.2 <-
                sub.muts %>%
                OrgEvents(sample.order=NULL, pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=2, allosome='merge', event.type='gene', run.type=NULL, debug) %>%
                .$sample %>%
                as.character %>%
                rev %>%
                unique
         }

         sub.cnas %>%
         OrgEvents(sample.order=sample.order.1, pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=1, allosome='merge', event.type='band', run.type=NULL, debug) %>%
         PlotVariants( output.file = str_c(sub.group.prefix, 'cna_heatmap/cna_heatmap_recurrent-muts_', sub.group$extension, '.pdf'),
                       clonal      = FALSE,
                       pathogenic  = FALSE,
                       ccf         = FALSE,
                       loh         = FALSE,
                       event.type  = 'band',
                       width       = (length(unique(.$sample))/2) + 4.5,
                       height      = (length(unique(.$gene))/6) + 3,
                       text.size   = 12 )

      }


      #--------------------------
      H2('combined variant plots')
      #--------------------------

      if(variant.plots != FALSE) {

          H3('variants cascade plot')

          sub.variants %>%
          OrgEvents(sample.order=sample.order.1, pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=1, allosome='merge', event.type='span', run.type=NULL, debug) %>%
          PlotVariants( output.file = str_c(sub.group.prefix, 'variant_heatmap/variant_heatmap_type_', sub.group$extension, '.pdf'),
                        clonal      = FALSE,
                        pathogenic  = FALSE,
                        ccf         = FALSE,
                        loh         = FALSE,
                        event.type  = 'span',
                        width       = (length(unique(.$sample))/2) + 4.5,
                        height      = (length(unique(.$gene))/6) + 3,
                        text.size   = 12 )


          H3('VARIANTS [recurrent] cascade plot')

         if( sub.variants %>% select(sample, span) %>% group_by(sample) %>% unique %>% ungroup %>% .$span %>% duplicated %>% any ) {
             sub.variants %>%
             OrgEvents(sample.order=sample.order.2, pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=2, allosome='merge', event.type='span', run.type=NULL, debug) %>%
             PlotVariants( output.file = str_c(sub.group.prefix, 'variant_heatmap/variant_heatmap_type_recurrent_', sub.group$extension, '.pdf'),
                           clonal      = FALSE,
                           pathogenic  = FALSE,
                           ccf         = FALSE,
                           loh         = FALSE,
                           event.type  = 'span',
                           width       = (length(unique(.$sample))/2) + 4.5,
                           height      = (length(unique(.$gene))/6) + 3,
                           text.size   = 14 )
         }
      }
   }

   # // -- cascade plots

   #------------------
   H2('venn diagrams')
   #------------------

   # // 

   if(run.venn != FALSE) {

       MakeVenn <- function(sets, file.name) {

          venn <- Venn(Sets=sets)

        pdf(file.name, 10, 10)
             if(length(sets) == 3) {
                plot(venn, show=list(Faces=FALSE))
             } else if (length(sets) < 5 ) {
               suppressWarnings(plot(venn, type='ellipses', show=list(Faces=FALSE)))
             } else {
                suppressWarnings(plot(venn, type='AWFE', show=list(Faces=FALSE)))
             }
        dev.off()

       }

       MakeVennSquare <- function(sets, file.name) {

          venn <- Venn(Sets=sets)

        pdf(file.name, 10, 10)
             if(length(sets) < 5) {
                suppressWarnings(plot(venn, type='squares', show=list(Faces=FALSE)))
             } else {
                suppressWarnings(plot(venn, type='AWFE', show=list(Faces=FALSE)))
             }
        dev.off()

       }

       if(length(sub.group.samples) < 8) {

          venn.list.base <- sub.group.samples %>% list.map(filter(sub.muts, sample == .) %>% .$sub.id %>% unique)
          venn.list.base.effect <- sub.muts %>% split(.$effect) %>% list.map(.$sub.id %>% unique)
          venn.list.base.cancer.gene <- list(cancer.gene=sub.muts %>% filter(k.l.c==TRUE) %>% .$sub.id %>% unique)

          venn.list.gene <- sub.group.samples %>% list.map(filter(sub.muts, sample == .) %>% .$gene %>% unique)
          venn.list.gene.effect <- sub.muts %>% split(.$effect) %>% list.map(.$gene %>% unique)
          venn.list.gene.cancer.gene <- list(cancer.gene=sub.muts %>% filter(k.l.c==TRUE) %>% .$gene %>% unique)

          if(length(sub.group.samples) > 1) {
             MakeVenn(sets=venn.list.base, file.name=str_c(sub.group.prefix, 'venn/muts_venn_base_', sub.group$extension, '.pdf'))
             MakeVenn(sets=venn.list.gene, file.name=str_c(sub.group.prefix, 'venn/muts_venn_gene_', sub.group$extension, '.pdf'))
          }

          if(any(duplicated(venn.list.base))) {
             MakeVennSquare(sets=c(venn.list.base, venn.list.base.effect['Frameshift In-Del']), file.name=str_c(sub.group.prefix, 'venn/muts_venn_base_effect_fs_indel_', sub.group$extension, '.pdf'))
             MakeVennSquare(sets=c(venn.list.base, venn.list.base.effect['Inframe In-Del']), file.name=str_c(sub.group.prefix, 'venn/muts_venn_base_effect_if_indel_', sub.group$extension, '.pdf'))
             MakeVennSquare(sets=c(venn.list.base, venn.list.base.effect['Missense SNV']), file.name=str_c(sub.group.prefix, 'venn/muts_venn_base_effect_mis_snv_', sub.group$extension, '.pdf'))
             MakeVennSquare(sets=c(venn.list.base, venn.list.base.effect['Silent']), file.name=str_c(sub.group.prefix, 'venn/muts_venn_base_effect_silent_', sub.group$extension, '.pdf'))
             MakeVennSquare(sets=c(venn.list.base, venn.list.base.effect['Splice site variant']), file.name=str_c(sub.group.prefix, 'venn/muts_venn_base_effect_splice_', sub.group$extension, '.pdf'))
             MakeVennSquare(sets=c(venn.list.base, venn.list.base.effect['Truncating SNV']), file.name=str_c(sub.group.prefix, 'venn/muts_venn_base_effect_trunc_snv_', sub.group$extension, '.pdf'))
             MakeVennSquare(sets=c(venn.list.base, venn.list.base.effect['Upstream, start/stop, or de novo modification']), file.name=str_c(sub.group.prefix, 'venn/muts_venn_base_effect_upstream_', sub.group$extension, '.pdf'))
             MakeVennSquare(sets=c(venn.list.base, venn.list.base.cancer.gene), file.name=str_c(sub.group.prefix, 'venn/muts_venn_base_effect_effect_cancer_gene_', sub.group$extension, '.pdf'))
          }

          if(any(duplicated(venn.list.gene))) {
             MakeVennSquare(sets=c(venn.list.gene, venn.list.gene.effect['Frameshift In-Del']), file.name=str_c(sub.group.prefix, 'venn/muts_venn_gene_effect_fs_indel_', sub.group$extension, '.pdf'))
             MakeVennSquare(sets=c(venn.list.gene, venn.list.gene.effect['Inframe In-Del']), file.name=str_c(sub.group.prefix, 'venn/muts_venn_gene_effect_if_indel_', sub.group$extension, '.pdf'))
             MakeVennSquare(sets=c(venn.list.gene, venn.list.gene.effect['Missense SNV']), file.name=str_c(sub.group.prefix, 'venn/muts_venn_gene_effect_mis_snv_', sub.group$extension, '.pdf'))
             MakeVennSquare(sets=c(venn.list.gene, venn.list.gene.effect['Silent']), file.name=str_c(sub.group.prefix, 'venn/muts_venn_gene_effect_silent_', sub.group$extension, '.pdf'))
             MakeVennSquare(sets=c(venn.list.gene, venn.list.gene.effect['Splice site variant']), file.name=str_c(sub.group.prefix, 'venn/muts_venn_gene_effect_splice_', sub.group$extension, '.pdf'))
             MakeVennSquare(sets=c(venn.list.gene, venn.list.gene.effect['Truncating SNV']), file.name=str_c(sub.group.prefix, 'venn/muts_venn_gene_effect_trunc_snv_', sub.group$extension, '.pdf'))
             MakeVennSquare(sets=c(venn.list.gene, venn.list.gene.effect['Upstream, start/stop, or de novo modification']), file.name=str_c(sub.group.prefix, 'venn/muts_venn_gene_effect_upstream_', sub.group$extension, '.pdf'))
             MakeVennSquare(sets=c(venn.list.gene, venn.list.gene.cancer.gene), file.name=str_c(sub.group.prefix, 'venn/muts_venn_gene_effect_cancer_gene_', sub.group$extension, '.pdf'))
          }
       }
    }

   # // -- venn diagrams

   #----------
   H2('TREES')
   #----------

   # // 

   if(run.trees != FALSE) {

      # sub muts tree
      sub.muts.tree <- tryCatch({
            sub.muts %>% MeltToTree(dist.method='hamming', clust.method='complete', sort.method='distance', span='gene', tree.samples=sub.group.samples)
      }, error = function(e) {
          Warn('muts MeltToTree encountered an error, skipping')
          NULL
      })

      # sub cnas tree
      sub.cnas.tree <- tryCatch({
            sub.cnas %>% MeltToTree(dist.method='hamming', clust.method='complete', sort.method='distance', span='band', tree.samples=sub.group.samples)
      }, error = function(e) {
          Warn('cnas MeltToTree encountered an error, skipping')
          NULL
      })

      # sub variants tree
      sub.variants.tree <- tryCatch({
            sub.variants %>% MeltToTree(dist.method='hamming', clust.method='complete', sort.method='distance', span='span', tree.samples=sub.group.samples)
      }, error = function(e) {
          Warn('variants MeltToTree encountered an error, skipping')
          NULL
      })

      # remove labels from tree
      if(tree.labels == FALSE){
          if(!is.null(sub.muts.tree)) { suppressWarnings(sub.muts.tree$dend %<>% set('labels', '')) }
          if(!is.null(sub.cnas.tree)) { suppressWarnings(sub.cnas.tree$dend %<>% set('labels', '')) }
          if(!is.null(sub.variants.tree)) { suppressWarnings(sub.variants.tree$dend %<>% set('labels', '')) }
      }

      #----------------------------------
      H2('sub group lineage calculation')
      #----------------------------------

      if(!is.null(sub.muts.tree)) {

         # determine possible tree sequences
         seq.patterns <- sub.muts.tree$event.matrix %>% rbind(., root=0) %>% phyDat(type='USER', levels=0:10/10)

         # tree hamming distance
         phylo.dist <- seq.patterns %>% dist.hamming %>% njs

         # determine most parsimonious arrangement
         phylo.parsimony <- seq.patterns %>% pratchet(start=phylo.dist) %>% acctran(data=seq.patterns)

         # root tree
         phylo.root <- phylo.parsimony %>% phytools::reroot(., node.number=which(.$tip.label == 'root'))


         # convert phylo to dentrogram using chronological estimation
         dend.root <-
             phylo.root %>%
             chronos(model='discrete') %>%
             as.dendrogram %>%
             invisible


         # tree aesthetics
         dend.ext <-
             dend.root %>%
             set('branches_lwd', 3) %>%
             set('labels_cex', 1) %>%
             set('branches_k_color', k=length(sub.group.samples)) %>%
             dendextend::ladderize # %>% prune('root')


         # default rectilinear tree
         pdf(str_c(sub.group.prefix, 'tree/muts_lineage_rect_', sub.group$extension, '.pdf'))
             plot(phylo.root)
         dev.off()


         # triangular tree
         pdf(str_c(sub.group.prefix, 'tree/muts_lineage_tri_', sub.group$extension, '.pdf'), width=8, height=6)
             par(mar=c(1,1,1,6))  # bottom, left, top, right
             plot(dend.ext, horiz=TRUE, type='triangle', edge.root=FALSE, axes=FALSE)
         dev.off()
      }
  }

  # // -- trees

  #-------------------------
  H2('TREES (experimental)')
  #-------------------------

   # // 

   if(run.experimental != FALSE) {

      #---------------------
      # distance tree ladder
      #---------------------

      pdf(str_c(sub.group.prefix, 'tree/genomic_distance_tree_ladder_', sub.group$extension, '.pdf'), 140, 38)
         par(mar=c(10, 10, 10, 10))  # bottom, left, top, right
         layout(matrix(c(1, 2), nrow=1), widths=c(10, 1))
         plot(sub.muts.tree$dend, cex.axis=6)

         colored_bars(colors               = sub.muts %>% select(pheno),
                      dend                 = sub.muts.tree$dend,
                      sort_by_labels_order = FALSE,
                      add                  = TRUE,
                      rowLabels            = sub.muts$sample,
                      y_scale              = 10,
                      cex.rowLabels        = 5.8)

         legend(x=3, y=3, legend=unique(sub.muts$sample), fill=unique(sub.muts$pheno), cex=4)
         legend(x=3, y=3, legend=unique(sub.muts$pheno), fill=unique(sub.muts$pheno), cex=4)
      dev.off()


      #--------------
      # weighted tree
      #--------------

      pdf(str_c(sub.group.prefix, 'tree/genomic_distance_tree_weighted_', sub.group$extension, '.pdf'),60,25)
         suppressWarnings(heatmap.2(
             sub.muts.tree$event.matrix,
             trace='none',
             Rowv=FALSE,
             hclustfun=function(x){hclust(x, 'ward.D2')},
             col=c('white','black'),
             reorderfun=function(d,w) rev(reorder(d,w)),
             distfun=function(x) as.dist(Hamming(x)),
             key=FALSE ))
      dev.off()


      #--------------
      # distance tree
      #--------------

      pdf(str_c(sub.group.prefix, 'tree/genomic_distance_tree_bw_', sub.group$extension, '.pdf'),60,25)
         heatmap.2(
             t(sub.muts.tree$event.matrix),
             trace='none',
             dendrogram='column',
             Colv=rev(sub.muts.tree$dend),
             col=c('white','black'),
             key=FALSE
         )
      dev.off()


      #----------------------------
      # CNAS LINEAGE MAP GENERATION
      #----------------------------

      lineage.matrix = cbind(sub.cnas.tree$event.matrix, Parental=FALSE) %>% t

      phylo.data <- phyDat(lineage.matrix, type='USER', levels=c(0, 1))
      phylo.hamming <- dist.hamming(phylo.data)
      phylo.tree <- njs(phylo.hamming)
      phylo.ratchet <- pratchet(phylo.data, start=phylo.tree) %>% acctran(phylo.data)
      lineage.dend <- root(phylo.ratchet, 'Parental')

      pdf(str_c(sub.group.prefix, 'tree/variant_lineage_cnas_', sub.group$extension, '.pdf'))
         plot(lineage.dend)
      dev.off()


      #-------------------------------
      # VARIANT LINEAGE MAP GENERATION
      #-------------------------------

      lineage.matrix = cbind(sub.variants.tree$event.matrix, Parental=FALSE) %>% t

      phylo.data <- phyDat(lineage.matrix, type="USER", levels=c(0, 1))
      phylo.hamming <- dist.hamming(phylo.data)
      phylo.tree <- njs(phylo.hamming)
      phylo.ratchet <- pratchet(phylo.data, start=phylo.tree) %>% acctran(phylo.data)
      lineage.dend <- root(phylo.ratchet, "Parental")

      pdf(str_c(sub.group.prefix, 'tree/variant_lineage_variants_', sub.group$extension, '.pdf'))
         plot(lineage.dend)
      dev.off()

  }

  # // -- trees (experimental)

   #-------------------
   H2("FISHER'S PLOTS")
   #-------------------

   # // 

   if(run.fishers.plots != FALSE & !is.na(sub.group$b)) {

      if(length(sub.group.samples) > 2) {

         #---------------------------
         # Fisher's exact CN plotting
         #---------------------------

         samples.a <- sub.sets[sub.group$a] %>% unlist %>% KeyMod(keys, debug)
         samples.b <- sub.sets[sub.group$b] %>% unlist %>% KeyMod(keys, debug)

         # copy number plotting
         gene.cn.a <- gene.cn[c('gene', 'chrom', 'start', 'end', samples.a)]
         gene.cn.b <- gene.cn[c('gene', 'chrom', 'start', 'end', samples.b)]

         # call fishers plots for copy number
         Fisher(plot.type       = 'copy number',
                gene.matrix.a   = gene.cn.a,
                gene.matrix.b   = gene.cn.b,
                plot.title.main = sub.group$comparison,
                plot.title.a    = sub.group$a,
                plot.title.b    = sub.group$b,
                allosome        = allosome,
                targets.file    = targets.file,
                suffix          = sub.group$extension,
                threshold.a     = FALSE,
                threshold.b     = FALSE,
                gene.names      = gene.names)

         #---------------------------------
         # Fisher's exact mutation plotting
         #---------------------------------

         # mutation plotting
         gene.muts.a <- sub.muts %>% select(sample, gene, chrom, pos, effect) %>% mutate(effect=1) %>% filter(sample %in% samples.a) %>% spread(sample, effect, fill=0)
         gene.muts.b <- sub.muts %>% select(sample, gene, chrom, pos, effect) %>% mutate(effect=1) %>% filter(sample %in% samples.b) %>% spread(sample, effect, fill=0)

         # call fishers function for mutations
         Fisher(plot.type       = 'mutation',
                gene.matrix.a   = gene.muts.a,
                gene.matrix.b   = gene.muts.b,
                plot.title.main = sub.group$comparison,
                plot.title.a    = sub.group[1],
                plot.title.b    = sub.group[2],
                allosome        = allosome,
                targets.file    = targets.file,
                suffix          = sub.group$extension,
                gene.names      = gene.names)
      }

   }

   # // -- fisher's plots

}

