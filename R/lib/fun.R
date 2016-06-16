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
                   'band'='band', 'Band'='band',
                   'cancer.gene.census'='cancer.gene.census', 'Cancer.Gene.Census'='cancer.gene.census', 'Cancer Gene Census'='cancer.gene.census','Cancer5000.S.genes..Lawrence.et.al.'='cancer.gene.census', 'cancer_gene_census'='cancer.gene.census', 'cancer.gene.census'='cancer.gene.census',
                   'ccf'='ccf', 'CCF'='ccf', 'cancer_cell_frac'='ccf', 'Cancer cell fraction'='ccf', 'Cancer.cell.fraction'='ccf',
                   'chasm'='chasm', 'CHASM'='chasm', 'Breast_chasm_score'='chasm',
                   'chrom'='chrom','Chrom'='chrom','Chromosome'='chrom','CHROM'='chrom',
                   'Clonal'='clonal', 'clonal'='clonal',
                   'cn'='cn', 'CN'='cn',
                   'ci95.low'='ci95.low', 'ccf_CI95_low'='ci95.low',
                   'dbNSFP_PROVEAN_pred'='provean', 'provean'='provean',
                   'NORMAL.DP'='depth.n', 'depth.n'='depth.n',
                   'TUMOR.DP'='depth.t', 'depth.t'='depth.t',
                   'Effect'='variant', 'Variant_Classification'='variant', 'ANN....EFFECT'='variant', 'ANN[*].EFFECT'='variant',
                   'effect'='effect',
                   'ExAC_AF'='ex.af', 'ex.af'='ex.af',
                   'end'='end', 'stop'='end', 'End_position'='end',
                   'fathmm'='fathmm', 'FATHMM'='fathmm', 'fathmm_pred'='fathmm',
                   'gene'='gene', 'Gene'='gene', 'Hugo_Symbol'='gene', 'GENE'='gene', 'hgnc'='gene', 'ANN....GENE'='gene', 'ANN[*].GENE'='gene', 'Gene.symbol'='gene', 'Hugo_Symbol'=='gene',
                   'haploinsufficient'='haplo', 'hap_insuf'='haplo', 'haplo'='haplo',
                   'ANN[*].HGVS_C'='hgvs.c',
                   'ANN[*].HGVS_P'='hgvs.p',
                   'ANN[*].IMPACT'='impact', 'impact'='impact',
                   'kandoth'='kandoth', '127 significantly mutated genes (Kandoth et al)'='kandoth', '127 significantly.mutated genes.(Kandoth et al)'='kandoth','X127.significantly.mutated.genes..Kandoth.et.al.'='kandoth',
                   'lawrence'='lawrence', 'Cancer5000-S genes (Lawrence et al)'='lawrence', 'Cancer5000-S genes.(Lawrence et al)'='lawrence',
                   'loh'='loh', 'LOH'='loh', 'Loss.of.heterozygocity.(LOH)'='loh', 'Loss of heterozygocity (LOH)'='loh','Loss.of.heterozygocity..LOH.'='loh',
                   'NORMAL_SAMPLE'='normal', 'normal'='normal',
                   'NORMAL_MAF'='maf.n', 'maf.n'='maf.n',
                   'maf.t'='maf.t', 'Mutant allele fraction', 'maf.t', 'Mutant.allele.fraction'='maf.t', 'TUMOR_MAF'='maf.t',
                   'mut.taster'='mut.taster', 'dbNSFP_MutationTaster_pred'='mut.taster',
                   'pathogenic'='pathogenic', 'Pathogenic'='pathogenic', 'pathogenicity'='pathogenic', 'Pathogenicity'='pathogenic',
                   'pheno'='pheno', 'pheno.bar'='pheno',
                   'pos'='pos','POS'='pos', 'Position'='pos', 'position'='position',
                   'pr.somatic.clonal'='pr.somatic.clonal', 'Pr_somatic_clonal'='pr.somatic.clonal',
                   'provean'='provean', 'Provean'='provean',
                   'purity'='purity',
                   'ref'='ref', 'REF'='ref', 'Reference.allele'='ref',
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
        col.names <- c( 'alt','band','cancer.gene.census','ccf','chasm','chrom','clonal','cn','ci95.low','effect', 'variant',
                        'end','fathmm','gene','haploinsufficient','kandoth','lawrence','loh','maf','mut.taster',
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
    events %<>% TypeCol( c('sample', 'gene', 'effect', 'variant', 'ref',  'alt',  'cancer.gene.census', 'kandoth', 'lawrence', 'haplo', 'fathmm', 'chasm'),
                         c('char',   'char', 'char',   'char',    'char', 'char', 'logic',       'logic',   'logic',    'logic', 'char',   'num') )

    if('variant' %in% colnames(events) & !'effect' %in% colnames(events)) {
        events %<>% mutate(effect=variant)
    }

    if('effect' %in% colnames(events)){

        # clip effect names to first if pipe or '&' seperated
        events %<>% rowwise %>% mutate(effect=effect %>% str_split('\\|') %>% .[[1]] %>% head(1)) %>% ungroup
        events %<>% rowwise %>% mutate(effect=effect %>% str_split('&') %>% .[[1]] %>% head(1)) %>% ungroup

        # rename variant classifications
        events %<>%
            tbl_df %>%
            { events <- .
                if(!all(is.na(events$effect))){ filter(events,!is.na(effect)) }
                events
            } %>%
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

    return(events)
}


#---------------------
# CNA heatmap function
#---------------------

PlotCNHeatmap <- function(gene.cn, file.name, sample.names=NULL, threshold=FALSE) {

    if(all(class(gene.cn) == 'character')) {
      gene.cn <- read.delim(gene.cn, sep='\t', stringsAsFactors=FALSE) %>% tbl_df
    }

    if(is.null(sample.names) & threshold==TRUE) {
        sample.names <- gene.cn %>% select(matches('threshold')) %>% names %>% sort
    } else if(is.null(sample.names)) {
        sample.names <- gene.cn %>% names %>% list.filter(! . %in% c('hgnc','gene','chrom','start','mid','end','band'))
    }

    # convert X & Y to 23
    gene.cn %<>% mutate(chrom=as.numeric(ifelse(chrom == 'X' | chrom == 'Y', '23', chrom)))

    # segment spans
    chr.rle <- gene.cn$chrom %>% rle
    chr.sep <- chr.rle$lengths %>% cumsum

    # midpoints
    chr.mid <- c(0, chr.sep[-length(chr.sep)]) + chr.rle$lengths/2

    # remove annotation cols
     gene.cn %<>% select(one_of(rev(sample.names)))

    pdf(file.name, width=24, height=2+length(sample.names)/3)

        par(mar=c(14, 14, 1, 1), oma=c(1, 1, 1, 1))  # bottom, left, top, right

        image(as.matrix(gene.cn), col=c('#cf3a3d', '#dc9493', '#FFFFFF', '#7996ba', '#2a4b94'), xaxt='n', yaxt='n', zlim=c(-2, 2))

        for (i in seq(-1, max(((2*(ncol(gene.cn)-1))+1),1), 2)) {
            abline(h=i/(2*(ncol(gene.cn)-1)), col="white", lwd=2)
        }

        for (i in (chr.sep*2)-1) {
            abline(v=i/((max(chr.sep)-1)*2), lwd=1.5, col='grey', lty="dashed")
        }

         axis( 1,
               at       = chr.mid/(max(chr.sep)-1),
               label    = chr.rle$values,
               cex.axis = 1.3,
               tick     = FALSE )

         axis( 2,
               at       = if(ncol(gene.cn)==1){0.5}else{seq(0, 1, 1/max((ncol(gene.cn)-1),1))},
               label    = if(threshold==TRUE){ sub('_LRR_threshold$', '', colnames(gene.cn)) }else{ colnames(gene.cn) },
               las      = 2,
               cex.axis = 1,
               tick     = FALSE )

        box()

         legend( 'bottom',
                 inset  = c(0, -0.15),
                 legend = c('Amplification', 'Gain', 'Loss', 'Homozygous deletion'),
                 fill   = c('#2a4b94', '#7996ba', '#dc9493', '#cf3a3d'),
                 xpd    = TRUE,
                 ncol   = 2,
                 cex    = 1.1 )
     dev.off()

}
