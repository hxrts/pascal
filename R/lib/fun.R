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
            mutate(effect=gsub('^\\s+|\\s+$', '', effect)) %>%
            unique %>%
            mutate(effect=
                ifelse(effect%in%c('STOP_GAINED','Nonsense_Mutation','stop_gained&splice_region_variant','stop_gained','Nonsense_Mutation','Stop_Codon_Ins','nonsense','truncating snv','Truncating snv','Truncating snv','Truncating SNV'),'Truncating SNV',
                ifelse(effect%in%c('FRAME_SHIFT','FRAME_SHIFT','Frame_Shift_Del','Frame_Shift_Ins','frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant','frameshift_variant','frameshift_variant&stop_gained','frameshift_variant&splice_region_variant','frameshift_variant&splice_acceptor_variant&splice_region_variant&splice_region_variant&intron_variant','Frame_Shift_Del','Frame_Shift_Ins','frame_shift_del','frame_shift_ins','frameshift indel','Frameshift indel','Frameshift In-Del','frameshift_variant'),'Frameshift In-Del',
                ifelse(effect%in%c('NON_SYNONYMOUS_CODING','STOP_LOST','Missense_mutation','missense_variant','missense_variant&splice_region_variant','missense_variant|missense_variant','Missense_Mutation','missense','missense snv','Missense snv','Missense SNV'),'Missense SNV',
                ifelse(effect%in%c('CODON_CHANGE_PLUS_CODON_DELETION','CODON_DELETION','CODON_INSERTION','In_Frame_Ins','In_Frame_Del','disruptive_inframe_deletion','disruptive_inframe_insertion','inframe_deletion','inframe_insertion','disruptive_inframe_deletion&splice_region_variant','inframe_deletion&splice_region_variant','In_Frame_Del','In_Frame_Ins','in_frame_del','in_frame_ins','inframe indel','Inframe indel','Inframe In-Del'),'Inframe In-Del',
                ifelse(effect%in%c('splice_donor_variant','splice_region_variant','splice_acceptor_variant','SPLICE_SITE_DONOR','SPLICE_SITE_ACCEPTOR','SPLICE_SITE_REGION','Splice_Site','splice_donor_variant&intron_variant','splice_acceptor_variant&intron_variant','splicing','splice_donor_variant&splice_region_variant&intron_variant','splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&splice_region_variant&intron_variant','Splice_Site','splice','splice site variant','Splice site variant','missense_variant & splice_region_variant'),'Splice site variant',
                ifelse(effect%in%c('STOP_LOST','START_LOST','START_GAINED','UTR_5_PRIME','start_lost','stop_lost',"5'UTR","5'Flank",'De_novo_Start_InFrame','De_novo_Start_OutOfFrame','Stop_Codon_Del','Start_Codon_SNP','Start_Codon_Ins','Start_Codon_Del','Nonstop_Mutation','nonstop','upstream, start/stop, or de novo modification','Upstream, start/stop, or de novo modification'),'Upstream, start/stop, or de novo modification',
                ifelse(effect%in%c('synonymous_variant','intron_variant','splice_region_variant&synonymous_variant','splice_region_variant&synonymous_variant','non_coding_exon_variant','upstream_gene_variant','downstream_gene_variant','intron_variant','frameshift_variant&splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant','non_coding_exon_variant|synonymous_variant','SYNONYMOUS_CODING','synonymous_variant|synonymous_variant','splice_region_variant&synonymous_variant|splice_region_variant&non_coding_exon_variant','splice_acceptor_variant & intron_variant','intragenic_variant',"3'UTR",'IGR','lincRNA','RNA','Intron','silent','intron_exon','silent','Silent','intron_variant & missense_variant'),'Silent',
                ifelse(effect%in%c('Amplification','amplification','amp','2'),'Amplification',
                ifelse(effect%in%c('Gain','gain','1'),'Gain',
                ifelse(effect%in%c('Loss','loss','-1'),'Loss',
                ifelse(effect%in%c('Deletion','deletion','del','-2'),'Deletion',
                ifelse(is.na(effect), NA, effect)))))))))))))

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


#----------------------------------
# prepare melted table for plotting
#----------------------------------

OrgEvents <- function(events, sample.order=NULL, pheno.order=NULL, sub.sets.pheno=NULL, recurrence=1, allosome='merge', event.type='gene', run.type=NULL, debug=FALSE) {

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
    col.names <- c('sample','gene','chrom','effect','pheno','band','span','pos','maf','ccf','loh','cancer.gene.census','clonality','pathogenic', 'k.l.c')

    if(!'ccf' %in% names(events)) {
      events %<>% mutate(ccf=1)
    }

    # add dummy column names
    events %<>% DummyCols(col.names, debug)

    # remove rows which can't be plotted
    events %<>% filter(!is.na(pheno)) %>% filter(!is.na(effect) | !is.na(ccf))

    if(is.null(sub.sets.pheno)) {
        sub.sets.pheno <- events %>% select(sample, pheno) %>% unique %>% spread(pheno, pheno)
    }

    if(is.null(pheno.order)) {
        pheno.order <- events$pheno %>% unique %>% sort
    }

    # clip gene names to first if pipe seperated
    events %<>% rowwise %>% mutate(gene=gene %>% str_split('\\|') %>% .[[1]] %>% head(1)) %>% ungroup

    # effect prescedence
    events %<>%
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

        events %<>%
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
        events <-
            events %>%
            # remove spans with lower prescedence
            group_by(sample, span) %>%
            arrange(span, precedence) %>%
            top_n(1, precedence) %>%
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
            arrange(desc(n.span), sample, desc(ccf), precedence, span)
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
        events <- invisible(full_join(events, events.fill))
    }

    if(is.null(sample.order)) {
        sample.order <- events$sample %>% unique %>% sort
    }

    # plot aesthetics
    events <-
        events %>%
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
        events %<>% filter(!is.na(gene))
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

PlotVariants <- function(events, output.file, clonal=FALSE, cancer.gene.census=FALSE, pathogenic=FALSE, ccf=FALSE, loh=TRUE, width=7, height=7, text.size=6, event.type='gene'){

    # set graphics device
    options(device=pdf)

    # rename for plot output
    events %<>% plyr::rename(replace=c(sample='Sample', gene='Gene', band='Band', span='Span', variant='Variant', effect='Effect', cancer.gene.census='Carcinogenic', pathogenic='Pathogenic', clonal='clonal', ccf='CCF', cn='CN'), warn_duplicated=FALSE)

    # plot aesthetic definitions
    palette  <- c( 'Truncating SNV'='#C84DDD',
                    'Frameshift In-Del'='#C17200',
                    'Missense SNV'='#00A5A8',
                    'Inframe In-Del'='#E44988',
                    'Splice site variant'='#008AE9',
                    'Upstream, start/stop, or de novo modification'='#749000',
                    'Silent'='#666666',
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
    if(ccf==TRUE) {
        hp <- hp +
        geom_tile(data=events, aes(fill=CCF, drop=FALSE), colour='white') +
        scale_fill_manual(breaks=names(palette), values=palette, na.value='gray90', drop=FALSE)
    } else {  # draw tiles and color
        hp <- hp +
        geom_tile(data=events, aes(fill=Effect, drop=FALSE), colour='white') +
        scale_fill_manual(breaks=names(palette), values=palette, na.value='gray90', drop=FALSE)
    }

    # loh slashes
    if(loh==TRUE & !all(is.na(events$loh))) {
        hp <- hp + geom_segment(data=ggplot_build(hp)$data[[1]][which(events$loh=='LOH'), ],
                         aes(x=xmin, xend=xmax, y=ymin, yend=ymax),
                         color='white',
                         size=1)
    }

    if(clonal==TRUE & !all(is.na(events$clonality))) {
        hp <- hp +
        geom_tile(data=events %>% filter(!is.na(Effect) & !is.na(clonal)), aes(colour=clonal), size=1, fill=NA) +
        scale_color_manual(values='#DCA43E')
    }

    # if(pathogenic==TRUE & !all(is.na(events$pathogenic))) {
    #     hp <- hp +
    #     geom_point(data=events, aes(shape=loh, stroke=1.5), size=2, colour='#ff0015') +
    #     scale_shape_manual(values=c(`LOH`=3), guide=guide_legend(colour = '#ff0015'))
    # }

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
           axis.text.x         = element_text(angle=90, vjust=0.5, hjust=1, margin=margin(0,10,0,0)),
           axis.text.y         = element_text(face='italic', hjust=1),
           axis.text           = element_text(size=text.size),
           axis.ticks.x        = element_blank(),
           axis.ticks.y        = element_blank(),
           legend.key          = element_rect(colour='black', fill=NULL, size=1),
           legend.key.size     = unit(1.8, 'lines'),
           legend.text         = element_text(size=text.size),
           strip.background    = element_rect(fill='white'),
           strip.text.x        = element_text(colour='white', size=text.size*1.2),
           panel.background    = element_rect(fill=NA, color=NA),
           panel.border        = element_rect(fill=NA, colour='black', size=2),
           plot.margin         = unit(c(1,1,1,1), 'pt') )

    # build grob object
    hpg <- suppressWarnings(ggplotGrob(hp))

    # # count number of samples in each group
    # plot.lengths <- events %>% split(.$pheno) %>% map(~ .x$Sample %>% unique %>% length) %>% unlist

    # # get the column indexcorresponding to the panels.
    # panelI <- hpg$layout$l[grepl('panel', hpg$layout$name)]

    # # replace the default panel widths with relative heights.
    # hpg$widths <- grid:::unit.list(hpg$widths)
    # hpg$widths[panelI] <- lapply(plot.lengths, unit, 'null')

    # # add extra width between panels
    # for(gap in 1:(length(panelI)-1)){
    #     hpg$widths[panelI[gap]+1]=list(unit(0.3, 'cm'))
    # }

    # draw plot
    pdf(output.file, width, height, bg='white')
        grid.draw(hpg)
    dev.off()
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
