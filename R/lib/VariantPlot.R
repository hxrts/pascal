#!/usr/bin/env Rscript

#----------
# LIBRARIES
#----------

pacman::p_load( dplyr,readr,tidyr,magrittr,purrr,stringr,rlist, # base
				crayon,colorspace,RColorBrewer, # coloring
				ggplot2,grid,gridExtra,gplots, # plot layout
				dendextend,dendextendRcpp,dynamicTreeCut,gclus, # dentrogram
				digest ) # hashing

#-------------
# INPUT PARAMS
#-------------

file.name      = 'analysis/plot.variants.pdf'
include.silent = TRUE
dist.method    = 'hamming'   # 'hamming', euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'     | see ?dist for details [hamming method implemented manually]
clust.method   = 'complete' # 'complete', ward', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid' | see ?hclust for details
exclude.values = c('', '.', 'Normal', 'Not Performed', 'Performed but Not Available', FALSE)
tree.sort      = 'ladderize' #'distance'  # ladder
tree.labels    = TRUE
pheno          = NULL
color.seed     = 3


#----------
# FUNCTIONS
#----------

# Hamming distance function
Hamming <- function(event.matrix) {
	D <- (1 - event.matrix) %*% t(event.matrix)
	D + t(D)
}


# wrap Hamming function, for easy specification
DistExtra <- function(event.matrix,dist.method){
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


# character hashing
hash <- function(x) {
	hstr <- digest(x, algo='xxhash32')
	as.numeric(paste0('0x', hstr)) %% 80
}


# clean muts nomenclature
FormatMuts <- function(muts) {

	# lookup table
	name.key <- c( 'sample'='sample', 'Sample'='sample', 'Tumor_Sample_Barcode'='sample', 'TUMOR_SAMPLE'='sample',
				   'gene'='gene', 'Gene'='gene', 'Hugo_Symbol'='gene', 'GENE'='gene', 'hgnc'='gene', 'ANN....GENE'='gene', 'ANN[*].GENE'='gene',
				   'effect'='effect', 'Effect'='effect', 'Variant_Classification'='effect', 'ANN....EFFECT'='effect', 'ANN[*].EFFECT'='effect',
				   'ccf'='ccf', 'CCF'='ccf', 'cancer_cell_frac'='ccf',
				   'loh'='loh', 'LOH'='loh',
				   'pheno'='pheno', 'pheno.bar'='pheno',
				   'chrom'='chrom','Chrom'='chrom','Chromosome'='chrom','CHROM'='chrom',
				   'start'='start', 'loc.start'='start',
				   'end'='end', 'loc.end'='end', 'end'='stop',
				   'POS'='pos', 'pos'='pos',
				   'pathogenic'='pathogenic', 'Pathogenic'='pathogenic',
				   'band'='band', 'Band'='band',
				   'clonality'='clonality', 'Clonality'='clonality' )

	# rename columns
	names(muts) <- name.key[names(muts)]

	# df typing
	muts %<>% as.data.frame

	DummyCol <- function(muts, col.name) {
		# add dummy ccf column if absent
		if(!col.name %in% colnames(muts)) {
			message(green(str_c('adding dummy column: ',col.name)))
			muts.names <- colnames(muts) %>% replace(which(is.na(.)),'empty') %>% list.filter(.!='add.col')
			muts$add.col <- NA
			muts %<>% setNames(c(muts.names,col.name))
		}
		return(muts)
	}

	muts %<>% DummyCol('gene')
	muts %<>% DummyCol('ccf')
	muts %<>% DummyCol('loh')
	muts %<>% DummyCol('pheno')
	muts %<>% DummyCol('chrom')
	muts %<>% DummyCol('band')
	muts %<>% DummyCol('start')
	muts %<>% DummyCol('end')
	muts %<>% DummyCol('pos')
	muts %<>% DummyCol('pathogenic')
	muts %<>% DummyCol('clonality')

	# select columns
	muts <- muts[c('sample','gene','effect','ccf','loh','pheno','chrom','start','end','pos','band','pathogenic','clonality')]

	# column typing
	muts %<>% mutate(sample=as.character(sample), gene=as.character(gene), effect=as.character(effect))

	# rename variant classifications & remove unmutated LOH if present
	muts %<>%
		tbl_df %>%
		filter(!is.na(effect)) %>%
		unique %>%
		mutate(effect=
			ifelse(effect%in%c('STOP_GAINED','Nonsense_Mutation','stop_gained&splice_region_variant','stop_gained','Nonsense_Mutation','Stop_Codon_Ins','nonsense','truncating snv','Truncating snv','Truncating snv','Truncating SNV'),'Truncating SNV',
			ifelse(effect%in%c('FRAME_SHIFT','FRAME_SHIFT','Frame_Shift_Del','Frame_Shift_Ins','frameshift_variant','frameshift_variant&stop_gained','frameshift_variant&splice_region_variant','frameshift_variant&splice_acceptor_variant&splice_region_variant&splice_region_variant&intron_variant','Frame_Shift_Del','Frame_Shift_Ins','frame_shift_del','frame_shift_ins','frameshift indel','Frameshift indel','Frameshift In-Del'),'Frameshift In-Del',
			ifelse(effect%in%c('NON_SYNONYMOUS_CODING','STOP_LOST','Missense_Mutation','missense_variant','missense_variant&splice_region_variant','missense_variant|missense_variant','Missense_Mutation','missense','missense snv','Missense snv','Missense SNV'),'Missense SNV',
			ifelse(effect%in%c('CODON_CHANGE_PLUS_CODON_DELETION','CODON_DELETION','CODON_INSERTION','In_Frame_Ins','In_Frame_Del','disruptive_inframe_deletion','disruptive_inframe_insertion','inframe_deletion','inframe_insertion','disruptive_inframe_deletion&splice_region_variant','inframe_deletion&splice_region_variant','In_Frame_Del','In_Frame_Ins','in_frame_del','in_frame_ins','inframe indel','Inframe indel','Inframe In-Del'),'Inframe In-Del',
			ifelse(effect%in%c('SPLICE_SITE_DONOR','SPLICE_SITE_ACCEPTOR','SPLICE_SITE_REGION','Splice_Site','splice_donor_variant&intron_variant','splice_acceptor_variant&intron_variant','splicing','splice_donor_variant&splice_region_variant&intron_variant','splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&splice_region_variant&intron_variant','Splice_Site','splice','splice site variant','Splice site variant'),'Splice site variant',
			ifelse(effect%in%c('STOP_LOST','START_LOST','START_GAINED','UTR_5_PRIME','start_lost','stop_lost',"5'UTR","5'Flank",'De_novo_Start_InFrame','De_novo_Start_OutOfFrame','Stop_Codon_Del','Start_Codon_SNP','Start_Codon_Ins','Start_Codon_Del','Nonstop_Mutation','nonstop','upstream, start/stop, or de novo modification','Upstream, start/stop, or de novo modification'),'Upstream, start/stop, or de novo modification',
			ifelse(effect%in%c('synonymous_variant','splice_region_variant&synonymous_variant','non_coding_exon_variant','upstream_gene_variant','downstream_gene_variant','intron_variant','frameshift_variant&splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant','non_coding_exon_variant|synonymous_variant','SYNONYMOUS_CODING','synonymous_variant|synonymous_variant','splice_region_variant&synonymous_variant|splice_region_variant&non_coding_exon_variant','intragenic_variant',"3'UTR",'IGR','lincRNA','RNA','Intron','silent','intron_exon','silent','Silent'),'Silent',
			ifelse(effect%in%c('Amplification','amplification','amp','2'),'Amplification',
			ifelse(effect%in%c('Deletion','deletion','del','-2'),'Deletion',
			ifelse(is.na(effect),NA,
		NA))))))))))) %>%
		mutate(loh=ifelse(loh%in%c('LOH','loh','Loss of heterozygosity') & !is.na(effect),'LOH',NA))

	# warn on unknown effects
	if(muts %>% filter(is.na(effect)) %>% nrow > 0){ message(yellow('warning: some variants not accounted for')) }

	return(muts)
}


# prepare muts table for plotting
OrgMuts <- function(muts, dend) {

# sort: default, ordinal
	DummyCol <- function(muts, col.name) {
		# add dummy ccf column if absent
		if(!col.name %in% colnames(muts)) {
			message(green(str_c('adding dummy column: ',col.name)))
			muts.names <- colnames(muts) %>% replace(which(is.na(.)),'empty') %>% list.filter(.!='add.col')
			muts$add.col <- NA
			muts %<>% setNames(c(muts.names,col.name))
		}
		return(muts)
	}

	muts %<>% DummyCol('start')
	muts %<>% DummyCol('end')
	muts %<>% DummyCol('chrom')
	muts %<>% DummyCol('pos')
	muts %<>% DummyCol('cn')

	muts %>%
	unique %>%
	filter(!is.na(effect)) %>%
	# variant importance for plot overlay
	mutate(precedence=
		ifelse(effect=='Deletion',1,
		ifelse(effect=='Amplification',2,
		ifelse(effect=='Truncating SNV',3,
		ifelse(effect=='Frameshift In-Del',4,
		ifelse(effect=='Missense SNV',5,
		ifelse(effect=='Inframe In-Del',6,
		ifelse(effect=='Splice site variant',7,
		ifelse(effect=='Upstreop, or de novo modification',8,
		ifelse(effect=='Silent',9,
		NA)))))))))) %>%
	# remove genes with lower prescedence
	group_by(sample,gene) %>%
	arrange(gene,precedence) %>%
	slice(rank(precedence, ties.method="first")==1) %>%
	# number of variants per gene
	unique %>%
	group_by(gene) %>%
	mutate(n.gene=n()) %>%
	ungroup %>%
	# number of variants per band
	unique %>%
	group_by(band) %>%
	mutate(n.band=n()) %>%
	ungroup %>%
	# define plot gene order
	arrange(desc(n.gene),sample,desc(ccf),precedence,gene) %>%
	#p(50)
	mutate(order=row_number()) %>%
	ungroup %>%
	unique %>%
	# fill empty tiles
	group_by(pheno) %>%
	spread(gene,effect) %>%
	gather(gene,effect,-matches(paste(c('^n.gene$|^precedence$|^n.band$|^order',names(muts),'$'),collapse='$|^'))) %>%
	filter(!is.na(effect)) %>%
	ungroup %>%
	select(-precedence,-n.gene) %>%
	mutate(order=ifelse(is.na(effect),-Inf,order)) %>%
	unique %>%
	# fix plot ordering
	arrange(!is.na(effect),desc(order)) %>%
	mutate(gene=factor(gene,levels=filter(.,!is.na(effect)) %>% .$gene %>% unique)) %>%
	mutate(chrom.num=as.numeric(ifelse(chrom=='X'|chrom=='Y',23,chrom))) %>%
	arrange(n.band, sample, desc(chrom), desc(start)) %>%
	rowwise %>%
	mutate(band=str_c('chr',chrom,': ',band)) %>%
	ungroup %>%
	mutate(band=factor(band, levels=filter(.,!is.na(effect)) %>% .$band %>% unique)) %>%
	select(-order) %>%
	mutate(sample=factor(sample,levels=labels(dend))) %>%
	mutate(clonality=ifelse(clonality=='Clonal',clonality,NA))

}


# main plotting function
PlotVariants <- function(muts, file.name, event.recode=NULL, ccf=FALSE, loh=TRUE, cn=FALSE, width=20, height=20, text.size=18){

	muts %<>% rename(Sample=sample, Gene=gene, Effect=effect,CCF=ccf, Band=band, CN=cn)

	# plot aesthetic definitions
	palette  <- c(
		`Truncating SNV`='#C84DDD',
		`Frameshift In-Del`='#C17200',
		`Missense SNV`='#00A5A8',
		`Inframe indel`='#E44988',
		`Splice site variant`='#008AE9',
		`Upstream, start/stop, or de novo modification`='#749000',
		`Silent`='#666666',
		`Amplification`='#333399',
		`Deletion`='#e60000')

	geometry <- c(`LOH`=3)

	# main plot params
	if(cn!=TRUE) {
		hp <- ggplot(muts,aes(Sample,Gene))
	} else {
		hp <- ggplot(muts,aes(Sample,Band))
	}

	if(ccf==TRUE & cn!=TRUE ){  # CCF coloring
		hp <- hp + geom_tile(data=muts, aes(fill=CCF), colour='white') +
		scale_fill_gradient(low='#b3d9ff', high='#000066') +
		geom_tile(data=muts %>% filter(!is.na(clonality)), aes(color=clonality, stroke=1.5), size=2, fill=NA)
	} else {  # draw tiles and color
		hp <- hp +
		geom_tile(data=muts,aes(fill=Effect),colour='white') +
		scale_fill_manual(breaks=names(palette), values=palette, na.value='gray90')
	}

	# add loh as +
	if(loh==TRUE & cn!=TRUE) {
		hp <- hp +
		geom_point(data=muts, aes(shape=loh, stroke=1.5), size=2) +
		scale_shape_manual(values=geometry, guide=guide_legend(colour = 'white'))
	}

	hp <- hp +
	# specify legend
	guides(colour='white') +
	guides(colour = guide_legend(override.aes=list(alpha = 1,fill=NA))) + 

	# tile groups
	facet_wrap(~pheno, nrow=1, scales='free_x') +
	scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +

	# theme params
	theme(	legend.title		= element_blank(),
			panel.grid.major	= element_blank(),
			panel.grid.minor	= element_blank(),
			text				= element_text(size=18),
			axis.title.x		= element_blank(), #element_text(vjust=7),
			axis.title.y		= element_blank(), #element_text(vjust=1),
			axis.text.x			= element_text(angle=90, vjust=0.5, hjust=1),
			axis.text.y			= element_text(face='italic', hjust=1),
			axis.text			= element_text(size=text.size),
			axis.ticks.x		= element_blank(),
			axis.ticks.y		= element_blank(),
			legend.key			= element_rect(colour='white', fill=NULL, size=0.2),
			legend.key.size		= unit(1.4, 'lines'),
			legend.text			= element_text(size=text.size*1.2),
			strip.background	= element_rect(fill='white'),
			strip.text.x		= element_text(colour='white', size=text.size*1.2),
			panel.background	= element_rect(fill=NA, color=NA),
			panel.border		= element_rect(fill=NA, colour='black', size=2),
			plot.margin			= unit(c(0,0,0,0), 'in')	)

	# build Grob object
	hpg <- ggplotGrob(hp)

	# count number of samples in each group
	plot.lengths <- muts %>% split(.$pheno) %>% map(~ .x$Sample %>% unique %>% length) %>% unlist

	# get the column indexcorresponding to the panels.
	panelI <- hpg$layout$l[grepl('panel', hpg$layout$name)]

	# replace the default panel widths with relative heights.
	hpg$widths <- grid:::unit.list(hpg$widths)
	hpg$widths[panelI] <- lapply(plot.lengths, unit, 'null')

	# add extra width between panels
	for(gap in 1:(length(panelI)-1)){
		hpg$widths[panelI[gap]+1]=list(unit(0.8, 'cm'))
	}

	# change group facet colors
	names.grobs <- grid.ls(grid.force(hpg),print=FALSE)$name
	strips <- names.grobs[which(grepl('strip.background',names.grobs))]
	if(is.null(event.recode)) {
		cols <- colorRampPalette(brewer.pal(8,'Dark2'))(length(strips))  # default coloring
	} else { 
		cols <- event.recode %>% select(pheno,pheno.color) %>% unique %>% arrange(pheno) %>% .$pheno.color
	}

	for(strip in 1:length(strips)){
		hpg=editGrob(grid.force(hpg), gPath(strips[strip]), gp=gpar(fill=cols[strip]))
	}

	# draw plot
	pdf(file.name,width,height,bg='white')
		grid.draw(hpg)
	dev.off()
}


#----------------
# data processing
#----------------

# random color selection
set.seed(color.seed)
color.palette <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#color.palette <- colorRampPalette(c('#2d4be0', '#20e6ae', '#ccb625', '#969696'))(100)
color.palette %<>% sample(.,length(.))


# format muts with function above
muts %<>% FormatMuts


# remove silent mutations if present
if(include.silent!=TRUE){
	muts %<>% filter(effect!='silent')
}


# recode table values as colors
event.recode <-
	muts %>%
	select(-sample,-gene,-effect,-ccf,-loh) %>%
	rename(pheno.color=pheno) %>%
	mutate_each(funs(as.character)) %>%
	map(~ { seed = hash(.x)
			values = .x %>% ifelse(. %in% c('.', NA), '', .)
			colors = color.palette[(1:length(unique(values)))+seed]
			colors[unique(values) == ''] <- 'grey'
			setNames(colors, unique(values))[values]
		}) %>%
	bind_cols(muts[c('gene', 'sample', 'pheno')], .) %>%
	mutate(exists=1) %>%
	select(sample,gene,exists,pheno,pheno.color,everything())  # table expansion & rearrange for column naming


# reshape melted data into wide-form for distance calculation
event.matrix <-
	event.recode %>%
	select(sample,gene,exists) %>%
	spread(sample,exists,fill=0) %>%
	data.frame(.,row.names=1,stringsAsFactors=FALSE,check.names=FALSE) %>%
	as.matrix %>%
	t


# compute distance matrix using method specified above
dist <- DistExtra(event.matrix,dist.method)

# cluster using method specified above
hc <- dist %>% hclust(method=clust.method)

# construct dendrogram
dend <- 
	hc %>%
	as.dendrogram %>% 
	set('branches_lwd', 10)

dend %<>% reorder(dist) # reorder rows using distance matrix min

# rotate tree to order using method specified above
if(tree.sort == 'distance') {
	dend %<>% rotate_DendSer(ser_weight=dist(x))
} else {
	dend %<>% ladderize
}

# # adaptive branch pruning, order by ladderized tree layout
# clusters <-
#   cutreeDynamic(
#           hc,
#           minClusterSize = 3,
#           distM = as.matrix(dist.matrix),
#           method = 'hybrid',
#           deepSplit = 4,              # clustering sensitivity [1-4]
#           maxCoreScatter = NULL,      # max scatter of the core for a branch to be a cluster given as absolute heights [0-1]
#           minGap = NULL,              # min cluster gap given as fraction of the difference between ‘cutHeight’ and the 5th percentile of joining heights [0-1]
#           maxAbsCoreScatter = NULL,   # max scatter of the core for a branch to be a cluster given as absolute heights
#           minAbsGap = NULL            # min cluster gap given as absolute height difference
#       ) %>%
#   .[order.dendrogram(dend)]

# # cluster coloring
# clusters_numbers <- unique(clusters.o) - (0 %in% clusters.o)
# n_clusters <- length(clusters_numbers)

# cluster palette
#palette <- colorRampPalette(brewer.pal(8,'Dark2'))(n_clusters) 

























# if(tree.labels == FALSE){
#     dend %<>% set('labels', '')
# }

dist.dend <- reorder(dend,dist.matrix)

pdf('tree.test.ladder2.pdf',140,38)
	par(mar = c(65,40,10,20))   # bottom, left, top, right
	layout(matrix(c(1,2),nrow=1), widths=c(10,1))
	plot(dend, cex.axis=6)
	colored_bars( colors = pheno,
				  dend = dend,
				  sort_by_labels_order = FALSE,
				  add = TRUE,
				  rowLabels = pheno,
				  y_scale = 10,
				  cex.rowLabels = 5.8 )
	legend.num = 0
	event.recode %>% select(-sample, -gene, -exists) %>% map(~ {
		legend.num <<- legend.num + 1
		print(legend.num)
		legend(x=3, y=legend.num*5.5, legend=unique(names(.x)), fill=unique(.x), cex=4)
	})
dev.off()

pdf('tree.charlotte.test.pdf',60,25)
	heatmap.2(
		event.matrix,
		trace='none',
		Rowv=FALSE,
		hclustfun=function(x){hclust(x, 'ward.D2')},
		col=c('white','black'),
		reorderfun=function(d,w) rev(reorder(d,w)),
		distfun=function(x) as.dist(Hamming(x)),
		key=FALSE
			)
dev.off()

pdf('heat.all.test2.pdf',60,25)
	heatmap.2(
		t(event.matrix),
		trace='none',
		dendrogram='column',
		Colv=rev(dend),
		col=c('white','black'),
		key=FALSE
	)
dev.off()






#----------------
# TREE GENERATION
#----------------

plot_file          = 'analysis/tree.pdf'
#colnames(muts)[1]  = 'Sample.ID'
colnames(muts)[7]  = 'Gene'
colnames(muts)[8]  = 'AA'
colnames(muts)[10] = 'Effect'
colnames(muts)[26] = 'CCF'
muts               = muts[!is.na(muts$ccf),]
tumor              = unique(muts$sample)


muts %<>% as.data.frame

sample_names   = as.list(sort(unique(muts$sample)))
mutation_genes = unique(muts$Gene)
rownames(muts) = 1:nrow(muts)
TCGA=FALSE



library(phangorn)

sample_names<-muts$sample %>% unique
mutation_genes<-muts$gene %>% unique

# Make a matrix of "blank" values the with nrow= #mutations and ncol=#samples
mutation_heatmap <- matrix(0, nrow=sum(unlist(lapply(sample_names, length))), ncol=sum(unlist(lapply(mutation_genes, length))))
rownames(mutation_heatmap) <- unlist(sample_names)
colnames(mutation_heatmap) <- unlist(mutation_genes)


# Make sure the sample and mutations are both in the list of gene mutations and gene samples
if (!TCGA) { smallmaf <- muts[which(muts$gene %in% unlist(mutation_genes) & muts$sample %in% unlist(sample_names)),]
}else { 
	muts$id <- unlist(lapply(muts$sample, function(x){substr(x, 1, 12)}))
	print(head(muts$id))
	print(head(unlist(sample_names)))
	smallmaf <- muts[which(muts$Hugo_Symbol %in% unlist(mutation_genes) & muts$id %in% unlist(sample_names)),] 
}

# For each row read the Effect and create the type based on which category it fits in
for (i in 1:nrow(smallmaf)) {
	if(!TCGA) { type = smallmaf$effect[i] } else { type = smallmaf$effect[i] }
	print(paste(i,type,sep="_"))

	if (!TCGA) { 
	if (mutation_heatmap[which(rownames(mutation_heatmap)==smallmaf$sample[i],), which(colnames(mutation_heatmap)==smallmaf$gene[i])] < type) {
		mutation_heatmap[which(rownames(mutation_heatmap)==smallmaf$sample[i],), which(colnames(mutation_heatmap)==smallmaf$gene[i])] <- type}
	} else { mutation_heatmap[which(rownames(mutation_heatmap)==smallmaf$id[i]), which(colnames(mutation_heatmap)==smallmaf$Hugo[i])] <- type }
}

smalltab = t(mutation_heatmap)
smalltab = smalltab>0

smalltab = cbind(smalltab, FALSE)
colnames(smalltab)[ncol(smalltab)] = "Parental"

smalltab = cbind(smalltab, FALSE)
colnames(smalltab)[ncol(smalltab)] = "Test"

#smalltab2 <- smalltab %>% as_data_frame %>% mutate(Parental=ifelse(LPSDD2==TRUE & LPSWD==TRUE,TRUE,FALSE)) %>% as.matrix
#rownames(smalltab2) <- row.names(smalltab)

pd <- phyDat(t(smalltab), type="USER", levels=c(FALSE, TRUE))
dm <- dist.hamming(pd)
tree <- njs(dm)
treeRatchet <- pratchet(pd, start=tree)
treeRatchet <- acctran(treeRatchet, pd)
rt <- root(treeRatchet, "Parental")

pdf("analysis/phylo_tree.pdf")
	plot(rt)
dev.off()





muts <- read.delim(mut.file ,sep='\t', stringsAsFactors=FALSE) %>% tbl_df



