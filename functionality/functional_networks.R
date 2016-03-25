#!/usr/bin/env Rscript

# create necessary directories
system("mkdir functional_networks &>/dev/null")

#-----
# init
#-----

#suppressMessages(pacman::p_load(GO.db))

# read metadata files & load base packages
source("pascal/lib/init.R")
source("pascal/lib/muts.R")


# http://software.broadinstitute.org/gsea/msigdb/collections.jsp [retreived 3.1.16]
# db.names <- as.list(c(
# 					hallmark    = "/ifs/e63data/reis-filho/reference/MSigDB/v5.1/gmt/h.all.v5.1.symbols.gmt",		# hallmark
# 					curated     = "/ifs/e63data/reis-filho/reference/MSigDB/v5.1/gmt/c2.all.v5.1.symbols.gmt",		# curated
# 					mined       = "/ifs/e63data/reis-filho/reference/MSigDB/v5.1/gmt/c4.all.v5.1.symbols.gmt",		# computationally defined via cancer gene micro array
# 					ontologic   = "/ifs/e63data/reis-filho/reference/MSigDB/v5.1/gmt/c5.all.v5.1.symbols.gmt",		# ontollogically grouped (not necessarily coexpressed)
# 					oncogenic	= "/ifs/e63data/reis-filho/reference/MSigDB/v5.1/gmt/c6.all.v5.1.symbols.gmt",		# oncogenic signatures
# 					immunologic = "/ifs/e63data/reis-filho/reference/MSigDB/v5.1/gmt/c7.all.v5.1.symbols.gmt"		# immunologic signatures
# 				))

# # read databases from file
# dbs <-
# 	db.names %>%
# 	map(~ {
# 			db <- read_tsv(.x,col_names=FALSE)
# 			db
# 		})

# db.tsv.names <- as.list(c(
# 					hallmark    = "/ifs/e63data/reis-filho/reference/MSigDB/v5.1/tsv/MSigDB.h-all.v5.1.tsv",		# hallmark
# 					curated     = "/ifs/e63data/reis-filho/reference/MSigDB/v5.1/tsv/MSigDB.c2-all.v5.1.tsv",		# curated
# 					mined       = "/ifs/e63data/reis-filho/reference/MSigDB/v5.1/tsv/MSigDB.c4-all.v5.1.tsv",		# computationally defined via cancer gene micro array
# 					ontologic   = "/ifs/e63data/reis-filho/reference/MSigDB/v5.1/tsv/MSigDB.c5-all.v5.1.tsv",		# ontollogically grouped (not necessarily coexpressed)
# 					oncogenic	= "/ifs/e63data/reis-filho/reference/MSigDB/v5.1/tsv/MSigDB.c6-all.v5.1.tsv",		# oncogenic signatures
# 					immunologic = "/ifs/e63data/reis-filho/reference/MSigDB/v5.1/tsv/MSigDB.c7-all.v5.1.tsv"		# immunologic signatures
# 				))

# # read databases from file
# dbs <-
# 	db.tsv.names %>%
# 	map( ~ {
# 			read_tsv(.x) %>%
# 			gather(gene,exists,-paths) %>%
# 			filter(exists==TRUE) %>%
# 			arrange(gene) %>%
# 			select(-exists)
# 		})

# # save database binary
# save(data=dbs,file="/ifs/e63data/reis-filho/reference/MSigDB/v5.1/dbs.RData")

n.genes <- muts.matched.filter[[1]] %>% select(gene,chrom,pos,alt,effect) %>% distinct %>% nrow
list.name <-"small_cell"

# load database binary
load("/ifs/e63data/reis-filho/reference/MSigDB/v5.1/dbs.RData",verbose=TRUE)

meta.db <-
	dbs %>%
	bind_rows(.id="db") %>%
	mutate(total.n=n()) %>%
	group_by(paths) %>%
	mutate(base.n=n()) %>%
	ungroup %>%
	inner_join(muts.matched.filter[[1]] %>% select(gene,tumor,chrom,pos,alt,effect)) %>%
	distinct %>%
	group_by(paths) %>%
	mutate(match.n = gene %>% unique %>% length) %>%
	mutate(p.hyper = phyper(match.n, base.n, total.n-base.n, n.genes, lower.tail=FALSE)) %>%
	ungroup %>%
	arrange(p.hyper)

directed.map <-
	meta.db %>%
	select(paths,gene) %>%
	distinct %>%
	group_by(paths) %>%
	filter(n()>1) %>%
	ungroup %>%
	split(.$paths) %>%
	map(~ expand.grid(source=.x$gene,target=.x$gene,stringsAsFactors=FALSE)) %>%
	bind_rows(.id="pathway") %>%
	select(source,target) %>%
	mutate(value=as.integer(1))

graph.nodes <-
	meta.db %>% select(name=gene,group=paths,size=1) %>% distinct


bar <- renderForceNetwork({
forceNetwork(Links = directed.map, Nodes = graph.nodes,
            Source = "source", Target = "target",
            Value = "value", NodeID = "name",
            Group = "group", opacity = input$opacity)
})


foo <- forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name",
              Group = "group", opacity = 0.4,
              colourScale = "d3.scale.category20b()")


simpleNetwork(networkData) %>%
saveNetwork(file = 'Net1.html')


meta.db %>% mutate(group.index=group_indices_(.,.dots="paths")) %>% glimpse


# gene-to-gene visualization
# http://cpdb.molgen.mpg.de/

# pathway analysis post
# http://www.gettinggeneticsdone.com/2012/03/pathway-analysis-for-high-throughput.html

# gene set enrichment
# http://software.broadinstitute.org/gsea/index.jsp

# Database for Annotation, Visualization and Integrated Discovery
# https://david.ncifcrf.gov/

# Reactome Pathway browser
# http://www.reactome.org/


#---------------------------------------------------


# IMPACT genes

impact.clinical <- read.delim("/home/limr/share/cbioportal/data-repos/msk-impact/msk-impact/data_clinical.txt",stringsAsFactors=FALSE,sep="\t",skip=5) %>%
tbl_df

impact.muts <- read.delim("/home/limr/share/cbioportal/data-repos/msk-impact/msk-impact/data_mutations_extended.txt",stringsAsFactors=FALSE,sep="\t",skip=5) %>%
tbl_df

impact.cna <- read.delim("/home/limr/share/cbioportal/data-repos/msk-impact/msk-impact/data_CNA.txt",stringsAsFactors=FALSE,sep="\t",skip=5) %>%
tbl_df


#genes.table <- read_tsv("/ifs/e63data/reis-filho/data/myoep/cell_lines_rnaseq/recurrent_pathways/Differentially_Expressed_Genes_AME_cell_lines.tsv")[,1:14]

# lapply(1:14, function(num){

# 	genes <- genes.table[,num] %>% unlist %>% list.filter(!is.na(.))
# 	list.name <- names(genes.table)[num] %>% strsplit(" ") %>% unlist %>% paste(sep="_",collapse="_")

# 	lapply(db.names,function(db.name) GetPathway(db.name,genes,list.name))

# })


#---------------------------------------------------

library(paxtoolsr)

pc.id.map <-
	muts.all$gene %>%
	unique %>%
	sort %>%
	data.frame(gene=.,stringsAsFactors=FALSE) %>%
	rowwise %>%
	mutate(pc.id=idMapping(gene,verbose=TRUE) %>% ifelse(is.null(.),"",.) %>% unlist)

pc.id.map %>%
	filter(pc.id!="")

muts %>%
left_join(pc.id.map,by="gene") %>%


pc.id.map$pc.id %>%
lapply(. %>%
	traverse(path = "ProteinReference/entityFeature:ModificationFeature") %>%
	xpathSApply("//value/text()", xmlValue)
	)













































