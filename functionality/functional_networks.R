#!/usr/bin/env Rscript

# create necessary directories
system("mkdir functional_networks &>/dev/null")

#-----
# init
#-----

# run-specific libraries
library(plyr)
#suppressMessages(pacman::p_load(GO.db))

# read metadata files & load base packages
source("pascal/lib/init.R")

# mutation plotting function
PlotVariants <- function(muts,filename){
	palette  <- c(`truncating snv`="#00A5A8",`frameshift indel`="#008AE9",`missense snv`="#C84DDD",`inframe indel`="#E44988",`splice site variant`="#C17200",`start/stop codon change`="#749000") #,`silent`="#00A24B")
	geometry <- c(`truncating snv`=20,`frameshift indel`=2,`missense snv`=3,`inframe indel`=4,`splice site variant`=5,`start/stop codon change`=6) #,`silent`=1)
	pdf(filename,12,16)
	print(
	ggplot(muts,aes(sample,gene)) +
	geom_tile(aes(fill=ifelse(is.na(variant),variant,NA)),colour="black") +

	geom_tile(data=muts %>% filter(variant%in%c("truncating snv","frameshift indel","splice site variant","stop codon change")),aes(fill=variant),colour="black") +
	scale_fill_manual(values=palette,na.value="gray90") +

	geom_point(data=muts %>% filter(variant%in%c("missense snv","inframe indel","silent")),aes(shape=variant),color="black") +
	scale_shape_manual(values=geometry,na.value=NA) +

	theme(
		panel.border=element_blank(),
		axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
		axis.text=element_text(size=10),
		legend.title=element_blank())
	)
	dev.off()
}

#-------------------
# load sufam results
#-------------------

muts.txt <- read_tsv("recurrent_mutations/sufam/all_sufam.txt") %>%
	select(sample,chrom,pos,ref=val_ref,alt=val_alt,cov,maf=val_maf)

muts.vcf <- read_tsv("recurrent_mutations/sufam/all_mutations.vcf") %>%
	select(chrom=`#CHROM`,pos=POS,gene=`ANN[*].GENE`,alt=ALT,effect=`ANN[*].EFFECT`)

#---------------
# join vcf & txt
#---------------

muts.all <- muts.vcf %>% full_join(muts.txt, by=c("chrom","pos","alt")) %>%
	rowwise() %>%
	mutate(gene=str_split(gene,"\\|") %>% unlist %>% head(1)) %>%
	mutate(effect=str_split(effect,"\\|") %>% unlist %>% head(1)) %>%
	ungroup() %>%
	select(sample,gene,chrom,pos,ref,alt,effect,cov,maf)

#--------------------
# pair tumor / normal
#--------------------

muts<-
	muts.all %>%
	filter(sample %in% samples$tumor) %>%
	rename(tumor=sample,cov.t=cov,maf.t=maf) %>%
	left_join(samples,by=c("tumor")) %>%
	left_join(muts %>%
		filter(sample %in% samples$normal) %>%
		select(-ref) %>%
		rename(normal=sample,cov.n=cov,maf.n=maf),
	by=c("normal","gene","chrom","pos","alt","effect")) %>%
	select(tumor,normal,gene,chrom,pos,ref,alt,cov.t,cov.n,maf.t,maf.n,effect) %>%
	filter(maf.t>0.01 & cov.t>5 & cov.n>5 & maf.n==0) %>%
	mutate(effect=
		ifelse(effect%in%c("STOP_GAINED","Nonsense_Mutation","stop_gained&splice_region_variant","stop_gained"),"truncating snv",
		ifelse(effect%in%c("FRAME_SHIFT","FRAME_SHIFT","Frame_Shift_Del","Frame_Shift_Ins","frameshift_variant","frameshift_variant&stop_gained","frameshift_variant&splice_region_variant","frameshift_variant&splice_acceptor_variant&splice_region_variant&splice_region_variant&intron_variant"),"frameshift indel",
		ifelse(effect%in%c("NON_SYNONYMOUS_CODING","STOP_LOST","Missense_Mutation","missense_variant","missense_variant&splice_region_variant","missense_variant|missense_variant"),"missense snv",
		ifelse(effect%in%c("CODON_CHANGE_PLUS_CODON_DELETION","CODON_DELETION","CODON_INSERTION","In_Frame_Ins","In_Frame_Del","disruptive_inframe_deletion","disruptive_inframe_insertion","inframe_deletion","inframe_insertion","disruptive_inframe_deletion&splice_region_variant","inframe_deletion&splice_region_variant"),"inframe indel",
		ifelse(effect%in%c("SPLICE_SITE_DONOR","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_REGION","Splice_Site","splice_donor_variant&intron_variant","splice_acceptor_variant&intron_variant","splicing","splice_donor_variant&splice_region_variant&intron_variant","splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&splice_region_variant&intron_variant"),"splice site variant",
		ifelse(effect%in%c("STOP_LOST","START_LOST","START_GAINED","UTR_5_PRIME","start_lost","stop_lost"),"start/stop codon change",
		#ifelse(effect%in%c("Amplification","Homozygous Deletion"),X #"CNA",
		ifelse(effect%in%c("synonymous_variant","splice_region_variant&synonymous_variant","non_coding_exon_variant","upstream_gene_variant","downstream_gene_variant","intron_variant","frameshift_variant&splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant","non_coding_exon_variant|synonymous_variant","SYNONYMOUS_CODING","synonymous_variant|synonymous_variant","splice_region_variant&synonymous_variant|splice_region_variant&non_coding_exon_variant","intragenic_variant"),"silent", # synonymous/noncoding/up/downstream/intragenic
		NA)))))))) %>%
	filter(effect!="silent") %>%
	arrange(-maf.t)

# muts %>% 
# select(tumor,gene,effect) %>%
# arrange(tumor) %>%
# spread(tumor,effect,fill=NA) %>%
# select(gene,samples$tumor %>% sort %>% .[c(1,length(.))] ) %>%
# distinct

# # melt MSK mutations & plot
# msk.list <- msk.muts %>%
# 	gather(sample,variant,ACC10A:MP3MET) %>%
# 	arrange(sample,gene,variant)

# # plot MSK mutations
# PlotVariants(msk.list,"msk.muts.suf.pdf")




GetPathway <- function(db.name,genes,list.name) {

	db <-  suppressWarnings(read_tsv(db.name,col_names=FALSE))
	db.id <- db.name %>% strsplit("/") %>% unlist %>% tail(1) %>% substr(1,nchar(.)-4)
	matches <- db %>% apply(1, . %>% list.filter(. %in% genes) %>% unname)
	file <- str_c("functional_networks/",list.name,"_pathways.tsv")
	
	cat(list.name %+% " -- " %+% db.id %+% "\n")

	sink(file,append=TRUE)
		cat("\n# " %+% db.id %+% "\n")
	sink()

	if(matches %>% length > 0){
		path <-
			matches %>%
			setNames(db$X1) %>%
			list.filter(length(.)>0) %>%
			list.sort(-length(.))

		path.table <- ldply(path,rbind) %>% filter(!is.na(`2`))

		write.table(path.table,file,quote=FALSE,row.names=FALSE,na="",col.names=FALSE,sep="\t",append=TRUE)

	}else{

		sink(file,append=TRUE)
			cat("NO MATCHES\n")
		sink()
	}

}

impact <- read_tsv("/home/limr/share/cbioportal/data-repos/msk-impact/msk-impact")
db.names <- list.files("/ifs/e63data/reis-filho/reference/MSigDB/v5.1/gmt",full.names=TRUE,pattern=".gmt")

genes <- muts$gene %>% unique %>% sort
list.name="small_cell"
lapply(db.names,function(db.name) GetPathway(db.name,genes,list.name))



#genes.table <- read_tsv("/ifs/e63data/reis-filho/data/myoep/cell_lines_rnaseq/recurrent_pathways/Differentially_Expressed_Genes_AME_cell_lines.tsv")[,1:14]

# lapply(1:14, function(num){

# 	genes <- genes.table[,num] %>% unlist %>% list.filter(!is.na(.))
# 	list.name <- names(genes.table)[num] %>% strsplit(" ") %>% unlist %>% paste(sep="_",collapse="_")

# 	lapply(db.names,function(db.name) GetPathway(db.name,genes,list.name))

# })











