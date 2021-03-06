#!/usr/bin/env Rscript

#---------------
# base libraries
#---------------

suppressMessages(pacman::p_load(dplyr,readr,tidyr,magrittr,purrr,stringr,rlist,ggplot2,crayon))

#-----------------
# path contingency
#-----------------

wd <- getwd()

if(!"samples.txt" %in% list.files()){
	setwd("..")
	if(!"samples.txt" %in% list.files()){
		setwd("..")
		if(!"samples.txt" %in% list.files()){
			setwd("..")
		}
	}
}

#-----------------
# load raw results
#-----------------

muts.raw <- read.delim("summary/tsv/mutation_summary.tsv",sep="\t",stringsAsFactors=FALSE) %>%
	tbl_df %>%
	select(tumor=`TUMOR_SAMPLE`,normal=`NORMAL_SAMPLE`,chrom=`CHROM`,pos=`POS`,ref=`REF`,alt=`ALT`,gene=`ANN....GENE`,effect=`ANN....EFFECT`,tumor.maf=`TUMOR_MAF`,normal.maf=`NORMAL_MAF`,tumor.dp=`TUMOR.DP`,normal.dp=`NORMAL.DP`)

#-------------------
# load sufam results
#-------------------

muts.txt <-
	read.delim("recurrent_mutations/sufam/all_sufam.txt",sep="\t",stringsAsFactors=FALSE) %>%
	tbl_df %>%
	select(sample,chrom,pos,ref=val_ref,alt=val_alt,cov,maf=val_maf)

muts.vcf <-
	read.delim("recurrent_mutations/sufam/all_mutations.vcf",sep="\t",stringsAsFactors=FALSE) %>%
	tbl_df %>%
	select(chrom=`X.CHROM`,pos=POS,gene=`ANN....GENE`,alt=ALT,effect=`ANN....EFFECT`)

#---------------
# join vcf & txt
#---------------

muts.all <-
	muts.vcf %>% full_join(muts.txt, by=c("chrom","pos","alt")) %>%
	rowwise() %>%
	mutate(gene=str_split(gene,"\\|") %>% unlist %>% head(1)) %>%
	mutate(effect=str_split(effect,"\\|") %>% unlist %>% head(1)) %>%
	ungroup() %>%
	select(sample,gene,chrom,pos,ref,alt,effect,cov,maf)

#--------------------
# pair tumor / normal
#--------------------

muts.tn <-
	muts.all %>%
	filter(sample %in% samples$tumor) %>%
	rename(tumor=sample,cov.t=cov,maf.t=maf) %>%
	left_join(samples,by=c("tumor")) %>%
	left_join(muts.all %>%
		filter(sample %in% samples$normal) %>%
		select(-ref) %>%
		rename(normal=sample,cov.n=cov,maf.n=maf),
	by=c("normal","gene","chrom","pos","alt","effect")) %>%
	select(tumor,normal,gene,chrom,pos,ref,alt,cov.t,cov.n,maf.t,maf.n,effect) %>%
	mutate(effect=
	ifelse(effect%in%c('STOP_GAINED','Nonsense_Mutation','stop_gained&splice_region_variant','stop_gained','Nonsense_Mutation','Stop_Codon_Ins','nonsense'),'truncating snv',
	ifelse(effect%in%c('FRAME_SHIFT','FRAME_SHIFT','Frame_Shift_Del','Frame_Shift_Ins','frameshift_variant','frameshift_variant&stop_gained','frameshift_variant&splice_region_variant','frameshift_variant&splice_acceptor_variant&splice_region_variant&splice_region_variant&intron_variant','Frame_Shift_Del','Frame_Shift_Ins','frame_shift_del','frame_shift_ins'),'frameshift indel',
	ifelse(effect%in%c('NON_SYNONYMOUS_CODING','STOP_LOST','Missense_Mutation','missense_variant','missense_variant&splice_region_variant','missense_variant|missense_variant','Missense_Mutation','missense'),'missense snv',
	ifelse(effect%in%c('CODON_CHANGE_PLUS_CODON_DELETION','CODON_DELETION','CODON_INSERTION','In_Frame_Ins','In_Frame_Del','disruptive_inframe_deletion','disruptive_inframe_insertion','inframe_deletion','inframe_insertion','disruptive_inframe_deletion&splice_region_variant','inframe_deletion&splice_region_variant','In_Frame_Del','In_Frame_Ins','in_frame_del','in_frame_ins'),'inframe indel',
	ifelse(effect%in%c('SPLICE_SITE_DONOR','SPLICE_SITE_ACCEPTOR','SPLICE_SITE_REGION','Splice_Site','splice_donor_variant&intron_variant','splice_acceptor_variant&intron_variant','splicing','splice_donor_variant&splice_region_variant&intron_variant','splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&splice_region_variant&intron_variant','Splice_Site','splice'),'splice site variant',
	ifelse(effect%in%c('STOP_LOST','START_LOST','START_GAINED','UTR_5_PRIME','start_lost','stop_lost',"5'UTR","5'Flank",'De_novo_Start_InFrame','De_novo_Start_OutOfFrame','Stop_Codon_Del','Start_Codon_SNP','Start_Codon_Ins','Start_Codon_Del','Nonstop_Mutation','nonstop'),'upstream, start/stop, or de novo modification',
	ifelse(effect%in%c('synonymous_variant','splice_region_variant&synonymous_variant','non_coding_exon_variant','upstream_gene_variant','downstream_gene_variant','intron_variant','frameshift_variant&splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant','non_coding_exon_variant|synonymous_variant','SYNONYMOUS_CODING','synonymous_variant|synonymous_variant','splice_region_variant&synonymous_variant|splice_region_variant&non_coding_exon_variant','intragenic_variant',"3'UTR",'IGR','lincRNA','RNA','Intron','silent','intron_exon'),'silent', # synonymous/noncoding/up/downstream/intragenic
	NA))))))))

muts.tn.filter <-
	muts.tn %>%
	filter(effect!="silent") %>%
	filter(maf.t>0.005 & cov.t>1 & cov.n>1 & maf.n==0) %>%
	arrange(-maf.t) 

muts.matched <-
	subsets %>%
	map(~
		muts.tn %>%
		filter(tumor %in% .x) %>%
		select(tumor,gene,chrom,pos,alt,maf.n) %>%
		spread(tumor,maf.n) %>%
		mutate(row.sum=select(.,matches(subsets %>% paste(collapse="|"))) %>% rowSums) %>%
		filter(row.sum==0) %>%
		select(-row.sum) %>%
		gather(tumor,maf.n,-c(gene,chrom,pos,alt)) %>%
		left_join(muts.tn,by=c("gene","chrom","pos","alt","tumor","maf.n")) %>%
		arrange(chrom,pos,gene,alt,tumor)
	) %>%
	setNames(subsets %>% names)

muts.matched.filter <-
	muts.matched %>%
	map(~ .x %>%
		filter(effect!="silent") %>%
		filter(maf.t>0.005 & cov.t>1 & cov.n>1 & maf.n==0) %>%
		arrange(-maf.t)
	)

#---------------
# print metadata
#---------------

cat(green("\n-muts.txt\n")) ; 				muts.txt %>% print(n=0)
cat(green("\n-muts.vcf\n")) ; 				muts.vcf %>% print(n=0)
cat(green("\n-muts.all\n")) ; 				muts.all %>% print(n=0)
cat(green("\n-muts.tn\n")) ; 				muts.tn %>% print(n=0)
cat(green("\n-muts.tn.filter\n")) ; 		muts.tn.filter %>% print(n=0)
cat(green("\n-muts.matched\n")) ; 			muts.matched %>% map(~ .x %>% print(n=0))
cat(green("\n-muts.matched.filter\n")) ;	muts.matched.filter %>% map(~ .x %>% print(n=0))


#--------------------------------
# move back to original directory
#--------------------------------

setwd(wd)
