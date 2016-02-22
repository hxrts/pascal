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
}

#-------------------
# load sufam results
#-------------------

muts.txt <-
	read_tsv("recurrent_mutations/sufam/all_sufam.txt") %>%
	select(sample,chrom,pos,ref=val_ref,alt=val_alt,cov,maf=val_maf)

muts.vcf <-
	read_tsv("recurrent_mutations/sufam/all_mutations.vcf") %>%
	select(chrom=`#CHROM`,pos=POS,gene=`ANN[*].GENE`,alt=ALT,effect=`ANN[*].EFFECT`)

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

#---------------
# print metadata
#---------------

cat(green("\n-") %+% " muts.txt | muts.vcf | muts.all | muts\n\n") ; muts ; cat("\n")

#--------------------------------
# move back to original directory
#--------------------------------

setwd(wd)
