# read & format mutations
#muts.vcf <- read.delim("recurrent_mutations/sufam/all_sufam.txt",stringsAsFactors=FALSE,sep="\t") %>%
muts.vcf <- read.delim("sufam/all_validated_sufam.txt",stringsAsFactors=FALSE,sep="\t") %>%
	select(sample.name=sample,chrom,pos,alt=val_alt,cov,maf=val_maf) %>%
	tbl_df

#muts.suf <- read.delim("recurrent_mutations/sufam/all_mutations.vcf",stringsAsFactors=FALSE,sep="\t") %>%
muts.suf <- read.delim("sufam/all_mutations.vcf",stringsAsFactors=FALSE,sep="\t") %>%
	select(chrom=`X.CHROM`,pos=POS,gene=`ANN....GENE`,alt=ALT,effect=`ANN....EFFECT`) %>%
	tbl_df

muts <- muts.vcf %>% full_join(muts.suf, by=c("chrom","pos","alt")) %>%
	rowwise() %>%
	mutate(gene=str_split(gene,"\\|") %>% unlist %>% head(1)) %>%
	mutate(effect=str_split(effect,"\\|") %>% unlist %>% tail(n=1)) %>%
	ungroup() %>%
	mutate(effect=
		ifelse(effect%in%c("STOP_GAINED","Nonsense_Mutation","stop_gained&splice_region_variant","stop_gained"),	"truncating snv",
		ifelse(effect%in%c("FRAME_SHIFT","FRAME_SHIFT","Frame_Shift_Del","Frame_Shift_Ins","frameshift_variant","frameshift_variant&stop_gained","frameshift_variant&splice_region_variant","frameshift_variant&splice_acceptor_variant&splice_region_variant&splice_region_variant&intron_variant"),	"frameshift indel",
		ifelse(effect%in%c("NON_SYNONYMOUS_CODING","STOP_LOST","Missense_Mutation","missense_variant","missense_variant&splice_region_variant","missense_variant|missense_variant"),"missense snv",
		ifelse(effect%in%c("CODON_CHANGE_PLUS_CODON_DELETION","CODON_DELETION","CODON_INSERTION","In_Frame_Ins","In_Frame_Del","disruptive_inframe_deletion","disruptive_inframe_insertion","inframe_deletion","inframe_insertion","disruptive_inframe_deletion&splice_region_variant","inframe_deletion&splice_region_variant"),	"inframe indel",
		ifelse(effect%in%c("SPLICE_SITE_DONOR","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_REGION","Splice_Site","splice_donor_variant&intron_variant","splice_acceptor_variant&intron_variant","splicing","splice_donor_variant&splice_region_variant&intron_variant","splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&splice_region_variant&intron_variant","splice_region_variant&intron_variant","frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant"),	"splice site variant",
		ifelse(effect%in%c("STOP_LOST","START_LOST","START_GAINED","UTR_5_PRIME","start_lost","stop_lost"),		"start/stop codon change",
		ifelse(effect%in%c("synonymous_variant","splice_region_variant&synonymous_variant","non_coding_exon_variant","upstream_gene_variant","downstream_gene_variant","intron_variant","frameshift_variant&splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant","non_coding_exon_variant|synonymous_variant","SYNONYMOUS_CODING","synonymous_variant|synonymous_variant","splice_region_variant&synonymous_variant|splice_region_variant&non_coding_exon_variant","intragenic_variant","intergenic_region","3_prime_UTR_variant","5_prime_UTR_premature_start_codon_gain_variant","5_prime_UTR_variant","intergenic_region"),		"silent", # synonymous/noncoding/up/downstream/intragenic
		NA)))))))) %>%
	distinct %>%
	select(-c(alt)) %>%
	mutate(chrom=ifelse(chrom=="X",23,ifelse(chrom=="Y",23,chrom))) %>%
	mutate(chrom=as.numeric(chrom))