#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

suppressMessages(pacman::p_load(yaml,dplyr,readr,tidyr,magrittr,purrr,stringr,rlist,crayon))

sysprefix="umask 002 && unset PYTHONPATH && source /home/bermans/miniconda2/envs/pyclone/bin/activate /home/bermans/miniconda2/envs/pyclone >/dev/null 2>&1 && "

# create necessary directories
system("mkdir pyclone pyclone/config pyclone/events pyclone/priors pyclone/tables pyclone/plots &>/dev/null")

#------
# input
#------

cat(green("\n-") %+% " reading input\n")

patients <- read.delim("sample_sets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE) %>%
	setNames(c("patient",1:(ncol(.)-1))) %>%
	group_by(patient) %>%
	filter(!grepl("#",patient)) %>%
	summarise_each(funs(ifelse(.=="","NA",.)))

subsets <- read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE) %>%
	setNames(c("subset",1:(ncol(.)-1))) %>%
	group_by(subset) %>%
	filter(!grepl("#",subset)) %>%
	summarise_each(funs(ifelse(.=="","NA",.)))
		
sample.t <- subsets %>%
	select(-subset) %>%
	unlist %>%
	list.filter(.!="NA") %>%
	sort

sample.n <- sample.t %>%
	lapply(.,function(patient)
		patients[which(apply(patients,1,function(sample)contains(sample,patient))),] %>%
		list.filter(.!="NA") %>%
		unlist %>%
		tail(n=1) 
	) %>%
	unlist

sample.tn <- data.frame(normal=sample.n,tumor=sample.t)


# read & format mutations
muts.vcf <- read.delim("recurrent_mutations/sufam/all_sufam.txt",stringsAsFactors=FALSE,sep="\t") %>%
	select(sample.name=sample,chrom,pos,alt=val_alt,cov,maf=val_maf) %>%
	tbl_df

muts.suf <- read.delim("recurrent_mutations/sufam/all_mutations.vcf",stringsAsFactors=FALSE,sep="\t") %>%
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
		#ifelse(effect%in%c("Amplification","Homozygous Deletion"),X #"CNA",
		ifelse(effect%in%c("synonymous_variant","splice_region_variant&synonymous_variant","non_coding_exon_variant","upstream_gene_variant","downstream_gene_variant","intron_variant","frameshift_variant&splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant","non_coding_exon_variant|synonymous_variant","SYNONYMOUS_CODING","synonymous_variant|synonymous_variant","splice_region_variant&synonymous_variant|splice_region_variant&non_coding_exon_variant","intragenic_variant","intergenic_region","3_prime_UTR_variant","5_prime_UTR_premature_start_codon_gain_variant","5_prime_UTR_variant","intergenic_region"),		"silent", # synonymous/noncoding/up/downstream/intragenic
		NA)))))))) %>%
	distinct %>%
	select(-c(alt)) %>%
	mutate(chrom=ifelse(chrom=="X",23,ifelse(chrom=="Y",23,chrom))) %>%
	mutate(chrom=as.numeric(chrom))


segfiles<-list.files("facets",pattern="*cncf.txt")

setwd("pyclone")	# for some reason PyClone needs to be run from the root directory

#---------------------------------
# variables for yaml configuration
#---------------------------------

num_iters <- as.integer(50000)
base_measure_params <- list(alpha=as.integer(1),beta=as.integer(1))
concentration <- list(value=as.integer(1),prior=list(shape=1.0,rate=0.001))
density <- "pyclone_beta_binomial"
beta_binomial_precision_params <- list(value=as.integer(1000),prior=list(shape=1.0,rate=0.0001),proposal=list(precision=0.01))
working_dir <- getwd()

#------------------------------
# main loop over sample subsets
#------------------------------

for (subnum in 1:nrow(subsets)){

	# input processing
	line <- as.vector(subsets[subnum,])
	subsamples <- line[line!="NA"][-1]
	subname <- line[[1]][1]

	cat(blue("\n--------------------------------\n  PYCLONE beginning subset ",subname,"\n--------------------------------\n",sep=""))

	system(str_c("mkdir ",subname," &>/dev/null"))

	#----------------------------
	# run-specific yaml variables
	#----------------------------

	samples <- lapply(subsamples,function(sample)
				list(
					mutations_file=str_c("priors/",sample,".priors.yaml"),
					tumour_content=list(value=1.0),
					error_rate=0.001)
				) %>% setNames(subsamples)


	#----------------
	# write yaml file
	#----------------

	cat(green("\n-") %+% " building configuration file:\n  config/",subname,".config.yaml\n",sep="")
	sink(file=str_c("config/",subname,".config.yaml"))
		cat(as.yaml(list(
			num_iters=num_iters,
			base_measure_params=base_measure_params,
			concentration=concentration,
			density=density,
			beta_binomial_precision_params=beta_binomial_precision_params,
			working_dir=working_dir,
			trace_dir=subname,
			samples=samples)))
	sink()

	#-------------------
	# build event tables
	#-------------------

	subevents=list()

	for (samplenum in 1:length(subsamples)) {

		sample.t <- subsamples[samplenum] %>% unlist
		sample.n <- sample.tn[which(sample.tn$tumor==sample.t),"normal"] %>% as.character

		seg <- read.delim(str_c("../facets/",grep(str_c(sample.t,"_",sep=""),segfiles,value=TRUE))) %>%
			select(chrom,start=loc.start,end=loc.end,tcn.em,lcn.em) %>%
			filter(!is.na(tcn.em)&!is.na(lcn.em)) %>% # remove all rows with unassigned CN so midpoint assignment will find next closest segment
			rowwise %>%
			mutate(mid=(start+end)/2.0) %>%
			select(-c(start,end))

		# assign muts to their nearest CN segment
		submuts <- filter(muts,sample.name==sample.t) %>%
			mutate(id=str_c(chrom,pos,gene,effect,sep=":")) %>%
			rename_("cov.t"="cov") %>%
			left_join(seg,by="chrom") %>%
			group_by(id) %>%
			slice(which.min(abs(mid-pos))) %>%
			ungroup %>%
			mutate(minor=ifelse(lcn.em==0,1,0)) %>%
			mutate(major=ifelse(lcn.em==0,tcn.em-1,tcn.em)) %>%
			bind_cols(data.frame(cov.n=filter(muts,sample.name==sample.n)$cov)) %>%
			filter(maf>0.05 & cov.t>10)

		# create events table
		events <- data.frame(
			mutation_id=submuts$id,
			ref_counts=round(submuts$cov.n),
			var_counts=round(submuts$cov.t),
			normal_cn=rep(2,nrow(submuts)),
			minor_cn=submuts$lcn.em,
			major_cn=submuts$tcn.em-submuts$lcn.em
		)
		subevents<-c(subevents,list(events))
	}

	#-----------------------------------------------------------
	# remove events with ref=0 & var=0 depth accross all samples
	#-----------------------------------------------------------

	rmrows <- subevents %>% lapply(., function(t) which(t$ref_counts==0 & t$var_counts==0)) %>% unlist
	if(length(rmrows)>0){
		subevents <- lapply(subevents,function(t) t[-rmrows,])
	}

	#-----------------------------
	# build event & mutation files
	#-----------------------------

	for (samplenum in 1:length(subsamples)){

		sample <- subsamples[samplenum] %>% unlist

		cat(green("\n-") %+% " building input files for sample ",sample,":",sep="")

		cat("\n  events/",sample,".events.tsv",sep="")
		write.table(subevents[samplenum],file=str_c("events/",sample,".events.tsv"),row.names=FALSE,quote=FALSE,sep="\t")

		cat("\n  priors/",sample,".priors.yaml\n",sep="")
		system(str_c(sysprefix,"PyClone build_mutations_file --in_file events/",sample,".events.tsv --out_file priors/",sample,".priors.yaml"))

	}

	#-----------------
	# pyclone analysis
	#-----------------
	cat(green("\n-") %+% " running MCMC simulation:\n")
	system(str_c(sysprefix,"PyClone run_analysis --config_file config/",subname,".config.yaml"))

	#-------------
	# build tables
	#-------------
	cat(green("\n-") %+% " building analysis tables:\n  tables/",subname,".loci.tsv",sep="")
	system(str_c(sysprefix,"PyClone build_table --config_file config/",subname,".config.yaml --out_file tables/",subname,".loci.tsv --table_type loci"))
	cat("\n  tables/",subname,".cluster.tsv\n",sep="")
	system(str_c(sysprefix,"PyClone build_table --config_file config/",subname,".config.yaml --out_file tables/",subname,".cluster.tsv --table_type cluster"))

	#---------
	# plotting
	#---------
	cat(green("\n-") %+% " plotting results:\n  plots/",subname,".loci.pdf",sep="")
	system(str_c(sysprefix,"xvfb-run PyClone plot_loci --config_file config/",subname,".config.yaml --plot_file plots/",subname,".loci.pdf --plot_type density"))
	cat("\n  plots/",subname,".cluster.pdf\n",sep="")
	system(str_c(sysprefix,"xvfb-run PyClone plot_clusters --config_file config/",subname,".config.yaml --plot_file plots/",subname,".cluster.pdf --plot_type density"))

}

#-------------
# PostPy paths
#-------------

#"/ifs/e63data/reis-filho/usr/PostPy/interval_analyser.py"
#"/ifs/e63data/reis-filho/usr/PostPy/CI_filter.py"
#"/ifs/e63data/reis-filho/usr/PostPy/pyclone_files_updater.py"











