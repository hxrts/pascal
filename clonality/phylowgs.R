#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

# load base libraries
suppressMessages(pacman::p_load(dplyr,readr,tidyr,magrittr,purrr,stringr,rlist,crayon))

# anaconda prefix
sysprefix="umask 002 && unset PYTHONPATH && source /home/bermans/miniconda2/envs/phylowgs/bin/activate /home/bermans/miniconda2/envs/phylowgs >/dev/null 2>&1 && "
# phyloWGS path
phylopath="python2 /ifs/e63data/reis-filho/usr/phylowgs/"

wd<-getwd()
# create necessary directories
system(str_c("mkdir ",wd,"/phylowgs ",wd,"/phylowgs/step0 ",wd,"/phylowgs/step0 ",wd,"/phylowgs/step1 ",wd,"/phylowgs/step2 ",wd,"/phylowgs/step3 ",wd,"/phylowgs/results &>/dev/null"))


#--------------------
# data pre-processing
#--------------------

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

sample.tn[] <- data.frame(normal=as.character(sample.n),tumor=as.character(sample.t)) %>% lapply(as.character) %>% tbl_df


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

muts.tn <-	muts %>%
	inner_join(sample.tn,by=c("sample.name"="tumor")) %>%
	rename(tumor=sample.name) %>%
	left_join(muts,by=c("normal"="sample.name","chrom","pos","effect","gene")) %>%
	rename(cov.t=cov.x,maf.t=maf.x,cov.n=cov.y,maf.n=maf.y) %>%
	filter(cov.t>1 & maf.t>0)


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




muts %>%
	mutate(ID=".",QUAL="990.0",FILTER=".",INFO="SOMATIC",FORMAT="TD:ND:TR:NR") %>%
	mutate(FINAL.DP=replace(FINAL.DP,is.na(FINAL.DP),0),FINAL.MAF=replace(FINAL.MAF,is.na(FINAL.MAF),0)) %>% #filter(TUMOR_SAMPLE==sample) %>% head(82) %>% tail(1) %>% glimpse
	rowwise() %>%
	mutate(TDR=NORMAL.DP,TDA=TUMOR.DP,NDR=round(G.FINAL.DP*(1-G.FINAL.MAF)),NDA=round(G.FINAL.DP*G.FINAL.MAF),TR=FINAL.DP,NR=G.FINAL.DP) %>%
	select(-SAMPLE) %>%  
	rename(SAMPLE=TUMOR_SAMPLE) %>%
	ungroup() %>%
	mutate(SAMPLE_NAME=str_c(TDR,",",TDA,":",NDR,",",NDA,":",TR,":",NR)) %>%
	select(SAMPLE,CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE_NAME) -> vmuts

segfiles<-list.files("facets",pattern="*cncf.txt")


#------------------------------
# main loop over sample subsets
#------------------------------

subnum=1

for (subnum in 1:nrow(subsets)){

	# input processing
	line<-as.vector(subsets[subnum,])
	subsamples<-line[line!=""][-1]
	subname<-line[[1]][1]

	cat(blue("\n--------------------------------\n  PHYLOWGS beginning subset ",subname,"\n-------------------------------\n",sep=""))

	#----
	# CNV
	#----

	getSeg <- function(sample){
		sample %>%
			str_c("_") %>%
			grep(segfiles,value=TRUE) %>%
			str_c("facets/",.) %>%
			read_tsv() %>%
			rowwise() %>%
			mutate(SEGMID=(loc.end-loc.start)/2)
		}

	subsamples %>%
		lapply(. %>% getSeg()) ->
		segs

	getBreaks <- function(segs){
		lapply(segs,function(seg){
			seg %>%
				mutate(loc.start,loc.end=loc.end+1) %>%
				gather(ordinal,position,loc.start:loc.end) %>%
				select(chrom,position)
		}) -> cseg
		do.call(rbind,cseg) %>%
		arrange(chrom,position) %>%
		distinct()
	}

	breaks<-getBreaks(segs)

	getChromSpan<-function(chrombreak){
		data.frame(START=chrombreak$position[-nrow(chrombreak)],END=(chrombreak$position[-1])-1)
	}

	lapply(unique(breaks$chrom),function(CHROM) cbind(CHROM,getChromSpan(subset(getBreaks(segs),chrom==CHROM)))) %>% do.call(rbind,.) ->
		segspans

	segspans %>%
		rowwise() %>%
		mutate(MID=(START+END)/2) ->
		breakmids

	exSeg <- function(breakmids){
		breakmids %>%
			rowwise() %>%
			filter(segs)
	}

	xsegs<- lapply(segs,function(seg){
		breakmids %>% left_join(seg,by=c("CHROM"="chrom")) %>% ungroup() %>% group_by(MID) %>% slice(which.min(abs(SEGMID-MID)))
	})

	segtables <- lapply(xsegs, function(xseg){
		xseg %>%
			rowwise() %>%
			mutate(major_cn=tcn.em-lcn.em) %>%
			#mutate(CHROM=replace(CHROM,CHROM=="X",23)) %>%
			mutate(cf=round(cf,2)) %>%
			ungroup() %>%
			arrange(CHROM,START) %>%
			select_("Chromosome"="CHROM","Start_Position(bp)"="START","End_Position(bp)"="END","copy_number"="tcn.em","MinorCN"="lcn.em","MajorCN"="major_cn","Clonal_Frequency"="cf")
	})

	#------------------
	# loop over samples
	#------------------

	for (i in 1:length(subsamples)) {

		sample <- subsamples[i]

		cat(blue("\n",sample,"\n",sep=""))

		segname <- str_c(wd,"/phylowgs/step0/",sample,".seg.txt")
		write_tsv(segtables[[i]],segname)

		#----
		# SSM
		#----

		vmuts %>%
			filter(SAMPLE==sample) %>%
			select(-SAMPLE) %>%
			rename_("#CHROM"="CHROM") %>%
			rename_(.dots=setNames("SAMPLE_NAME",sample)) ->
			submuts

		vcfname <- str_c(wd,"/phylowgs/step0/",sample,".vcf")
		sink(vcfname)
			cat('##fileformat=VCFv4.1									
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
##FORMAT=<ID=TD,Number=.,Type=Integer,Description="Tumor allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=ND,Number=.,Type=Integer,Description="Normal allelic depths for the ref and alt alleles in the order listed">
##INFO=<ID=TR,Number=1,Type=Integer,Description="Approximate tumor read depth; some reads may have been filtered">
##INFO=<ID=NR,Number=1,Type=Integer,Description="Approximate normal read depth; some reads may have been filtered">\n')
		sink()
		submuts %>% write_tsv(path=vcfname,append=TRUE,col_names=TRUE)

		cellu<-muts %>% filter(TUMOR_SAMPLE==sample) %$% cancer.cell.frac %>% max

		str_c(phylopath,"parser/parse_cnvs.py -f titan -c ",cellu," --cnv-output ",wd,"/phylowgs/step1/",sample,".seg.txt ",wd,"/phylowgs/step0/",sample,".seg.txt") %>%
			system()

		str_c(phylopath,"parser/create_phylowgs_inputs.py --cnvs ",wd,"/phylowgs/step1/",sample,".seg.txt --output-cnvs ",wd,"/phylowgs/step2/",sample,".cnv.txt -v mutect_tcga ",wd,"/phylowgs/step0/",sample,".vcf --output-variants ",wd,"/phylowgs/step2/",sample,".ssm.txt") %>%
			system()
	}

	#------------
	# combine SSM
	#------------

	str_c(wd,"/phylowgs/step2/",subsamples,".ssm.txt") %>%
		lapply(function(ssm) read_tsv(ssm)) ->
		ssms
	
	ssms[[1]] %>%
		left_join(ssms[[2]],by="id") ->
		ssm2

	ssm2 %>%
		rowwise() %>%
		mutate(a=str_c(a.x,",",a.y),
			   d=str_c(d.x,",",d.y)) %>%
		select(id,gene=gene.x,a,d,mu_r=mu_r.x,mu_v=mu_v.y) %>%
		write_tsv(str_c(wd,"/phylowgs/step3/",subname,".ssm.txt"))

	#------------
	# combine CNV
	#------------

	str_c(wd,"/phylowgs/step2/",subsamples,".cnv.txt") %>%
		lapply(function(cnv) read_tsv(cnv)) ->
		cnvs
	
	cnvs[[1]] %>%
		left_join(cnvs[[2]],by="cnv") ->
		cnv2

	gluessm<-function(ssms.x,ssms.y){
		paste(
			c(ssms.x,ssms.y)[!is.na(c(ssms.x,ssms.y))],collapse=";"
			)}

	cnv2 %>%
		rowwise() %>%
		mutate(a=str_c(a.x,",",a.y),
			   d=str_c(d.x,",",d.y),
			   ssms=gluessm(ssms.x,ssms.y)) %>%
		select(cnv,a,d,ssms) %>%
		write_tsv(str_c(wd,"/phylowgs/step3/",subname,".cnv.txt"))

	#-------------
	# run phyloWGS
	#-------------

	cat(blue("\n-------------------------------\n  starting PHYLOWGS calculation ",subname,"\n-------------------------------\n",sep=""))

	system(str_c("mkdir phylowgs/results/",subname," &>/dev/null"))
	# setwd(str_c("phylowgs/results/",subname))
	# str_c(phylopath,"evolve.py -s 20 -i 20 ",wdorigin,"//phylowgs/step3/",subname,".ssm.txt ",wdorigin,"//phylowgs/step3/",subname,".cnv.txt") %>%
	# 	system()
	# setwd("../../..")

	str_c(phylopath,"evolve.py ",wd,"/phylowgs/step3/",subname,".ssm.txt ",wd,"/phylowgs/step3/",subname,".cnv.txt") %>%
		system()

}







