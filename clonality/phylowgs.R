#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

# load base libraries
suppressMessages(pacman::p_load(dplyr,readr,tidyr,magrittr,stringr,crayon))

# anaconda prefix
sysprefix="umask 002 && unset PYTHONPATH && source /home/bermans/miniconda2/envs/phylowgs/bin/activate /home/bermans/miniconda2/envs/phylowgs >/dev/null 2>&1 && "
# phyloWGS path
phylopath="python2 /ifs/e63data/reis-filho/usr/phylowgs/"

wd<-getwd()
print(wd)
# create necessary directories
system(str_c("mkdir ",wd,"/phylowgs ",wd,"/phylowgs/step0 ",wd,"/phylowgs/step0 ",wd,"/phylowgs/step1 ",wd,"/phylowgs/step2 ",wd,"/phylowgs/step3 ",wd,"/phylowgs/results &>/dev/null"))

#------
# input
#------

cat(green("\n-") %+% " reading input\n")
subsets<-read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)

rmrow<-grep("#",subsets[,1])
if(length(rmrow)>0){subsets<-subsets[-rmrow,]}

muts<-read_tsv("summary/muts.tsv")
#muts<-read.xlsx("summary/mutation_summary.xlsx",sheet="SNV_HIGH_MODERATE_SUMMARY",check.names=TRUE)
segfiles<-list.files("facets",pattern="*cncf.txt")

#--------------------
# data pre-processing
#--------------------

muts %>%
	mutate(ID=".",QUAL="15.0",FILTER=".",INFO="SOMATIC",FORMAT="TD:ND:TR:NR") %>%
	mutate(FINAL.DP=replace(FINAL.DP,is.na(FINAL.DP),0),FINAL.MAF=replace(FINAL.MAF,is.na(FINAL.MAF),0)) %>% #filter(TUMOR_SAMPLE==sample) %>% head(82) %>% tail(1) %>% glimpse
	rowwise() %>%
	mutate(TDR=NORMAL.DP,TDA=TUMOR.DP,NDR=round(G.FINAL.DP*(1-G.FINAL.MAF)),NDA=round(G.FINAL.DP*G.FINAL.MAF),TR=FINAL.DP,NR=G.FINAL.DP) %>%
	select(-SAMPLE) %>%  
	rename(SAMPLE=TUMOR_SAMPLE) %>%
	ungroup() %>%
	mutate(SAMPLE_NAME=str_c(TDR,",",TDA,":",NDR,",",NDA,":",TR,":",NR)) %>%
	select(SAMPLE,CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE_NAME) -> vmuts

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







