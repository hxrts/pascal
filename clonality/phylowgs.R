#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

#suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(crayon))

# create necessary directories
system("mkdir phylowgs &>/dev/null")

# phyloWGS path
phylopath="/ifs/e63data/reis-filho/usr/phylowgs"

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

##fileformat=VCFv4.1									
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
##FORMAT=<ID=TD,Number=.,Type=Integer,Description="Tumor allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=ND,Number=.,Type=Integer,Description="Normal allelic depths for the ref and alt alleles in the order listed">
##INFO=<ID=TR,Number=1,Type=Integer,Description="Approximate tumor read depth; some reads may have been filtered">
##INFO=<ID=NR,Number=1,Type=Integer,Description="Approximate normal read depth; some reads may have been filtered">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TCGA-BJ-A191-01A-11D-A13W-08
1	241152	.	C	T	7.312238	.	SOMATIC	TD:ND:TR:NR	15,3:31,2:18:33



muts %>%
	#mutate(ID=str_c(CHROM,POS,ANN....GENE,sep=":")) %>%
	mutate(ID=".",QUAL="15",FILTER=".",INFO="SOMATIC",FORMAT="TD:ND:TR:NR") %>%
	rowwise() %>%
	mutate(TDR=NORMAL.DP,TDA=TUMOR.DP,NDR=G.FINAL.DP*(1-G.FINAL.MAF),NDA=G.FINAL.DP*G.FINAL.MAF,TR=FINAL.DP,NR=G.FINAL.DP) %>%
	mutate(SAMPLE_NAME=str_c(TDR,",",TDA,":",NDR,",",NDA,":",TR,":",NR)) %>%
	select(SAMPLE=TUMOR_SAMPLE,CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE_NAME) -> vmuts

#------------------------------
# main loop over sample subsets
#------------------------------

for (subnum in 1:nrow(subsets)){

	# input processing
	line<-as.vector(subsets[subnum,])
	subsamples<-line[line!=""][-1]
	subname<-line[[1]][1]

	cat(blue("\n--------------------------------\n  PHYLOWGS beginning subset ",subname,"\n-------------------------------\n",sep=""))

	# per sample?

	submuts <- filter(vmuts,SAMPLE==subsample)
	colnames(foo)[c(2,11)] <- c("#CHROM","hiya")

	colnames(muts)<-c("SAMPLE","#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO",FORMAT,SAMPLE_NAME)



}


