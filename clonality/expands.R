#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

#suppressMessages(require(openxlsx))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(crayon))
suppressMessages(library(expands))

# create necessary directories
suppressMessages(system("mkdir expands &>/dev/null"))

subnum=samplenum=1

#------
# input
#------

cat(green("\n-") %+% " reading input\n")
subsets<-read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)

rmrow<-grep("#",subsets[,1])
if(length(rmrow)>0){subsets<-subsets[-rmrow,]}

muts<-read_tsv("excel/muts.tsv")
segfiles<-list.files("facets",pattern="*cncf.txt")

setwd("expands")

#------------------------------
# main loop over sample subsets
#------------------------------

for (subnum in 1:nrow(subsets)){

	# input processing
	line<-as.vector(subsets[subnum,])
	subsamples<-line[line!=""][-1]
	subname<-line[[1]][1]

	cat(blue("\n--------------------------\n  beginning subset ",subname,"\n--------------------------\n\n",sep=""))

	#-------------------
	# build event tables
	#-------------------

	for (samplenum in 1:length(subsamples)) {

		sample=subsamples[samplenum]

		# SNV
		muts %>%
			filter(TUMOR_SAMPLE==sample) %>%
			mutate(PN_B=0) %>%
			select(chr=CHROM,startpos=POS,AF_Tumor=cancer.cell.frac,PN_B) %>%
			mutate(chr=replace(chr,chr=="X",23)) ->
			SNV

		# CBS
		seg<-read_tsv(str_c("../facets/",grep(str_c(sample,"_"),segfiles,value=TRUE)))
		
		seg %>%
			mutate(TUMOR_SAMPLE=unlist(lapply(ID,function(x)strsplit(x,"_")))[1]) %>%
			filter(TUMOR_SAMPLE==sample) %>%
			select(chr=chrom,startpos=loc.start,endpos=loc.end,CN_Estimate=lcn.em) %>%
			mutate(CN_Estimate=replace(CN_Estimate,which(is.na(CN_Estimate)),1)) ->
			mseg

		exo<-runExPANdS(data.matrix(SNV),data.matrix(mseg),precision=0.01,snvF=str_c(sample,".",subname))

	}


}


