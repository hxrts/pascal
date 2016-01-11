#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

suppressMessages(require(openxlsx))
suppressMessages(require(plyr))
suppressMessages(require(readr))
suppressMessages(require(stringr))
suppressMessages(require(crayon))

# create necessary directories
suppressMessages(system("mkdir phylowgs &>/dev/null"))

# phyloWGS path
phylopath="/ifs/e63data/reis-filho/usr/phylowgs"

#------
# input
#------

cat(green("\n-") %+% " reading input\n")
subsets<-read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)

rmrow<-grep("#",subsets[,1])
if(length(rmrow)>0){subsets<-subsets[-rmrow,]}

muts<-read.xlsx("excel/mutation_summary.xlsx",sheet="SNV_HIGH_MODERATE_SUMMARY",check.names=TRUE)
segfiles<-list.files("facets",pattern="*cncf.txt")

#--------------------
# data pre-processing
#--------------------

# assign unique mutation IDS
muts$ID<-paste(muts$CHROM,muts$POS,muts$ANN....GENE,sep=":")

muts %>%
	mutate(ID=str_c(CHROM,POS,ANN....GENE,sep=":")) %>%
	select(CHROM=CHROM,POS=POS,ID=ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE,GENE=ANN....GENE,,)
	glimpse()


#------------------------------
# main loop over sample subsets
#------------------------------

for (subnum in 1:nrow(subsets)){

	# input processing
	line<-as.vector(subsets[subnum,])
	subsamples<-line[line!=""][-1]
	subname<-line[[1]][1]

	cat(blue("\n--------------------------------\n  PHYLOWGS beginning subset ",subname,"\n-------------------------------\n",sep=""))

	system(paste("mkdir",subname,"&>/dev/null"))

}


