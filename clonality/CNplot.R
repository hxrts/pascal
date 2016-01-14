#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

#suppressMessages(library(openxlsx))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
#suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(crayon))

# create necessary directories
system("mkdir CN &>/dev/null")

#------
# input
#------

cat(green("\n-") %+% " reading input\n")
subsets<-read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)

rmrow<-grep("#",subsets[,1])
if(length(rmrow)>0){subsets<-subsets[-rmrow,]}

segfiles<-list.files("facets",pattern=".cncf.txt")

printPlot<-function(CN,weights,xlabel,plottitle,pdfname){
	pdf(pdfname)
		# histogram with density plot
		print(ggplot(data.frame(CN), aes(x=CN,weight=weights))+
		geom_histogram(aes(y=..density..,weight=weights), binwidth=0.015, colour="black", fill="grey")+
		stat_density(aes(y=..density..,weight=weights),adjust=0.1,alpha=.5,colour="black",fill="darkblue")+
		scale_color_brewer(palette="Accent")+
		expand_limits(x = 0, y = 0)+
	  	theme_minimal(base_size = 10)+
	  	theme(plot.title=element_text(face="bold",size=20,vjust=0.5),
	  									axis.text=element_text(size=16),
	  									axis.title=element_text(size=16,face="bold"),
										axis.title.x=element_text(vjust=0.5),
										axis.title.y=element_text(angle=90,vjust=0.5),
										plot.margin = unit(c(1,1,1,1),"cm"))+
	  	labs(x=xlabel)+
		ggtitle(plottitle))
	dev.off()
}

#------------------------------
# main loop over sample subsets
#------------------------------

for (subnum in 1:nrow(subsets)){

	# input processing
	line<-as.vector(subsets[subnum,])
	subsamples<-line[line!=""][-1]
	subname<-line[[1]][1]

	cat(blue("\n--------------------------------------\n  CN plot - beginning subset ",subname,"\n-------------------------------------\n\n") %+% green("-") %+% " plotting MAF for samples: \n\n",sep="")

	system(str_c("mkdir CN/",subname," &>/dev/null"))

	#------------------
	# loop over samples
	#------------------

	for (sample in subsamples) {

		seg<-read_tsv(str_c("facets/",grep(str_c(sample,"_"),segfiles,value=TRUE)))

		cat("  ",sample,"\n",sep="")

		CN<-seg$cnlr.median
		weights<-(seg$loc.end-seg$loc.start)/sum(as.numeric(seg$loc.end-seg$loc.start))
		xlabel<-"Copy number log ratio"
		plottitle<-str_c(subname," CN log ratio distribution")
		pdfname<-str_c("CN/",subname,"/",sample,"_LOG_RATIO_CN_distribution.pdf")

		printPlot(CN,weights,xlabel,plottitle,pdfname)

		CN<-seg$tcn.em
		weights<-(seg$loc.end-seg$loc.start)/sum(as.numeric(seg$loc.end-seg$loc.start))
		plottitle<-str_c(subname," CN total copy number distribution")
		xlabel<-"Total copy number"
		pdfname<-str_c("CN/",subname,"/",sample,"_TOTAL_CN_distribution.pdf")

		printPlot(CN,weights,xlabel,plottitle,pdfname)

		CN<-seg$lcn.em
		weights<-(seg$loc.end-seg$loc.start)/sum(as.numeric(seg$loc.end-seg$loc.start))
		plottitle<-str_c(subname," CN lower copy number distribution")
		xlabel<-"Lower copy number"
		pdfname<-str_c("CN/",subname,"/",sample,"_LOWER_CN_distribution.pdf")

		printPlot(CN,weights,xlabel,plottitle,pdfname)

	}

}

cat(green("\n-") %+% " done\n")

