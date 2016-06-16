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
system("mkdir MAF &>/dev/null")

#------
# input
#------

cat(green("\n-") %+% " reading input\n")
subsets<-read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)

rmrow<-grep("#",subsets[,1])
if(length(rmrow)>0){subsets<-subsets[-rmrow,]}

muts<-read_tsv("summary/muts.tsv")

#------------------------------
# main loop over sample subsets
#------------------------------

for (subnum in 1:nrow(subsets)){

	# input processing
	line<-as.vector(subsets[subnum,])
	subsamples<-line[line!=""][-1]
	subname<-line[[1]][1]

	cat(blue("\n--------------------------------------\n  MAF plot - beginning subset ",subname,"\n-------------------------------------\n\n") %+% green("-") %+% " plotting MAF for samples: \n\n",sep="")

	system(str_c("mkdir MAF/",subname," &>/dev/null"))

	#------------------
	# loop over samples
	#------------------

	for (sample in subsamples) {

		submuts<-muts %>% filter(TUMOR_SAMPLE==sample)
		MAF<-submuts$FINAL.MAF[!is.na(submuts$FINAL.MAF)]

		cat("  ",sample,"\n",sep="")

		pdf(str_c("MAF/",subname,"/",sample,"_MAF_distribution.pdf"))

			# histogram with density plot
			print(ggplot(data.frame(MAF), aes(x=MAF))+
			geom_histogram(aes(y=..density..), binwidth=0.009, colour="darkgrey", fill="lightgrey")+
			stat_density(adjust=0.6,alpha=.5,colour="black",fill="darkblue")+
			scale_color_brewer(palette="Accent")+
			expand_limits(x = 0, y = 0)+
		  	theme_minimal(base_size = 10)+
		  	theme(plot.title=element_text(face="bold",size=20,vjust=0.5),
		  									axis.text=element_text(size=16),
		  									axis.title=element_text(size=16,face="bold"),
											axis.title.x=element_text(vjust=0.5),
											axis.title.y=element_text(angle=90,vjust=0.5),
											plot.margin = unit(c(1,1,1,1),"cm"))+
			ggtitle(str_c(sample," MAF distribution")))

		dev.off()

	}

}

cat(green("\n-") %+% " done\n")

