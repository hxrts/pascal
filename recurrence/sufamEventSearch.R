#!/usr/bin/env Rscript

#------
# input
#------

if(file.exists("events.vcf")==FALSE){print("* no events.vcf file")}
subsets<-read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)

#----------
# functions
#----------

# run sufam for each sample
runSufam<-function(samples){
	for (samplenum in 1:length(samples)){
		sample<-samples[samplenum]
		cmd<-paste("sufam /ifs/e63data/reis-filho/reference/human_g1k_v37.fa recurrent_mutations/events.vcf bam/",sample,".bam 2> log/",sample,"_sufamEventSearch.log > recurrent_mutations/",sample,"_sufamEventSearch.tsv",sep="")
		print(paste("running sufam on sample",sample,"against events.vcf with command",cmd))
		system(cmd)
	}
}

#----------
# main loop
#----------

for (subnum in 1:nrow(subsets)){

	line<-as.vector(subsets[subnum,])
	subset<-line[line!=""][-1]
	# samples to vector
	samples<-unlist(unique(subset))

	# execute functions
	runSufam(samples)

}

cat("* complete\n")




