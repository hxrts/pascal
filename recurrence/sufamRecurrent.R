#!/usr/bin/env Rscript

#------
# input
#------

subsets<-read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)
mutfile<-read.delim("recurrent_mutations/recurrent_mutations.tsv",sep="\t",stringsAsFactors=FALSE)

#---------------------
# function definitions
#---------------------

# create VCF from subsetted recurrent mutations file
buildVCF<-function(muts,name){
	muts$ID<-apply(muts[,c("CHROM","POS","REF","ALT","ANN....GENE_SPLIT")],1,function(x)paste(gsub("^\\s+|\\s+$","",x),sep="",collapse=":"))
	vcf<-muts[,c("CHROM","POS","ID","REF","ALT")]
	colnames(vcf)<-c("#CHROM","POS","ID","REF","ALT")
	write.table(vcf,file=paste("recurrent_mutations/",name,"_recurrent.vcf",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
}

# run sufam for each sample
runSufam<-function(muts,name,subset){
	for (samplenum in 1:length(subset)){
		sample<-subset[samplenum]
		cmd<-paste("sufam /ifs/e63data/reis-filho/reference/human_g1k_v37.fa recurrent_mutations/",name,"_recurrent.vcf bam/",sample,".bam 2> log/",sample,"_sufamRecurrent.log > recurrent_mutations/",sample,"_sufamRecurrent.tsv",sep="")
		print(sprintf("  sufam command: ",cmd,"\n"))
		system(cmd)
	}
})

#----------
# main loop
#----------

for (subnum in 1:nrow(subsets)){

	# input processing
	line<-as.vector(subsets[subnum,])
	subset<-line[line!=""][-1]
	name<-line[[1]][1]
	muts<-mutfile[which(mutfile$TUMOR_SAMPLE%in%subset),]

	print(sprintf("* running sample set ",name,"\n"))
	print(sprintf("* building VCF files\n")
	buildVCF(muts,name)
	print(sprintf("running sufam on sample: ",name,"\n"))
	runSufam(muts,name,subset)

}

cat("* complete\n")
