#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(crayon))
suppressMessages(library(clonevol))

# create necessary directories
system("mkdir clonevol &>/dev/null")

# for testing purposes
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

#------------------------------
# main loop over sample subsets
#------------------------------

for (subnum in 1:nrow(subsets)){

	# input processing
	line<-as.vector(subsets[subnum,])
	subsamples<-line[line!=""][-1]
	subname<-line[[1]][1]

	cat(blue("\n--------------------------\n  beginning subset ",subname,"\n--------------------------\n\n",sep=""))

	system(str_c("mkdir clonevol/",subname," &>/dev/null"))

	#-----------------------
	# clusterVAF preperation
	#-----------------------

	cat(green("\n-") %+% " clusterVAF preperation\n",sep="")

	# cluster | prim.vaf | met1.vaf | ... | metn.vaf
	loci<-read_tsv(str_c("pyclone/tables/",subname,".loci.tsv"))
	loci %>%
		spread(key=sample_id,value=variant_allele_frequency,fill=0) %>%
		#rename_(prim.vaf=samples[1],met1.vaf=samples[2]) %>%
		select(cluster=cluster_id,one_of(subsamples)) %>%
		arrange(cluster)->
		#write_tsv(mID,str_c("clonevol/clusterVAF/clusterVAF.",subname,".tsv"))
		CD

	CD[,-1]<-round(CD[,-1]*100)

	#-------------------
	# clonevol algorithm
	#-------------------

	cat(green("\n-") %+% " running clonevol algorithm\n",sep="")

	CM <- infer.clonal.models(variants=CD,
		cluster.col.name="cluster",
		vaf.col.names=subsamples,
		subclonal.test="bootstrap",
		subclonal.test.model="non-parametric",
		cluster.center="mean",
		num.boots=1000,
		founding.cluster=as.numeric(loci[which.max(loci$variant_allele_frequency),"cluster_id"]),
		min.cluster.vaf=0.01,
		p.value.cutoff=0.01,
		alpha=0.1,
		random.seed=63108)

	cat(green("\n-") %+% " clonevol plot 1\n",sep="")

	plot.clonal.models((CM$models)[1],
		matched=CM$matched,
 		variants=CD,
		box.plot=FALSE,
		out.format="pdf",
		overwrite.output=TRUE,
		scale.monoclonal.cell.frac=TRUE,
 		cell.frac.ci=TRUE,
		tree.node.shape="circle",
		tree.node.size=40,
		tree.node.text.size=0.65,
		width=8, height=5,
		out.dir=str_c("clonevol/",subname,"/evo"))

	cat(green("\n-") %+% " clonevol plot 2\n",sep="")

	var.to.highlight = filter(CD,cluster==as.numeric(loci[which.max(loci$variant_allele_frequency),"cluster_id"]))
	colnames(var.to.highlight) = c("cluster", "variant.name")
	plot.clonal.models((CM$models)[1],
		matched=CM$matched,
		variants=CD,
		box.plot=FALSE,
		out.format="pdf",
		overwrite.output=TRUE,
		scale.monoclonal.cell.frac=TRUE,
		cell.frac.ci=TRUE,
		variants.to.highlight=var.to.highlight,
		variant.color="blue",
		variant.angle=60,
		tree.node.shape="circle",
		tree.node.size=40,
		tree.node.text.size=0.65,
		width=8, height=5,
		out.dir=str_c("clonevol/",subname,"/evohilight"))

}
