#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

# load base libraries
suppressMessages(pacman::p_load(dplyr,readr,stringr,tidyr,magrittr,crayon,clonevol))

# create necessary directories
system("mkdir clonevol &>/dev/null")

# for testing purposes
subnum=2
samplenum=1

#------
# input
#------

cat(green("\n-") %+% " reading input\n")
subsets<-read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)

rmrow<-grep("#",subsets[,1])
if(length(rmrow)>0){subsets<-subsets[-rmrow,]}

muts<-read_tsv("summary/muts.tsv")
segfiles<-list.files("facets",pattern="*cncf.txt")

#------------------------------
# main loop over sample subsets
#------------------------------

for (subnum in 1:nrow(subsets)){

	# input processing
	line<-as.vector(subsets[subnum,])
	subsamples<-line[line!=""][-1]
	subname<-line[[1]][1]

	cat(blue("\n--------------------------------\n  CLONEVOL beginning subset ",subname,"\n--------------------------------\n\n",sep=""))

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

	modeltype="non-parametric"		#'normal', 'normal-truncated', 'beta', 'binomial', 'beta-binomial', or 'non-parametric'

	CM <- infer.clonal.models(variants=CD,
		cluster.col.name="cluster",
		vaf.col.names=subsamples,
		subclonal.test="bootstrap",
		#subclonal.test.model="non-parametric",
		subclonal.test.model=modeltype,
		cluster.center="mean",
		num.boots=1,
		founding.cluster=as.numeric(loci[which.max(loci$variant_allele_frequency),"cluster_id"]),
		min.cluster.vaf=0.03,
		p.value.cutoff=0.03,
		alpha=0.1,
		random.seed=1)


	var.to.highlight = filter(CD,cluster==as.numeric(loci[which.max(loci$variant_allele_frequency),"cluster_id"]))
	colnames(var.to.highlight) = c("cluster", "variant.name")

	try(
		plot.clonal.models((CM$models),
			matched=CM$matched,
			variants=CD,
			box.plot=TRUE,
			out.format="pdf",
			overwrite.output=TRUE,
			scale.monoclonal.cell.frac=TRUE,
			cell.frac.ci=TRUE,
			#variants.to.highlight=var.to.highlight,
			variant.color="blue",
			variant.angle=60,
			tree.node.shape="circle",
			tree.node.size=40,
			tree.node.text.size=0.65,
			width=8, height=5,
			out.dir=str_c("clonevol/",subname,"/",modeltype))
	)

	num.clusters <- length(unique(CD$cluster))
	pdf(str_c("clonevol/",subname,"/",modeltype,"/",subname,".violin.pdf"))
		variant.box.plot(CD,
		                 vaf.col.names=subsamples,
		                 variant.class.col.name=NULL,
		                 cluster.axis.name="",
		                 vaf.limits=70,
		                 violin=TRUE,
		                 box=TRUE,
		                 order.by.total.vaf=FALSE,
		                 jitter=FALSE,
		                 jitter.center.method="mean",
		                 jitter.center.size=0.5,
		                 jitter.center.color="darkgray",
		                 jitter.shape=1,
		                 jitter.color=get.clonevol.colors(num.clusters),
		                 jitter.size=2,
		                 jitter.alpha=1,
		                 #highlight="is.cancer.gene",
		                 #highlight.note.col.name="gene",
		                 #highlight.shape=19,
		                 display.plot=TRUE)
	dev.off()

}
