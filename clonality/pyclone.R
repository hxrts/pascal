#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

suppressMessages(require(yaml))
suppressMessages(require(openxlsx))
suppressMessages(require(plyr))
suppressMessages(require(dplyr))
suppressMessages(require(crayon))

# create necessary directories
suppressMessages(system("mkdir pyclone pyclone/config pyclone/events pyclone/priors pyclone/tables pyclone/plots &>/dev/null"))

#------
# input
#------

cat(green("\n-") %+% " reading input\n")
subsets<-read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)

rmrow<-grep("#",subsets[,1])
if(length(rmrow)>0){subsets<-subsets[-rmrow,]}

muts<-read.xlsx("excel/mutation_summary.xlsx",sheet="SNV_HIGH_MODERATE_SUMMARY",check.names=TRUE)
segfiles<-list.files("facets",pattern="*cncf.txt")

setwd("pyclone")	# for some reason PyClone needs to be run from the root directory

#--------------------
# data pre-processing
#--------------------

# assign unique mutation IDS
muts$ID<-paste(muts$CHROM,muts$POS,muts$ANN....GENE,sep=":")

#---------------------------------
# variables for yaml configuration
#---------------------------------

num_iters=as.integer(100000)
base_measure_params<-list(alpha=as.integer(1),beta=as.integer(1))
concentration<-list(value=as.integer(1),prior=list(shape=1.0,rate=0.001))
density<-"pyclone_beta_binomial"
beta_binomial_precision_params<-list(value=as.integer(1000),prior=list(shape=1.0,rate=0.0001),proposal=list(precision=0.01))
working_dir<-getwd()

#------------------------------
# main loop over sample subsets
#------------------------------

for (subnum in 1:nrow(subsets)){

	# input processing
	line<-as.vector(subsets[subnum,])
	subsamples<-line[line!=""][-1]
	subname<-line[[1]][1]

	cat(blue("\n--------------------------\n  beginning subset ",subname,"\n--------------------------\n",sep=""))

	system(paste("mkdir",subname,"&>/dev/null"))

	#----------------------------
	# run-specific yaml variables
	#----------------------------

	trace_dir=subname
	samples<-lapply(subsamples,function(sample)
				list(
					mutations_file=paste("priors/",sample,".priors.yaml",sep=""),
					tumour_content=list(value=1.0),
					error_rate=0.001)
			)
	names(samples)<-subsamples

	#----------------
	# write yaml file
	#----------------

	cat(green("\n-") %+% " building configuration file:\n  config/",subname,".config.yaml\n",sep="")
	sink(file=paste("config/",subname,".config.yaml",sep=""))
		cat(as.yaml(list(num_iters=num_iters,base_measure_params=base_measure_params,concentration=concentration,density=density,beta_binomial_precision_params=beta_binomial_precision_params,working_dir=working_dir,trace_dir=trace_dir,samples=samples)))
	sink()

	#-------------------
	# build event tables
	#-------------------

	subevents=list()
	for (samplenum in 1:length(subsamples)) {

		sample=subsamples[samplenum]

		submuts<-filter(muts,TUMOR_SAMPLE==sample)
		submuts[submuts$CHROM=="X","CHROM"]<-23

		seg<-read.delim(paste("../facets/",grep(paste(sample,"_",sep=""),segfiles,value=TRUE),sep=""))
		seg<-seg[!is.na(seg$tcn.em)&!is.na(seg$lcn.em),]	# remove all rows with unassigned CN so midpoint assignment will find next closest segment

		# assign muts to their nearest CN segment
		seg$mid<-rowMeans(select(seg,loc.start,loc.end))
		submutseg<-adply(submuts,1,transform,seg=slice(filter(seg,chrom==CHROM),which.min(abs(mid-POS))))

		events<-data.frame(		# variables for alternate PyClone run: variant_case, variant_freq, genotype
			mutation_id=submutseg$ID,
			ref_counts=round(submutseg$NORMAL.DP),
			var_counts=round(submutseg$TUMOR.DP),
			normal_cn=rep(2,nrow(submuts)),		#1+(submutseg$CHROM!="Y"),	# all chromosomes set to 2 except Y set to 1				!! must use more robust method here !!
			minor_cn=submutseg$seg.lcn.em,
			major_cn=submutseg$seg.tcn.em-submutseg$seg.lcn.em
		)
		subevents<-c(subevents,list(events))
	}

	#---------------------------------------
	# remove events with ref=0 & var=0 depth
	#---------------------------------------

	rmrows<-unlist(lapply(subevents,function(t) which(t$ref_counts==0&t$var_counts==0)))
	if(length(rmrows)>0){
		subevents<-lapply(subevents,function(t) t[-rmrows,])
	}

	#-----------------------------
	# build event & mutation files
	#-----------------------------

	for (samplenum in 1:length(subsamples)){

		sample<-subsamples[samplenum]

		cat(green("\n-") %+% " building input files for sample ",sample,":",sep="")

		cat("\n  events/",sample,".events.tsv",sep="")
		write.table(subevents[samplenum],file=paste("events/",sample,".events.tsv",sep=""),row.names=FALSE,quote=FALSE,sep="\t")

		cat("\n  priors/",sample,".priors.yaml\n",sep="")
		system(paste("source activate pyclone >/dev/null 2>&1 && PyClone build_mutations_file --in_file events/",sample,".events.tsv --out_file priors/",sample,".priors.yaml",sep=""))

	}

	#-----------------
	# pyclone analysis
	#-----------------
	cat(green("\n-") %+% " running MCMC simulation:\n")
	system(paste("source activate pyclone >/dev/null 2>&1 && PyClone run_analysis --config_file config/",subname,".config.yaml",sep=""))

	#-------------
	# build tables
	#-------------
	cat(green("\n-") %+% " building analysis tables:\n  tables/",subname,".loci.tsv",sep="")
	system(paste("source activate pyclone >/dev/null 2>&1 && PyClone build_table --config_file config/",subname,".config.yaml --out_file tables/",subname,".loci.tsv --table_type loci",sep=""))
	cat("\n  tables/",subname,".cluster.tsv\n",sep="")
	system(paste("source activate pyclone >/dev/null 2>&1 && PyClone build_table --config_file config/",subname,".config.yaml --out_file tables/",subname,".cluster.tsv --table_type cluster",sep=""))

	#---------
	# plotting
	#---------
	cat(green("\n-") %+% " plotting results:\n  plots/",subname,".loci.pdf",sep="")
	system(paste("source activate pyclone >/dev/null 2>&1 && PyClone plot_loci --config_file config/",subname,".config.yaml --plot_file plots/",subname,".loci.pdf --plot_type density",sep=""))
	cat("\n  plots/",subname,".cluster.pdf\n",sep="")
	system(paste("source activate pyclone >/dev/null 2>&1 && PyClone plot_clusters --config_file config/",subname,".config.yaml --plot_file plots/",subname,".cluster.pdf --plot_type density",sep=""))

	#system(paste("source activate pyclone >/dev/null 2>&1 && PyClone plot_cellular_frequencies --config_file config/",subname,".config.yaml --plot_file plots/",subname,".cluster.pdf --plot_type density",sep=""))
	#PyClone plot_cellular_frequencies config.yaml cellular_frequencies/ --burnin 1000

}

#-------
# PostPy
#-------

#"/ifs/e63data/reis-filho/usr/PostPy/interval_analyser.py"
#"/ifs/e63data/reis-filho/usr/PostPy/CI_filter.py"
#"/ifs/e63data/reis-filho/usr/PostPy/pyclone_files_updater.py"











