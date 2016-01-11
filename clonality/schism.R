#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

suppressMessages(library(yaml))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(crayon))

# create necessary directories
system("mkdir schism schism/mutID schism/mutRCNA &>/dev/null")

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

	cat(blue("\n--------------------------------\n  SCHISM beginning subset ",subname,"\n--------------------------------\n",sep=""))

	#-------------------
	# build mutID tables
	#-------------------

	loci<-read_tsv(str_c("pyclone/tables/",subname,".loci.tsv"))

	loci %>%
		select(mutationID=mutation_id,clusterID=cluster_id) ->
		mID
	write_tsv(mID,str_c("schism/mutID/mutID.",subname,".tsv"))

	#---------------------
	# build mutRCNA tables
	#---------------------

	loci %>%
		mutate(cellularity=round(cellular_prevalence,4),sd=round(cellular_prevalence_std,4)) %>%
		select(sampleID=sample_id,mutationID=mutation_id,cellularity,sd) ->
		mRC
	write_tsv(mRC,str_c("schism/mutRCNA/mutRCNA.",subname,".tsv"))

	system(str_c("mkdir schism/",subname," &>/dev/null"))

	#----------------
	# write yaml file
	#----------------

	cat(green("\n-") %+% " building configuration file:\n  schism/config/",subname,".config.yaml\n",sep="")

	sink(file=str_c("schism/",subname,".config.yaml"))
		cat(as.yaml(list(

			working_dir=str_c(getwd(),"/schism"),
			mutation_to_cluster_assignment=str_c("mutID/mutID.",subname,".tsv"),
			mutation_cellularity_input=str_c("mutRCNA/mutRCNA.",subname,".tsv"),
			output_prefix=str_c(subname,"/",subname),
			cellularity_estimation=str_c("mutRCNA/mutRCNA.",subname,".tsv"),

			hypothesis_test=list(
				test_level="mutations",
				significance_level=0.05,
				store_pvalues="True"),

			genetic_algorithm=list(
				instance_count=as.integer(10),
				generation_count=as.integer(50),
				generation_size=as.integer(100),
				random_object_fraction=0.2,
				mutation_probability=0.9,
				crossover_probability=0.25,
				fitness_coefficient=5.0,
				verbose="False")
		)))
	sink()

	cat(green("\n-") %+% " running SCHISM algorithm\n",sep="")

	system(str_c("runSchism analyze -c schism/",subname,".config.yaml"))

}
