#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

#suppressMessages(require(openxlsx))
suppressMessages(require(dplyr))
suppressMessages(require(readr))
suppressMessages(require(stringr))
suppressMessages(require(tidyr))

system("mkdir subseek subseek/clusters subseek/db subseek/results &>/dev/null")
# system("mkdir subseek subseek/segs subseek/clusters subseek/vcf subseek/db subseek/results &>/dev/null")

subseekpath="/ifs/e63data/reis-filho/usr/SubcloneSeeker/"

#setwd("subseek")

#------
# input
#------

cat("reading input\n")
subsets<-read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)
#muts<-read.xlsx("excel/mutation_summary.xlsx",sheet="SNV_HIGH_MODERATE_SUMMARY",check.names=TRUE)
#segfiles<-list.files("facets",pattern="*cncf.txt")
clusterfiles<-list.files("pyclone/tables",pattern="*cluster.tsv")

#--------------------
# data pre-processing
#--------------------

# assign unique mutation IDS
#muts$ID<-paste(muts$CHROM,muts$POS,muts$ANN....GENE,sep="-")

#------------------------------
# main loop over sample subsets
#------------------------------

for (subnum in 1:nrow(subsets)){

	# subsets input processing
	line<-as.vector(subsets[subnum,])
	subsamples<-line[line!=""][-1]
	subname<-line[[1]][1]

	#----------------------------
	# build cluster files
	#----------------------------

	clusters<-read.delim(paste("pyclone/tables/",subname,".cluster.tsv",sep=""))

	clusters %>% 
		select(sample_id,cluster_id,mean) %>%
		arrange(desc(mean)) %>%
		filter(!is.na(mean)) %>%
		spread(sample_id,mean) %>%
		select_(.dots=subsamples) %>%
		replace(is.na(.), 0) ->
		mclusters

		clustertsv=str_c("subseek/clusters/",subname,".cluster.tsv")
		write_tsv(mclusters,clustertsv)	

		clusterdb=str_c("subseek/db/",subname,".cluster.sqlite")
		cmd=str_c(subseekpath,"utils/cluster2db ",clustertsv," ",clusterdb)
		system(cmd)

}

	# for(sample in subsamples){

	# 	cat("building CCF file: samples/",sample,"\n",sep="")
	# 	submuts<-filter(muts,TUMOR_SAMPLE==sample)
	# 	submuts[submuts$CHROM=="X","CHROM"]<-23

	# 	seg<-read_tsv(paste("facets/",grep(paste(sample,"_",sep=""),segfiles,value=TRUE),sep=""))

	# 	# reshape segments table
	# 	rseg <- seg %>% mutate(ID=paste(sample,chrom,loc.start,loc.end,sep=":")) %>% select(ID=ID,Chrom=chrom,StartLoc=loc.start,EndLoc=loc.end,numMark=num.mark,segMean=cnlr.median.clust)

	# 	segtsv=str_c("subseek/segs/",sample,".seg.tsv")
	# 	write_tsv(rseg,segtsv)

	# 	segsq=str_c("db/",sample,".seg.sqlite")
	# 	cmd=str_c(subseekpath,"utils/segtxt2db ",segtsv," ",segsq)
	# 	system(cmd)

	# 	#less alltables/allTN.mutect.dp_ft.som_ad_ft.target_ft.pass.dbsnp.eff.cosmic.nsfp.gene_ann.cn_reg.clinvar.exac_nontcga.chasm.fathmm.tab.high_moderate.novel.txt

	# }

}










