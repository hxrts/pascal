#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

#suppressMessages(require(openxlsx))
suppressMessages(require(dplyr))
suppressMessages(require(readr))
suppressMessages(require(stringr))
suppressMessages(require(tidyr))
suppressMessages(require(crayon))

system("mkdir subseek subseek/input subseek/subclones subseek/clusters subseek/dot &>/dev/null")
# system("mkdir subseek subseek/segs subseek/clusters subseek/vcf subseek/db subseek/results &>/dev/null")

subseekpath="/ifs/e63data/reis-filho/usr/SubcloneSeeker/"

#setwd("subseek")

#------
# input
#------

cat("reading input\n")
subsets<-read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)

rmrow<-grep("#",subsets[,1])
if(length(rmrow)>0){subsets<-subsets[-rmrow,]}

#muts<-read.xlsx("excel/mutation_summary.xlsx",sheet="SNV_HIGH_MODERATE_SUMMARY",check.names=TRUE)
#segfiles<-list.files("facets",pattern="*cncf.txt")
clusterfiles<-list.files("pyclone/~archive/10000_steps/tables",pattern="*cluster.tsv")

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

	cat(blue("\n--------------------------\n  beginning subset ",subname,"\n--------------------------\n\n",sep=""))

	#----------------------------
	# build cluster files
	#----------------------------

	clusters<-read.delim(paste("pyclone/~archive/10000_steps/tables/",subname,".cluster.tsv",sep=""))

	clusters %>% 
		select(sample_id,cluster_id,mean) %>%
		arrange(desc(mean)) %>%
		filter(!is.na(mean)) %>%
		spread(sample_id,mean) %>%
		select_(.dots=subsamples) %>%
		replace(is.na(.), 0) %>%
		format(scientific = FALSE) ->
		mclusters

		# write cluster tsv files
		clustertsv=str_c("subseek/input/",subname,".cluster.tsv")
		write_tsv(mclusters,clustertsv,col_names=FALSE)	

		# cluster2db
		cmd=str_c(subseekpath,"utils/cluster2db ",clustertsv)
		system(cmd)

		# move primary db files to clusters dir
		priclust=str_c("subseek/clusters/",subname,".cluster.pri.sqlite")
		cmd=str_c("mv ",clustertsv,"-pri.sqlite ",priclust)
		system(cmd)

		# move relapse db files to clusters dir
		relclust=str_c("subseek/clusters/",subname,".cluster.rel.sqlite")
		cmd=str_c("mv ",clustertsv,"-rel.sqlite ",relclust)
		system(cmd)

		# run ssmain on primary
		prisub=str_c("subseek/subclones/",subname,".sub.pri.sqlite")
		cmd=str_c(subseekpath,"utils/ssmain ",priclust," ",prisub," >/dev/null 2>&1")
		system(cmd)

		# run ssmain on relapse
		relsub=str_c("subseek/subclones/",subname,".sub.rel.sqlite")
		cmd=str_c(subseekpath,"utils/ssmain ",relclust," ",relsub," >/dev/null 2>&1")
		system(cmd)

		# run treemerge on primary + relapse
		cat(green("\n-") %+% " run treemerge on primary + relapse\n")
		cmd=str_c(subseekpath,"utils/treemerge ",prisub," ",relsub)
		system(cmd)

		# treeprint list
		cat(green("\n-") %+% " treeprint list\n")
		cmd=str_c(subseekpath,"utils/treeprint -l ",prisub)
		system(cmd)

		# # treeprint struct
		# cat(green("\n-") %+% " treeprint struct\n")
		# clusternum=2
		# cmd=str_c(subseekpath,"utils/treeprint -r clusternum ",prisub)
		# system(cmd)

		# # create dot file
		# cat(green("\n-") %+% " create dot file\n")
		# dotfile=str_c("subseek/dot/",subname,".dot")
		# cmd=str_c(subseekpath,"utils/treeprint -r clusternum -g ",prisub," | tee ",dotfile)
		# system(cmd)

		# # graphviz graph from dot
		# cat(green("\n-") %+% " graphviz graph from dot\n")
		# cmd=str_c("dot -Tpng -O ",dotfile)
		# system(cmd)

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

#}










