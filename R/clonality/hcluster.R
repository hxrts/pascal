#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

# https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html

# load base libraries
suppressMessages(pacman::p_load(dplyr,readr,tidyr,magrittr,stringr,crayon,mclust))

# create necessary directories
system("mkdir SNPRelate hcluster &>/dev/null")

#------
# input
#------

cat(green("\n-") %+% " reading input\n")
subsets<-read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)

rmrow<-grep("#",subsets[,1])
if(length(rmrow)>0){subsets<-subsets[-rmrow,]}

muts<-read_tsv("summary/muts.tsv")

#subnum=1

#------------------------------
# main loop over sample subsets
#------------------------------

for (subnum in 1:nrow(subsets)){

	# input processing
	line<-as.vector(subsets[subnum,])
	subsamples<-line[line!=""][-1]
	subname<-line[[1]][1]

	cat(blue("\n----------------------------------\n  CLUSTER beginning subset ",subname,"\n----------------------------------\n",sep=""))

	submuts <- subset(muts,TUMOR_SAMPLE%in%subsamples)

	# submuts %>%
	# select(TUMOR_SAMPLE,REF,ALT,CHROM,POS,CCF=cancer.cell.frac) %>%
	# rowwise() %>%
	# mutate(GENO=if(length(CCF)>0){"1|1"}) %>%
	# spread(TUMOR_SAMPLE,GENO,fill=0) %>%
	# arrange(CHROM,POS) %>%
	# mutate(ID=".",QUAL="100",FILTER="PASS",INFO=".",FORMAT="GT",GERM="0|0") %>%
	# select(-CCF) %>%
	# select(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,GERM,everything()) %>%
	# rename_("#CHROM"="CHROM") -> vmuts

	submuts %>%
	select(TUMOR_SAMPLE,CHROM,POS,CCF=cancer.cell.frac) %>%
	rowwise() %>%
	#mutate(GENO=as.numeric(CCF>0)) %>%
	mutate(GENO=CCF) %>%
	spread(TUMOR_SAMPLE,GENO,fill=0) %>%
	select_(.dots=subsamples) %>% as.matrix -> idt
	#kmeans(4,iter.max=100)
	#hc<-hclust(dist(idt))
	pdf(str_c("hcluster/",subname,".hcluster.pdf"))
		plot(as.dendrogram(hclust(dist(idt))))
	dev.off()

	pdf(str_c("hcluster",subname,".mcluster.pdf"))
		#plot(mclustBIC(idt))
		plot(Mclust(idt))
		#plclust(idt)
	dev.off()

	# arrange(CHROM,POS) %>%
	# mutate(ID=".",QUAL="100",FILTER="PASS",INFO=".",FORMAT="GT",GERM="0|0") %>%
	# select(-CCF) %>%
	# select(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,GERM,everything()) %>%
	# rename_("#CHROM"="CHROM") -> vmuts


	# vcfname <- str_c("SNPRelate/vcf/",subname,".vcf")
	# sink(vcfname)
	# 	cat('##fileformat=VCFv4.0
	# 		##source=BCM:SNPTools:hapfuse
	# 		##reference=1000Genomes-NCBI37
	# 		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
	# sink()
	# vmuts %>% write_tsv(path=vcfname,append=TRUE,col_names=TRUE)

	# gdsname<-str_c("SNPRelate/gds/",subname,".gds")
	# snpgdsVCF2GDS(vcfname,gdsname)

	# gdata <- snpgdsOpen(gdsname)
	# snpgdsSummary(gdsname)

	# set.seed(1000)
	# # Try different LD thresholds for sensitivity analysis
	# #snpset <- snpgdsLDpruning(gdata, ld.threshold=0.2)
	# #snpset.id <- unlist(snpset)
	# #pca <- snpgdsPCA(gdata, snp.id=snpset.id, num.thread=2)
	# pca <- snpgdsPCA(gdata, num.thread=2)
	# pc.percent <- pca$varprop*100
	# #head(round(pc.percent, 2))

	# tab <- data.frame(sample.id = pca$sample.id,
 #    EV1 = pca$eigenvect[,1],    # the first eigenvector
 #    EV2 = pca$eigenvect[,2],    # the second eigenvector
 #    stringsAsFactors = FALSE)

	# pdfname<-str_c("SNPRelate/plot/",subname,".pdf")
	# pdf(pdfname)
	# 	plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
	# dev.off()


	# # ibs
	# ibs <- snpgdsIBS(gdata, num.thread=2)

	# #convert to distance
	# hc <- hclust(as.dist(1-ibs$ibs))
	 
	# #obtain the different populations
	# # pop <- read.gdsn(index.gdsn(genofile, path="sample.annot/pop.group"))
	# # pop <- unclass(factor(pop))
	 
	# #define some colours
	# mycolor <- rainbow(2)

	# #colorindex<-lapply(subsamples,function(sample) as.numeric(vmuts[,sample]!=0)) %>% do.call(cbind,.) %>% as.data.frame()  %>% rowSums()

	# pdfname<-str_c("SNPRelate/plot/",subname,".tree.pdf")
	# pdf(pdfname)
	# #plot the dendrogram
	# plot(hc)
	#  ColorDendrogram(hc,branchlength=0.06)
	# dev.off()
	#add a legend
	#legend(x = 260, y = 0.05, fill = my_colour, legend = levels(pop)
}

