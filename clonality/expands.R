#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

#suppressMessages(require(openxlsx))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(crayon))
suppressMessages(library(expands))
suppressMessages(library(RColorBrewer))
suppressMessages(library(colorspace))

# create necessary directories
suppressMessages(system("mkdir expands &>/dev/null"))

# loop defaults for development
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

setwd("expands")

for (subnum in 1:nrow(subsets)){

	# input processing
	line<-as.vector(subsets[subnum,])
	subsamples<-line[line!=""][-1]
	subname<-line[[1]][1]

	cat(blue("\n---------------------------------\n  EXPANDS beginning subset ",subname,"\n---------------------------------\n",sep=""))

	suppressMessages(system(str_c("mkdir ",subname," &>/dev/null")))
	setwd(subname)

	exolist=list()
	for (samplenum in 1:length(subsamples)) {

		sample=subsamples[samplenum]

		#-------------------
		# build event tables
		#-------------------

		cat(green("\n-") %+% " building event tables\n")

		# build SNV df
		muts %>%
			filter(TUMOR_SAMPLE==sample) %>%
			mutate(PN_B=0) %>%
			select(chr=CHROM,startpos=POS,AF_Tumor=cancer.cell.frac,PN_B) %>%
			mutate(chr=replace(chr,chr=="X",23)) ->
			SNV

		# build CBS df
		seg<-read_tsv(str_c("../../facets/",grep(str_c(sample,"_"),segfiles,value=TRUE)))
		seg %>%
			mutate(segmentLength=1+loc.end-loc.start) %>%
			mutate(TUMOR_SAMPLE=unlist(lapply(ID,function(x)strsplit(x,"_")))[1]) %>%
			filter(TUMOR_SAMPLE==sample) %>%
			select(chr=chrom,startpos=loc.start,endpos=loc.end,CN_Estimate=lcn.em,segmentLength) %>%
			mutate(CN_Estimate=replace(CN_Estimate,which(is.na(CN_Estimate)),1)) %>%
			mutate(CN_Estimate=round(CN_Estimate)) ->
			CBS

		CBS %>% write_tsv(str_c(sample,".",subname,".cbs"))

		#----------------------
		# run expands algorithm
		#----------------------

		exo<-runExPANdS(data.matrix(SNV),data.matrix(CBS),precision=0.01,snvF=str_c(sample,".",subname))
		exolist<-c(exolist,list(exo))

		#---------
		# plotting
		#---------

		cat(green("\n-") %+% " plotting sample " %+% sample %+% "\n")

		pdf(str_c(sample,".subpop.pdf"))
			plotSPs(exo$dm,sample,cex=1)
		dev.off()

		pdf(str_c(sample,".tree.pdf"))
			plot(exo$tree,cex=2.5)
		dev.off()

	}

	setwd("..")

	#----------------------
	# phylogeny render prep
	#----------------------

	# list CBS/SBS files
	CBSfiles<-str_c(subname,"/",list.files(subname,pattern="*.cbs"))
	SPSfiles<-str_c(subname,"/",list.files(subname,pattern="*.sps"))

	# build a sample group per patient to calculate combined phylogeny
	samplegroup=list(cbs=sort(CBSfiles,decreasing=TRUE),sps=sort(SPSfiles,decreasing=TRUE),labels=sort(subsamples,decreasing=TRUE))

	colorize<-function(){
	# color tree tips according to sample origin
	colmap=brewer.pal(length(samplegroup$labels),"Paired")
	colors<-rep(colmap[1],each=length(subphylo$tip.label))
	for (i in 1:length(samplegroup$labels)){
		ii=grep(samplegroup$labels[[i]],subphylo$tip.label)
		colors[ii]=colmap[i]}
	}

	#---------------------------------------
	# render tree meta-population resolution
	#---------------------------------------

	subphylo=buildMultiSamplePhylo(samplegroup,str_c(subname,"/",subname,".meta"),ambigSg=FALSE,plotF=0,spRes=FALSE)

	colorize()

	# plot rooted inter-sample phylogeny
	pdf(str_c(subname,"/",subname,".meta.rooted.pdf"))
		plot(subphylo,cex=1.4)
	dev.off()

	# plot unrooted inter-sample phylogeny
	pdf(str_c(subname,"/",subname,".meta.unrooted.pdf"))
		plot(subphylo, cex = 1.4, type="u")
	dev.off()

	#--------------------------------------
	# render tree sub-population resolution
	#--------------------------------------

	subphylo=buildMultiSamplePhylo(samplegroup,str_c(subname,"/",subname),ambigSg=F,plotF=0)

	colorize()

	# plot rooted inter-sample phylogeny
	pdf(str_c(subname,"/",subname,".sub.rooted.pdf"))
		plot(subphylo,cex=1.4)
	dev.off()

	# plot unrooted inter-sample phylogeny
	pdf(str_c(subname,"/",subname,".sub.unrooted.pdf"))
		plot(subphylo, cex = 1.4, type="u")
	dev.off()

	# tree<-root(drop.tip(subphylo,4),1)
	# cls <- list(c1=c("P10_SP_0.66"),
	#             c2=c("M10_SP_0.66"),
	#             c3=c("FM10_SP_0.83"))
	# tree <- groupOTU(tree, cls)

	# pdf(str_c(subname,".pdf"))
	# 	ggtree(tree, aes(color=group)) + geom_text(size=4.5,aes(label=label),vjust = -1) + scale_x_continuous(expand = c(.1, .1))
	# dev.off()

}


