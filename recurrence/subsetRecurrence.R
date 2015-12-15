#!/usr/bin/env Rscript

#------
# input
#------

subsets<-read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)
mutfile<-read.delim("recurrent_mutations/recurrent_mutations.tsv",sep="\t",stringsAsFactors=FALSE)

#---------------------
# function definitions
#---------------------

# create mutation count data.frame
buildD<-function(samples){

	Dl<-lapply(samples,function(sample)
		unlist(lapply(rev(names(sort(table(muts[which(muts$TUMOR_SAMPLE%in%samples),]$ANN....GENE_SPLIT	)))),function(gene)
			nrow(subset(muts[which(muts$TUMOR_SAMPLE%in%samples),],TUMOR_SAMPLE==sample&ANN....GENE_SPLIT==gene)))
			)
		)

	D<-data.frame(do.call(cbind,Dl))
	row.names(D)<-rev(names(sort(table(muts[which(muts$TUMOR_SAMPLE%in%samples),]$ANN....GENE_SPLIT))))
	colnames(D)<-samples
	return(D)
}

# draw heatmap
heatM<-function(D,file){
	colorRamp<-colorRampPalette(c("white","blue"))(max(D)+1)
	matrixpar=list(mfrow=c(2,1),mar=c(1,3,1,1),oma=c(2,2,2,2))  # mar/oma - bottom, left, top, right
	pdf(file)
		par(matrixpar)
		layout(matrix(1:2,ncol=2), widths=c(11,1), heights=c(1,1))

		# plot
		image(1:ncol(D),1:nrow(D),z=t(as.matrix(D[nrow(D):1,])),col=colorRamp,xlab="",ylab="",axes=FALSE)
		title(main=unlist(strsplit(file,".pdf")),line=0.5,cex.main=0.8)
		axis(1,at=ncol(D):1,labels=names(D),cex.axis=0.6,las=2)			# x-axis
		axis(2,at=nrow(D):1,labels=row.names(D),cex.axis=0.45,las=1)	# y-axis
		box("plot")

		# legend
		par(mar=c(24,1,1,0))  # mar/oma - bottom, left, top, right
		colorLevels<-min(D):max(D)
		image(1,colorLevels,matrix(data=colorLevels,ncol=length(colorLevels),nrow=1),col=colorRamp,xlab="",ylab="",xaxt="n",axes=FALSE)
		axis(2,at=min(D):max(D),labels=min(D):max(D),las=1,cex.axis=0.6)
		box("plot")

	dev.off()
}

#----------
# main loop
#----------

for (subnum in 1:nrow(subsets)){

	line<-as.vector(subsets[subnum,])
	subset<-line[line!=""][-1]
	name<-line[[1]][1]

	muts<-mutfile[which(mutfile$TUMOR_SAMPLE%in%subset),]

	write.table(muts,file=paste("recurrent_mutations/",name,"_recurrent.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)

	heat<-buildD(subset)
	heatM(heat,paste(name,"_recurrent.pdf",sep=""))

}


