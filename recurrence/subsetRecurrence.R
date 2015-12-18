#!/usr/bin/env Rscript

#------
# input
#------

subsets<-read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)
mutfile<-read.delim("recurrent_mutations/recurrent_mutations.tsv",sep="\t",stringsAsFactors=FALSE)

runSufam<-commandArgs(TRUE)[1]=="TRUE"

#---------------------
# function definitions
#---------------------

# create mutation count data.frame
buildD<-function(muts,samples,groupEventsByGene){

	if(groupEventsByGene==TRUE){
		Dl<-lapply(samples,function(sample)
			unlist(lapply(rev(names(sort(table(muts[which(muts$TUMOR_SAMPLE%in%samples),]$ANN....GENE_SPLIT)))),function(gene)
				nrow(subset(muts[which(muts$TUMOR_SAMPLE%in%samples),],TUMOR_SAMPLE==sample&ANN....GENE_SPLIT==gene)))
				)
			)
		D<-data.frame(do.call(cbind,Dl))
		row.names(D)<-rev(names(sort(table(muts[which(muts$TUMOR_SAMPLE%in%samples),]$ANN....GENE_SPLIT))))
	}
	if(groupEventsByGene==FALSE){
		Dl<-lapply(samples,function(sample)
			unlist(lapply(rev(names(sort(table(muts[which(muts$TUMOR_SAMPLE%in%samples),]$ID)))),function(id)
				nrow(subset(muts[which(muts$TUMOR_SAMPLE%in%samples),],TUMOR_SAMPLE==sample&ID==id)))
				)
			)
		D<-data.frame(do.call(cbind,Dl))
		row.names(D)<-rev(names(sort(table(muts[which(muts$TUMOR_SAMPLE%in%samples),]$ID))))
	}

	colnames(D)<-samples
	return(D)
}

# draw heatmap
heatM<-function(D,name,groupEventsByGene){

	if(groupEventsByGene==TRUE){
		pdfname<-paste("recurrent_mutations/",name,"_recurrent.pdf",sep="")
	}else{
		pdfname<-paste("recurrent_mutations/",name,"_eventsGroupedByGene_recurrent_test.pdf",sep="")
	}

	colorRamp<-colorRampPalette(c("white","blue"))(max(D)+1)
	matrixpar=list(mfrow=c(2,1),mar=c(6,15,3,1),oma=c(1,1,1,1))  # mar/oma - bottom, left, top, right

	pdf(pdfname,width=10,height=22)
		par(matrixpar)
		layout(matrix(1:2,ncol=2), widths=c(11,1), heights=c(1,1))

		# plot
		image(1:ncol(D),1:nrow(D),z=t(as.matrix(D[nrow(D):1,])),col=colorRamp,xlab="",ylab="",axes=FALSE)
		title(main=name,line=0.7,cex.main=2)
		axis(1,at=ncol(D):1,labels=rev(names(D)),cex.axis=1.5,las=2)	# x-axis
		axis(2,at=nrow(D):1,labels=row.names(D),cex.axis=1.1,las=1)	# y-axis
		grid(nx=ncol(D),ny=nrow(D),col="gray92",lty=1,lwd=par("lwd"),equilogs=TRUE)
		box("plot")

		# legend
		par(mar=c(99,3,3,0))  # mar/oma - bottom, left, top, right
		colorLevels<-min(D):max(D)
		image(1,colorLevels,matrix(data=colorLevels,ncol=length(colorLevels),nrow=1),col=colorRamp,xlab="",ylab="",xaxt="n",axes=FALSE)
		axis(2,at=min(D):max(D),labels=min(D):max(D),las=1,cex.axis=1.5)
		grid(nx=1,ny=length(colorLevels),col="gray92",lty=1,lwd=par("lwd"),equilogs=TRUE)
		box("plot")

	dev.off()
}

# draw heatmap
heatMS<-function(D,name){

	pdfname<-paste("recurrent_mutations/",name,"_sufamRecurrent.pdf",sep="")

	E<-as.data.frame(do.call(cbind,lapply(1:ncol(D),function(samplenum){
		sample<-names(D)[samplenum]
		sufam<-read.delim(paste("recurrent_mutations/",sample,"_sufamRecurrent.tsv",sep=""),sep="\t",stringsAsFactors=FALSE)

		Dpos<-as.data.frame(do.call(rbind,strsplit(row.names(D),":")))
		vf<-apply(Dpos,1,function(x) subset(sufam,chrom==x[2]&pos==x[3])[1,"val_maf"])
		return(vf)
	})))
	row.names(E)<-row.names(D)
	colnames(E)<-names(D)



	colorRamp<-colorRampPalette(c("white","blue"))(8)
	matrixpar=list(mfrow=c(2,1),mar=c(6,15,3,1),oma=c(1,1,1,1))  # mar/oma - bottom, left, top, right

	pdf(pdfname,width=10,height=22)
		par(matrixpar)
		layout(matrix(1:2,ncol=2), widths=c(11,1), heights=c(1,1))

		# plot
		image(1:ncol(E),1:nrow(E),z=t(as.matrix(E[nrow(E):1,])),col=colorRamp,xlab="",ylab="",axes=FALSE)
		title(main=name,line=0.7,cex.main=2)
		axis(1,at=ncol(E):1,labels=rev(names(E)),cex.axis=1.5,las=2)	# x-axis
		axis(2,at=nrow(E):1,labels=row.names(E),cex.axis=1.1,las=1)	# y-axis
		grid(nx=ncol(E),ny=nrow(E),col="gray92",lty=1,lwd=par("lwd"),equilogs=TRUE)
		box("plot")

		for(x in 1:ncol(D)){
			for(y in 1:nrow(D)){
				if(D[y,x]==1){
					rect(xleft=x-0.5,ybottom=nrow(D)-y+0.5,xright=x+0.5,ytop=nrow(D)-y+1.5,density=6)
				}
			}
		}

		# legend
		par(mar=c(99,3,3,0))  # mar/oma - bottom, left, top, right
		colorLevels<-1:8
		image(1,colorLevels,matrix(data=colorLevels,ncol=length(colorLevels),nrow=1),col=colorRamp,xlab="",ylab="",xaxt="n",axes=FALSE)
		axis(2,at=c(1,8),labels=c("0",round(max(E),2)),las=1,cex.axis=1)
		grid(nx=1,ny=length(colorLevels),col="gray92",lty=1,lwd=par("lwd"),equilogs=TRUE)
		box("plot")

	dev.off()


}


#----------
# main loop
#----------

for (subnum in 1:nrow(subsets)){

	# input processing
	line<-as.vector(subsets[subnum,])
	samples<-line[line!=""][-1]
	name<-line[[1]][1]
	muts<-mutfile[which(mutfile$TUMOR_SAMPLE%in%samples),]

	# write out tables
	write.table(muts,file=paste("recurrent_mutations/",name,"_recurrent.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)

	# create unique mut IDs
	muts$ID<-apply(muts[,c("CHROM","POS","ANN....GENE_SPLIT")],1,function(x)paste(strip(x[c(3,1,2)]),collapse=":"))

	# group events by gene
	heat<-buildD(muts,samples,groupEventsByGene=TRUE)
	heatM(heat,name,groupEventsByGene=TRUE)

	# stratify events within the same gene
	heat<-buildD(muts,samples,groupEventsByGene=FALSE)
	heatM(heat,name,groupEventsByGene=FALSE)

	# sufam variant count coloring
	if(runSufam==TRUE){
		heat<-buildD(muts,samples,groupEventsByGene=FALSE)
		heatMS(heat,name)
	}

}

cat("* complete\n")