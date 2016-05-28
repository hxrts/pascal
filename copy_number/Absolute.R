#---------------
# load libraries
#---------------

library(ABSOLUTE)
library(parallel)
library(xlsx)
#library(snow)

#---------
# run info
#---------

primary.disease="hsarc_ecd1"      # default = NULL
platform="Illumina_WES"           # possible values = "SNP_250K_STY", "SNP_6.0", "Illumina_WES"
lohcut=0.2
ncores=2

#------
# setup
#------

# rewrite seg files with new proper header
dir.create("segs",showWarnings=FALSE)
segfiles<-list.files("../",pattern="*.cncf.txt")
#segfiles<-list.files("segs",pattern="*.cncf.txt")
newHeader<-function(x){
  seg<-read.table(paste("../",x,sep=""),stringsAsFactors=FALSE,header=TRUE,sep="\t")
  #seg<-read.table(x,stringsAsFactors=FALSE,header=TRUE,sep="\t")
  colnames(seg)<-c("ID","Chromosome","Start","End","seg","Num_Probes","nhet","Segment_Mean","mafR","segclust","cnlr.median.clust","mafR.clust","cf","tcn","lcn","cf.em","tcn.em","lcn.em")
  write.table(seg,file=paste("segs/h.",x,sep=""),sep="\t",quote=FALSE,row.names=FALSE,append=FALSE)
  return(seg)
}
segs<-do.call(rbind,lapply(segfiles,function(x) newHeader(x)))
segs$Tumor_Sample<-unlist(lapply(strsplit(segs$ID,"_"),function(x)x[[1]]))

# strip sample suffixes for naming
segsbase<-unlist(lapply(segfiles,function(x)strsplit(x,"_")[[1]][1]))

# build boilerplate params table & write out file
init_params<-data.frame(
  include=rep("Y",length(segfiles)),
  patient=rep(".",length(segfiles)),
  sample=rep(".",length(segfiles)),
  seg.dat.fn=paste("segs/h.",segfiles,sep=""),
  sigma.p=rep("0.01",length(segfiles)),
  max.sigma.h=rep("0.05",length(segfiles)),
  min.ploidy=rep("1.8",length(segfiles)),
  max.ploidy=rep("8",length(segfiles)),
  primary.disease=rep(if(is.null(primary.disease)){system("basename $(cd ../../..; pwd)")}else{primary.disease},length(segfiles)),
  platform=rep(platform,length(segfiles)),
  sample.name=paste(segfiles,".facets",sep=""),
  results.dir=segsbase,
  max.as.seg.count=rep("800",length(segfiles)),
  copy_num_type=rep("total",length(segfiles)),
  max.neg.genome=rep("0.05",length(segfiles)),
  max.non.clonal=rep("0.5",length(segfiles)),
  maf.fn=paste("maf/",segsbase,".maf",sep=""),
  min.mut.af=rep("0",length(segfiles)),
  output.fn.base=paste(segsbase,".facets",sep="")
)

write.table(init_params,file="absolute_params.txt",quote=FALSE,sep="\t") 


#read in SNV & indel events, format tables & combine
snv_hm<-read.xlsx("../../excel/mutation_summary.xlsx",sheetName="SNV_HIGH_MODERATE_SUMMARY",header=TRUE)
snv_sy<-read.xlsx("../../excel/mutation_summary.xlsx",sheetName="SNV_SYNONYMOUS_SUMMARY",header=TRUE)
indel_hm<-read.xlsx("../../excel/mutation_summary.xlsx",sheetName="INDEL_HIGH_MODERATE_SUMMARY",header=TRUE)
indel_hm<-cbind(indel_hm[,1:14],fathmm_pred=rep(".",nrow(indel_hm)),Breast_chasm_score=rep(".",nrow(indel_hm)),indel_hm[,15:21])

eventt<-rbind(snv_hm,snv_sy,indel_hm)

segmatch<-lapply(1:nrow(eventt),function(x)subset(segs,Tumor_Sample==eventt$TUMOR_SAMPLE[[x]]&Chromosome==eventt$CHROM[[x]]))
segmatch<-do.call(rbind,lapply(1:nrow(eventt),function(x)segmatch[[x]][which.min(abs(((segmatch[[x]]$Start+segmatch[[x]]$End)/2)-eventt$POS[[x]])),]))

eventloh<-rep(".",nrow(eventt))
eventloh[segmatch$mafR>lohcut]<-"LOH"

eventt<-data.frame(
  Tumor_Sample_Barcode=eventt$TUMOR_SAMPLE,
  Hugo_Symbol=eventt$ANN....GENE,
  Effect=eventt$ANN....EFFECT,
  Chromosome=as.numeric(replace(replace(as.character(eventt$CHROM),which(as.character(eventt$CHROM)=="X"),"23"),which(replace(as.character(eventt$CHROM),which(as.character(eventt$CHROM)=="X"),"23")=="Y"),"24")),
  Start_position=eventt$POS,
  Reference_Allele=eventt$REF,
  Alternate_Allele=eventt$ALT,
  Mutant_Allele.fraction=eventt$TUMOR.DP/(eventt$NORMAL.DP+eventt$TUMOR.DP),
  Loss.of.Heterozigosity=eventloh,
  Depth.at.mutation..x.=eventt$TUMOR.DP+eventt$NORMAL.DP,
  t_alt_count=eventt$TUMOR.DP,
  t_ref_count=eventt$NORMAL.DP,
  MutationTaster=eventt$dbNSFP_MutationTaster_pred,
  #CHASM..breast.=eventt$,
  #Provean=eventt$,
  Cancer.Gene.Census=eventt$cancer_gene_census,
  X127.Significantly.Mutated.Genes..Kandoth.et.al.=eventt$kandoth,
  Cancer5000.S.Genes..Lawrence.et.al.=eventt$lawrence,
  dbSNP_Val_Status=eventt$dbNSFP_ExAC_Adj_AF
  )

  # NORMAL_SAMPLE
  # ANN....HGVS_P
  # ANN....HGVS_C
  # TUMOR_MAF
  # NORMAL_MAF
  # TUMOR.DP
  # NORMAL.DP
  # dbNSFP_ExAC_Adj_AF
  # fathmm_pred
  # Breast_chasm_score
  # hap_insuf
  # ANN....IMPACT


#lapply(1:nrow(mutlist),function(x) subset(eventt,TUMOR_SAMPLE==mutlist[x,"Tumor.sample"] & CHROM==mutlist[x,"CHROM"] & )


# break event table into per-sample files
dir.create("maf",showWarnings=FALSE)
lapply(unique(eventt$Tumor_Sample_Barcode),function(x){
    samp<-subset(eventt,Tumor_Sample_Barcode==x)
    write.table(samp,paste("maf/",x,".maf",sep=""),sep="\t",quote=FALSE,row.names=FALSE,append=FALSE)
  }
)

run_params <- read.table("absolute_params.txt",stringsAsFactors=FALSE,header=TRUE,sep="\t")
run_params <- run_params[which(run_params$include=="Y"), , drop=FALSE]

#-----------------------------
# ABSOLUTE solution generation
#-----------------------------

runAbsolute<-function(x) {
  RunAbsolute(seg.dat.fn=paste("segs/h.",x$seg.dat.fn,sep=""),
    sigma.p=as.numeric(x$sigma.p),
    max.sigma.h=as.numeric(x$max.sigma.h),
    min.ploidy=as.numeric(x$min.ploidy),
    max.ploidy=as.numeric(x$max.ploidy),
    primary.disease=x$primary.disease,
    platform=as.character(x$platform),
    sample.name=x$results.dir,
    results.dir=x$results.dir,
    max.as.seg.count=as.numeric(x$max.as.seg.count),
    copy_num_type=x$copy_num_type,
    max.neg.genome=as.numeric(x$max.neg.genome),
    max.non.clonal=as.numeric(x$max.non.clonal), 
    maf.fn=paste("maf/",x$maf.fn,sep=""),
    min.mut.af=as.numeric(x$min.mut.af),
    output.fn.base=x$output.fn.base, 
    verbose=TRUE)
}

mclapply(1:nrow(run_params),function(x) runAbsolute(run_params[x,]),mc.cores=ncores)

lapply(1:nrow(run_params),function(x) runAbsolute(run_params[x,]))


#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------


###################################################
###################################################

#-------
# STEP 2
#-------
# Go were the output of the the first step is (you can do this using the "setwd("your output path")" command in R) ###


library(ABSOLUTE)

results1.dir="absolute_out/step_1"
setwd(results1.dir)

#sids<-unlist(lapply(strsplit(absolute.files,"/"),function(x) x[[1]]))

obj.name         = "hsarc_ecd1"
absolute.files   = list.files(path=".",recursive=-TRUE,pattern="*.RData")
indv.results.dir = "../step_2" 
copy_num_type    = "total"



lapply(1:length(absolute.files),function(x){
  load(absolute.files[[x]])
  seg.dat$sample.name<-list.files()[[x]]
  save(seg.dat,file=absolute.files[[x]])
  })



CreateReviewObject(obj.name         = obj.name, 
                   absolute.files   = absolute.files, 
                   indv.results.dir = indv.results.dir, 
                   copy_num_type    = copy_num_type, 
                   plot.modes       = TRUE, 
                   verbose          = TRUE)

# add "solution" as first column

#-------
# STEP 3
#-------

calls.path       = file.path("../step_2", paste(obj.name, ".PP-calls_tab.txt", sep=""))
modes.path       = file.path("../step_2", paste(obj.name, ".PP-modes.data.RData", sep=""))
output.path      = file.path("../step_3", "Output")

ExtractReviewedResults(calls.path,"hxrts",modes.path,output.path,"absolute","total",verbose=TRUE)

# will create abs maf files with 3 cols: cancer cell frac, pr_somatic_clonal, 95ci_low & make a 4th col. if pr_somatic_clonal is >= 0.5 and ci95low is >= 0.9 then mark as 'clonal' else mark as 'subclonal'

res3t_path="../step_3/Output/reviewed/absolute.hxrts.ABSOLUTE.table.txt"

res3t<-read.delim(res3t_path,header=TRUE,sep="\t",stringsAsFactors=FALSE)
#res3t<-res3t[,2:12]
#row.names(result3_table)<-result3_table$sample

#result3_table_wmuts<-result3_table$sample[result3_table$sample%in%unlist(lapply(strsplit(list.files("absolute_results3/Output/reviewed/SEG_MAF/",pattern="MAF"),"_"),function(x)x[[1]]))]


#samples <- read.xlsx("../use_samples.xlsx", "pheno", stringsAsFactors=F)
#samples<-samples$SAMPLE

clonality<-lapply(list.files("absolute_out/step_1"),function(x){
  file=paste("absolute_out/step_3/Output/reviewed/SEG_MAF/",x,".facets_ABS_MAF.txt",sep="")
  res3s<-read.delim(file,header=TRUE,sep="\t",stringsAsFactors=FALSE)
  cloneb<-res3s$Pr_somatic_clonal>=0.5|res3s$ccf_CI95_low>=0.9
  clone<-rep("subclonal",length(cloneb))
  clone[cloneb]<-"clonal"
  prob<-res3s$Pr_somatic_clonal
  pr95<-res3s$ccf_CI95_low
  sample<-res3s$sample
  ccf<-res3s$cancer_cell_frac
  hugo_symbol<-res3s$Hugo_Symbol
  data.frame(id=paste(sample,hugo_symbol,sep="_"),sample=sample,hugo_symbol=hugo_symbol,clone=clone,prob=prob,pr95=pr95,ccf=ccf)
})

clonality<-do.call(rbind,clonality)

##unlist(lapply(clonality,function(x) x$sample[[1]]))

# result3_table[result3_table_wmuts,"clonality"]<-clone
# result3_table[result3_table_wmuts,"prob"]<-clonality$prob
# result3_table[result3_table_wmuts,"lb"]<-clonality$lb

# write results out to clonality.txt
write.table(data.frame(sample=clonality$sample,gene=clonality$hugo_symbol,clone=clonality$clone,prob=clonality$prob,pr95=clonality$pr95,ccf=clonality$ccf),"clonality.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)



#######################################################################################################
###################################### MATRIX PLOTTING ################################################
#######################################################################################################

suptable <- read.xlsx("../supplementary_table_4.xlsx", "Sheet1", stringsAsFactors=F,header=TRUE,startRow=2)
lum <- read.xlsx("../20150921_pheno.xlsx", "pheno", stringsAsFactors=F)
LumA<-sort(subset(lum,SAMPLE.TYPE=="LumA")$SAMPLE)
LumB<-sort(subset(lum,SAMPLE.TYPE=="LumB")$SAMPLE)



clonet<-data.frame(matrix(0, nrow=length(unique(clonality$hugo_symbol)), ncol=length(c(LumA,LumB))))
row.names(clonet)<-names(rev(sort(table(clonality$hugo_symbol))))
samplenv<-paste("MB-T",substr(c(LumA,LumB),4,6),sep="") #   [!c(LumA,LumB)%in%c("MBT41B","MBT65","MBT68","MBT82","MBT92","MBT06","MBT38","MBT62","MBT63")],4,6),sep="") # missing 9 samples here
colnames(clonet)<-samplenv
colnames(clonet)<-substr(samplenv,1,6)

clonev<-rep(0,length(clonality$ccf))
clonev[which(clonality$ccf>0&clonality$ccf<=0.2)]<-1
clonev[which(clonality$ccf>0.2&clonev<=0.4)]<-2
clonev[which(clonality$ccf>0.4&clonality$ccf<=0.6)]<-3
clonev[which(clonality$ccf>0.6&clonality$ccf<=0.8)]<-4
clonev[which(clonality$ccf>0.8&clonality$ccf<=1)]<-5


lapply(1:nrow(clonality),function(x) clonet[which(row.names(clonet)==clonality$hugo_symbol[x]),which(colnames(clonet)==clonality$sample[x])]<<-clonev[x])


# fillm<-lapply(1:nrow(clonality),function(x) data.frame(i=which(row.names(clonet)==clonality$hugo_symbol[x]), j=which(colnames(clonet)==clonality$sample[x]), k=clonality[x,"ccf"]))
# fillm<-do.call(rbind,fillm)

# clonet[cbind(fillm$i, fillm$j)]<-fillm$k



#-------------------------
# amp matrix plot function
#-------------------------

plotMatrix<-function(D,subtype=NA,xgroup,plotpos=1,ptitle=NA,...){

  # plot variable setup

  min=0
  max=5
  yLabels<-rownames(D)
  xLabels<-colnames(D)
  D<-as.matrix(D)

  # check for additional function arguments
  if(length(list(...))){
    Lst<-list(...)
    if( !is.null(Lst$zlim) ){
       min<-Lst$zlim[1]
       max<-Lst$zlim[2]
    }
    if(!is.null(Lst$yLabels)){
       yLabels<-c(Lst$yLabels)
    }
    if(!is.null(Lst$xLabels)){
       xLabels<-c(Lst$xLabels)
    }
    if(!is.null(Lst$title)){
       title<-Lst$title
    }
  }

  # check for null values
  if(is.null(xLabels)){
    xLabels<-c(1:ncol(D))
  }
  if(is.null(yLabels)){
    yLabels<-c(1:nrow(D))
  }

  # build color palette
  colorRamp<-c("#DEDEDE",colorRampPalette(c("#9C8FDB","#310170"))(5))

  # reverse Y axis
  reverse<-nrow(D):1
  yLabels<-yLabels[reverse]
  D<-D[reverse,]

  # data Map

  image(1:length(xLabels),1:length(yLabels),z=t(D),col=colorRamp,xlab="",ylab="",axes=FALSE,zlim=c(min,max))


  # add plot title if desired
  # if(!is.na(ptitle)){
  mtext("Male BC cancer mutation cell fractions",side=3,at=25,cex=6,line=2,padj=0,las=1)


  # add sample labels to left
  axis(2,at=1:length(yLabels),labels=yLabels,las=1,cex.axis=2.6)
  axis(1,at=1:length(xLabels),labels=xLabels,las=2,cex.axis=2.6)

  # render grid lines
  grid(nx=ncol(D),ny=nrow(D),col="white",lty=1,lwd=4,equilogs=TRUE)

  # add chromosome demarcations

  mtext(xgroup,side=1,at=c(6.5,32),cex=4,line=2,padj=4.8,las=1)

  clones<-subset(clonality,clone=="clonal")
  lapply(1:nrow(clones),function(x){
    yy<-which(yLabels==clones$hugo_symbol[x])
  xx<-which(xLabels==clones$sample[x])
    rect(xx-0.5,yy-0.5,xx+0.5,yy+0.5,fill=NULL,border="goldenrod1",lwd=10)
  })
  abline(v=17.5,col="black",lty=1,lwd=10)

  #rect(xleft, ybottom, xright, ytop, density = NULL, angle = 45, col = NA, border = NULL, lty = par("lty"), lwd = par("lwd"))

  # add box around each plot
  box("plot",lwd=3)
}

#--------------------------------
# amp matrix legend plot function
#--------------------------------

plotScale<-function(pmar){
  par(mar=pmar)  # mar/oma - bottom, left, top, right   MALE
  colorRamp<-c("#FFFFFF","#FFFFFF","#DEDEDE",colorRampPalette(c("#9C8FDB","#310170"))(5))
  colorLevels <- 1:8
  image(1,colorLevels,matrix(data=colorLevels, ncol=length(colorLevels),nrow=1),col=colorRamp,xlab="",ylab="",xaxt="n",axes=FALSE)
  mtext("Cancer cell fraction",side=3,at=-0.1,cex=3,line=3,padj=0,las=1)
  axis(2,labels=c("Clonal mutation","","0%",">0%-20%",">20%-40%",">40%-60%",">60%-80%",">80%-100%"),at=1:8,las=1,cex.axis=2.6,tic=FALSE)
  grid(nx=1,ny=8,col="white",lty=1,lwd=4,equilogs=TRUE)
  rect(0.65,0.65,1.35,1.35,border="goldenrod1",lwd=10)
}
 
#---------------------
# amp matrix execution
#---------------------

matrixpar=list(mfrow=c(1,2),mar=c(0,0,2,16),oma=c(17,14,7,3))  # mar/oma - bottom, left, top, right
pdf("ccf_matrix.pdf",width=36, height=60)
  par(matrixpar)
  layout(matrix(1:2,ncol=2), widths=c(15,1), heights=c(30,1))
  plotMatrix(clonet,xgroup=c("Luminal A-like","Luminal B-like"),plotpos=0)
  plotScale(pmar=c(250,3,2,2))
dev.off()


fintable<-clonality[,c("sample","hugo_symbol","ccf","prob","pr95","clone")]
fintable$group<-c("LumA",rep("",37),"LumB",rep("",153))
fintable<-fintable[,c("group","sample","hugo_symbol","ccf","prob","pr95","clone")]
fintable$prob[which(fintable$prob<0.005)]<-0
fintable$pr95[which(fintable$pr95<0.005)]<-0
fintable$ccf[which(fintable$ccf<0.005)]<-0
fintable$sample<-as.character(fintable$sample)
fintable$sample[unlist(lapply(unique(fintable$sample),function(x)which(x==fintable$sample)[-1]))]<-rep("",length(unlist(lapply(unique(fintable$sample),function(x)which(x==fintable$sample)[-1]))))
colnames(fintable)<-c("Group","Sample","Gene","Cancer Cell Fraction (CCF) (ABSOLUTE)","P(Clonal)","Lower bound of 95% confidence interval","Clonal/Subclonal mutation")
stargazer(fintable, type = "html",summary=FALSE,digits.extra=0,column.labels=colnames(fintable),digits=2,rownames=FALSE,style="all")

