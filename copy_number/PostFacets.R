#!/usr/bin/env Rscript

# parser <- OptionParser(usage = '%prog [geneCN file] [output file] [max sequential gene amp/del]')

# arguments <- parse_args(parser, positional_arguments = TRUE)
# opt <- arguments$options

# if (length(arguments$args) != 2) {
# 	#cat('Need one geneCN file and one output plot file\n')
# 	print_help(parser)
# 	stop()
# } else {
# 	geneCN <- arguments$args[1]
# 	outFile <- arguments$args[2]
# }

suppressMessages(pacman::p_load(dplyr,readr,stringr,tidyr,broom,purrr,magrittr,rlist,crayon,colorspace,ggplot2,grid,gridExtra,RColorBrewer))

cncf.path <- '../../facets.clean/'
cnv.file <- '../../facets.clean/geneCN.txt'
plot.file <- 'geneCN.pdf'

cnv.matrix <- read.delim(cnv.file, sep='\t', header=T, stringsAsFactors=F) %>% tbl_df
#plot_heatmap(cnv.matrix,plot.file)


cnv.matrix %<>%
	mutate(chrom=
		ifelse(chrom=='X',23,
		ifelse(chrom=='Y',23,chrom))) %>%
	mutate(chrom=as.integer(chrom))

chrom <- cnv.matrix$chrom

cnv.matrix %<>% arrange(chrom,start,end)

cnv.matrix[is.na(cnv.matrix)] <- 3	# replace NA with numeric for use in rle function
breaks <- c(0,cumsum(table(chrom)[chrom %>% table %>% names %>% as.numeric %>% order])) # start-1 == end

for (sample.name in names(cnv.matrix) %>% list.filter(!. %in% c('chrom','start','end','hgnc','band'))){

	cat('\n * sample:',sample.name)

	#loop over chromosomes in junction table
	for(chromosome in 1:length(cnv.matrix$chrom %>% unique)){

		#cat('\n   chr',formatC(chromosome, width=2, flag='0'),': ',sep='')
		start <- breaks[chromosome]+1
		end <- breaks[chromosome+1]
		chbit <- cnv.matrix[start:end,sample.name] %>% unlist %>% rle

		# replace an entirely empty chromosome with calls of 0
		if((chbit$lengths %>% length)==1){
			cnv.matrix[start:end,sample.name] <- 0
		}else{
			# extend chromosome start calls from nearest integer on same chromosome
			if(cnv.matrix[start,sample.name]==3){
				#cat('start..')
				NAbottom <- min(which(chbit$values==3))
				Ibottom <- start + cumsum(chbit$lengths)[NAbottom]
				cnv.matrix[start:(Ibottom-1),sample.name] <- cnv.matrix[Ibottom,sample.name]
			}
			# extend chromosome end calls from nearest integer on same chromosome
			if(cnv.matrix[end,sample.name]==3){
				#cat('end..')
				NAtop <- min(which(rev(chbit$values==3)))
				Itop <- end - cumsum(rev(chbit$lengths))[NAtop]
				cnv.matrix[(Itop+1):end,sample.name] <- cnv.matrix[Itop,sample.name]
			}
		}

		nav <- which(chbit$values[-c(1,length(chbit$values))]==3)+1
		if(length(nav>0)){
			#cat('gaps..')
			# find nearest neighbors of inner-chromosome gaps
			gstarts <- cumsum(chbit$lengths)[nav-1]+start
			glengths <- chbit$lengths[nav]

			nrows <- 
				gstarts %>%
				map2(glengths, ~ .x:(.x+.y-1)) %>%
				unlist %>%
				unname

			nfills <-
					cnv.matrix %>%
					select(chrom,start,end,hgnc,get(sample.name)) %>%
					tbl_df %>%
					slice(nrows) %>%
					mutate(mid=rowMeans(.[,c('start','end')])) %>%
					rowwise %>%
					inner_join(
						cnv.matrix %>%
						filter(chrom==chromosome) %>%
						mutate(qmid=rowMeans(.[,c('start','end')])) %>% 
						select(chrom,qmid,fill=get(sample.name)),
						by='chrom'
						) %>%
					ungroup %>%
					group_by(chrom,start,end,hgnc) %>%
					filter(fill!=3) %>%
					slice(which.min(abs(mid-qmid))) %>%
					.$fill

			# write neighbor calls to master table
			if(length(nfills)>0){
				cnv.matrix[nrows,sample.name] <- nfills
			}else{
				cnv.matrix[nrows,sample.name] <- NA
			}
		}
	} # loop over chromosomes

	#cat('\n')
}

cnv.matrix$chrom <- chrom
cnv.matrix %<>% mutate(mid=(start+end)/2) %>% select(chrom,start,mid,end,everything())

# pre.cn.calls <- lapply(foo,function(sample.name) cnv.matrix[,sample.name] %>% unlist %>% unname %>% rle )


OverCall <- function(sample.name){

	sample.rle <- cnv.matrix[,sample.name] %>% unlist %>% unname %>% rle

	mark.amp <- intersect(which(sample.rle$lengths > 100),which(sample.rle$values == 2))
	loc.amp <- mapply(function(x,y) x:y ,c(0,cumsum(sample.rle$lengths))[mark.amp]+1 , cumsum(sample.rle$lengths)[mark.amp]) %>% unlist %>% as.vector

	cncf <-
		read.delim(str_c(cncf.path,'/',sample.name %>% strsplit('_') %>% unlist %>% head(2) %>% paste(collapse='_'),'.cncf.txt'),sep='\t',stringsAsFactors=FALSE) %>%
		tbl_df %>%
		mutate(loc.mid=(loc.start+loc.end)/2) # add genome-wide stats

	logic.amp <-
		cnv.matrix[loc.amp,c('chrom','mid')] %>%
		add_rownames(var='index') %>%
		mutate(index=as.numeric(index)) %>%
		group_by(index) %>%
		left_join(cncf,by='chrom') %>%
		slice(which.min(abs(mid-loc.mid))) %>%
		filter(cnlr.median>1.1) %>%
		.$index

	loc.amp <- loc.amp[logic.amp]
	cnv.matrix[loc.amp,sample.name] <<- 1

	# summary stats, need to be corrected for log odds ratio removal
	# data.frame(
	# 	amp.starts=sample.rle$lengths[mark.amp],
	# 	amp.lengths=c(0,cumsum(sample.rle$lengths))[mark.amp]+1
	# )

	mark.del <- intersect(which(sample.rle$lengths > 100),which(sample.rle$values == -2))
	loc.del <- mapply(function(x,y) x:y ,c(0,cumsum(sample.rle$lengths))[mark.del]+1 , cumsum(sample.rle$lengths)[mark.del]) %>% unlist %>% as.vector

	logic.del <-
		cnv.matrix[loc.del,c('chrom','mid')] %>%
		add_rownames(var='index') %>%
		mutate(index=as.numeric(index)) %>%
		group_by(index) %>%
		left_join(cncf,by='chrom') %>%
		slice(which.min(abs(mid-loc.mid))) %>%
		filter(cnlr.median<0.5) %>%
		.$index

	loc.del <- loc.del[logic.del]
	cnv.matrix[loc.del,sample.name] <<- -1

	# summary stats, need to be corrected for log odds ratio removal
	# data.frame(
	# 	del.start=c(0,cumsum(sample.rle$lengths))[mark.del]+1,
	# 	del.length=sample.rle$lengths[mark.del]
	# )

}


names(cnv.matrix) %>%
list.filter(!. %in% c('chrom','start','mid','end','hgnc','band')) %>% list.map(.,.) %>%
map(~ OverCall(.x))

# post.cn.calls <- lapply(foo,function(sample.name) cnv.matrix[,sample.name] %>% unlist %>% unname %>% rle )

write_tsv(cnv.matrix,'../../facets.clean/geneCN.fill.txt')

cn.map <-
	cnv.matrix %>%
	select(chrom,start,end,matches('_threshold')) %>%
	gather(sample,cn,-c(chrom,start,end)) %>%
	group_by(sample,chrom) %>% arrange(sample,chrom,start) %>%
	mutate(chrom.max=max(end)) %>%
	mutate(chrom.min=min(start)) %>%
	mutate(lag.end=lag(end)) %>%
	mutate(lag.end=ifelse(is.na(lag.end),end,lag.end)) %>%
	filter(cn!=lag(cn,default=0)) %>%
	mutate(lead.lag.end=lead(lag.end)) %>%
	mutate(lead.lag.end=ifelse(is.na(lead.lag.end),chrom.max,lead.lag.end)) %>%
	ungroup %>% rowwise %>% arrange(chrom,start) %>%
	mutate(seg.length=as.numeric(lead.lag.end-start+1)) %>%
	ungroup %>% group_by(sample) %>% arrange(sample,chrom,start) %>%
	mutate(seg.sum=cumsum(seg.length)-seg.length) %>%
	mutate(genome.length=max(seg.sum)) %>%
	ungroup %>% rowwise %>%
	mutate(seg.frac=seg.length/genome.length) %>%
	ungroup %>% group_by(sample) %>% arrange(sample,chrom,start) %>%
	mutate(scale.start=(lag(seg.sum,default=0,order_by=sample)+1)/genome.length,scale.end=seg.sum/genome.length) %>%
	ungroup %>% arrange(sample,chrom,start) %>%
	mutate(sample.max = abs(group_indices_(., .dots='sample')-n())) %>%
	mutate(sample.min = sample.max-1) %>%
	mutate(cn=
		ifelse(cn==-2,'deletion',
		ifelse(cn==-1,'loss',
		ifelse(cn==0,'no change',
		ifelse(cn==1,'gain',
		ifelse(cn==2,'amplification',
			NA)))))) %>%
	mutate(cn=factor(cn,levels=c('amplification','gain','no change','loss','deletion')))

pdf('cn.map.uncorrected.pdf',width=14,height=,6,bg='white')
	palette  <- c(`deletion`='#D13838',`loss`='#DD9492',`no change`='#EAEAEA',`gain`='#7895BC',`amplification`='#284996')
	cn.plot <- ggplot() + 
	geom_rect(cn.map,mapping=aes(xmin=scale.start,xmax=scale.end,ymin=sample.min,ymax=sample.max,fill=cn)) +
	scale_fill_manual(values=palette,na.value='gray90') +
	theme(	legend.title 		= element_blank(),
			panel.grid.major 	= element_blank(),
			panel.grid.minor    = element_blank(),
			text 				= element_text(size=15),
			axis.title.x		= element_text(vjust=15),
			axis.title.y		= element_text(vjust=1),
			axis.text.y  		= element_text(face='italic'),
			axis.text    		= element_text(size=15),
			legend.key   		= element_rect(colour='white',fill=NULL,size=0.2),
			legend.key.size 	= unit(1.4, 'lines'),
			legend.text			= element_text(size=15),
			panel.background 	= element_rect(fill='white',color='white')
			#axis.line 			= element_line(colour = 'black'),
			#panel.border = element_rect(fill = NA, colour = 'black')
		)
	grid.draw(cn.plot)
dev.off()

cat('done\n')
