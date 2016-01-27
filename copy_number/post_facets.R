#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

# load base libraries
suppressMessages(pacman::p_load(dplyr,readr,stringr,tidyr,purrr,magrittr,crayon))

# read 
gcn <- "facets/geneCN.txt" %>% read.delim(stringsAsFactors=FALSE,header=TRUE)

chrom <- suppressWarnings(as.numeric(gcn$chrom))
chrom[which(is.na(chrom))] <- 23
gcn[is.na(gcn)] <- 3	# replace NA with numeric for use in rle function
breaks <- c(0,cumsum(table(chrom)[chrom %>% table %>% names %>% as.numeric %>% order]))	# start-1 == end

#-----
# main
#-----

for (samplename in grep("threshold",names(gcn),value=TRUE)){

	cat(blue("\n*") %+% " sample: " %+% samplename)

	#loop over chromosomes in junction table
	for(chromosome in 1:23){

		cat("\n   chr",formatC(chromosome, width=2, flag="0"),": ",sep="")
		start<-breaks[chromosome]+1
		end<-breaks[chromosome+1]
		chbit <- rle(gcn[start:end,samplename])

		# replace an entirely empty chromosome with calls of 0
		if((chbit$lengths %>% length)==1){
			gcn[start:end,samplename] <- 0
		}else{
			# extend chromosome start calls from nearest integer on same chromosome
			if(gcn[start,samplename]==3){
				cat("start..")
				NAbottom <- min(which(chbit$values==3))
				Ibottom <- start + cumsum(chbit$lengths)[NAbottom]
				gcn[start:(Ibottom-1),samplename] <- gcn[Ibottom,samplename]
			}
			# extend chromosome end calls from nearest integer on same chromosome
			if(gcn[end,samplename]==3){
				cat("end..")
				NAtop <- min(which(rev(chbit$values==3)))
				Itop <- end - cumsum(rev(chbit$lengths))[NAtop]
				gcn[(Itop+1):end,samplename] <- gcn[Itop,samplename]
			}
		}

		nav <- which(chbit$values[-c(1,length(chbit$values))]==3)+1
		if(length(nav>0)){
			cat("gaps..")
			# find nearest neighbors of inner-chromosome gaps
			gstarts <- cumsum(chbit$lengths)[nav-1]+start
			glengths <- chbit$lengths[nav]
			map2(gstarts,glengths, ~ .x:(.x+.y-1)) %>%
				unlist ->
				nrows
			nrows %>%
				slice(gcn,.) %>%
				mutate(mid=rowMeans(.[,c("start","end")])) %>%
				select(chrom,mid,hgnc) %>%
				rowwise %>%
				inner_join(gcn %>%
					mutate(qmid=rowMeans(.[,c("start","end")])) %>% 
					select(chrom,qmid,get(samplename)),by="chrom") %>%
				ungroup %>%
				filter_(str_c(samplename,"!=3")) %>%
				group_by(hgnc) %>%
				arrange_(str_c("desc(",samplename,")")) %>%
				slice(which.min(abs(mid-qmid))) %>%
				ungroup %>%
				select_(samplename) %>%
				unlist ->
				nfills
			# write neighbor calls to master table
			if(length(nfills)>0){
				gcn[nrows,samplename] <- nfills
			}else{
				gcn[nrows,samplename] <- NA
			}
		}
	}
	cat("\n")
}

#-----------
# finalizing
#-----------

# write updated table
write_tsv(gcn,"facets/geneCN_fill.txt")

cat(green("\n  [done]\n\n"))

