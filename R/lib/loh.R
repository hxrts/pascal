#!/usr/bin/env Rscript

#---------------
# base libraries
#---------------

suppressMessages(pacman::p_load(dplyr,readr,tidyr,magrittr,purrr,stringr,rlist,ggplot2,crayon))

#-----------------
# path contingency
#-----------------

wd <- getwd()

if(!"samples.txt" %in% list.files()){
	setwd("..")
}

#----------------
# read cncf files
#----------------

samples %>%
mutate(cncf=str_c("facets/",tumor,"_",normal,".cncf.txt")) %$%
cncf %>%
map(~
	read.delim(.x,sep="\t",stringsAsFactors=FALSE) %>%
	tbl_df %>%
	select(tcn.em,lcn.em) %>%
	mutate(loh=ifelse(tcn.em==0,TRUE,))
	)



















#--------------------------------
# move back to original directory
#--------------------------------

setwd(wd)