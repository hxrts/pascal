#!/usr/bin/env Rscript

# create necessary directories
system("mkdir functional_networks &>/dev/null")

#-----
# init
#-----

source("pascal/lib/init.R")

#----------
# libraries
#----------

suppressMessages(pacman::p_load(GO.db))

muts <-
	read.delim("excel/tsv/mutation_summary.tsv",stringsAsFactors=FALSE,sep="\t") %>%
	tbl_df
