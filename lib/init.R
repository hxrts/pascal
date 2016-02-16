#---------------
# base libraries
#---------------

suppressMessages(pacman::p_load(dplyr,readr,tidyr,magrittr,purrr,stringr,rlist,crayon))

#------
# input
#------

# read subsets file
subsets <-
	if("subsets.txt" %>% file.exists){
		read.delim("subsets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE) %>%
		setNames(c("subset",1:(ncol(.)-1))) %>%
		group_by(subset) %>%
		filter(!grepl("#",subset)) %>%
		summarise_each(funs(ifelse(.=="","NA",.))) %>%
		t %>%
		as.data.frame(stringsAsFactors=FALSE) %>%
		setNames(.[1,]) %>%
		slice(-1) %>%
		tbl_df
}

# read samples file
samples <-
	read.delim("sample_sets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE) %>%
	setNames(c("patient",1:(ncol(.)-1))) %>%
	group_by(patient) %>%
	filter(!grepl("#",patient)) %>%
	summarise_each(funs(ifelse(.=="","NA",.))) %>%
	setNames(c("tumor","normal"))

# add all samples as subset if no subsets file exists
if("subsets" %>% exists){
	subsets <-
		samples$tumor %>%
		as.data.frame %>%
		setNames(getwd() %>% strsplit("/") %>% unlist %>% tail(1))
}

# filter samples as those in subsets
samples %<>% filter(tumor %in% unlist(subsets))


#---------------
# print metadata
#---------------

cat(green("\n-") %+% " samples\n\n") ; print(samples %>% as.data.frame,row.names=FALSE,right=FALSE) ; cat("\n")
cat(green("\n-") %+% " subsets\n\n") ; print(subsets %>% as.data.frame,row.names=FALSE,right=FALSE) ; cat("\n")
