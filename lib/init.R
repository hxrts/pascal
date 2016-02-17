#---------------
# base libraries
#---------------

suppressMessages(pacman::p_load(dplyr,readr,tidyr,magrittr,purrr,stringr,rlist,crayon))

#------
# input
#------

# read samples file
samples <-
	read.delim("sample_sets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE) %>%
	filter(!grepl("#",V1)) %>%
	(function(sets){
		tumor <-
			apply(sets,1, function(row) row %>% list.filter (.!="") %>% head(-1)) %>% unlist
		normal <-
			apply(sets,1, function(row) row %>% list.filter (.!="") %>% tail(1)) %>%
			rep(apply(sets,1,function(row) row %>% list.filter(.!="") %>% length-1))
		data.frame(normal,tumor,stringsAsFactors=FALSE) %>%
		tbl_df
	})

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
	} else {	# add all samples as subset if no subsets file exists
		samples$tumor %>%
		as.data.frame %>%
		setNames(getwd() %>% strsplit("/") %>% unlist %>% tail(1))
	}


# filter samples as those in subsets
samples %<>% filter(tumor %in% unlist(subsets))


#---------------
# print metadata
#---------------

cat(green("\n-") %+% " samples\n\n") ; print(samples %>% as.data.frame,row.names=FALSE,right=FALSE)
cat(green("\n-") %+% " subsets\n\n") ; print(subsets %>% as.data.frame,row.names=FALSE,right=FALSE) ; cat("\n")
