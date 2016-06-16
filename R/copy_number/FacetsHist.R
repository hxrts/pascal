
# facets copy number histogram plots

#----------
# libraries
#----------

pacman::p_load( dplyr,lazyeval,readr,tidyr,magrittr,purrr,stringr,rlist,  # base
                crayon,colorspace,RColorBrewer,                           # coloring
                ggplot2,grid,gridExtra,gplots,scales )                    # plotting


#-----------
# parameters
#-----------

absolute.dir    = 'absolute'
cncf.dir        = 'facets/cncf'
hist.dir        = 'facets/cn_hist'


#----------
# functions
#----------

# make directories
MakeDirs <- function(dir.list, debug=TRUE) {
    dir.list %>% list.map({
        dir <- .
        if(debug){ message(blue(str_c('- making directory: ', dir))) }
        dir.create(dir, recursive=TRUE, showWarnings=FALSE)
    }) %>% invisible
}


# read cncf files
ReadCncf <- function(cncf.file, cncf.dir) {

   # read cncf files
   cncf.file %>%
   str_c(cncf.dir, ., sep='/') %>%
   read.delim(stringsAsFactors=FALSE, sep='\t') %>%
   tbl_df

}


# plot histogram
PlotHist <- function(cncf.h, sample.name, hist.dir) {

   message(blue(str_c('- plotting: ', sample.name)))

   hp <-
      ggplot(cncf.h, aes(cn, ..density.., fill=cn.type)) +
      geom_histogram(aes(weight=weight), binwidth=0.2) +
      stat_density(bw=0.1, alpha=0.2, color='black') +
      xlab('copy number') +
      ylab('density') +
      scale_x_continuous(breaks=pretty_breaks()) +
      facet_wrap(~ cn.type, nrow=2)

   pdf(str_c(hist.dir, '/', sample.name, '.pdf'), 5, 5)
      plot(hp)
   dev.off()

}


#-----
# main
#-----

MakeDirs(hist.dir)

# process cncfs
cncfs <-
   list.files(cncf.dir, pattern='*.cncf.txt') %>%
   list.map(ReadCncf(., cncf.dir)) %>%
   bind_rows %>%
   mutate(sample=ID) %>%
   separate(sample, into=c('sample', 'normal'), sep='_')

# extract copy number information
cncfs.cn <-
    cncfs %>%
    mutate(weight=loc.end - loc.start) %>%
    select(ID, weight, tcn.em, lcn.em) %>%
    gather(cn.type, cn, tcn.em:lcn.em) %>%
    filter(!is.na(cn)) %>%
    split(.$ID)

# build histograms
map2(cncfs.cn, names(cncfs.cn), ~ { PlotHist(.x, .y, hist.dir) }) %>% invisible

message(green('[ done ]'))
