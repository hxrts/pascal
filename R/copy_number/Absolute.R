
# format variants from mutation summary as ABSOLUTE input

#---------------
# LOAD LIBRARIES
#---------------

pacman::p_load(dplyr, lazyeval, readr, tidyr, magrittr, purrr, stringr, rlist, crayon, openxlsx)

source('~/pascal/R/lib/fun.R')


#-----------
# PARAMETERS
#-----------

build.maf      = TRUE
build.cncf     = FALSE

absolute.dir   = 'absolute'
platform       = 'Illumina_WES'  # possible values: SNP_250K_STY | SNP_6.0 | Illumina_WES
muts.file      = 'summary/mutation_summary.xlsx'

allosome       = 'merge'
keys           = NULL

#----------
# FUNCTIONS
#----------

CncfHeader <- function(cncf.file, cncf.dir){

   # read cncf files
   cncf.file %>%
   str_c(cncf.dir, ., sep='/') %>%
   read.delim(stringsAsFactors=FALSE, sep='\t') %>%
   tbl_df
}

CncfWrite <- function(cncf, absolute.dir) {

   out.file <- str_c(absolute.dir, 'segment', str_c(cncf.file %>% str_split('.cncf.') %>% list.map(.[[1]]) %>% unlist, '.cncf.txt'), sep='/')

   # write copies using correct column headers
   cncf %>%
   rename(Chromosome=chrom, Start=loc.start, End=loc.end, Num_Probes=num.mark, cncfment_Mean=cnlr.median) %>%
   write_tsv(out.file)

}


#-----
# MAFS
#-----

muts <- ReadMuts(muts.file)

# format muts
muts %<>%
    FormatEvents(keys=keys, allosome=allosome) %>%
    select(sample, chrom, pos, gene, effect, everything())

if(build.maf == TRUE) {

   message(green('  [ building maf files ]'))

   # create mafs directory
   MakeDirs(str_c(absolute.dir, '/maf'))

   # process mafs
   muts.maf <-
      muts %>%
      mutate(dbSNP_Val_Status='validated') %>%
      arrange(sample, chrom, pos) %>%
      mutate(t_ref_count = round(depth.t * (1-maf.t))) %>%
      mutate(t_alt_count = round(depth.t * maf.t)) %>%
      mutate(Tumor_Sample_Barcode = str_c(sample, normal, sep='_')) %>%
      select(Tumor_Sample_Barcode,
             Hugo_Symbol            = gene,
             t_ref_count,
             t_alt_count,
             dbSNP_Val_Status,
             Chromosome             = chrom,
             Start_position         = pos) %>%
      arrange(Chromosome, Start_position) %>%
      filter(t_alt_count>0) %>%
      unique %>%
      split(.$Tumor_Sample_Barcode)

   map2(.x=muts.maf, .y=names(muts.maf), ~ {
      write_tsv(.x, str_c(absolute.dir, '/maf/', .y, '.maf.txt'))
   }) %>% invisible

}


#-----
# SEGS
#-----

if(build.cncf == TRUE) {

   message(green('  [ building cncf files ]'))

   # create segs directory
   MakeDirs(str_c(absolute.dir, '/segment'))

   # process cncfs
   cncfs <-
      list.files('facets', pattern='*.cncf.txt') %>%
      list.map(CncfHeader(., 'facets/cncf')) %>%
      list.map(CncfWrite(., absolute.dir)) %>%
      bind_rows %>%
      mutate(sample=ID) %>%
      separate(sample, into=c('sample', 'normal'), sep='_')

      CncfWrite(cncfs, absolute.dir)

}

message(green('  [ done ]'))
