#!/usr/bin/env Rscript

#-------
# params
#-------

suf.conda = '/home/debruiji/share/usr/anaconda-envs/sufam-0.4.2'


#------
# input
#------

MakeDirs('summary/sufam/log')

vcfs <-
	samples %>%
	rename(sample=tumor) %>%
	split(.$name) %>%
	map2(., names(.), ~ {

		vcf.file <- str_c('summary/sufam/', .y, '.vcf')

		vcf <-
			.x %>%
			left_join(muts) %>%
			arrange(sample, chrom, pos, ref, alt) %>%
			mutate(chrom=as.character(ifelse(chrom==23, 'X', chrom))) %>%
			select(`#CHROM`=chrom, POS=pos, ID=mut.id, REF=ref, ALT=alt)

		write_tsv(vcf, vcf.file)
		print(vcf)

		return(vcf.file)

		} )

vcfs %>%
map2(., names(.), ~ {

		conda <- str_c('unset PYTHONPATH && source ', suf.conda, '/bin/activate ', suf.conda, ' && ')
		cmd.m <- str_c(conda, 'sufam /ifs/e63data/reis-filho/reference/human_g1k_v37.fa ', .x, ' bam/', .y, 'M.bam 2> summary/sufam/log/', .y,'M_sufam.log > summary/sufam/', .y,'M_sufam.tsv')
		cmd.p <- str_c(conda, 'sufam /ifs/e63data/reis-filho/reference/human_g1k_v37.fa ', .x, ' bam/', .y, 'P.bam 2> summary/sufam/log/', .y,'P_sufam.log > summary/sufam/', .y,'P_sufam.tsv')

		message(str_c('running sufam on sample ', .y, 'M'))
		message(cmd.m)
		system(cmd.m)

		message('\n+\n')

		message(str_c('running sufam on sample ', .y, 'P'))
		message(cmd.p)
		system(cmd.p)
		message('\n-----\n')

} ) %>% invisible

sufam <-
	list.files('summary/sufam', pattern='_sufam.tsv$', full.names=TRUE) %>%
	list.map(read_tab(.) %>%
		rowwise %>%
		TypeCol(c('chrom', 'X..1', 'X..2', 'most_common_indel', 'most_common_indel_type', 'val_al_type', 'val_ref', 'val_alt', 'most_common_al'),
				c('char',  'char', 'char', 'char',              'char',                   'char',        'char',    'char',    'char')) %>%
		ungroup
		) %>%
	bind_rows %>%
	separate(sample, into = c('dir', 'sample'), sep = "\\/") %>%
	separate(sample, into = c('sample', 'ext'), sep = "\\.") %>%
	select(-dir, -ext) %>%
	mutate(val_alt=ifelse(val_alt=='TRUE','T',val_alt)) %>%
	arrange(sample, chrom, pos, ref, val_alt) %>%
	rowwise %>%
	mutate(alt_init=substr(val_alt,1,1)) %>%
	mutate(ref_init=substr(ref,1,1)) %>%
	ungroup

write_tsv(sufam, 'summary/sufam/sufam_all.tsv')


message(green(' [complete]'))


muts <-
	read_tab('summary/mutation_summary.tsv') %>%
	FormatEvents %>%
	mutate(chrom=ifelse(chrom==23, 'X', chrom))

muts %<>%
	rowwise %>%
	mutate(alt_init=substr(alt,1,1)) %>%
	mutate(ref_init=substr(ref,1,1)) %>%
	ungroup



muts.sufam <-
	sufam %>%
	FormatEvents %>%
	left_join(muts %>% mutate(chrom=as.numeric(ifelse(chrom=='X',23,chrom))) %>% select(-sample, -alt.count.t, -ref.count.t, -ref2, -alt2, -maf.t, -depth.t, -pr.somatic.clonal, -ci95.low, -tcn, -lcn, -ccf, -ref, -alt, -loh, -clonal), by=c('chrom','pos','ref_init','alt_init')) %>%
	mutate(alt=val_alt) %>%
	select(-val_alt) %>%
	#rename(maf.t=val_maf, alt.depth.t=val_al_count) %>%
	rename(maf.t=most_common_maf, alt.depth.t=most_common_count) %>%
	#TypeCol(c('lcn', 'tcn'), c('num', 'num')) %>%
	select(sample, chrom, pos, gene, effect, everything()) %>%
	mutate(provean=NA) %>%
	mutate(depth.t=cov) %>%
	mutate(hotspot=NA) %>%
	rowwise %>%
	mutate(patient=as.numeric(substr(sample,9,11))) %>%
	mutate(loc=substr(sample,12,12)) %>%
	ungroup %>%
	unique %>%
	#filter(maf > 0) %>%
	mutate(mut.id=str_c(sample, ':', row_number())) %>%
	select(mut.id, sample, case, chrom, pos, gene, effect, everything())



sufam %>% write_tsv('sufam_all.tsv')

muts.sufam %>% write_tsv('muts_sufam.tsv')


muts.extra <- data_frame(sample=c('JRFbcmet013M','JRFbcmet018M'),chrom=c(19,6),pos=c(15271691,157099165), gene=c('',''))

cmds <-
	muts.extra %>%
	rowwise %>%
	mutate(cmd=str_c('samtools view -hb bam/', sample, '.bam ',chrom,':',pos,'-',pos,' > summary/bams/',sample,'_',chrom,'_',pos,'.bam && samtools index summary/bams/',sample,'_',chrom,'_',pos,'.bam')) %>%
	.$cmd




cmds %>%
list.map(system(.))




