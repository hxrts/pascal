# boinformatics functions [requires Rsamtools]
bam_coverage <- function (bamfile) {

	    # read in bam file
	    bam <- scanBam(bamfile)[[1]]

    # filter reads without match position
    ind <- !is.na(bam$pos)

        # remove non-matches
        bam <- lapply(bam, function(x) x[ind])
        ranges <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))

	    # names of bam data frame: qname, flag, rname, strand, pos, qwidth, mapq, cigar, mrnm, mpos, isize, seq, qual
	    # construc: genomic ranges object containing all reads
	    ranges <- GRanges( seqnames=Rle(bam$rname), ranges=ranges, strand=Rle(bam$strand), flag=bam$flag, readid=bam$rname )

	    # returns coverage for each reference sequence (chromosome) in bam
	    return (mean(coverage(ranges)))
}

