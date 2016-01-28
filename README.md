# Pascal
This package hooks into the ```jrflab/modules``` sequencing pipeline for tailored / dynamic analysis.

### Setup

From root project directory:

1. initialize with: ```sh pascal/init.sh```
2. build subsets file by editing the ```subsets.txt``` file as follows (space delimited):

```
subsetName1 sample1 sample2
subsetName2 sample1 sample3 sample4
```

## subsetRecurrence
Builds a prevelance-ordered, subsetted heatmap + table of recurrent mutations from the master recurrent mutation table.

### Prerequisites

* requires: ```make recurrent_mutations```
* working installation of [sufam](https://github.com/inodb/sufam)

### Run
```Rscript pascal/recurrence/subsetRecurrence.R```

Generates heatmap .pdf and table files in recurrent_mutations directory within project root.

## sufamRecurrent
Performs base-truth recurrent mutation lookups in the .bam files of all samples within subsets.

### Prerequisites

* requires: ```make recurrent_mutations```

### Run
```Rscript pascal/recurrence/sufamRecurrent.R```

Generates log file in log directory, vcf and sufam output .tsv in recurrent_mutations directory

## sufamEventSearch

Use Sufam to arbitrarily search for base-truth genomic events in a sample subset by specifying a vcf file.

### Prerequisites

* add genomic positions of interest to events.vcf, column headings required ```#CHROM	POS	ID	REF	ALT```
* will run on all samples found in subsets.txt

### Run

```Rscript pascal/recurrence/sufamEventSearch.R```

## PyClone

### Prerequisites

* make mutation_summary
* make facets
* will run on all samples found in subsets.txt

### Run

```Rscript pascal/clonality/pyclone.R```

### PhyloWGS

### EXPANDS

### SCHISM
