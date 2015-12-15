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

Generates log file in log directory, vcf and sufam output .tsv in recurrent_mutations directory.
