# Pascal
This package hooks into the ```jrflab/modules``` sequencing pipeline for tailored / dynamic analysis.

## subsetRecurrence
[in progress]
Builds a prevelance-ordered subsetted heatmap + table of recurrent mutations from the master recurrent mutation table.

### Prerequisites
From root project directory:

1. requires: ```make recurrent_mutations```
2. initialize with: ```sh pascal/init.sh```
3. build subsets file by editing the ```subsets.txt``` file as follows (space delimited):

```
subsetName1 sample1 sample2
subsetName2 sample1 sample3 sample4
```

# Run
```Rscript pascal/recurrence/subsetRecurrence.R```

Generates heatmap .pdf and table files in ```recurrent_mutations``` directory within project root.
