# Pascal
This package hooks into the ```jrflab/modules``` sequencing pipeline for tailored / dynamic analysis.

## subsetRecurrence
Builds a prevelance-ordered heatmap of recurrent mutations for a subset of sample from the master recurrent mutation table.

### Prerequisites
From root project directory:

1. requires: ```make recurrent_mutations```
2. initialize with: ```sh pascal/init.sh```
3. build subsets file by editing the ```subsets.txt``` file as follows (space delimited):

```
subsetName1 sample1 sample2
subsetName2 sample1 sample3 sample4
```

