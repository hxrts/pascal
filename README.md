# Pascal
This package hooks into the ```jrflab/modules``` sequencing pipeline for more customized and dynamic analysis.

## subsetRecurrence
[component in-progress]

### Prerequisites
From root project directory:
1. requires: ```make recurrent_mutations```
2. initialize with: ```sh pascal/init.sh```
3. build subsets file by editing the ```subsets.txt``` file as follows (space delimited):
```
subsetName1 sample1 sample2
subsetName2 sample1 sample3 sample4
```

