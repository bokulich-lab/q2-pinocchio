# q2-long-reads-qc

### QIIME 2 Plugin for quality control of long read sequences using the LongQC tool



### Installation
#### Step 1:
* Create and activate a conda environment:
```shell
mamba create -n longReadQC -c conda-forge -c bioconda -c https://packages.qiime2.org/qiime2/2024.2/shotgun/passed/ -c defaults q2cli q2-types minimap2 bs4

conda activate longReadQC

```

#### Step 2:
* Install necessary libraries:
```shell
make install_dependencies
make all
```

