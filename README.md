# q2-LongReadQC

### QIIME 2 Plugin for quality control of long read sequences using the LongQC tool



### Installation
#### Step 1:
* Create and activate a conda environment:
```shell
mamba create -y -n longReadQC \   
   -c https://packages.qiime2.org/qiime2/2024.2/shotgun/passed/ \
   -c conda-forge -c bioconda -c defaults \
   q2cl

conda activate longReadQC

```
