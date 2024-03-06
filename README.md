# q2-long-reads-qc

### QIIME 2 Plugin for quality control of long read sequences using the LongQC tool



#### Step 1: Create q2-long-reads-qc environment
```shell
mamba create -n q2-long-reads-qc -c conda-forge -c bioconda -c https://packages.qiime2.org/qiime2/2024.2/shotgun/passed/ -c defaults q2cli q2-types minimap2 bs4
```

#### Step 2: Activate q2-long-reads-qc environment
```shell
conda activate long-reads-qc
```

#### Step 3: Installing python package
```shell
make dev
```

