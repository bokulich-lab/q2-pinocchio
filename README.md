# q2-long-reads-qc

### QIIME 2 Plugin for quality control of long read sequences using the LongQC tool



##### Step 1: Create q2-long-reads-qc environment
```shell
mamba create -n q2-long-reads-qc -c conda-forge -c bioconda -c https://packages.qiime2.org/qiime2/2024.2/shotgun/passed/ -c defaults q2cli q2-types minimap2 bs4
```

##### Step 2: Activate q2-long-reads-qc environment
```shell
conda activate long-reads-qc
```

##### Step 3: Installing python package
```shell
make dev
```

##### Step 4: Download input data
[Download](https://polybox.ethz.ch/index.php/s/Y81jl4JAtPjuKH6)

##### Step 5: Execution

* minimap2-build 
```shell
qiime long-reads-qc minimap2-build --i-sequences reference.qza --o-database database.qza --verbose
```

* filter-reads
  - keep mapped
  ```shell
  qiime long-reads-qc filter-reads --i-minimap2-index database.qza --i-sequences reads.qza --o-filtered-sequences mapped.qza --verbose
  ```
  - keep unmapped
  ```shell
  qiime long-reads-qc filter-reads --i-minimap2-index database.qza --i-sequences reads.qza --p-exclude-mapped True --o-filtered-sequences unmapped.qza --verbose
  ```
  - keep mapped reads with mapping percentage >= 85%
  ```shell
  qiime long-reads-qc filter-reads --i-minimap2-index database.qza --i-sequences reads.qza --p-min-per-identity 0.85 --o-filtered-sequences mapped_over_85.qza --verbose
  ```
