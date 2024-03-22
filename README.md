# q2-minimap2

## QIIME 2 Plugin for quality control and taxonomic classification of long read sequences using Minimap2



### Step 1: Create q2-minimap2 environment
```shell
mamba create -n q2-minimap2 -c conda-forge -c bioconda -c https://packages.qiime2.org/qiime2/2024.2/shotgun/passed/ -c defaults q2cli q2-types q2-feature-classifier minimap2 bs4 samtools
```

### Step 2: Activate q2-minimap2 environment
```shell
conda activate q2-minimap2
```

### Step 3: Installing python package
```shell
make dev
```

### Step 4: Download input data
[Download](https://polybox.ethz.ch/index.php/s/Y81jl4JAtPjuKH6)

### Step 5: Execution

* build-index
  - Build Minimap2 index database
  ```shell
  qiime minimap2 build-index --i-sequences reference.qza --o-database database.qza --verbose
  ```

<br>

* minimap2-search
  - Generate both hits and no hits for each query. Keep a maximum of one hit per query (primary).
  ```shell
  qiime minimap2 minimap2-search --i-query-reads fasta_reads.qza --i-index-database database.qza --o-search-results paf.qza
  ```

  - Generate only hits for each query. Keep a maximum of one hit per query (primary mappings).
  ```shell
  qiime minimap2 minimap2-search --i-query-reads fasta_reads.qza --i-index-database database.qza --o-search-results paf_only_hits.qza
  ```

  - Generate only hits for each query. Keep a maximum of 10 hits per query.
  ```shell
  qiime minimap2 minimap2-search --i-query-reads fasta_reads.qza --i-index-database database.qza --p-maxaccepts 10 --p-output-no-hits False --o-search-results paf_only_hits_up_to_10_per_query.qza
  ```

<br>

* filter-reads
  - Keep mapped
  ```shell
  qiime minimap2 filter-reads --i-query-reads reads.qza --i-index-database database.qza --o-filtered-query-reads mapped.qza --verbose
  ```
  - Keep unmapped
  ```shell
  qiime minimap2 filter-reads --i-query-reads reads.qza --i-index-database database.qza --p-keep "unmapped" --o-filtered-query-reads unmapped.qza --verbose
  ```
  - Keep mapped reads with mapping percentage >= 85%
  ```shell
  qiime minimap2 filter-reads --i-query-reads reads.qza --i-index-database database.qza --p-min-per-identity 0.85  --o-filtered-query-reads mapped_over_85p_id.qza --verbose
  ```
  
  - Using the reference sequences instead of the index database
  ```shell
  qiime minimap2 filter-reads --i-query-reads reads.qza --i-reference-reads reference.qza --o-filtered-query-reads mapped.qza --verbose
  ```

<br>

* classify-consensus-minimap2
  - Assign taxonomy to query sequences using Minimap2
  ```shell
  qiime minimap2 classify-consensus-minimap2 --i-index-database classification_input/index.qza --i-query classification_input/n1K_initial_reads_SILVA132.fna.qza --i-reference-taxonomy classification_input/raw_taxonomy.qza --p-num-threads 6 --output-dir outDir --verbose
  ```


