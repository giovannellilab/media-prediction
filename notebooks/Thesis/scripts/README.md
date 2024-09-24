# Assigning gene clusters for a given set of taxon IDs


## 1. Getting protein annotations for the given taxa

### 1.1. Downloading protein annotations

[NCBI's datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/getting_started/) is used for querying each taxon ID at a time (to avoid IP ban).

```bash
bash ncbi-download-taxon.sh <TAXON_ID_LIST.txt>
```

### 1.2. Extracting actual annotations and merging them together

The results from the previous step are extracted from their respective folders.
Addition of taxon ID to the sequences header is performed (right after the ">") in order to allow for further identification in the downstream analysis.

```bash
bash ncbi-preprocess-data.sh <OUTPUT_DIR>
```


## 2. Gene clustering using [MMseqs2](https://github.com/soedinglab/mmseqs2) and [UniRef90](https://www.uniprot.org/help/uniref)

### 2.1. Preparing the target (UniRef90) and query (NCBI protein annotations) databases

```bash
# Target database: UniRef90
mmseqs databases UniRef90 /path/to/UniRef90 tmp

# Query database: NCBI protein annotations
mmseqs createdb all-proteins.faa queryDB
```

### 2.2. Search query sequences in UniRef90

```bash
mmseqs search queryDB /path/to/UniRef90 resultDB tmp
```

### 2.3. Format output results

```bash
mmseqs convertalis queryDB /path/to/UniRef90 resultDB resultDB.m8
```
