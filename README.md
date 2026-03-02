# Table of Contents
- [Getting started](#getting-started)
- [Installation](#installation)
  - [Obtaining BINge](#obtaining-binge)
  - [Prerequisites](#prerequisites)
  - [Installation suggestions](#installation-suggestions)
- [How to use BINge](#how-to-use-binge)
  - [Running an analysis](#running-an-analysis)
  - [Rerunning an analysis](#rerunning-an-analysis)
- [How to cite](#how-to-cite)

# Getting started
```
# Obtain BINge
conda create -n binge bioconda::binge
conda activate binge

# Initialise a directory with query files and target genomes for clustering
BINge initialise -d /location/of/working/directory \
    -i /location/of/genome1.gff3,/location/of/genome1.fasta \ # .gff3 file optional for 'pre-seeding'
    --ix /location/of/transcriptome1.fasta \ # mRNA or CDS
    --ig /location/of/annotation1.gff3,/location/of/annotation1.fasta \ # extracts mRNA/CDS/protein sequences
    --threads 8

# Cluster all sequences within the initialised directory
BINge cluster -d /location/of/working/directory \
    --threads 8

## Further steps for filtering or DGE analysis follow:

# BLAST sequences against a database
BINge_post blast -d /location/of/working/directory \
    -t /location/of/database.fasta \ # ideally a UniRef## sequence file
    -s protein \ # or 'nucleotide' depending on the database molecule type
    --threads 8

# Salmon quantification of sequences
BINge_post salmon -d /location/of/working/directory \
    -r /location/of/reads -s .fq.gz \ # or '.fastq.gz' or '.fq' etc.
    --threads 8

# Filter sequence clusters that are low-quality
BINge_post filter -d /location/of/working/directory \
    --useBLAST \ # optional, only if you ran 'BINge_post blast'
    --useSalmon \ # optional, only if you ran 'BINge_post salmon'
    --readLength 150 # only used if '--useSalmon' is also provided

# Pick representative sequences for each cluster
BINge_post representatives -d /location/of/working/directory \
    --useGFF3 \ # optional, prioritises sequences in a GFF3 given to -i or --ig
    --useBLAST --useSalmon # optional, only if you ran the BINge_post functions

# Generate R script to streamline data loading and subsequent DESeq2 DGE analysis
BINge_post dge -d /location/of/working/directory

# Generate an annotation table for downstream investigations
BINge_post annotate -d /location/of/working/directory \
    -id /location/of/idmapping_selected.tab \ # UniProtKB download file
    -io /location/of/go.obo # gene ontology in OBO format
```

# Installation
It is recommended that you set up an Anaconda or Miniconda environment to run BINge. A conda package is provided through bioconda to simplify the installation of BINge and its prerequisites. To ensure package version compatibility, obtain BINge in its own environment by doing:

```
conda create -n binge bioconda::binge
```

See the [Installing BINge wiki page](https://github.com/zkstewart/BINge/wiki/Installing-BINge) for more details.

# How to use BINge
On the command line, you can always ask BINge to provide help information for its different functions by doing:

```
BINge initialise -h
BINge cluster -h
BINge_post filter -h
... etc ...
```

For the 'cluster' function only, many parameters are hidden since their default values are what most people should use. If you want to see those, try:

`BINge cluster --help-long`

Otherwise, refer to the [BINge wiki](https://github.com/zkstewart/BINge/wiki) for more detailed information on the program. The wiki also contains a tutorial guiding you through a full replication of one part of the BINge study which can be easily applied to your own dataset.

Some short but important details concerning BINge follow below.

## Running an analysis
When running BINge, each analysis will be given its own hash code to identify it e.g., `run_379df7a68ec7f39d3f02`. If you are running multiple BINge analyses using the same input files (in other words, you are just tweaking optional clustering parameters like the `--identity` value), you may want to specify which run to use with the `--analysis` parameter. This is not relevant to salmon or BLAST as they use the input sequences themselves which are not subject to any optional parameters.

However, if you do not specify a `--analysis` value, all BINge functions will default to using the most recently generated results from any prerequisite steps. In practice, if you are running BINge without tweaking any optional parameters during the clustering process, you do not need to think about run identifiers or hash codes. Just use BINge and it will handle that all for you!

If you are not sure what run you want to be using, try using `BINge view --analysis ...` substituting the `...` with the folder/run name of any outputs you are curious about. It will give you a complete overview of what files and parameters pertain to that analysis.

## Rerunning an analysis
If you try to rerun a previously completed analysis, you might encounter an error indicating that the output file already exists. BINge does not like to overwrite existing results, so you should move, delete, or rename any results before rerunning something performed previously.

# How to cite
A publication is hopefully forthcoming which can be referred to when using this program. Until then, you can link to this repository.
