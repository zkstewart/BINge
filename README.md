# Table of Contents
- [Getting started](#getting-started)
- [Introduction](#introduction)
  - [Reason for BINge existing](#reason-for-binge-existing)
  - [What does BINge actually do?](#what-does-binge-actually-do)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Installation suggestions](#installation-suggestions)
- [How to cluster](#how-to-cluster)
  - ['initialise' a working directory](#initialise-a-working-directory)
  - ['cluster' a working directory](#cluster-a-working-directory)
  - [Interpreting results](#interpreting-results)
- [Downstream processing of clustering results](#downstream-processing-of-clustering-results)
  - ['blast' a working directory](#blast-a-working-directory)
  - [use 'salmon' on a working directory](#use-salmon-on-a-working-directory)
  - ['filter' clustering results in a working directory](#filter-clustering-results-in-a-working-directory)
  - [pick 'representatives' from clustering results in a working directory](#pick-representatives-from-clustering-results-in-a-working-directory)
  - ['dge' analysis of a working directory](#dge-analysis-of-a-working-directory)
  - ['annotate' a working directory](#annotate-a-working-directory)
- [How to cite](#how-to-cite)

# Getting started
```
# Download this repository [making sure you have prerequisite Python packages]
git clone https://github.com/zkstewart/BINge.git
cd BINge
git submodule update --init --recursive

# Initialise a directory with query files and target genomes for clustering
python BINge.py initialise -d /location/of/working/directory \
    -ix /location/of/transcriptome1.fasta \ # mRNA or CDS
    -ig /location/of/annotation1.gff3,/location/of/annotation1.fasta \ # extracts mRNA/CDS/protein sequences automatically
    -t /location/of/genome1.gff3,/location/of/genome1.fasta \ # .gff3 file optional for 'pre-seeding'
    --threads 8

# Cluster all sequences within the initialised directory
python BINge.py cluster -d /location/of/working/directory \
    --threads 8

## Further steps for filtering or DGE analysis follow:

# BLAST sequences against a database
python BINge_post.py blast -d /location/of/working/directory \
    -t /location/of/database.fasta \ # ideally a UniRef## sequence file
    -s protein \ # or 'nucleotide' depending on the database molecule type
    --threads 8

# Salmon quantification of sequences
python BINge_post.py salmon -d /location/of/working/directory \
    -r /location/of/reads -s .fq.gz \ # or '.fastq.gz' or '.fq' etc.
    --threads 8

# Filter sequence clusters that are low-quality
python BINge_post filter -d /location/of/working/directory \
    --useBLAST \ # optional, only if you ran 'BINge_post.py blast'
    --useSalmon \ # optional, only if you ran 'BINge_post.py salmon'
    --readLength 150 # only used if '--useSalmon' is also provided

# Generate R script to streamline data loading and subsequent DESeq2 DGE analysis
python BINge_post dge -d /location/of/working/directory

# Pick representative sequences for each cluster
python BINge_post representatives -d /location/of/working/directory \
    --useGFF3 \ # optional, prioritises sequences in a GFF3 given to -ig or -t
    --useBLAST --useSalmon # optional, only if you ran the BINge_post.py functions

# Generate an annotation table for downstream investigations
python BINge_post.py annotate -d /location/of/working/directory \
    -id /location/of/idmapping_selected.tab \ # UniProtKB download file
    -io /location/of/go.obo # gene ontology in OBO format
```

# Introduction
## Reason for BINge existing
BINge (**Bin** **G**enes for **E**xpression analyses) is a sequence clustering script which aims to group alternative isoforms, orthologs, and paralogs into sequence clusters.
Its algorithm can do this accurately while receiving multiple species' sequences as input, a feat that currently existing clustering programs cannot achieve.
In doing so, BINge enables differential gene expression analyses to take place that compare one or many species in a way that minimises biases associated with:
- Data mappability (e.g., reads from species A not aligning optimally to sequence models derived from species B), and
- Ortholog sequence length (e.g., a gene in species A may be longer than the ortholog in species B which means, at an equal sequencing depth and amount of expression, we might expect
more reads to be sequenced from species A's ortholog)

The solution to these biases is to not just map to a single species' reference models - we should map to all species' sequences. This way, we allow the reads from species A to align
optimally to its own sequences, and so on. Then we can tally the number of alignments made to all the sequences in a cluster after they've been quantified taking sequence length bias into account. 

## What does BINge actually do?
BINge automates the alignment of the input sequences against the genome sequences using GMAP. A **_bin_** will be positioned along the genome where these sequences align, or where a gene has been annotated onto the
genome in a GFF3 file. This **_bin_** loosely corresponds to an exon, and BINge will remember any sequences which align over that _binned_ position.

After this, BINge will look at each input sequence to see how many times other sequences were placed into the same bin with it. Fundamentally, this is measuring how many exons the input sequence shares in
common with any other sequence. When two sequences share many exon alignments in common with each other, they will be clustered together; this is the 'clustering by co-occurrence' algorithm that BINge employs.

Because BINge aligns each input sequence against many positions in each genome, but only looks at the best and near-best alignments, we can expect isoforms, orthologs, and paralog sequences to co-occur with each other
in bins. Hence, these sequences get clustered together.

# Installation
## Prerequisites
BINge is written in Python and requires a modern version 3 [note: development occurred using >= 3.9, and I'm not sure how far backwards its compatibility goes].

BINge and its associated utility scripts makes use of several other Python packages which include:
- networkx (https://networkx.org/)
- sci-kit-learn (https://scikit-learn.org/)
- intervaltree (https://github.com/chaimleib/intervaltree)
- goatools (https://github.com/tanghaibao/goatools)
- biopython (https://biopython.org/)
- pyfaidx (https://pypi.org/project/pyfaidx/)
- numpy (https://numpy.org/)

It also calls upon external programs including:
- GMAP (http://research-pub.gene.com/gmap/)
- MMseqs2 (https://github.com/soedinglab/MMseqs2)
- salmon (https://combine-lab.github.io/salmon)

And optionally:
- CD-HIT (https://sites.google.com/view/cd-hit)

BINge has been developed on Linux and within Windows Subsystem for Linux (WSL). If you want to run BINge in WSL, you must ensure that developer mode is enabled so that symbolic links are functional in Windows.
You can do that by going to Settings > For developers > Toggle 'Developer Mode' On.

## Installation suggestions
If you are interested in running BINge, you should ideally set up an Anaconda or Miniconda environment containing a recent Python 3 version. All Python packages listed above are available through community channels including conda-forge.

The external programs can be installed through the bioconda channel, or you may opt to install them yourself by following any instructions detailed at the provided websites.

# How to cluster
On the command line, you can always ask BINge.py to provide help information for its submodules by doing:

```
python /location/of/BINge.py initialise -h
python /location/of/BINge.py cluster -h
python /location/of/BINge_post.py filter -h
... etc ...
```

For the 'cluster' submodule only, many parameters are hidden since their default values are what most people should use. If you want to see those, try:

`python /location/of/BINge.py cluster --help-long`

Otherwise, BINge tries to keep the clustering process simple, and minimally requires you to ...

## 'initialise' a working directory
The first step of a BINge clustering analysis is to 'initialise' a directory using the correspondingly named submodule in BINge.py. Sequences can be provided in a variety of formats. 

- Using `-ix` you can provide FASTA sequence files such as *de novo* assembled transcriptome(s) or reference genome mRNA or CDS models.
  - You can optionally provide three inputs here separated by commas with a format like `mrna.fasta,cds.fasta,protein.fasta` if you have predicted ORFs yourself. It is *essential* that you provide the three files in that order - mRNA, then CDS, then protein.
- Using `-ig` you can provide pairs of genome GFF3 and FASTA annotation files from which protein-coding mRNA sequences will be extracted for you.
- Using `-t` you can provide pairs of GFF3 and FASTA annotation files, which will provide the genome(s) as target for clustering and will also result in their sequences being extracted and included in the clustering process.
  - If you are providing pairs of files to `-t` it becomes optional to include files with `-ix` or `-ig`. If you only provide the genome target to `-t` then you must provide at least one input to either `-ix` or `-ig`.

`-t` must receive at least one genome FASTA file to serve as a target for clustering. Multiple can be given, and if you provide a paired GFF3 annotation you will enable the clustering algorithm to be *pre-seeded* with *bins* to guide the clustering process.

If you are inputting a GFF3 file for a microbial organism without alternative splicing, you must specify the `--microbial` argument so that BINge knows to interpret the GFF3 following a hierarchy of gene->CDS, rather than gene->mRNA->exon/CDS.
You should also specify the NCBI translation table number with `--translation` so that the correct codon table can be used for protein translations.

## 'cluster' a working directory
After a directory has been 'initialise'd, all the prerequisite files for clustering will exist in the working directory. You can make use of the default (and recommended) BINge parameters by leaving them blank.
Hence, clustering can be performed as simply as:

`python BINge.py cluster -d /location/of/working/directory --threads 8`

If you specified `--microbial` during 'initialise' you should do so again.

## Interpreting results
BINge will run everything within the specified `-d` working directory, and produce a result file with the name `BINge_clustering_result.tsv` within a run folder identified by
a hash of its parameters. An example directory structure might be like:

`/location/of/working/directory/analysis/run_55f3403b67b7371cb2eb/BINge_clustering_result.tsv`

The output file is tab separated into three columns with a format like:

```
#BINge clustering information file
cluster_num  sequence_id  cluster_type
0  sequence1  binned
0  sequence2  binned
1  sequence3  binned
1  sequence4  binned
2  sequence5  binned
...
10001  sequenceN  unbinned
10002  sequenceN+1  unbinned
...
```

Hence, it contains one comment line at the start, followed by a header row, then by the actual results rows. Sequences with the same cluster number (first column) belong to the same
cluster group, with the sequence's ID (second column) followed by whether the cluster was binned or unbinned (third column). 

# Downstream processing of clustering results
`BINge_post.py` offers several utilities to make use of a clustering result and/or convert it into more useful formats.
Several of these utilities rely on or benefit from running a MMseqs2 search and/or mapping reads to sequences using salmon.

## 'blast' a working directory
The `BINge_post.py blast` module enables one to BLAST (using it as a verb, not as the program name!) the sequences that were clustered to a database of sequences.
These results can be used for later filtering of clusters, or output of an annotation table. This can be done like in the following example:

`python BINge_post.py blast -d /location/of/working/directory -t database.fasta -s <protein / nucleotide> --threads 8`

- The file given to `-t` can be any sequence file, but for compatibility with later annotation table generation, you should provide a UniRef database e.g., UniRef90.
- The value given to `-s` should indicate the molecule type of the `-t` sequence.
- The `--threads` value should be whatever your system can spare for quicker operation (for this and any other BINge process).

You will find the resulting file at `/location/of/working/directory/blast/MMseqs2_results.tsv`

## use 'salmon' on a working directory
The `BINge_post.py salmon` module enables one to use salmon to map reads against the sequences that were clustered. These results can be used for later filtering of
clusters, or for use in DGE analysis. A simple example is as follows:

`python BINge_post.py salmon -d /location/of/working/directory -r /location/of/reads -s .fq.gz --threads 8`

- One or more directories can be provided to `-r` which will be searched for all files ending in the suffix given to `-s`.
- For paired reads, you must ensure that `1` and `2` immediately precedes the suffix given to `-s`.
  - For example, if you indicate `-s .fq.gz` then we would expect to find files like sampleA_**1.fq.gz** and sampleB_**2.fq.gz**.
- If you are using single end reads instead, you should specify `--singleEnd`.

You will find the resulting files at `/location/of/working/directory/salmon` in subdirectories with names corresponding to the prefix of your read files.

## 'filter' clustering results in a working directory
The raw results of BINge clustering contain both *binned* and *unbinned* cluster sequences. In many cases you might want to exclude unbinned clusters since they did not align well against any of the input genomes.
Or, if your reference genome(s) are incomplete, you might instead just want to filter these *unbinned* clusters to remove spurious or low-quality ones.
You can do that like:

`python BINge_post.py filter -d /location/of/working/directory --useBLAST --useSalmon --readLength 150`

- The options `--useBLAST` and `--useSalmon` are available if you've run their corresponding modules beforehand. Doing so allows you to filter out clusters that fail to obtain a good BLAST result,
or which fail to obtain much read coverage.
- You can use `--useGFF3` if you'd like to retain clusters that contain a sequence given in a GFF3 file during working directory initialisation.
  - The logic here is that a reference annotation sequence is likely to be genuine, whereas sequences given in a transcriptome may be spurious.
- The `--readLength` option should be specified if you `--useSalmon` since the read length is used to derive some automatic coverage cutoffs.
- Alternatively, you might specify `--justDropUnbinned` to remove any *unbinned* sequence clusters and ignore any other fancy filtering.
  - This is *probably ideal* if you are confident that your genome reference(s) are highly complete, in which case *unbinned* sequences are likely to be spurious. This is *not ideal* if unbinned sequences may correspond to orphan genes that don't exist in your reference genome(s)!

After running 'filter', any downstream use of 'representatives' or 'dge' will automatically preference this result. You can override the automatic preferencing of `BINge_post.py` at any point
by using the `--analysis` argument to specifically indicate which run folder you want to use. Otherwise, we default to 1) using filtered results over raw results, and 2) using the most recent
results.

You will find the resulting file `BINge_clustering_result.filtered.tsv` at `/location/of/working/directory/filter/run_<hash>` or `/location/of/working/directory/filter/most_recent` if you have not run any other 'filter' analyses.

## pick 'representatives' from clustering results in a working directory
Since each cluster represents a putative locus, you may want just one sequence from each cluster which you can take as being *representative* of all the other sequences in that cluster.
This could enable completeness analysis of clustering results using BUSCO/compleasm, or functional annotation of sequence clusters. You can do this like:

`python BINge_post.py representatives -d /location/of/working/directory --useBLAST --useSalmon`

- The options here are like that seen when using 'filter'.
- Instead, we can `--useBLAST` evidence to select a representative from a cluster which has the best BLAST result to your chosen database.
- We can also `--useSalmon` evidence to select a representative from a cluster based on its read depth, selecting for one which has the most alignments in your dataset.
- You can use `--useGFF3` to select representatives that show up in a reference annotation. These sequences may have identifiers which make it easier to relate them to other datasets.

You will find the resulting files at `/location/of/working/directory/representatives/run_<hash>` or `/location/of/working/directory/representatives/most_recent` if you have not run any other 'representatives' analyses.
Three files should exist, starting with the prefix `BINge_clustering_representatives` and ending with a file suffix where `.aa` contains protein sequences, `.cds` contains CDS, and `.mrna` contains the mRNA transcripts.

## 'dge' analysis of a working directory
After running 'salmon', BINge will help you to make use of these read counts alongside the most recent 'filter' analysis or, if no filtering has been done, the most recent raw clustering result.
This is done like:

`python BINge_post.py dge -d /location/of/working/directory`

No optional parameters are needed for this submodule. It will create outputs at `/location/of/working/directory/dge/run_<hash>` or `/location/of/working/directory/dge/most_recent` if you have no run any other 'dge' analyses.
The contents of this directory include:

- A `tx2gene.tsv` file which enables summarisation of transcript quantification values to the locus level as predicted by BINge clustering.
- A `salmon_qc_statistics.tsv` file which contains three columns to indicate 1) the sample name, 2) the percentage of reads which mapped to your clustered sequences, and 3) the percentage of read mappings that remain after filtering of clusters.
  - You would likely report this file as a Supplementary Table in any manuscript to establish that your read mapping was successful, and that the filtering was not too strict.
- A `BINge_DGE.R` file which provides the initial script setup for a DGE analysis in DESeq2, making use of your salmon mapping results and the `tx2gene.tsv` file.
- A `salmon_samples.txt` which is used by the R script for locating your salmon alignment files.

## 'annotate' a working directory
If you've run 'blast' against a UniRef database and also picked 'representatives', then the `BINge_post.py annotate` submodule will allow you to generate an annotation table relating clusters to their best BLAST hits, gene names, and gene ontologies.
To do this, you should run something like:

`python BINge_post.py annotate -d /location/of/working/directory -id /location/of/idmapping_selected.tab -io /location/of/go.obo`

As indicated, you must provide the `idmapping_selected.tab` file as found from [the UniProtKB FTP site](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/) as well as
a `go.obo` file as found from [the Gene Ontology downloads site](https://geneontology.org/docs/download-ontology/).

You can also provide the `--largeTable` option if you'd like the resulting table file to contain the full outfmt6 details rather than just the percentage identity, E-value, and bitscore of hits.

You will find the resulting file `BINge_annotation.tsv` at `/location/of/working/directory/annotate/run_<hash>` or `/location/of/working/directory/annotate/most_recent` if you have not run any other 'annotate' analyses.
The file contains a header to help with interpretation. Note however that the `Best_mapped_GOs_+_parents` column contains the same values as `Best_mapped_GOs` in addition to all their parent/ancestor terms. You should use
the `Best_mapped_GOs_+_parents` GO terms in any subsequent enrichment analyses.

# How to cite
A publication is hopefully forthcoming which can be referred to when using this program. Until then, you can link to this repository.
