# Table of Contents
- [Getting started](#getting-started)
- [Introduction](#introduction)
  - [Reason for BINge existing](#reason-for-binge-existing)
  - [What does BINge actually do?](#what-does-binge-actually-do)
- [Installation](#installation)
- [How to use](#how-to-use)
  - [Input sequences to cluster](#input-sequences-to-cluster)
  - [Input genomes to use as targets](#input-genomes-to-use-as-targets)
  - [Configurable parameters](#configurable-parameters)
  - [Interpreting results](#interpreting-results)
- [Filtering results](#filtering-results)
- [Picking representative sequences](#picking-representative-sequences)
- [Using results for differential gene expression](#using-results-for-differential-gene-expression)
- [Other utilities provided with BINge](#other-utilities-provided-with-binge)
  - [generate_annotation_table.py](#generate_annotation_tablepy)
  - [tabulate_salmon_qc.py](#tabulate_salmon_qcpy)
- [A typical analysis pipeline](#a-typical-analysis-pipeline)
- [How to cite](#how-to-cite)

# Getting started
```
# Download this repository [making sure you have prerequisite Python packages]
git clone https://github.com/zkstewart/BINge.git

# Cluster transcriptome sequences against a reference genome
python BINge/BINge.py -i /location/of/transcriptome1.fasta -g /location/of/genome1.fasta -o /location/to/write/outputs
# Cluster an assortment of input sequences against a reference genome
python BINge/BINge.py -i /location/of/transcriptome1.fasta /location/of/transcriptome2.fasta /location/of/genome1.fasta,/location/of/genome1.gff3 /location/of/genemodels.fasta \
    -g /location/of/genome1.fasta \
    -o /location/to/write/outputs
# Cluster sequences against multiple reference genomes with or without annotations
python BINge/BINge.py -i /location/of/transcriptome1.fasta \
    -g /location/of/genome1.fasta /location/of/genome2.fasta,/location/of/genome2.gff3 \
    -o /location/to/write/outputs

# Other parameters to include during clustering [for MMseqs2 OR CD-HIT clustering of unbinned sequences]
python BINge/BINge.py -i <...> -g <...> -o <...> --threads 8 --clusterer mmseqs-cascade --mmseqs /location/of/mmseqs/bin
python BINge/BINge.py -i <...> -g <...> -o <...> --threads 8 --clusterer cd-hit --cdhit /location/of/cd-hit/bin
```

# Introduction
## Reason for BINge existing
BINge (**Bin** **G**enes for **E**xpression analyses) is a sequence clustering program which aims to group alternative isoforms, orthologs, and paralogs into sequence clusters.
Its algorithm is capable of doing this accurately while receiving multiple species' sequences as input, a feat that currently existing clustering programs cannot achieve.
In doing so, BINge enables differential gene expression analyses to take place that compare one or many species in a way that minimises biases associated with:
- Data mappability (e.g., reads from species A not aligning optimally to sequence models derived from species B), and
- Ortholog sequence length (e.g., a gene in species A may be longer than the ortholog in species B which means, at an equal sequencing depth and amount of expression, we might expect
more reads to be sequenced from species A's ortholog)

The solution to these biases is to not just map to a single species' reference models - we should map to all species' sequences. This way, we allow the reads from species A to align
optimally to its own sequences, and so on. Then we can tally the amount of alignments made to all the sequences in a cluster after they've been quantified taking sequence length bias into account. 

## What does BINge actually do?
BINge automates the alignment of the input sequences against the genome sequences using GMAP. A **_bin_** will be positioned along the genome where these sequences align, or where a gene has been annotated onto the
genome in a GFF3 file. This **_bin_** loosely corresponds to an exon, and BINge will remember any sequences which align over that _binned_ position.

After this, BINge will look at each input sequence to see how many times other sequences were placed into the same bin with it. Fundamentally, this is measuring how many exons the input sequence shares in
common with any other sequence. When two sequences share many exon alignments in common with each other, they will be clustered together; this is the 'clustering by co-occurrence' algorithm that BINge employs.

Because BINge aligns each input sequence against many positions in each genome, but only looks at the best and near-best alignments, we can expect isoforms, orthologs, and paralog sequences to co-occur with each other
in bins. Hence, these sequences get clustered together.

# Installation
BINge is written in Python and requires a modern version 3 [note: development occured using >= 3.9, and I'm not sure how far backwards its compatibility goes].

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

And one of either:
- MMseqs2 (https://github.com/soedinglab/MMseqs2)
- CD-HIT (https://sites.google.com/view/cd-hit)

With all dependencies accounted for, BINge can be obtained by cloning this repository.

`git clone https://github.com/zkstewart/BINge.git`

BINge has been developed on Linux and within Windows Subsystem for Linux (WSL). If you want to run BINge in WSL, you must ensure that developer mode is enabled so that symbolic links are functional in Windows.
You can do that by going to Settings > For developers > Toggle 'Developer Mode' On.

# How to use
On the command line, you can always ask BINge to provide help information by doing:

`python /location/of/BINge.py -h`

This provides a short form of the help information showing only the parameters that you should specify or may want to configure. Many other parameters
are hidden since their default values are what most people should use. If you want to see those, try:

`python /location/of/BINge.py --help-long`

Otherwise, BINge tries to keep the clustering process simple, and minimally requires you to ...

## Input sequences to cluster
The sequences you want to cluster can be in the form of FASTA sequence files such as *de novo* assembled transcriptome(s) or as reference genome CDS models. Alternatively, you may use pairs of genome FASTA and GFF3 annotation files from which CDS' will be extracted.

On the command line, providing one or more inputs might look something like:

`python /location/of/BINge.py -i input1.fasta input2.fasta <...>`

Whereas providing genome FASTA and GFF3 annotation files for extraction is done like:

`python /location/of/BINge.py -i input1.fasta,input1.gff3 input2.fasta,input2.gff3 <...>`

And you can mix these inputs too, like:

`python /location/of/BINge.py -i input1.fasta input2.fasta,input2.gff3 <...>`

## Input genomes to use as targets
The genomes you want to use while clustering must be provided as FASTA files with or without a reference GFF3 annotation to use when positioning bins along the genome.

On the command line this would be similar to how we provided input sequences, with a mixed input looking like:

`python /location/of/BINge.py -i <...> -g genome1.fasta genome2.fasta,genome2.gff3 <...>`

## Configurable parameters
Many of BINge's operations can be run in parallel, and you can use the `--threads` argument to indicate how many CPUs you want to run any parallelisable steps with. This extends
to GMAP alignment and external clustering.

For the external clustering of any sequences that do not align against your reference genomes, you need to tell BINge where to find either MMseqs2 or CD-HIT. You do that through
the `--clusterer` argument (to tell it whether you're using MMseqs2 or CD-HIT) and either the `--mmseqs` argument (to point it to folder storing the `mmseqs` executable file) or
the `--cdhit` argument (to point it to the folder storing the `cd-hit-est` executable file).

In terms of clustering behaviour, most parameters are tuned via objective evaluation to be set at appropriate values for your use. But, there are three parameters of particular note
that you *might* want to modify. These are:

- `--gmapIdentity` which sets a minimum value for the best GMAP alignment of any sequence to be considered.
  - In practice, changing this value actually has very little impact since most recognised orthologs maintain high sequence similarity. The default of 0.95 works in most
instances, but you could lower this to 0.90 if you're clustering very distantly related organisms.
- `--clusterVoteThreshold` which sets the proportion of bins two sequences must co-occur within to cluster together. Based on objective evaluation, you should leave this at
0.66 or, if you want marginally tighter clustering (e.g., fewer clusters) you can lower this to 0.5. **You probably shouldn't deviate from these recommendations.**
- `--microbial` can be specified if you're clustering a microbial organism since it changes how any of your input GFF3s will be interpreted. If you're not using any
GFF3 files for your microbial organism, this parameter is not needed.

## Interpreting results
BINge will run everything within the specified `-o` output directory, and produce a result file that looks like `BINge_clustering_result.<hash>.tsv`. The hash uniquely identifies
your run parameters, and you can consult the `param_cache.json` file if you want to see what files you used as inputs and the GMAP identity and vote threshold parameters.

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
cluster group, with the sequence's ID shown followed by whether the cluster was binned or unbinned. In many cases you might want to excluded unbinned clusters since they did not
align well against any of the input genomes.

# Filtering results
The raw results of BINge clustering contain both *binned* and *unbinned* cluster sequences. `BINge_filter.py` facilitates the filtering of these results to either:

1. Remove *unbinned* clusters to focus analysis only on gene sequences present in a genome assembly; this intrinsically removes contaminant or spurious transcript assemblies
that might be present in your *de novo* transcriptomes.
2. Or, you can retain only the *unbinned* transcripts that appear to be real based on BLAST and/or read alignment evidence; this could be desirable if you know the genome
assembly is incomplete or not assembled at a chromosome level.

See example code below.

```
# Obtain a sequence FASTA file for use prior to cluster filtering and representative picking
## Note that this is a convenience script which merely concatenates the input sequence files together
python BINge/utilities/extract_clustered_sequences.py -i /location/to/write/outputs/BINge_clustering_result.<hash>.tsv \
    -d /location/to/write/outputs -o BINge_raw_sequences.nucl

# BLAST/MMseqs2 query against sequence database [e.g., RefSeq or UniRef90]
<refer to BLAST or MMsesqs2 manuals to produce outfmt6 results, making use of the BINge_raw_sequences.nucl file as input>

# Salmon align reads to sequences
<likewise, refer to Salmon manual, mapping each sample's reads to the BINge_raw_sequences.nucl file>

# Filter cluster file to drop unbinned clusters
python BINge/BINge_filter.py -i /location/to/write/outputs/BINge_clustering_result.<hash>.tsv \
    -f /location/to/write/outputs -o BINge_clustering_result.filtered.tsv -s cds --be_tolerant \
    --justDropUnbinned
# Or, filter cluster file to only drop low-quality unbinned clusters
python BINge/BINge_filter.py -i /location/to/write/outputs/BINge_clustering_result.<hash>.tsv \
    -f /location/to/write/outputs -o BINge_clustering_result.filtered.tsv -s cds --be_tolerant \
    --annot /location/of/genome1.gff3 --length <minimum protein size> \
    --blast /location/of/blast_or_mmseqs2/results.outfmt6 --evalue <E-value threshold> \
    --salmon /location/of/salmon/results/sample1/quant.sf /location/of/salmon/results/sample2/quant.sf <...> \
    --require1x --readLength <length in bp of your reads>
```

# Picking representative sequences
Sequence clusters can contain many members, but for interpretation of results you probably want to have one sequence you can use as the representative of that cluster; in essence,
we can say that all the sequences in a cluster are isoforms, orthologs, paralogs, or even just fragments of that chosen representative. You can hence functionally annotate the sequence
chosen as representative and just consider that (and not the other cluster members) for further interpretation and analysis.

`BINge_representatives.py` will handle this process for you. It will receive the same files used for filtering but uses that information to instead pick out the best representative from
each cluster. In this context, _best_ refers to the sequence which has the best score without ties for any of the metrics in the following order of priority:

1. If the sequence is from a reference annotation, we'd prefer to use it; otherwise, we'll choose the sequence with the ...
2. Highest bitscore from BLAST alignment
3. Most read alignments
4. Longest length

The resulting FASTA file will have sequence identifiers with a format like `>cluster-N representative=sequenceID`.

See example code below.

```
python BINge/BINge_representatives.py -i BINge_clustering_result.filtered.tsv -f /location/to/write/outputs -o BINge_filtered_sequences.nucl \
  --annot /location/of/genome1.gff3 \
  --blast /location/of/blast_or_mmseqs2/results.outfmt6 --evalue <E-value threshold> \
  --salmon /location/of/salmon/results/sample1/quant.sf /location/of/salmon/results/sample2/quant.sf <...>
```

# Using results for differential gene expression
The `utilities/tx2gene.py` script will reformat a BINge cluster file for use when loading Salmon read alignments into R, as exampled below.

`python BINge/utilities/tx2gene.py -i BINge_clustering_result.filtered.tsv -o BINge_clustering_result.filtered.tx2gene.tsv`

Then in R, you can follow the instructions of [the tximport vignette](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html) and do something like seen in the example code below.

```
# Locate all salmon quant files
files <- file.path(dir, "salmon", samples$run, "quant.sf")
# Load in tx2gene file
tx2gene <- read_delim(Tx2GENE_FILE, delim="\t")
# Obtain read counts
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
# Produce a DESeq2 dataset for analysis
dds = DESeqDataSetFromTximport(txi = txi.salmon, colData = <...>, design = ~ <...>)
```

# Other utilities provided with BINge
There are several additional scripts to enable common use cases with BINge's clustering results. The most useful ones are...

## generate_annotation_table.py
This script will receive inputs including:

- The BINge representatives FASTA file
- The FASTA file used as a BLAST database/target
- The output from BLASTing your representatives FASTA to the database
- The `idmapping_selected.tab` file available from the UniProt FTP (https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz)
- A `go.obo` file available from the Gene Ontology website (https://geneontology.org/docs/download-ontology/ and https://purl.obolibrary.org/obo/go.obo)

You can specify the E-value cut-off to employ and how many results you want to tabulate for each sequence. The result will be a tab-separated values (TSV) file listing the best hits
for each sequence, including the taxon of the hit, its E-value and bit score, as well as their GO annotations if any exist.

The `Best_mapped_GOs_+_parents` column indicates the GO annotations most suitable for GOseq analysis (for example). As the column header indicates, this contains not just the GO
annotations associated with the sequence as found in the `idmapping_selected.tab` file, it will also present all parents(/ancestors) of those terms as obtained from the `go.obo` file.

## tabulate_salmon_qc.py
If you've run Salmon with stderr redirected to individual files per sample, this script will receive those stderr output files to produce a table of quality control statistics which indicate:

1. What percentage of your reads mapped against the entire set of clustered sequences, and
2. What percentage of your read mappings remain after filtering of your clustering results

This information should likely be reported as a supplement in any manuscript making use of BINge to show what percentage of reads were filtered out as a result of running `BINge_filter.py`.

# A typical analysis pipeline
This section will give a brief runthrough of how BINge might be used for a DGE experiment. Step-by-step, you might:

1. *De novo* assemble your RNAseq data using your pipeline of choice; I would perform a multiple assembly using Trinity, SOAPdenovo-Trans and others combined using EvidentialGene.
2. Obtain any reference genome FASTA files for the species you are investigating with or without their GFF3 annotations. If you do not have any reference genome for your species being investigated, find closely related species' genomes instead.
3. Provide your *de novo* transcriptome assembly alongside any reference gene models you obtained as the sequences to be clustered (`-i`). Indicate any reference genomes you obtained as targets to be clustered against (`-g`).
4. Run `BINge.py` and note the file name of your output cluster file (i.e., the `BINge_clustering_result.<hash>.tsv` file found within the folder specified by `-o`).
5. Follow the example code in [Filtering results](#filtering-results) where `utilities/extract_clustered_sequences.py` is used to generate a single FASTA file for your gene sequences, then BLAST/MMseqs2 query these sequences against a database such as RefSeq or UniRef90.
6. Salmon align your RNAseq reads to these sequences.
7. Decide if you want to include unbinned clusters in your analysis.
    - If your genome assemblies are high-quality and complete, you *probably* should **drop** unbinned clusters (especially if you've used *de novo* transcript assemblies as input).
    - If you do not have a genome for your studied species or the genome assemblies are not high-quality, you *probably* should **not drop** unbinned clusters since they may contain relevant genes.
8. Filter results (see [Filtering results](#filtering-results)) to either drop unbinned clusters, or just filter to remove low-quality sequence clusters.
9. Pick representatives (see [Picking representative sequences](#picking-representative-sequences)) from the **filtered** cluster file using `BINge_representatives.py` with the BLAST and Salmon output files as evidence.
10. Load read counts into R (see [Using results for differential gene expression](#using-results-for-differential-gene-expression)) by use of `utilities/tx2gene.py` and the `tximport` package.

Some extra steps that may be relevant to you include:

1. Use `utilities/tabulate_salmon_qc.py` if you've run Salmon with per-sample stderr output to obtain a table indicating what percentage of reads aligned against your sequences.
2. Run `utilities/generate_annotation_table.py` to generate functional annotations and obtain putative gene names for each of your sequence clusters. The gene ontologies generated in this file can be obtained for use with GOseq analysis.

# How to cite
A publication is hopefully forthcoming which can be referred to when using this program. Until then, you can link to this repository.
