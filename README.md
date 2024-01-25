# BINge
Clusters _de novo_ transcripts and gene models using a reference genome to provide gene-level counts

## Overview
BINge (**Bin** **G**enes for **E**xpression analyses) is a sequence clustering program which aims to group alternative isoforms, orthologs, and paralogs into sequence clusters.
Its algorithm is capable of doing this accurately while receiving multiple species' sequences as input, a feat that currently existing clustering programs cannot achieve.
In doing so, BINge enables differential gene expression analyses to take place that compare one or many species in a way that minimises biases associated with:
- Data mappability (e.g., reads from species A not aligning optimally to sequence models derived from species B), and
- Ortholog sequence length (e.g., a gene in species A may be longer than the ortholog in species B which means, at an equal sequencing depth and amount of expression, we might expect
more reads to be sequenced from species A's ortholog)

The solution to these biases is to not just map to a single species' reference models - we should map to all species' sequences. This way, we allow the reads from species A to align
optimally to its own sequences, and so on. Then we can tally the amount of alignments made to all the sequences in a cluster after they've been quantified taking sequence length bias into account. 

## Data required to run
Inputs can be categorised as being:
1. The input sequences to be clustered, which can consist of reference gene models with the possible addition of _de novo_ assembled transcripts.
2. One or more genome sequences to use while clustering, which can be provided with their GFF3 annotations (to assist the clustering process) or just as plain FASTA files which will be aligned against.

## How does BINge work?
BINge automates the alignment of the input sequences against the genome sequences. A **_bin_** will be positioned along the genome where these sequences align, or where a gene has been annotated onto the
genome in a GFF3 file. This **_bin_** loosely corresponds to an exon, and BINge will remember any sequences which align over that _binned_ position.

After this, BINge will look at each input sequence to see how many times other sequences were placed into the same bin with it. Fundamentally, this is measuring how many exons the input sequence shares in
common with any other sequence. When two sequences share many exon alignments in common with each other, they will be clustered together; this is the 'clustering by co-occurrence' algorithm that BINge employs.

Because BINge aligns each input sequence against many positions in each genome, but only looks at the best and near-best alignments, we can expect isoforms, orthologs, and paralog sequences to co-occur with each other
in bins. Hence, these sequences get clustered together!

## Installation
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

It has been developed on Linux and within Windows Subsystem for Linux (WSL). If you want to run BINge in WSL, you must ensure that developer mode is enabled so that symbolic links are functional in Windows.
You can do that by going to Settings > For developers > Toggle 'Developer Mode' On.

## How to use
On the command line, you can always ask BINge to provide help information by doing:

`python /location/of/BINge.py -h`

This provides a short form of the help information showing only the parameters that you should specify or may want to configure. Many other parameters
are hidden since their default values are what most people should use. If you want to see those, try:

`python /location/of/BINge.py --help-long`

Otherwise, BINge tries to keep the clustering process simple, and minimally requires you to ...

### Input sequences to cluster
The sequences you want to cluster can be in the form of FASTA sequence files or as pairs of genome FASTA and GFF3 annotation files (from which coding DNA sequences will be extracted).

On the command line, providing one or more FASTA files might look something like:

`python /location/of/BINge.py -i input1.fasta input2.fasta <...>`

Whereas providing genome FASTA and GFF3 annotation files for extraction is done like:

`python /location/of/BINge.py -i input1.fasta,input1.gff3 input2.fasta,input2.gff3 <...>`

And you can mix these inputs too, like:

`python /location/of/BINge.py -i input1.fasta input2.fasta,input2.gff3 <...>`

### Input genomes to use as targets
The genomes you want to use while clustering must be provided as FASTA files with or without a reference GFF3 annotation to use when positioning bins along the genome.

On the command line this would be similar to how we provided input sequences, with a mixed input looking like:

`python /location/of/BINge.py -i <...> -g genome1.fasta genome2.fasta,genome2.gff3 <...>`

### Configurable parameters
Most of BINge's parameters are tuned via objective evaluation to be set at appropriate values for your use.

Many of BINge's operations can be run in parallel, and you can use the `--threads` argument to indicate how many CPUs you want to run any parallelisable steps with. This extends
to GMAP alignment and external clustering.

For the external clustering of any sequences that do not align against your reference genomes, you need to tell BINge where to find either MMseqs2 or CD-HIT. You do that through
the `--clusterer` argument (to tell it whether you're using MMseqs2 or CD-HIT) and either the `--mmseqs` argument (to point it to folder storing the `mmseqs` executable file) or
the `--cdhit` argument (to point it to the folder storing the `cd-hit-est` executable file).

In terms of clustering behaviour, there are only three parameters that are of particular note and that you might want to modify. These are:

- `--gmapIdentity` which sets a minimum value for the best GMAP alignment of any sequence to be considered.
  - In practice, changing this value actually has very little impact since most recognised orthologs maintain high sequence similarity. The default of 0.95 works in most
instances, but you could lower this to 0.90 if you're clustering very distantly related organisms.
- `--clusterVoteThreshold` which sets the proportion of bins two sequences must co-occur within to cluster together. Based on objective evaluation, you should leave this at
0.66 or, if you want marginally tighter clustering (e.g., fewer clusters) you can lower this to 0.5. You probably shouldn't deviate from these recommendations.
- `--microbial` can be specified if you're clustering a microbial organism since it changes how any of your input GFF3s will be interpreted. If you're not using any
GFF3 files for your microbial organism, this parameter is not needed.

## Interpreting results
BINge will run everything within the specified output directory, and produce a file that will look something like `BINge_clustering_result.<hash>.tsv`. The hash uniquely identifies
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

## Filtering results
If you want to analyse your unbinned sequences, you may first want to filter them for quality. `BINge_filter.py` helps to do this by receiving your BINge cluster file, the sequences
you used as input for clustering, as well as a BLAST (or MMseqs2) outfmt6 results file and/or salmon alignment files. The filtering then can be used to remove short sequences, unbinned
clusters which do not have much read alignments, or those which do not receive significant E-value scores from alignment against a database of your choice.

## Picking representative sequences
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

## Using results for differential gene expression
The `BINge_counter.py` script will allow you to take BINge's clustering results, alongside salmon read alignments to the sequences you used as input for clustering, and produce files
which will enable DGE analysis. It generates three outputs (`.abundance`, `.counts`, and `.length`) and provides instructions for how to load the outputs into R as a DESeq2 dataset for analysis.

## Other utilities provided with BINge
There are several additional scripts to enable common use cases with BINge's clustering results. The most useful ones are...

### generate_annotation_table.py
This script will receive inputs including:

- The BINge representatives FASTA file
- The FASTA file used as a BLAST database/target
- The output from BLASTing your representatives FASTA to the database
- The `idmapping_selected.tab` file available from the UniProt FTP (https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz)
- A `go.obo` file available from the Gene Ontology website (https://geneontology.org/docs/download-ontology/ and https://purl.obolibrary.org/obo/go.obo)

You can specify the E-value cut-off to employ and how many results you want to tabulate for each sequence. Otherwise, it will use the input files to indicate the best hits for each sequence
alongside their GO annotations which can be useful for GOseq analysis as an example.

### tabulate_salmon_qc.py
This script will receive the stderr output files from salmon and produce quality control statistics to indicate:

1. What percentage of your reads mapped against the entire set of clustered sequences, and
2. What percentage of your read mappings remain after filtering of your clustering results

## A typical analysis pipeline
This section will give a brief runthrough of how BINge might be used for a DGE experiment. Step-by-step, you might:

1. _De novo_ assemble your RNAseq data using your pipeline of choice; I would perform a multiple assembly using Trinity, SOAPdenovo-Trans and others combined using EvidentialGene.
2. Obtain any reference genome FASTA files for the species you are investigating with or without their GFF3 annotations. If you did not have any reference genome for your species being investigated, find closely related species' genomes instead.
3. Provide your _de novo_ transcriptome assembly alongside any reference gene models you obtained as the sequences to be clustered (`-i`). Indicate any reference genomes you obtained as targets to be clustered against (`-g`).
4. Run `BINge.py` and note the file name of your output cluster file (found wherever you specified `-o`).
5. Decide if you want to include unbinned clusters in your analysis. If you do not, remove them from the output clustering file using the utility program `<to_be_coded.py>`.
6. BLAST (/ MMseqs2) query your sequences that were clustered against a database containing relevant sequences.
7. Salmon align your RNAseq reads to the sequences that were clustered.
8. If you chose to analyse unbinned clusters, you should filter your clustering result using `BINge_filter.py` using your BLAST and Salmon output files; otherwise, this step can be omitted.
9. Pick out representatives for each cluster using `BINge_representatives.py` using the same BLAST and Salmon output files.
10. Run the utility script `generate_annotation_table.py` to generate functional annotations and obtain putative gene names for each of your sequence clusters.
11. Generate a counts table for DGE analysis using `BINge_counter.py` using the Salmon output files you've already generated.
12. Run `tabulate_salmon_qc.py` using the stderr outputs from Salmon alignment to obtain a table indicating what percentage of reads aligned against your sequences. If you removed unbinned clusters, that's all you'll find of
interest. If you did not remove unbinned clusters, a second statistic will show what percentage of sequence alignments were removed during filtration.

## How to cite
A publication is hopefully forthcoming which can be referred to when using this program.
