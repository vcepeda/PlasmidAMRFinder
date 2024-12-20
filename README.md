PlasmidAMRFinder
==========
PlasmidAMRFinder is a comprehensive pipeline designed to identify circular plasmids in bacterial genome assemblies (or complete genomes) and detect antimicrobial resistance genes. The pipeline supports both long and short reads, making it a versatile tool for plasmid identification and antimicrobial resistance characterization.

### Context
The COVID-19 pandemic underscored the critical need for robust genomic surveillance systems. Inspired by efforts to monitor SARS-CoV-2 variants, we established a central data hub coupled with an automated bioinformatics workflow for the collection, storage, and analysis of bacterial sequencing data from clinical patients. This project is aimed at the early detection and identification of emerging pathogens and outbreaks, and the monitoring of antimicrobial resistance (AMR) in bacterial pathogens. This initiative is part of the Genomic Pathogen Surveillance in Germany (GenSurv and GenSurv+) initiatives.

A central data hub was developed to collect and store short- and long-read bacterial sequencing data. The submission server was designed to efficiently handle data and associated metadata. Key features include customizable configurations, Docker-based deployment, and support for persistent file and database storage. Data collected by the hub is processed through a comprehensive Nextflow pipeline that automates genome assembly, annotation, plasmid identification, and AMR gene identification.

### Note:
This pipeline is under active development and currently in early-mid development stages. The pipeline will be extended to automate genome assembly and annotation. Contributions and suggestions are welcome to help improve and extend the functionality of PlasmidAMRFinder.

Pipeline Steps
------------

- Long Reads Analysis:
  - Gene prediction with Glimmer3
  - Identification of antimicrobial resistance genes using the CARD Database RGI
  - Long read alignment against assembly
  - Coverage analysis with Mosdepth
  - GC Content and GC Skew calculation
  - Identification of reads that overlap the gap in the plasmid, indicating circular reads

- Short Reads Analysis:
  - Quality control using FastQC
  - Trimming of adapters and low-quality sequences using Trimmomatic
  - Alignment of short reads to the assembly using BWA
  - Plasmid detection using PlasmidFinder
  - Coverage analysis with Mosdepth
  - GC Content and GC Skew calculation
  - Identification of antimicrobial resistance genes using the CARD Database RGI

Requirements
------------
- Linux or Mac OS
- Java 8.x
- Docker or Singularity container application or Conda package manager

Installation
------------

- Install Nextflow

```
curl -s https://get.nextflow.io | bash
```
This creates the Nextflow executable in the current directory.

- Download Pipeline
You can either get the latest version by cloning this repository:

```
git clone https://github.com/vcepeda/PlasmidAMRFinder
```
or download one of the releases.

- Download Dependencies
All the dependencies for this pipeline can be downloaded in a [docker (https://docs.docker.com/install/) container.

```
docker pull caspargross/plasmident
```
- Alternative dependency installations:
  - Singularity Container (docs/alternative_installation.md#singularity_container)
  - Use conda environment  (no docker)(docs/alternative_installation.md#conda_environment)

Run Application
---------------

The pipeline requires an input file with a sample ID (string) and paths for the assembly file in .fasta format and long reads in .fastq or .fastq.gz. The paths can either be absolute or relative to the launch directory. In normal configuration (with docker), it is not possible to follow symbolic links.

The file must be tab-separated with three columns:

```
sample_id	assembly_fasta	longread_fq
```
The file must not have a header line and start directly with the data. Here is an example file:

```
myid1	/path/to/assembly1.fasta	/path/to/reads1.fastq.gz
myid2	/path/to/assembly2.fasta	/path/to/reads2.fastq.gz
```
To include short reads analysis, add paths for the short reads in .fastq or .fastq.gz format. The file should now be tab-separated with four columns:

```
sample_id	assembly_fasta	longread_fq	shortread_fq
```
Example:


```
myid1	/path/to/assembly1.fasta	/path/to/reads1.fastq.gz	/path/to/shortreads1.fastq.gz
myid2	/path/to/assembly2.fasta	/path/to/reads2.fastq.gz	/path/to/shortreads2.fastq.gz
```
The pipeline is started with the following command:

```
nextflow run plasmident --input read_locations.tsv
```
For short reads analysis, use:

```
nextflow run plasmident --input read_locations.tsv --shortreads
```
There are other run profiles for specific environments.

### Optional Run Parameters
--outDir Path of output folder
--seqPadding Number of bases added at contig edges to improve long read alignment [Default: 1000]
--covWindow Moving window size for coverage and GC content calculation [Default: 50]
--max_cpu Number of threads used per process [Default: 4]
--max_memory Maximum amount of memory available
--targetCov Large read files are subsampled to this target coverage to speed up the process [Default: 50]

Results
-------

The pipeline creates the following output folders:

  - alignment: Contains the long read alignment (full genome)
  - coverage: Long read coverage for the whole input genome (compressed bedfile)
  - gc: Windowed GC content (full genome)
  - genes: Predicted gene locations
  - plasmids: Nucleotide sequences for all confirmed plasmids in separate FASTA files
  - resistances: GFF file with locations of identified antimicrobial resistance genes
Additionally, for short reads:

  - short_alignment: Contains the short read alignment (full genome)
  - short_coverage: Short read coverage for the whole input genome (compressed bedfile)
Additional File:
  - sampleID_summary.csv: Tabular text file with contig lengths, plasmid status, and identified antimicrobial resistance genes.
