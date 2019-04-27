[TOC levels=1-3]: #

# Table of Contents
- [Name](#name)
- [Purpose](#purpose)
- [Data downloading](#data-downloading)
    - [Sequencing data](##sequencing-data)
    - [Reference genome data](##reference-genome-data)
- [Quality control and trimming](#quality-control-and-trimming)
- [Methylation analysis](#methylation-analysis)
    - [Genome indexing](##genome-indexing)
    - [Read alignment](##read-alignment)
    - [Aligned reads deduplication](##aligned-reads-deduplication)
    - [Methylation information extracting](##methylation-information-extracting)
- [Downstream analysis](#downstream-analysis)
- [Reference](#reference)
- [Author](#author)

# Name
5hmC finding procedure - bioinfo procedures of finding the 5hmC DMLs in mouse samples based on the ACE-seq data coming from [pushlished paper](https://www.nature.com/articles/nbt.4204).

# Purpose
* Learn the methylation data analysis procedure based on software tool [`bismark`](https://github.com/FelixKrueger/Bismark) and [`DSS`](http://bioconductor.org/packages/release/bioc/html/DSS.html) package.
* Get familiar with the Linux sever workspace and the basic bioinfo data analysis protocol.

# Data downloading
* The sequencing data comes from the NBT paper and the GEO accession number is [GSE116016](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116016)
* The Mus musculus reference genome data is downloaded from [Ensembl](https://www.ensembl.org/info/data/ftp/index.html)
* All data are downloaded via FTP address by a lightweight downloading tool [aria2](https://aria2.github.io/)
* Detailed information about data downloading ==> [Blog](https://www.jianshu.com/u/3fcc93cd84c1)
* Data used in this sample procedure:
    - Wild type: SRR7368841 & SRR7368842
    - TetTKO type: SRR7368845
    
## Sequencing data
Open browser, vist [EBI-ENA search page](https://www.ebi.ac.uk/ena) and search for the GEO accession. On the [result page](https://www.ebi.ac.uk/ena/data/view/PRJNA476795), we could get the expected FTP address. Use `aria2` to download the data and we could get `.fastq.gz` files.

```bash
mkdir -p $HOME/NBT_repeat/data/seq_data/WT_mESC_rep1
mkdir -p $HOME/NBT_repeat/data/seq_data/TetTKO_mESC_rep1
cd $HOME/NBT_repeat/data/seq_data/

aria2c -d ./WT_mESC_rep1/ -Z ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR736/001/SRR7368841/SRR7368841.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR736/002/SRR7368842/SRR7368842.fastq.gz
aria2c -d ./TetTKO_mESC_rep1/ -Z ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR736/005/SRR7368845/SRR7368845.fastq.gz
```
## Reference genome data 
Open browser, visit [ensembl download page](https://www.ensembl.org/info/data/ftp/index.html), and find the DNA `.fasta` file. On the [result page](ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/), we could get the expected FTP address.

```bash
mkdir -p $HOME/NBT_repeat/data/genome_data/
cd $HOME/NBT_repeat/data/genome_data/

# Basic command line download
aria2c -d ./ -Z ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.1.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.2.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.3.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.4.fa.gz
aria2c -d ./ -Z ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.5.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.6.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.7.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.8.fa.gz
aria2c -d ./ -Z ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.9.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.10.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.11.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.12.fa.gz
aria2c -d ./ -Z ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.13.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.14.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.15.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.16.fa.gz
aria2c -d ./ -Z ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.17.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.18.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.19.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.X.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.Y.fa.gz 

# Alternative option: script download
sh $HOME/Scripts/shell/genome_data_download.sh $HOME/NBT_repeat/data/genome_data/
```
# Quality control and trimming
Performing some quality control is highly recommended for all high throughput sequencing to tell straight whether the dataset is of good quality or whether there
were any fundamental problems with either the library or the sequencing itself.

Here, we use [FastQc](www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Trim Galore](www.bioinformatics.babraham.ac.uk/projects/trim_galore/) to do the quality control and adapter trimming respectively. Trimming reports could be found in sub-directory [`trimming-reports/`](trimming-reports/).

```bash
cd $HOME/NBT_repeat/data/seq_data/

# quality control to see the quality of the raw seq-data
fastqc --threads 3 ./WT_mESC_rep1/*.fastq.gz ./TetTKO_mESC_rep1/*.fastq.gz

# quality and adapter trimming, followd by fastqc operation
trim_galore -o ./WT_mESC_rep1/trimmed_data/ --fastqc ./WT_mESC_rep1/*.fastq.gz
trim_galore -o ./TetTKO_mESC_rep1/trimmed_data/ --fastqc ./TetTKO_mESC_rep1/*.fastq.gz
```

`FastQc` used options:
* `-t/--threads <int>`: Specifies the number of files which can be processed simultaneously.

`Trim Galore` used options:
* `-o/--output_dir <DIR>`: If specified all output will be written to this directory instead of the current directory. If the directory doesn't exist it will be created for you.  
* `--fastqc`: Run FastQC in the default mode on the FastQ file once trimming is complete.

# Methylation analysis
After the quality control and adapter trimming, BS-seq data analysis protocal based on [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) could be followed to do the methylation analysis for the ACE-seq data as well.

## Genome indexing 
Before alignments can be carried out, the genome of interest needs to be bisulfite converted *in-silico* and indexed. Based on `Bismark`, we can use both `Bowtie` and `Bowtie2` to do the indexing work.

```bash
bismark_genome_preparation --bowtie2 $HOME/NBT_repeat/data/genome_data/
```

`bismark_genome_preparation` used options:
* `--bowtie2`: This will create bisulfite indexes for Bowtie 2. (Default: ON).

## Read alignment
The core of the methylation data analysis procedure is to align the sequencing reads to the reference genome, and it is assumed that all data have been quality and adapter trimmed. The alignment reports could be found in sub-directory [`bismark-alignment-reports/`](bismark-alignment-reports/).

```bash
genome_path="$HOME/NBT_repeat/data/genome_data/"
cd $HOME/NBT_repeat/data/seq_data/

# read alignment
bismark -o ./WT_mESC_rep1/bismark_result/ --parallel 4 --genome_folder ${genome_path} ./WT_mESC_rep1/*.fq.gz
bismark -o ./TetTKO_mEC_rep1/bismark_result/ --parallel 4 --genoem_folder ${genome_path} ./TetTKO_mESC_rep1/*fq.gz

# merge the two WT_mESC_rep1 result .bam file
samtools cat -o SRX4241790_trimmed_bismark_bt2.bam ./WT_mESC_rep1/bismark_result/*.bam
```

`bismark` used options:
* `-o/--output_dir <dir>`: Write all output files into this directory. By default the output files will be written into the same folder as the input file(s). If the specified folder does not exist, Bismark will attempt to create it first.
* `--parallel <int>`: (May also be `--multicore <int>`) Sets the number of parallel instances of Bismark to be run concurrently. This forks the Bismark alignment step very early on so that each individual Spawn of Bismark processes only every n-th sequence (n being set by --parallel).
* `--genome_folder`: **Arguments**, The path to the folder containing the unmodified reference genome as well as the subfolders created by the `bismark_genome_preparation` script.

`samtools` used command and options:
* `cat`: **Command** concatenate BAMs
    - `-o`: output file format

## Aligned reads deduplication

## Methylation information extracting


# Downstream analysis

# Reference
* [Aria2 Manual](https://aria2.github.io/manual/en/html/index.html)

# Author 
Jihong Tang &lt;njutangjihong@gmail.com&gt;