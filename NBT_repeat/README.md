[TOC levels=1-3]: #

# Table of Contents
- [Name](#name)
- [Purpose](#purpose)
- [Data downloading](#data-downloading)
    - [Sequencing data](#sequencing-data)
    - [Reference genome data](#reference-genome-data)
- [Quality control and trimming](#quality-control-and-trimming)
- [Methylation analysis](#methylation-analysis)
    - [Genome indexing](#genome-indexing)
    - [Read alignment](#read-alignment)
    - [Aligned reads deduplication](#aligned-reads-deduplication)
    - [Methylation information extracting](#methylation-information-extracting)
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
Mammalian genomes are so huge that it is rather unlikely to encounter several genuinely independent fragments which align to the very same genomic position. It is much more likely that such reads are a result of PCR amplification. For large genomes, removing duplicate
reads is therefore a valid route to take.

The deduplication step could be finished through the `deduplicate_bismark` script, and the deduplication report could be found in sub-directory [`deduplication-reports/`](##deduplication-reports/).

```bash
cd $HOME/NBT_repeat/data/seq_data/
mkdir -p ./WT_mESC_rep1/deduplicated_result/
mkdir -p ./TetTKO_mESC_rep1/deduplicated_result/

# aligned reads deduplication
deduplicate_bismark --bam --output_dir ./WT_mESC_rep1/deduplicated_result/ ./WT_mESC_rep1/bismark_result/SRX4241790_trimmed_bismark_bt2.bam

deduplicate_bismark --bam --output_dir ./TetTKO_mESC_rep1/deduplicated_result/ ./TetTKO_mESC_rep1/bismark_result/*.bam
```

`deduplicate_bismark` used options:
* `--bam`: The output will be written out in BAM format instead of the default SAM format.
* `--output_dir`: Output directory, either relative or absolute. Output is written to the current directory if not specified explicitly.

## Methylation information extracting
The `bimark_methylation_extractor` script in the `Bismark` package could extract the methylation information from the alignment result files and act as the endpoint of the `Bismark` package. In addition, the methylation information could be easily transfered into other format facilitating the downstream analysis in this script.

```bash
genome_path= "$HOME/NBT_repeat/data/genome_data/"
cd $HOME/NBT_repeat/data/seq_data/

# methylation information extracting
bismark_methylation_extractor --single-end --gzip --parallel 4 --bedGraph \
--cytosine_report --genome_folder ${genome_path} \
-o ./WT_mESC_rep1/deduplicated_result/ ./WT_mESC_rep1/deduplicated_result/*.bam

bismark_methylation_extractor --single-end --gzip --parallel 4 --bedGraph \
--cytosine_report --genome_folder ${genome_path} \
-o ./TetTKO_mESC_rep1/deduplicated_result/ ./TetTKO_mESC_rep1/deduplicated_result/*.bam
```

`bismark_methylation_extractor` used options:
* `--single-end`: Input file(s) are Bismark result file(s) generated from single-end read data. If neither -s nor -p is set the type of experiment will be determined automatically.
* `--gzip`: The methylation extractor files (CpG_OT_..., CpG_OB_... etc) will be written out in a GZIP compressed form to save disk space. This option is also passed on to the genome-wide cytosine report.
* `--parallel <int>`: May also be `--multicore <int>`. Sets the number of cores to be used for the methylation extraction process.
* `--bedGraph`: After finishing the methylation extraction, the methylation output is written into a sorted bedGraph file that reports the position of a given cytosine and its methylation state.
* `--cytosine_report`: After the conversion to bedGraph has completed, the option `--cytosine_report` produces a genome-wide methylation report for all cytosines in the genome. By default, the output uses 1-based chromosome coordinates (zero-based start coords are optional) and reports CpG context only (all cytosine context is optional). The output considers all Cs on both forward and reverse strands and reports their position, strand, trinucleotide content and methylation state (counts are 0 if not covered).
* `--genome_folder <path>`: Enter the genome folder you wish to use to extract sequences from (full path only).
* `-o/ --output`: Allows specification of a different output directory (absolute or relative path). If not specified explicitly, the output will be written to the current directory.

# Downstream analysis
Based on the methylation information result got from the methylation analysis procedure, we could make some basic downstream analysis work including finding specific locus and detecting differential methylation loci (DML) or differential methylation regions (DMR). Here we use R package [`DSS`](http://bioconductor.org/packages/release/bioc/html/DSS.html) to make the differential methylation analysis.

## Input data preparation
`DSS` requires data from each BS-seq like experiment to be summarized into following information for each CG position: chromosome number, genomic coordinate, total number of reads, and number of reads showing methylation. 

The required input data could be transfered from the bismark result `.cov` file since the count file contain following columns: chr, start, end, methylation%,
count methylated, count unmethylated.

```bash
mkdir -p $HOME/NBT_repeat/R_analysis/WT_data
mkdir -p $HOME/NBT_repeat/R_analysis/TetTKO_data
cd $HOME/NBT_repeat/R_analysis/
# store the file path information
file_WT_path="$HOME/NBT_repeat/data/seq_data/WT_mESC_rep1/deduplicated_result/SRX4241790_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"
file_TetTKO_path="$HOME/NBT_repeat/data/seq_data/TetTKO_mESC_rep1/deduplicated_result/SRR7368845_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"
# copy the result file to the R analysis folder
cp file_WT_path ./WT_data/
cp file_TetTKO_path ./TetTKO_data/
# unzip
gunzip -d ./WT_data/SRX4241790_trimmed_bismark_bt2.deduplicated.bismark.cov.gz
gunzip -d ./TetTKO_data/SRR7368845_trimmed_bismark_bt2.deduplicated.bismark.cov.gz
# transfer the .cov file to .txt file
cp ./WT_data/SRX4241790_trimmed_bismark_bt2.deduplicated.bismark.cov ./WT_data/SRX4241790_methylation_result.txt
cp ./TetTKO_data/SRR7368845_trimmed_bismark_bt2.deduplicated.bismark.cov ./TetTKO_data/SRR7368845_methylation_result.txt

# basic command line environment procedure 
R
```
```r
library(tidyr)
library(dplyr)

file_names <- c("./WT_data/SRX4241790_methylation_result.txt", "./TetTKO_data/SRR7368845_methylation_result.txt")

func_read_file <- function(file_name){
	dir_vec <- strsplit(file_name, split = "/")[[1]]
	len <- length(dir_vec)
	file_prefix = substring(dir_vec[len], 0, nchar(dir_vec[len]) - 4)
	file_save_path = substring(file_name, 0, nchar(file_name) - nchar(dir_vec[len]))
	print(paste("File", file_name, "is being importing and this may take a while..."), sep = "")
	rawdata_df <- read.table(file_name, header = F, stringsAsFactors = F)
	print("Importing file is finished!")
	colnames(rawdata_df) <- c("chr", "start", "end", "methyl%", "methyled", "unmethyled")
	write.table(rawdata_df, paste(file_save_path, file_prefix, "_transfered.txt", sep = ""), row.names = F )
}

lapply(file_names, func_read_file)

q()
```

```bash
# Rscript procedure
Rscript $HOME/Scrpits/R/bismark_result_transfer.R ./WT_data/SRX4241790_methylation_result.txt ./TetTKO_data/SRR7368845_methylation_result.txt
```

## DML/DMR detection
After the input data preparation, we could use `DSS` package to find DMLs or DMRs. The detailed result files of both DML and DMR detection could be found in the sub-directory [`methylation-differential-analysis-reports/`](methylation-differential-analysis-reports/).

```bash

```
# Reference
* [Aria2 Manual](https://aria2.github.io/manual/en/html/index.html)
* [Bismark User Guide](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html)
* [DSS package Manual](http://bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.pdf)

# Author 
Jihong Tang &lt;njutangjihong@gmail.com&gt;