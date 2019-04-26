[TOC levels=1-3]: #

# Table of Contents
- [Name](#name)
- [Purpose](#purpose)
- [Data downloading](#data-downloading)
    - [Secquencing data](##sequencing-data)
    - [Reference genome data](##reference-genome-data)
- [Quality control](#quality-control)
- [Methylation analysis](#methylation-analysis)
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

## Scequencing data
Open browser, vist [EBI-ENA search page](https://www.ebi.ac.uk/ena) and search for the GEO accession. On the [result page](https://www.ebi.ac.uk/ena/data/view/PRJNA476795), we could get the expected FTP address.

```bash
mkdir -p $HOME/NBT_repeat/data/seq_data/WT_mESC_rep1
mkdir -p $HOME/NBT_repeat/data/seq_data/TetTKO_mESC_rep1
cd $HOME/NBT_repeat/data/seq_data/

aria2c -d ./WT_mESC_rep1/ -Z ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR736/001/SRR7368841/SRR7368841.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR736/002/SRR7368842/SRR7368842.fastq.gz

aria2c -d ./TetTKO_mESC_rep1 -Z ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR736/005/SRR7368845/SRR7368845.fastq.gz
```
## Reference genome data 
Open browser, visit [ensembl download page](https://www.ensembl.org/info/data/ftp/index.html), and find the DNA `.fasta` file. On the [result page](ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/), we could get the expected FTP address.

```bash
mkdir -p $HOME/NBT_repeat/data/genome_data/
cd $HOME/NBT_repeat/data/genome_data/

#Basic command line download
aria2c -d ./ -Z ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.1.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.2.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.3.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.4.fa.gz

aria2c -d ./ -Z ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.5.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.6.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.7.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.8.fa.gz

aria2c -d ./ -Z ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.9.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.10.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.11.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.12.fa.gz

aria2c -d ./ -Z ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.13.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.14.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.15.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.16.fa.gz

aria2c -d ./ -Z ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.17.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.18.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.19.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.X.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.Y.fa.gz 

# Alternative option: script download
sh $HOME/Scripts/shell/genome_data_download.sh $HOME/NBT_repeat/data/genome_data/
```
# Quality control

# Methylation analysis

# Downstream analysis

# Reference

# Author 
Jihong Tang &lt;njutangjihong@gmail.com&gt;