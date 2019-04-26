#!/bin/bash/
# This simple script is used to download the mus musculus reference genome data from ensembl site to a given directory.
# Argument: [1] the given directory to store the downloaded data files
# Author: Jihong Tang
# Date: April 23, 2019
 
prefix='ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.'
suffix='.fa.gz'

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y
do 
    ftpadd="${prefix}${chr}${suffix}"
   # echo $ftpadd
    aria2c -d $1 -Z $ftpadd
done

