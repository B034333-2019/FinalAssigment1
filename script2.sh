#!/bin/sh

#VARIABLES
Assignment_directory="/localdisk/home/s1239445/Assignment1"
sample_file_dir="/localdisk/data/BPSM/Assignment1/fastq";
fastqc_sample_file="/localdisk/data/BPSM/Assignment1/fastq/fqfiles";
fastqc_output_dir="/localdisk/home/s1239445/Assignment/Assignment1/FastQC_Output/"
tbb_genome_raw_data="/localdisk/data/BPSM/Assignment1/Tbb_genome/Tb927_genome.fasta.gz"
tbb_genome_unzipped="/localdisk/home/s1239445/Assignment/Assignment1/TBB_Database/Tb927_genome.fasta"
tbb_genome_database="/localdisk/home/s1239445/Assignment/Assignment1/TBB_Database/"

#PART THREE - CREATE GENOME DATABASE

#unzip TBB fasta genome file 
#gzip -d "/localdisk/home/s1239445/Assignment/Assignment1/TBB_Database/Tb927_genome.fasta"

#build Bowtie2 index database using bowtie2-build
bowtie2-build -f $tbb_genome_unzipped $tbb_genome_database 
