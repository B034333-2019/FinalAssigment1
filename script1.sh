#!/bin/sh

#VARIABLES
Assignment_directory="/localdisk/home/s1239445/Assignment1"
sample_file_dir="/localdisk/data/BPSM/Assignment1/fastq";
fastqc_sample_file="/localdisk/data/BPSM/Assignment1/fastq/fqfiles";
fastqc_output_dir="/localdisk/home/s1239445/Assignment/Assignment1/FastQC_Output/"

#PART ONE - FASTQC

#checking to see if can access sample file
#more $fastqc_sample_file

#moving to sample directory 
cd $sample_file_dir

#Run fastqc into $fastqc_output_dir with 8 threads to speed up simultaneous file creation, leaving output files unzipped
#fastqc -o $fastqc_output_dir -t 8 --extract *.fq.gz 

#PART TWO - SUMMARY OF FASTQC OUTPUT

#concatenate all summary.txt files into one merged file advising on numbers and quality of fastqc output
cat $fastqc_output_dir/*_L8_*_fastqc/summary.txt >> $fastqc_output_dir/fast_qc_overall_quality.txt
#display on-screen the quality check
print $fastqc_output_dir/fast_qc_overall_quality.txt

