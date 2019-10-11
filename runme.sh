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
#echo -e "Fastqc running..."
#fastqc -o $fastqc_output_dir -t 8 --extract *.fq.gz 

#PART TWO - SUMMARY OF FASTQC OUTPUT

#concatenate all summary.txt files into one merged file advising on numbers and quality of fastqc output
#cat $fastqc_output_dir/*_L8_*_fastqc/summary.txt >> $fastqc_output_dir/fast_qc_overall_quality.txt
#display on-screen the quality check
#print $fastqc_output_dir/fast_qc_overall_quality.txt


#VARIABLES
tbb_genome_raw_data="/localdisk/data/BPSM/Assignment1/Tbb_genome/Tb927_genome.fasta.gz"
tbb_genome_unzipped="/localdisk/home/s1239445/Assignment/Assignment1/TBB_Database/Tb927_genome.fasta"
tbb_genome_database="/localdisk/home/s1239445/Assignment/Assignment1/TBB_Database/"

#PART THREE - CREATE GENOME DATABASE

#unzip TBB fasta genome file 
#echo -e "unzipping TTB raw data using gzip..."
#gzip -d "/localdisk/home/s1239445/Assignment/Assignment1/TBB_Database/Tb927_genome.fasta"

#build Bowtie2 index database using bowtie2-build
#echo e- "building indexed database using bowtie2-build..."
#bowtie2-build -f $tbb_genome_unzipped $tbb_genome_database 

#PART FOUR - ALIGN PAIRED SEQUENCES WITH INDEX

#VARIABLES
sam_bam_output_dir="/localdisk/home/s1239445/Assignment/Assignment1/SAM_BAM_Output"

#create output SAM files
#cd $sam_bam_output_dir
#touch 216_slender.sam 218_slender.sam 219_slender.sam 220_stumpy.sam 221_stumpy.sam 222_stumpy.sam sam_bam_output_216.sam sam_bam_output_218.sam sam_bam_output_219.sam sam_bam_output_220.sam sam_bam_output_221.sam sam_bam_output_222.sam

#run bowtie2
#echo -e- "bowtie2 running, aligning read pairs..."
#bowtie2 -x $tbb_genome_database -1 $sample_file_dir/216_L8_1.fq.gz -2 $sample_file_dir/216_L8_2.fq.gz -S $sam_bam_output_dir/sam_bam_output_216.sam
#bowtie2 -x $tbb_genome_database -1 $sample_file_dir/218_L8_1.fq.gz -2 $sample_file_dir/218_L8_2.fq.gz -S $sam_bam_output_dir/sam_bam_output_218.sam
#bowtie2 -x $tbb_genome_database -1 $sample_file_dir/219_L8_1.fq.gz -2 $sample_file_dir/219_L8_2.fq.gz -S $sam_bam_output_dir/sam_bam_output_219.sam
#bowtie2 -x $tbb_genome_database -1 $sample_file_dir/220_L8_1.fq.gz -2 $sample_file_dir/220_L8_2.fq.gz -S $sam_bam_output_dir/sam_bam_output_220.sam
#bowtie2 -x $tbb_genome_database -1 $sample_file_dir/221_L8_1.fq.gz -2 $sample_file_dir/221_L8_2.fq.gz -S $sam_bam_output_dir/sam_bam_output_221.sam
#bowtie2 -x $tbb_genome_database -1 $sample_file_dir/222_L8_1.fq.gz -2 $sample_file_dir/222_L8_2.fq.gz -S $sam_bam_output_dir/sam_bam_output_222.sam

#convert SAM to BAM format
#echo -e "samtools converting SAM to BAM files"
#samtools view -S -b $sam_bam_output_dir/sam_bam_output_216.sam > $sam_bam_output_dir/bam_output_216.bam
#samtools view -S -b $sam_bam_output_dir/sam_bam_output_218.sam > $sam_bam_output_dir/bam_output_218.bam
#samtools view -S -b $sam_bam_output_dir/sam_bam_output_219.sam > $sam_bam_output_dir/bam_output_219.bam
#samtools view -S -b $sam_bam_output_dir/sam_bam_output_220.sam > $sam_bam_output_dir/bam_output_220.bam
#samtools view -S -b $sam_bam_output_dir/sam_bam_output_221.sam > $sam_bam_output_dir/bam_output_221.bam
#samtools view -S -b $sam_bam_output_dir/sam_bam_output_222.sam > $sam_bam_output_dir/bam_output_222.bam

#sort BAM files
#echo -e "samtools sorting BAM files"
#samtools sort $sam_bam_output_dir/bam_output_216.bam -o $sam_bam_output_dir/bam_output_216_sorted.bam
#samtools sort $sam_bam_output_dir/bam_output_218.bam -o $sam_bam_output_dir/bam_output_218_sorted.bam
#samtools sort $sam_bam_output_dir/bam_output_219.bam -o $sam_bam_output_dir/bam_output_219_sorted.bam
#samtools sort $sam_bam_output_dir/bam_output_220.bam -o $sam_bam_output_dir/bam_output_220_sorted.bam
#samtools sort $sam_bam_output_dir/bam_output_221.bam -o $sam_bam_output_dir/bam_output_221_sorted.bam
#samtools sort $sam_bam_output_dir/bam_output_222.bam -o $sam_bam_output_dir/bam_output_222_sorted.bam

#index BAM files
echo -e "samtools indexing sorted BAM files"
samtools index $sam_bam_output_dir/bam_output_216_sorted.bam 
samtools index $sam_bam_output_dir/bam_output_218_sorted.bam
samtools index $sam_bam_output_dir/bam_output_219_sorted.bam 
samtools index $sam_bam_output_dir/bam_output_220_sorted.bam 
samtools index $sam_bam_output_dir/bam_output_221_sorted.bam 
samtools index $sam_bam_output_dir/bam_output_222_sorted.bam 

#PART FIVE - GENERATE COUNTS DATA

#VARIABLES
tbb_bedfile="/localdisk/data/BPSM/Assignment1/Tbbgenes.bed"
bedfile_output_dir="/localdisk/home/s1239445/Assignment/Assignment1/Bedtools_Output/"
slender_output="/localdisk/home/s1239445/Assignment/Assignment1/Bedtools_Output/final_output_slender.txt"
stumpy_output="/localdisk/home/s1239445/Assignment/Assignment1/Bedtools_Output/final_output_stumpy.txt"

#run bedtools
echo -e "bedtools multicov running..."

#output for slender reads
bedtools multicov -bams $sam_bam_output_dir/bam_output_216_sorted.bam $sam_bam_output_dir/bam_output_218_sorted.bam $sam_bam_output_dir/bam_output_219_sorted.bam -bed $tbb_bedfile  > $bedfile_output_dir/final_output_slender.txt	

#bedtools multicov -bams $sam_bam_output_dir/bam_output_216_sorted.bam -bed $tbb_bedfile > $bedfile_output_dir/final_output_216.txt

#bedtools multicov -bams $sam_bam_output_dir/bam_output_218_sorted.bam -bed $tbb_bedfile > $bedfile_output_dir/final_output_218.txt
#bedtools multicov -bams $sam_bam_output_dir/bam_output_219_sorted.bam -bed $tbb_bedfile > $bedfile_output_dir/final_output_219.txt

#output for stumpy reads
bedtools multicov -bams $sam_bam_output_dir/bam_output_220_sorted.bam $sam_bam_output_dir/bam_output_221_sorted.bam $sam_bam_output_dir/bam_output_222_sorted.bam  -bed $tbb_bedfile > $bedfile_output_dir/final_output_stumpy.txt

#bedtools multicov -bams $sam_bam_output_dir/bam_output_221_sorted.bam -bed $tbb_bedfile > $bedfile_output_dir/final_output_221.txt
#bedtools multicov -bams $sam_bam_output_dir/bam_output_222_sorted.bam -bed $tbb_bedfile > $bedfile_output_dir/final_output_222.txt

echo -e "bedtools done..."
