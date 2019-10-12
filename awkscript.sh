#VARIABLES
final_output_dir="/localdisk/home/s1239445/Assignment/Assignment1/Final_Output"
tbb_bedfile="/localdisk/data/BPSM/Assignment1/Tbbgenes.bed"
bedfile_output_dir="/localdisk/home/s1239445/Assignment/Assignment1/Bedtools_Output/"

#Summing bedtools multicov output to give average of three slender sequences
awk '{sum = $7 + $8 + $9; avg = sum / 3; print $4"\t"avg}' Bedtools_Output/final_output_slender.txt > $final_output_dir/slender_mean_final_output.txt

#Summing bedtools multicov output to give average of three stumpy sequences
awk '{sum = $7 + $8 + $9; avg = sum / 3;  print $4"\t"avg}' Bedtools_Output/final_output_stumpy.txt > $final_output_dir/stumpy_mean_final_output.txt
 
echo - "Final Gene names and Average Counts for slender and stumpy sequences can be found in: $final_output_dir"
