# inspector_clusto
Combination of bash and python. 

Input files: 
# limma-based csv file. 
# transcriptome gff3 file. 
# antismash gbk file. 

Uses a sliding window based approach to search for closely physically located genes (of given upregulation/expression) within the genome and allocates
to clusters. 
Includes code to determine whether said clusters appear in the commonly used Antismash-software output. 
