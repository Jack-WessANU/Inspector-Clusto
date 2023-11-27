## Firstly, we have to change the master table from a comma separated file to space separated text file

mkdir clusto_after_R

cp input_data/penx_master_table_inplanta_vs_invitro.csv clusto_after_R

cd clusto_after_R

sed 's/\,/\t/g' penx_master_table_inplanta_vs_invitro.csv > penx_master_table_inplanta_vs_invitro.txt

sed -i '1d' penx_master_table_inplanta_vs_invitro.txt

## We then need to generate a list of gene ids from the master table and cast it against the antismash output to determine which were found by both antismash and clusto

awk '{print $1}' penx_master_table_inplanta_vs_invitro.txt > gene_ids.txt

while read LINE; do

    if  grep "$LINE" ../input_data/all_antismash.gbk;
    then
        echo $LINE "Y" >> gene_ids_antismash.txt
    else
        echo $LINE "N" >> gene_ids_antismash.txt
    fi

done<gene_ids.txt

## Then i needed to process teh clusto output so it can be combined with the antismash data and then join the master table
## Pull out single line-clusters into multiple lines

grep "%" 500.txt > clusto_summary.txt

## Remove [number] from Rstudio output

awk '{$1=""}1' clusto_summary.txt | awk '{$1=$1}1' > clusto_summary2.txt

## Separate out into clusters on one line each. 

cat clusto_summary2.txt | tr -s "%" "\n" > clusto_summary3.txt

## Code to combine values in columns into 1 column for ease of use

awk '{print $1, $2"-"$3"-"$4"-"$5}' clusto_summary3.txt > clusto_summary4.txt

## Code to combine the clusto and antismash

awk 'NR==FNR {h[$1] = $2; next} {print $1, $2"-"h[$1]}' clusto_summary4.txt gene_ids_antismash.txt > gene_ids_clusto_antismash.txt
