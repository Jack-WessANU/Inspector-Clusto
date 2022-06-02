##!/bin/bash

mkdir above_4

cp input_data/inplantavsinvitro.csv above_4

cp input_data/Penicillium_sp._X.gff3 above_4

cd above_4

## Convert our csv file to a tab-delimited one. 

sed 's/\,/\t/g' inplantavsinvitro.csv > inplantavsinvitro.txt

## Write an if loop that if the LFC (column 2) is above a specified value, then pull the ID (column 1) to a new txt file for grepping against gff3 file. 

while read LINE; do

awk '{ minimumLFC=4 ;

if ($2 >= minimumLFC)

    print $1 > "selective_geneids.txt"

}'

done<inplantavsinvitro.txt

## Use this new selective gene id file to modify a gff3 to only contain genes regulated above or below the specified spot. 

grep -f selective_geneids.txt Penicillium_sp._X.gff3 > Penicillium_sp._X_above_four.gff3.txt

## Then grep it down to only mRNA 

grep "mRNA" Penicillium_sp._X_above_four.gff3.txt > Penicillium_sp._X_above_four_onlymrna.gff3.txt

## This step will make sure that the files contain all of the functional annotation, before hand it would cut off between hypothetical and protein as there was a space!
## If the text is too unwieldy on the page, simply delete this line of code, it won't have an effect on the rest of the script. 

 awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9$10$11$12$13$14$15$16$17$18$19$20 }' Penicillium_sp._X_above_four_onlymrna.gff3.txt > Penicillium_sp._X_above_four_onlymrnacomb.gff3.txt


## Then we need to separate them into contigs. 

mkdir passed_to_awk

## for loop to do the separating, BUG; works but cant get the $i in output file name to be in middle of the file name where I want it to be. 

for i in {1..8};

do

    grep 'tig0000000'$i Penicillium_sp._X_above_four_onlymrnacomb.gff3.txt > passed_to_awk/Penicillium_sp._X_above_four_tig$i

done

## Step X. Modify gff3 to remove un-needed data and calculate mid-point of gene. 

for FILE in passed_to_awk/*;

do

    awk -i inplace '{print $9, ($4+$5)/2}' $FILE

done 

## The final selected-gene-location files will be in the passed_to_awk directory for use. 

cd passed_to_awk

## Then use this to run the window code and output to a new file 

python3 ~/inspector_clusto/inspector_package/scripts/clusto_window_looping.py > ../mega_clusters.txt

cd ..

## Now we need to separate the cluster mega into contig-clusters. 
## Make a new directory for the new files

mkdir separated_clusters

for i in {1..8}
do
    grep 'tig'$i mega_clusters.txt > separated_clusters/'clustersontig'$i
done

## We can then move into the new directory

cd separated_clusters

## I then add '.txt' to the end of each file as I had to remove it before. 

find . -type f -exec bash -c 'mv "$0" "$0.txt"' {} \;

python3 ~/inspector_clusto/inspector_package/scripts/clusto_grab.py > jim.txt


## Then create a list of files to use for grepping the next file. 

for i in {0..10000000..30000};
do
    echo $i >> cluster_starts.txt
done
