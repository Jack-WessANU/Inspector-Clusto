##!/bin/bash
echo "RUNNING CLUSTER-FINDER INSPECTOR_CLUSTO"

echo "BUILDING SUBFOLDERS"

mkdir above_4_below_negative_4

cp input_data/inplantavsinvitro.csv above_4_below_negative_4

cp input_data/Penicillium_sp._X.gff3 above_4_below_negative_4

cp input_data/penx_master_table_inplantavsinvitro.csv above_4_below_negative_4

cd above_4_below_negative_4

## Convert our csv file to a tab-delimited one

sed 's/\,/\t/g' inplantavsinvitro.csv > inplantavsinvitro.txt

sed 's/\,/\t/g' penx_master_table_inplantavsinvitro.csv > penx_master_table_inplantavsinvitro.txt

## Write an if loop that if the LFC (column 2) is above a specified value, then pull the ID (column 1) to a new txt file for grepping against gff3 file

echo "CUTTING DOWN GFF3 TO JUST IMPORTANT PARTS + DIVIDING GFF3 UP INTO RESPECTIVE CONTIGS"

while read LINE; do

awk '{ minimumLFC=4 ;
if ($2 >= 4 || $2 <= -4)
    print $1 > "selective_geneids.txt"
}'

done<inplantavsinvitro.txt

## Use this new selective gene id file to modify a gff3 to only contain genes regulated above or below the specified spot

grep -f selective_geneids.txt Penicillium_sp._X.gff3 > Penicillium_sp._X_above_four_below_negative_four.gff3.txt

## Then grep it down to only mRNA

grep "mRNA" Penicillium_sp._X_above_four_below_negative_four.gff3.txt > Penicillium_sp._X_above_four_below_negative_fouronlymrna.gff3.txt

## This step will make sure that the files contain all of the functional annotation, before hand it would cut off between hypothetical and protein as there was a space!
## If the text is too unwieldy on the page, simply delete this line of code, it won't have an effect on the rest of the script

 awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9$10$11$12$13$14$15$16$17$18$19$20 }' Penicillium_sp._X_above_four_below_negative_fouronlymrna.gff3.txt > Penicillium_sp._X_above_four_below_negative_fouronlymrnacomb.gff3.txt


## Then we need to separate them into contigs

mkdir passed_to_awk

## for loop to do the separating, BUG; works but cant get the $i in output file name to be in middle of the file name where I want it to be.
echo "separating out tigs and cutting out the useless stuff"

for i in {1..7};

do

    grep 'tig0000000'$i Penicillium_sp._X_above_four_below_negative_fouronlymrnacomb.gff3.txt > passed_to_awk/Penicillium_sp._X_above_four_below_negative_fourtig$i

done

