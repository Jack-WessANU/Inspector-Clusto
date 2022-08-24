##!/bin/bash
echo "RUNNING CLUSTER-FINDER INSPECTOR_CLUSTO"

echo "BUILDING SUBFOLDERS"

mkdir backbones

cp input_data/inplantavsinvitro.csv backbones

cp input_data/Penicillium_sp._X.gff3 backbones

cp input_data/penx_master_table_inplantavsinvitro.csv backbones

cd backbones

## Convert our csv file to a tab-delimited one

sed 's/\,/\t/g' inplantavsinvitro.csv > inplantavsinvitro.txt

sed 's/\,/\t/g' penx_master_table_inplantavsinvitro.csv > penx_master_table_inplantavsinvitro.txt

## Then grep it down to only mRNA

grep "mRNA" Penicillium_sp._X.gff3 > Penicillium_sp._X_onlymrna.gff3

## This step will make sure that the files contain all of the functional annotation, before hand it would cut off between hypothetical and protein as there was a space!
## If the text is too unwieldy on the page, simply delete this line of code, it won't have an effect on the rest of the script

awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9$10$11$12$13$14$15$16$17$18$19$20 }' Penicillium_sp._X_onlymrna.gff3 > Penicillium_sp._X_onlymrna_modded.gff3

## Then we need to separate them into contigs

mkdir passed_to_awk

## for loop to do the separating, BUG; works but cant get the $i in output file name to be in middle of the file name where I want it to be.
echo "separating out tigs and cutting out the useless stuff"

for i in {1..7};

do

    grep 'tig0000000'$i Penicillium_sp._X_onlymrna_modded.gff3 > passed_to_awk/Penicillium_sp._X_onlymrna_modded.gff3.tig$i

done

for FILE in passed_to_awk/*;

do

    awk -i inplace '{print $9, $4, $5}' $FILE
    sed -i 's/\;/\t/g' $FILE
    awk -i inplace '{print $1, $NF-1, $NF1}' $FILE
    sed -i 's/\=/\t/g' $FILE
    sed -i 's/\ /\t/g' $FILE
    awk -i inplace '{print $2, $(NF-1), $NF}' $FILE
    awk -i inplace 'FNR==NR { a[$1]=$2; next } $1 in a { print $1, a[$1], $2, $3 }' penx_master_table_inplantavsinvitro.txt $FILE 

done



