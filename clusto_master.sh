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

grep "mRNA" Penicillium_sp._X.gff3 > Penicillium_sp._X_above_four_below_negative_fouronlymrna.gff3.txt

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

## Step X. Modify gff3 to remove un-needed data and calculate mid-point of gene

for FILE in passed_to_awk/*;

do

    awk -i inplace '{print $9, ($4+$5)/2}' $FILE

done

## The final selected-gene-location files will be in the passed_to_awk directory for use

cd passed_to_awk

## Then use this to run the window code and output to a new file

echo "SLIDING THE WINDOW FROM ATTACHED PYTHON-SCRIPT OVER MODIFIED GFF3"

python3 ~/inspector_clusto/inspector_package/scripts/clusto_window_looping.py > ../mega_clusters.txt

cd ..

## Now we need to separate the cluster mega into contig-clusters
## Make a new directory for the new files

mkdir separated_clusters

for i in {1..8}
do
    grep 'tig'$i mega_clusters.txt > separated_clusters/'clustersontig'$i
done

## We can then move into the new directory

cd separated_clusters

## I then add '.txt' to the end of each file as I had to remove it before

find . -type f -exec bash -c 'mv "$0" "$0.txt"' {} \;

#3
## This first loop is annotated to explain their purposes

echo "WINDOW HAS BEEN RUN, SEPARATING + ANNOTATING WINDOW-OUTPUT, REMOVING EMPTY AND SINGLE-GENE WINDOW-FILES + TAKING INFORMATION FROM ANTISMASH GBK TO COMPARE
    ANTISMASH TO INSPECTOR_CLUSTO"

for c in {1..8}
do
    mkdir $c
    mkdir tidied$c
    mv clustersontig$c.txt $c
    cd $c

    for i in {0..10000000..30000} ## Make files for each possible window
    do
        grep ' '$i' ' clustersontig* > $i ## Pull out genes to individual files that exist within those windows
    done

    for FILE in *; do         awk -i inplace '{print $2}' $FILE; done
    for FILE in *; do         sed -i 's/\;/\t/g' $FILE; done
    for FILE in *; do         awk -i inplace '{print $1}' $FILE; done
    for FILE in *; do         sed -i 's/.*=//' $FILE; done

    for FILE in *
    do
        grep -f $FILE ~/inspector_clusto/inspector_package/above_4_below_negative_4/penx_master_table_inplantavsinvitro.txt > ../tidied$c/$FILE.txt
    done

    cd ..
    cd tidied$c

    find . -type f -empty -print -delete > /dev/null ## Remove all empty and therefore redundant files (empty windows) 
    
    for filename in ./*.txt; do
            if [ "$( wc -l <"$filename" )" -eq 1 ]; then
                rm -f "$filename"
            fi
    done

    find -name '*.txt' -execdir bash -c \ 'mv -v "$0" "${0%.txt}_$(wc -l < "$0").txt"' {} \; > /dev/null
    
    rm clustersontig$c*
    
    for file in ./*;do 

    while read LINE; do 
        gene_names=$(awk '{print $1}' $file)
        if  grep "$gene_names" ~/inspector_clusto/inspector_package/input_data/all_antismash.gbk; 
        then
            awk '{print "Y\t"$0}' $file > $file.inspector_output
        else
            awk '{print "N\t"$0}' $file > $file.inspector_output
        fi

    done<$file

    done


    for file in ./*.inspector_output
    do
        sort -o $file $file
        echo " "
        echo " "
        echo " "
        echo "###############################"$file"###############################"
        awk 'NR==1{first = $11} END{print $11 - first,"bp",NR,"genes"}' $file
        awk '{print $10,$2,"LFC",$3,"AveExpr",$4,"antismash?",$1}' $file 
    
    done
    cd ..
done
