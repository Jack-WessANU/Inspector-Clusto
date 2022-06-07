##!/bin/bash

## This first loop is annotated to explain their purposes

mkdir 1 ## Make a directory for output
mv clustersontig1.txt 1 ## Move the cluster file into it 
cd 1
for i in {0..10000000..30000} ## Make files for each possible window
do
    grep ' '$i' ' clustersontig1.txt > $i ## Pull out genes to individual files that exist within those windows
done
find . -type f -empty -print -delete ## Remove all empty and therefore redundant files (empty windows)
cd ..
            ## Repeat for all contigs
mkdir 2
mv clustersontig2.txt 2
cd 2
for i in {0..10000000..30000}
do
    grep ' '$i' ' clustersontig2.txt > $i
done
find . -type f -empty -print -delete
cd ..

mkdir 3
mv clustersontig3.txt 3
cd 3
for i in {0..10000000..30000}
do
    grep ' '$i' ' clustersontig3.txt > $i
done
find . -type f -empty -print -delete
cd ..

mkdir 4
mv clustersontig4.txt 4
cd 4
for i in {0..10000000..30000}
do
    grep ' '$i' ' clustersontig4.txt > $i
done
find . -type f -empty -print -delete
cd ..

mkdir 5
mv clustersontig5.txt 5
cd 5
for i in {0..10000000..30000}
do
    grep ' '$i' ' clustersontig5.txt > $i
done
find . -type f -empty -print -delete
cd ..

mkdir 6
mv clustersontig6.txt 6
cd 6
for i in {0..10000000..30000}
do
    grep ' '$i' ' clustersontig6.txt > $i
done
find . -type f -empty -print -delete
cd ..

mkdir 7
mv clustersontig7.txt 7
cd 7
for i in {0..10000000..30000}
do
    grep ' '$i' ' clustersontig7.txt > $i
done
find . -type f -empty -print -delete
cd ..

mkdir 8
mv clustersontig8.txt 8
cd 8
for i in {0..10000000..30000}
do
    grep ' '$i' ' clustersontig8.txt > $i
done
find . -type f -empty -print -delete
cd ..

## This is just to show a distance between the empty files being deleted and the cluster count of
## the following for loop. 

echo   '        ###SEPARATION OF EMPTY AND CONTAINING###
        #########################################
        #########################################
        #########################################'


## This will print the number of genes in each window to the standard input, which is then copied
## to excel for investigation
for i in {1..8}
do
    for file in $i/*
    do
        wc -l $file 
    done
done
