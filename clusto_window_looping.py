## Python code to slide over a modified gff3 file and pull out gene names as they appear in the sliding window (of a given nucleotide size)
import os

path = '/home/jack.a/inspector_clusto/inspector_package/above_4_below_negative_4/passed_to_awk'

folder = os.fsencode(path)

for file in os.listdir(folder):
    with open(file, "r") as a_file:
        for line in a_file:
                columns = line.split()
                gene_loc = float(columns[1]) ## Create the element gene_loc, the location of each gene.
                gene_id = columns[0] ## Create the element gene_id, the name and functional annoation of each gene.
                for cstart in range(0, 10000000, 15000):
                        if cstart < gene_loc < cstart+20000:
                                print(file, gene_id, cstart, cstart+20000)
