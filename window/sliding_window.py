## Python code to create a 90kbp sliding window that increments in 30kbp steps over a modified gff3 file to identify closely located genes. 

## First off we must open the modified gff3 file to access the elements within. 
with open("Penicillium_sp._X_above_four_tig1_onlymrna.gff3_gene_locations.txt", "r") as a_file:
        for line in a_file:
                columns = line.split()
                gene_loc = float(columns[1]) ## Create the element gene_loc, the location of each gene. 
                gene_id = columns[0] ## Create the element gene_id, the name and functional annoation of each gene. 
                for cstart in range(0, 10000000, 30000): ## This range will 
                        if cstart < gene_loc < cstart+90000:
                                print("the gene", gene_id, "lies in the window starting at", cstart, "bp and ending at", cstart+90000, "bp")

## This will print the standard output to the screen, which can then be collected by bash to manipulate etc. 
