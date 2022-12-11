sink("500.txt") ## This code will print the entire standard output to a file called "GAP.txt"

## Loading in modded versions of the gff3 files for each PenX contig into dataframes
tig_1_modded.gff3 <- read.csv("Penicillium_sp._X_onlymrna_modded.gff3.tig1", sep = " ", header = FALSE)
tig_2_modded.gff3 <- read.csv("Penicillium_sp._X_onlymrna_modded.gff3.tig2", sep = " ", header = FALSE)
tig_3_modded.gff3 <- read.csv("Penicillium_sp._X_onlymrna_modded.gff3.tig3", sep = " ", header = FALSE)
tig_4_modded.gff3 <- read.csv("Penicillium_sp._X_onlymrna_modded.gff3.tig4", sep = " ", header = FALSE)
tig_5_modded.gff3 <- read.csv("Penicillium_sp._X_onlymrna_modded.gff3.tig5", sep = " ", header = FALSE)
tig_6_modded.gff3 <- read.csv("Penicillium_sp._X_onlymrna_modded.gff3.tig6", sep = " ", header = FALSE)
tig_7_modded.gff3 <- read.csv("Penicillium_sp._X_onlymrna_modded.gff3.tig7", sep = " ", header = FALSE)

## Adding an extra column to each dataframe with the number of its column
tig_1_modded.gff3["contig"] <- "1" 
tig_2_modded.gff3["contig"] <- "2" 
tig_3_modded.gff3["contig"] <- "3" 
tig_4_modded.gff3["contig"] <- "4" 
tig_5_modded.gff3["contig"] <- "5" 
tig_6_modded.gff3["contig"] <- "6" 
tig_7_modded.gff3["contig"] <- "7" 

## Adding column names to dataframes
colnames(tig_1_modded.gff3) <- c("gene_id","LFC","bp_start","bp_end","contig")
colnames(tig_2_modded.gff3) <- c("gene_id","LFC","bp_start","bp_end","contig")
colnames(tig_3_modded.gff3) <- c("gene_id","LFC","bp_start","bp_end","contig")
colnames(tig_4_modded.gff3) <- c("gene_id","LFC","bp_start","bp_end","contig")
colnames(tig_5_modded.gff3) <- c("gene_id","LFC","bp_start","bp_end","contig")
colnames(tig_6_modded.gff3) <- c("gene_id","LFC","bp_start","bp_end","contig")
colnames(tig_7_modded.gff3) <- c("gene_id","LFC","bp_start","bp_end","contig")

## Combining all 7 contig-dataframes into 1 super contig-dataframe
penx_transcriptome_modded.gff3 <- rbind(tig_1_modded.gff3,tig_2_modded.gff3,tig_3_modded.gff3,tig_4_modded.gff3,tig_5_modded.gff3,tig_6_modded.gff3,tig_7_modded.gff3)

clus_count <- 0 ## Create variable clus_count that will increase by one every time a new cluster is defined, basically counting the number of clusters in genome
list_of_gene_count_all_clusters <- list() ## Create a list that will be appended to with the number of genes in each cluster
for (i in 1:7){
  master_mix <- penx_transcriptome_modded.gff3[penx_transcriptome_modded.gff3$contig==i,] ## Will go through the mega-dataframe contig by contig based on column and complete as below
  for (row in 1:nrow(master_mix)){
    gap <- 750 ## Variable "gap" denotes the size of the search area between genes
    LFC_curr <- master_mix[row, "LFC"] ## This variable is used to find backbone genes, any above or below the values below are considered backbones and used as centres of clusters
    if(LFC_curr >= 4 | LFC_curr <= -4) {
      genes_in_this_cluster <- list()
      
      num_r <- 1
      num_f <- 1
      first_r <- row-1
      first_n <- row+1
      pos_f <- row+2
      pos_r <- row-2
      pos_c <- row
      pos_c_b <- row-1
      pos_c_f <- row+1
      clus_count <- clus_count+1
      up_reg <- 0
      down_reg <- 0
      
      gene_id_curr <- master_mix[row, "gene_id"]
      LFC_curr <- master_mix[row, "LFC"]
      bp_start_curr <- master_mix[row, "bp_start"]
      bp_end_curr <- master_mix[row, "bp_end"]
      
      gene_id_first_prev <- master_mix[first_r, "gene_id"]
      LFC_first_prev <- master_mix[first_r, "LFC"]
      bp_start_first_prev <- master_mix[first_r, "bp_start"]
      bp_end_first_prev <- master_mix[first_r, "bp_end"]
      
      gene_id_first_next <- master_mix[first_n, "gene_id"]
      LFC_first_next <- master_mix[first_n, "LFC"]
      bp_start_first_next <- master_mix[first_n, "bp_start"]
      bp_end_first_next <- master_mix[first_n, "bp_end"]
      
      gene_id_prev <- master_mix[pos_r, "gene_id"]
      gene_id_next <- master_mix[pos_f, "gene_id"]
      bp_start_next <- master_mix[pos_f,"bp_start"]
      bp_start_prev <- master_mix[pos_r,"bp_start"]
      bp_end_next <- master_mix[pos_f, "bp_end"]
      bp_end_prev <- master_mix[pos_r,"bp_end"]
      LFC_next <- master_mix[pos_f, "LFC"]
      LFC_prev <- master_mix[pos_r, "LFC"]
      bp_start_curr_p <- master_mix[pos_c_b, "bp_start"]
      bp_end_curr_p <- master_mix[pos_c_b, "bp_end"]
      bp_start_curr_n <- master_mix[pos_c_f, "bp_start"]
      bp_end_curr_n <- master_mix[pos_c_f, "bp_end"]
      LFC_curr_p <- master_mix[pos_c_b, "LFC"]
      LFC_curr_n <- master_mix[pos_c_f, "LFC"]
      gene_id_curr_p <- master_mix[pos_c_b, "gene_id"]
      gene_id_curr_n <- master_mix[pos_c_f, "gene_id"]
      cat('\n')
      genes_in_this_cluster <- append(genes_in_this_cluster, list(gene_id_curr))
      print(paste(gene_id_curr, LFC_curr, bp_start_curr, bp_end_curr, "BACKBONE_GENE", "#", clus_count, ":)"), quote=FALSE)
      if (LFC_curr >= 2){
        up_reg <- up_reg+1
      } else if (LFC_curr <= -2){
        down_reg <- down_reg+1
      }
      try(if ((bp_start_curr - bp_end_first_prev)<gap){
        print(paste(gene_id_first_prev, LFC_first_prev, bp_start_first_prev, bp_end_first_prev, "previous_gene #", num_r, (bp_start_curr - bp_end_first_prev)), quote=FALSE)
        genes_in_this_cluster <- append(genes_in_this_cluster, list(gene_id_first_prev))
        if (LFC_first_prev >= 2){
          up_reg <- up_reg+1
        } else if (LFC_first_prev <= -2){
          down_reg <- down_reg+1
        }
        num_r <- num_r+1
        gene_id_curr <- master_mix[row, "gene_id"]
        LFC_curr <- master_mix[row, "LFC"]
        bp_start_curr <- master_mix[row, "bp_start"]
        bp_end_curr <- master_mix[row, "bp_end"]
        gene_id_prev <- master_mix[pos_r, "gene_id"]
        gene_id_next <- master_mix[pos_f, "gene_id"]
        bp_start_next <- master_mix[pos_f,"bp_start"]
        bp_start_prev <- master_mix[pos_r,"bp_start"]
        bp_end_next <- master_mix[pos_f, "bp_end"]
        bp_end_prev <- master_mix[pos_r,"bp_end"]
        LFC_next <- master_mix[pos_f, "LFC"]
        LFC_prev <- master_mix[pos_r, "LFC"]
        bp_start_curr_p <- master_mix[pos_c_b, "bp_start"]
        bp_end_curr_p <- master_mix[pos_c_b, "bp_end"]
        bp_start_curr_n <- master_mix[pos_c_f, "bp_start"]
        bp_end_curr_n <- master_mix[pos_c_f, "bp_end"]
        LFC_curr_p <- master_mix[pos_c_b, "LFC"]
        LFC_curr_n <- master_mix[pos_c_f, "LFC"]
        gene_id_curr_p <- master_mix[pos_c_b, "gene_id"]
        gene_id_curr_n <- master_mix[pos_c_f, "gene_id"]
        
        try(while ((bp_start_curr_p - bp_end_prev)<gap){
          print(paste(gene_id_prev, LFC_prev, bp_start_prev, bp_end_prev, "previous_gene #", num_r, (bp_start_curr_p - bp_end_prev)), quote=FALSE)
          genes_in_this_cluster <- append(genes_in_this_cluster, list(gene_id_prev))
          if (LFC_prev >= 2){
            up_reg <- up_reg+1
          } else if (LFC_prev <= -2){
            down_reg <- down_reg+1
          }
          num_r <- num_r+1
          pos_c_b <- pos_c_b-1
          pos_r <- pos_r-1
          gene_id_curr <- master_mix[row, "gene_id"]
          LFC_curr <- master_mix[row, "LFC"]
          bp_start_curr <- master_mix[row, "bp_start"]
          bp_end_curr <- master_mix[row, "bp_end"]
          gene_id_prev <- master_mix[pos_r, "gene_id"]
          gene_id_next <- master_mix[pos_f, "gene_id"]
          bp_start_next <- master_mix[pos_f,"bp_start"]
          bp_start_prev <- master_mix[pos_r,"bp_start"]
          bp_end_next <- master_mix[pos_f, "bp_end"]
          bp_end_prev <- master_mix[pos_r,"bp_end"]
          LFC_next <- master_mix[pos_f, "LFC"]
          LFC_prev <- master_mix[pos_r, "LFC"]
          bp_start_curr_p <- master_mix[pos_c_b, "bp_start"]
          bp_end_curr_p <- master_mix[pos_c_b, "bp_end"]
          bp_start_curr_n <- master_mix[pos_c_f, "bp_start"]
          bp_end_curr_n <- master_mix[pos_c_f, "bp_end"]
          LFC_curr_p <- master_mix[pos_c_b, "LFC"]
          LFC_curr_n <- master_mix[pos_c_f, "LFC"]
          gene_id_curr_p <- master_mix[pos_c_b, "gene_id"]
          gene_id_curr_n <- master_mix[pos_c_f, "gene_id"]
        })
      })
      try(if ((bp_start_first_next - bp_end_curr)<gap){
        print(paste(gene_id_first_next, LFC_first_next, bp_start_first_next, bp_end_first_next, "next_gene #", num_f, (bp_start_first_next - bp_end_curr)), quote=FALSE)
        genes_in_this_cluster <- append(genes_in_this_cluster, list(gene_id_first_next))
        if (LFC_first_next >= 2){
          up_reg <- up_reg+1
        } else if (LFC_first_next <= -2){
          down_reg <- down_reg+1
        }
        num_f <- num_f+1
        gene_id_curr <- master_mix[row, "gene_id"]
        LFC_curr <- master_mix[row, "LFC"]
        bp_start_curr <- master_mix[row, "bp_start"]
        bp_end_curr <- master_mix[row, "bp_end"]
        gene_id_prev <- master_mix[pos_r, "gene_id"]
        gene_id_next <- master_mix[pos_f, "gene_id"]
        bp_start_next <- master_mix[pos_f,"bp_start"]
        bp_start_prev <- master_mix[pos_r,"bp_start"]
        bp_end_next <- master_mix[pos_f, "bp_end"]
        bp_end_prev <- master_mix[pos_r,"bp_end"]
        LFC_next <- master_mix[pos_f, "LFC"]
        LFC_prev <- master_mix[pos_r, "LFC"]
        bp_start_curr_p <- master_mix[pos_c_b, "bp_start"]
        bp_end_curr_p <- master_mix[pos_c_b, "bp_end"]
        bp_start_curr_n <- master_mix[pos_c_f, "bp_start"]
        bp_end_curr_n <- master_mix[pos_c_f, "bp_end"]
        LFC_curr_p <- master_mix[pos_c_b, "LFC"]
        LFC_curr_n <- master_mix[pos_c_f, "LFC"]
        gene_id_curr_p <- master_mix[pos_c_b, "gene_id"]
        gene_id_curr_n <- master_mix[pos_c_f, "gene_id"]
        try(while ((bp_start_next - bp_end_curr_n)<gap){
          print(paste(gene_id_next, LFC_next, bp_start_next, bp_end_next, "next_gene #", num_f, (bp_start_next - bp_end_curr_n)), quote=FALSE)
          genes_in_this_cluster <- append(genes_in_this_cluster, list(gene_id_next))
          if (LFC_next >= 2){
            up_reg <- up_reg+1
          } else if (LFC_next <= -2){
            down_reg <- down_reg+1
          }
          num_f <- num_f+1
          pos_c_f <- pos_c_f+1
          pos_f <- pos_f+1
          gene_id_curr <- master_mix[row, "gene_id"]
          LFC_curr <- master_mix[row, "LFC"]
          bp_start_curr <- master_mix[row, "bp_start"]
          bp_end_curr <- master_mix[row, "bp_end"]
          gene_id_prev <- master_mix[pos_r, "gene_id"]
          gene_id_next <- master_mix[pos_f, "gene_id"]
          bp_start_next <- master_mix[pos_f,"bp_start"]
          bp_start_prev <- master_mix[pos_r,"bp_start"]
          bp_end_next <- master_mix[pos_f, "bp_end"]
          bp_end_prev <- master_mix[pos_r,"bp_end"]
          LFC_next <- master_mix[pos_f, "LFC"]
          LFC_prev <- master_mix[pos_r, "LFC"]
          bp_start_curr_p <- master_mix[pos_c_b, "bp_start"]
          bp_end_curr_p <- master_mix[pos_c_b, "bp_end"]
          bp_start_curr_n <- master_mix[pos_c_f, "bp_start"]
          bp_end_curr_n <- master_mix[pos_c_f, "bp_end"]
          LFC_curr_p <- master_mix[pos_c_b, "LFC"]
          LFC_curr_n <- master_mix[pos_c_f, "LFC"]
          gene_id_curr_p <- master_mix[pos_c_b, "gene_id"]
          gene_id_curr_n <- master_mix[pos_c_f, "gene_id"]
          
        })
      })
      ## Below prints out the information on the genes in the cluster: gene_ids, cluster_id, number of genes in cluster, number upregulated, number down regulated. "%" is used
      ## later in bash to help separate out elements
      print(paste("%", genes_in_this_cluster, clus_count, num_f+num_r-1, up_reg, down_reg,"%"), quote=FALSE) 
      
      ## Counter that keeps track of the total number of genes in the current cluster
      number_of_genes_in_this_cluster <- (num_f+num_r-1)
      list_of_gene_count_all_clusters <- append(list_of_gene_count_all_clusters, list(number_of_genes_in_this_cluster))
    }
  }
}
