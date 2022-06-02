for i in {0..10000000..30000};
do
    echo $i >> cluster_starts.txt
done

import os 

path = '/home/jack.a/inspector_clusto/inspector_package/above_4/separated_clusters'

folder = os.fsencode(path)
for file in os.listdir(folder):
    with open(file, "r") as a_file:
        columns = line.split()
        cluster_begin = float(columns[2])
        gene_id = columns[1]
        for clusterint in range(0,10000000, 30000):
            if cluster_begin == clusterint:
                print(file, gene_id, clusterint)


for file in X: 

for i in (0,10000000, 30000):
do
    grep i file i > i.txt
    wc -l i.txt
done

done
