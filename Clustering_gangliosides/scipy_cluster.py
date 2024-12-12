import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster, linkage
from collections import Counter

# v0.1 Martin Calvelo (ITQB NOVA) martin.calvelo(at)gmail.com
# Script for calculating cluster size of gangliosides in membrane. First, with gromacs, the COM of each ganglioside molecule was calculated:
# gmx traj -s topol.tpr -f traj.xtc -com -n index_clusters.ndx -ox A.xvg
# There is an index group per ganglioside molecule, created with index_cluster.py

def read_xvg( xvg ):
    with open(xvg, 'r') as f:
        lines = f.readlines()
    data = [line.split() for line in lines if not line.startswith(('#', '@'))]
    data = np.array(data, dtype=float)
    x = data[:, 1]  
    y = data[:, 2]
    
    return np.column_stack((x, y))

def clusters( coords_frame, cut_off ):
    global Results, num_mol
    # Calcute distance_matrix
    dists = pdist(coords_frame, metric='euclidean')
    # group
    links = linkage(dists, method='single')
    # Determinate to which cluster belong
    label = fcluster(links, t=cut_off, criterion='distance')
    # Count size clusters
    count_size = Counter(np.bincount(label)); keys=list(count_size.keys())
    #Output in list
    list_clusters = []
    for i in range(1,num_mol+1):
        if i in keys:
            list_clusters.append(count_size[i])
        else:
            list_clusters.append(0)
    #Ensure sum equal to number of GM1
    check=0
    for c in range(len(list_clusters)):
        check=check+(list_clusters[c]*(c+1))
    if check != 24:
        print(check, frame, "ERROR!!!")

    #Store results in dataframe
    for i in range(len(list_clusters)):
        name_key=str(i+1)
        Results [ name_key ] = list_clusters[i]

files=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X"]
Results = { "Frame": {}, "1": {}, "2": {}, "3": {}, "4": {}, "5": {}, "6": {}, "7": {}, "8": {}, "9": {}, "10": {}, "11": {}, "12": {}, "13": {},
            "14": {}, "15": {}, "16": {}, "17": {}, "18": {}, "19": {}, "20": {}, "21": {}, "22": {}, "23": {}, "24": {} }
coords = []
for xvgfile in files:
    name_file=str(xvgfile)+".xvg"
    coord = read_xvg(name_file)
    coords.append(coord)
coords = np.array(coords)
num_mol=coords.shape[0]
frames=coords.shape[1]

# Dataframe for results
df_results = pd.DataFrame()

for frame in range(frames):
    Results[ "Frame" ] = frame
    coords_frame = coords[:, frame, :]
    clusters( coords_frame, 1.5 )
    to_df = pd.DataFrame(Results, index=[0])
    df_results = pd.concat([df_results, to_df], ignore_index=True)
df_results.to_csv('Cluster_vs_time.csv', index=False)

