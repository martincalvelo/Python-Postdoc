import argparse
import pandas as pd
import MDAnalysis as mda
import numpy as np


#**************************************************************************************************************************************************************************
#Arguments:
desc="""
Analysis of the trajectories for the ABlina Database (Xunta de Galicia Project 2023).
Results stored in a csv file for later Data Analysis (ABlinaresults.csv)
python analysis_ABlina_database.py -pdb ref.pdb -xtc traj.xtc
March 2024, by Martin Calvelo martin.calvelo@gmail.com
"""

parser = argparse.ArgumentParser(description = desc )
parser.add_argument( "-pdb", "--pdb", type = str, default = 'check.pdb',
    help = """PDB File from the simulation.\n
    Default: %(default)s """ )
parser.add_argument( "-xtc", "--xtc", type = str, default = 'check.xtc',
    help = """XTC File with the trajectory of the system.\n
    Default: %(default)s """ )
args = parser.parse_args()

#**************************************************************************************************************************************************************************
#Functions:

def clusters( frame ):
    ''' Calculate number of peptide clusters (also most repeated and biggest but not important so far). Mostly based on mdanalysis library
    '''
    global Results
    Groups = []
    #Iterate by beta amyloid. For each of them, check which ones of the others are close to it.
    #Stored as a list of tuples. 
    for i,pep in enumerate(peptides.segments):
        sel1="protein and segid Pep"+str(i+1); sel2="protein and not segid Pep"+str(i+1)
        pep1=u.select_atoms(sel1); pep2=u.select_atoms(sel2)
        if len(mda.lib.NeighborSearch.AtomNeighborSearch(pep2, frame.dimensions).search(pep1, 6, level='S')) > 0:
            Groups.append(tuple(set(mda.lib.NeighborSearch.AtomNeighborSearch(pep2, frame.dimensions).search(pep1, 6, level='S').segids)))
        else:
            Groups.append(tuple(pep1.segments.segids))
    Cluster = []
    #Using intersection, we clusterize (check the elements in common in the tuples)
    for g1 in Groups:
        for g2 in Groups*2:
            if set( g1 ).intersection( set( g2 ) ):
                g1 = tuple( set( tuple( g1 ) + tuple( g2 ) ) )
        Cluster.append( tuple( sorted( g1 ) ) )
    Cluster = tuple( set( Cluster ) )
    #Order results for having most repeated and number of cluster
    list_clusters =([ np.sum([ 1 for j in range(len(Cluster)) if len(Cluster[j])==i]) for i in range(1,len(peptides.segments)+1) ])
    if list_clusters.count(max(list_clusters)) > 1:
        most_repeated=[ ]
        for i,number in enumerate( list_clusters ):
            if number == max(list_clusters):
                most_repeated.append(i+1)
    else:
        most_repeated=list_clusters.index(max(list_clusters))+1
    #With this, we have the cluster with the highest number of elements
    last_cluster=next(( i for i in reversed(list_clusters) if i != 0 ), None)
    position_last_cluster=len(list_clusters) - list_clusters[::-1].index(last_cluster)
    Results[ "Num_clusters" ] = len(Cluster)
    return

def aggregation( frame ):
    global Results
    ''' Calculate percentage of peptide aggregation. Follows same stragy than for clusterization
    '''
    aggregated=0
    for i,pep in enumerate(peptides.segments):
        sel1="protein and segid Pep"+str(i+1); sel2="protein and not segid Pep"+str(i+1)
        pep1=u.select_atoms(sel1); pep2=u.select_atoms(sel2)
        if len(mda.lib.NeighborSearch.AtomNeighborSearch(pep2, frame.dimensions).search(pep1, 6, level='S')) > 0:
            aggregated+=1
    #Beta amyloids that are INSIDE a cluster, so they are not isolated
    aggregation=round(((aggregated)/(len(peptides.segments)))*100, 3)
    Results[ "Aggregation" ] = aggregation


def depostion( frame ):
    global Results
    '''Calculate percentage of peptides deposited on the membrane (deposition). Analogous to aggregation
    '''
    contacts = []
    lipids=u.select_atoms("resname POPC POPE POPS CHOL DPSM GM1 GM3 GB3") #Type of lipids simulated so far
    in_membrane=len(mda.lib.NeighborSearch.AtomNeighborSearch(peptides, frame.dimensions).search(lipids, 6, level='S').segids)
    aggregation=round(((in_membrane)/(len(peptides.segments)))*100, 3)
    #Number of Beta amiloyds deposited in the lipid membrane as percentage and raw number
    Results[ "Deposition" ] = aggregation
    Results[ "AB_in_membrane" ] = in_membrane
    return

def calculate_distance(point1, point2):
    return np.sqrt(np.sum((point1 - point2)**2))

def distances_mem( frame ):
    global Results
    '''Calculate the average and minimum distance between the membrane and peptides
    '''
    Bilayer = u.select_atoms('resname POPC POPE POPS CHOL DPSM GM1 GM3 GB3').center_of_mass()
    distance = []
    for i,pep in enumerate(peptides.segments):
        sel1="protein and segid Pep"+str(i+1); pep1=u.select_atoms(sel1).center_of_mass()
        distance.append(calculate_distance(Bilayer, pep1))
    Results[ "Av_distance" ] = round(np.mean(distance), 3)
    Results[ "Min_distance" ] = round(min(distance), 3)
    return

def contacts_mol_lipid( frame ):
    global Results
    '''Number of X lipid molecules that are in contact with at least one peptide
    '''
    list_lipids=['POPC','POPE','POPS','CHOL','DPSM','GM1','GM3','GB3']
    total=0
    for i in list_lipids:
        name_key="Contacts_m"+str(i)
        sel1="resname "+str(i); lipid=u.select_atoms(sel1)
        if i == 'GM3':
            contacts=len(mda.lib.NeighborSearch.AtomNeighborSearch(lipid, frame.dimensions).search(peptides, 6, level='S'))
        else:
            contacts=len(mda.lib.NeighborSearch.AtomNeighborSearch(lipid, frame.dimensions).search(peptides, 6, level='R'))
        Results[ name_key ] = contacts
        total+=contacts
    Results[ "Contacts_mtotal" ] = total
    return

def contacts_atom_lipid( frame ):
    global Results
    '''Number of X lipid atoms (CG level) that are in contact with at least one peptide. Analogous to contacts_mol_lipid, but in this case we identify the exact atom that are in contact.
    '''
    list_lipids=['POPC','POPE','POPS','CHOL','DPSM','GM1','GM3','GB3']
    total=0
    for i in list_lipids:
        name_key="Contacts_a"+str(i)
        sel1="resname "+str(i); lipid=u.select_atoms(sel1)
        contacts=len(mda.lib.NeighborSearch.AtomNeighborSearch(lipid, frame.dimensions).search(peptides, 6, level='A'))
        Results[ name_key ] = contacts
        total+=contacts
    Results[ "Contacts_atotal" ] = total
    return

#**************************************************************************************************************************************************************************

if __name__ == "__main__" :
    # Load the structure and the trajectory
    u = mda.Universe( args.pdb, args.xtc )

   #Select bilayer, all possible lipids simulated so far (Jan 2024)
    Bilayer = u.select_atoms('resname POPC POPE POPS CHOL DPSM GM1 GM3 GB3')

    #Select AB42 peptides
    peptides = u.select_atoms('protein')
    number=0
    #Next loop is for assing a segid number (pdb format) different to each beta amiloyd. Very useful for the analysis, simplyfies a lot 
    for i, atom in enumerate( peptides ):
        if atom.name == 'BB' and atom.resid == 1:
            number+=1
            interval=atom.index+96; selection="index "+str(atom.index)+" to "+str(atom.index+96); name_segment="Pep"+str(number)
            core_atoms = u.select_atoms(selection); core_segment = u.add_Segment(segid=name_segment)
            core_atoms.residues.segments = core_segment

    # Dictionary for storing the results
    Results = { "Time" : {},
                "Contacts_mPOPC" : {},
                "Contacts_mPOPE" : {},
                "Contacts_mPOPS" : {},
                "Contacts_mCHOL" : {},
                "Contacts_mDPSM" : {},
                "Contacts_mGM3" : {},
                "Contacts_mGM1" : {},
                "Contacts_mGB3" : {},
                "Contacts_mtotal" : {},
                "Contacts_aPOPC" : {},
                "Contacts_aPOPE" : {},
                "Contacts_aPOPS" : {},
                "Contacts_aCHOL" : {},
                "Contacts_aDPSM" : {},
                "Contacts_aGM3" : {},
                "Contacts_aGM1" : {},
                "Contacts_aGB3" : {},
                "Contacts_atotal" : {},
                "Aggregation" : {},
                "Num_clusters" : {},
                "Av_distance" : {},
                "Min_distance" : {},
                "AB_in_membrane" : {},
                "Deposition" : {} }

    # Dataframe for results
    df_results = pd.DataFrame()
    #Process trajectory
    for id_f, frame in enumerate(u.trajectory[:]):
        Results[ "Time" ] = [ frame.time/1000 ]
        #Compute
        contacts_mol_lipid(frame)
        contacts_atom_lipid(frame)
        clusters(frame)
        aggregation(frame)
        depostion(frame)
        distances_mem(frame)
        to_df = pd.DataFrame(Results)
        df_results = pd.concat([df_results, to_df], ignore_index=True)

    #Print summary, also check if the scripts works
    print("Aggregation: "+str(df_results['Aggregation'].mean()))
    print("Deposition: "+str(df_results['Deposition'].mean()))

    #Create a logfile with a summary of the results. So far, the most important is the analysis of the peptide aggregation and deposition on the membrane
    with open('analysis.log', 'w') as fout:
        fout.write("Aggregation: "+str(df_results['Aggregation'].mean())+"\n")
        fout.write("Deposition: "+str(df_results['Deposition'].mean())+"\n")

    #Print the csv file
    df_results.to_csv('ABlinaresults.csv', index=False)

