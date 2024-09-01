import argparse
import pandas as pd
import MDAnalysis as mda
import numpy as np

#**************************************************************************************************************************************************************************
#Arguments:
desc="""
Analysis of the trajectories for the ABlina Database (Xunta de Galicia Project 2023).
Contacs between gangliosides and beta amyloids, split by ganglioside CG atom
python contacts_GM.py -pdb ref.pdb -xtc traj.xtc -lip GM1
March 2024, by Martin Calvelo martin.calvelo@gmail.com
"""

parser = argparse.ArgumentParser(description = desc )
parser.add_argument( "-pdb", "--pdb", type = str, default = 'check.pdb',
    help = """PDB File from the simulation.\n
    Default: %(default)s """ )
parser.add_argument( "-xtc", "--xtc", type = str, default = 'check.xtc',
    help = """XTC File with the trajectory of the system.\n
    Default: %(default)s """ )
parser.add_argument( "-lip", "--lipid", type = str, default = 'GM1',
    help = """Ganglioside of interest (GM1 or GM3 so far).\n
    Default: %(default)s """)
args = parser.parse_args()

#**************************************************************************************************************************************************************************
#Functions:

def contacts_atom_ganglioside( frame, lipid ):
    global Results
    '''Number of contacts beteen ganglioside CG atoms and peptides
    '''
    #Dictionary depending on lipid type
    if lipid == "GM1":
        GM = { "glc" : "name A2 B2 C2 V2",
               "gal-int" : "name A3 B3 C3 V3",
               "nana" : "name A4 B4 C4 D4 E4 V4",
               "galnac" : "name A5 B5 C5 D5 V5",
               "gal-ext" : "name A6 B6 C6 V6",
               "tails" : "name AM1 AM2 T1A C2A C3A C1B C2B C3B C4B",
             }
    elif lipid == "GM3":
        GM = { "glc" : "resid 2",
                "gal-int" : "resid 3",
                "nana" : "resid 4",
                "galnac" : "name A5 B5 C5 D5 V5",
                "gal-ext" : "name A6 B6 C6 V6",
                "tails" : "resid 1",
              }
    
    list_atoms=list(GM.keys())

    total=0
    for i in list_atoms:
        name_key="Contacts_"+str(i)
        sel1="resname "+str(lipid)+" and "+str(GM[i]); sel_atoms=u.select_atoms(sel1)
        contacts=len(mda.lib.NeighborSearch.AtomNeighborSearch(sel_atoms, frame.dimensions).search(peptides, 6, level='A'))
        Results[ name_key ] = contacts
        total+=contacts
    Results[ "Contacts_total" ] = total
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
    for i, atom in enumerate( peptides ):
        if atom.name == 'BB' and atom.resid == 1:
            number+=1
            interval=atom.index+96; selection="index "+str(atom.index)+" to "+str(atom.index+96); name_segment="Pep"+str(number)
            core_atoms = u.select_atoms(selection); core_segment = u.add_Segment(segid=name_segment)
            core_atoms.residues.segments = core_segment

    # Dictionary for the results
    Results = { "Time" : {},
                "Contacts_glc" : {},
                "Contacts_gal-int" : {},
                "Contacts_nana" : {},
                "Contacts_galnac" : {},
                "Contacts_gal-ext" : {},
                "Contacts_tails" : {} ,
                "Contacts_total" : {} }

    # Dataframe for results
    df_results = pd.DataFrame()
    #Process trajectory
    for id_f, frame in enumerate(u.trajectory[:]):
        Results[ "Time" ] = [ frame.time/1000 ] 
        #Compute
        contacts_atom_ganglioside(frame, args.lipid)
        print(Results)
        to_df = pd.DataFrame(Results)
        df_results = pd.concat([df_results, to_df], ignore_index=True)

    #Store results in a csv
    df_results.to_csv('atoms-ganglioside-contacts.csv', index=False)


