import argparse
import pandas as pd
import MDAnalysis as mda
import numpy as np

#**************************************************************************************************************************************************************************
#Arguments:
desc="""
Create csv with top info of a trajectory included in ABlina database. Useful for data analysis of ABlina.
Results stored in a csv file for later Data Analysis (info_top.csv)
python get_top_info.py -pep 1ba4 -mem PC -pdb ref.pdb
March 2024, by Martin Calvelo martin.calvelo@gmail.com
"""

parser = argparse.ArgumentParser(description = desc )
parser = argparse.ArgumentParser(description =
    '''Create csv with top info of a trajectory included in ABlina database. Useful for data analysis of ABlina''')
parser.add_argument( "-pdb", "--pdb", type = str, default = 'ref.pdb',
    help = """PDB File from the simulation.\n
    Default: %(default)s """ )
parser.add_argument( "-mem", "--mem", type = str, default = 'PC',
    help = """Membrane name.\n
    Default: %(default)s """ )
parser.add_argument( "-pep", "--pep", type = str, default = '1ba4',
    help = """PDB code of the peptide. In ablina different pdb codes, choose one \n
    Default: %(default)s """ )
args = parser.parse_args()

#**************************************************************************************************************************************************************************
#Functions
def get_lipids():
    global Results
    total_lipids=0
    list_lipids=['POPC','POPE','POPS','CHOL','DPSM','GM1','GM3','GB3']
    for lip in list_lipids:
        name_key="#_"+str(lip)
        sel1="resname "+str(lip); lipid=u.select_atoms(sel1)
        if lip == 'GM3': #Needed because "mistake" in GM3 topology. Each GM3 molecule is split in 4 residues. 
            Results[ name_key ] = int(len(lipid.residues)/4)
            total_lipids=total_lipids+int(len(lipid.residues)/4)
        else:
            Results[ name_key ] = len(lipid.residues)
            total_lipids=total_lipids+int(len(lipid.residues))
    Results[ "#_total_lipids" ] = total_lipids


#**************************************************************************************************************************************************************************

if __name__ == "__main__" :

    # Dictionary for the results
    Results = { "membrane" : {},
                "PDB ID" : {},
                "#_peptides" : {},
                "#_POPC" : {},
                "#_POPE" : {},
                "#_POPS" : {},
                "#_CHOL" : {},
                "#_DPSM" : {},
                "#_GM3" : {},
                "#_GM1" : {},
                "#_GB3" : {},
                "#_total_lipids" : {} }
    Results[ "membrane" ] = ( args.mem )
    Results[ "PDB ID" ] = ( args.pep ).upper()

    # Load the structure and the trajectory
    u = mda.Universe( args.pdb )

    #Count lipids
    get_lipids()

    #Count AB42 peptides
    peptides = u.select_atoms('protein')
    Results[ "#_peptides" ] = int(len(peptides)/96)

    #Store results in a csv file
    df_results = pd.DataFrame(Results, index=[0])
    df_results.to_csv('info_top.csv', index=False)

