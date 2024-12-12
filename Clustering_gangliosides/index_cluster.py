import argparse
import MDAnalysis as mda

#**************************************************************************************************************************************************************************
#Arguments:
parser = argparse.ArgumentParser(description =
    ''' Index for clustering''')
parser.add_argument( "-pdb", "--pdb", type = str, default = '1ms.pdb',
    help = """PDB File from the simulation.\n
    Default: %(default)s """ )
args = parser.parse_args()

#**************************************************************************************************************************************************************************
#Functions:
def get_centers(molecule, frame):
    '''Calculate COM of every GM1 '''
    global Results
    for mol in molecule.segments:
        name_key=str(mol.segid)
        sel1="resname GM1 and segid "+str(mol.segid); mol1=u.select_atoms(sel1) #If only sugars remove CER1
        center=mol1.center_of_mass(); center_output=list(center[:2])
        Results[ name_key ] = center_output
    return 

#**************************************************************************************************************************************************************************

if __name__ == "__main__" :
    # Load the structure and the trajectory
    u = mda.Universe( args.pdb )

    # Select GM1
    gm1=u.select_atoms("resname GM1")

    for mol in gm1.segments:
        sel1="resname GM1 and segid "+str(mol.segid); mol1=u.select_atoms(sel1)
        with mda.selections.gromacs.SelectionWriter('index_clusters.ndx', mode='a') as ndx:
            ndx.write(mol1, name=str(mol.segid))

