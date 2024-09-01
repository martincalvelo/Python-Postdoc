#!/usr/bin/env python

import optparse
import fileinput

desc="""
Create topology for cyclic peptides (CPs) at Coarse Grained (CG) level using Martini 2 Force Field, following same strategy than https://pubs.acs.org/doi/full/10.1021/jp9064196
Starts from CG topolgy of linear CP, obtained using martinize.py (Available at Martini webpage)
python create_topology_CPs.py -f Protein.itp -m MOL -o CP.itp
January 2017, by Martin Calvelo martin.calvelo@gmail.com
"""
def main() :
    parser = optparse.OptionParser(description=desc, version='%prog version 1.0')
    parser.add_option('-f', '--infile', help='input top file1)', action='store')
    parser.add_option('-o', '--outfile', help='output itp file', action='store')
    parser.add_option('-m', '--molname', help='name MOL', action='store')
    parser.set_defaults(infile='topol.top', outfile='CP', molname='MOL')
    options, arguments = parser.parse_args()
#**************************************************************************************************

#Read inputs 
    molname=str(options.molname)
    infile=open(str(options.infile), "r")
    fout=open(str(options.outfile), "w")

#1st read of topology. Read index number
    read_atoms=0
    bb=[]
    for line in infile :
        if line.strip() :
            fields=line.split()
            if (len(fields)>1):
                if fields[1]=="atoms":
                    read_atoms=1
                elif fields[1]=="bonds":
            	    read_atoms=0
                if read_atoms == 1 and len(fields) > 3:
                    if fields[4]=="BB":
                        bb.append([fields[0], fields[2]])
    infile.close()

#2nd read of topology. Read everything, split by topology (GROMACS format)
    atoms=[]; bonds=[]; constraints=[]; angles=[]; dihedrals=[]
    #If this variable is 0, don't read. If is 1, read. One for each section of gromacs topology
    read_at=0; read_bonds=0; read_constraints=0; read_angles=0; read_dih=0
    infile=open(infile_name, "r")
    for line in infile :
        if line.strip() :
            fields=line.split()
            if (len(fields)>1):
                #When the name of each section is found, variable turns 1. With this we will read the info in the next loop
                if fields[0]=="[" and fields[1]=="atoms" and fields[2]=="]":
                    read_at=1
                elif fields[0]==";" and fields[1]=="Sidechain" and fields[2]=="bonds":
                    read_bonds=1
                elif fields[1]=="constraints":
                    read_constraints=1; read_bonds=0
                elif fields[0]==";" and fields[1]=="Backbone-sidechain" and fields[2]=="angles":
                    read_angles=1
                elif  fields[0]=="[" and fields[1]=="dihedrals" and fields[2]=="]":
                    read_angles = 0; read_dih=1
                elif fields[0]=="[":
                    read_at=0; read_constraints=0
                #Read each section if variable is 1
                if read_at > 0 and fields[0] != "[":
                    atoms.append(fields)
                elif read_bonds > 0 and fields[0] != "[":
                    bonds.append(fields)
                elif read_constraints > 0 and fields [0] != "[":
                    constraints.append(fields)
                elif read_angles > 0 and fields[0] != "[":
                    angles.append(fields)
                elif read_dih > 0 and fields[0] != "[":
                    dihedrals.append(fields)
    infile.close()

#Some proteins have exclusions or virtual sites, but not all of them. Read again, just in case
    infile=open(infile_name, "r")
    exclusions=[]
    read_exclusions=0
    for line in infile :
        if line.strip() :
            fields=line.split()
            if (len(fields)>1):
                if fields[0]=="[" and fields[1]=="exclusions" and fields[2]=="]":
                    read_exclusions=1
                elif fields[0]=="[":
                    read_exclusions=0
                if read_exclusions==1 and fields[0] != "[": 
                    exclusions.append(fields)
    infile.close()
    infile=open(infile_name, "r")
    virtual_sites=[]
    read_virtual=0
    for line in infile :
        if line.strip() :
            fields=line.split()
            if (len(fields)>1):
                if fields[0]=="[" and fields[1]=="virtual_sites2" and fields[2]=="]":
                    read_virtual=1
                elif fields[0]=="[":
                    read_virtual=0
                if read_virtual==1 and fields[0] != "[":
                    virtual_sites.append(fields)
    infile.close()
    
# Create list with bonds between CG-atoms of backbone, adding an extra one for closing CP. Distances and Force constants in Tarek's paper
    bonds_BB=[]
    for i in range(len(bb)-1):
    	bonds_BB.append([bb[i][0], bb[i+1][0], "1", "0.38000", "6275", "; Not in martinize"])
    bonds_BB.append([bb[-1][0], bb[0][0], "1", "0.38000", "6275", "; Not in martinize"]) #Close CP

# Now, angles
    angles_BB=[]
    for i in range(len(bb)-2):
        angles_BB.append([bb[i][0], bb[i+1][0], bb[i+2][0], "1", "135", "627", "; Not in martinize"])
    angles_BB.append([bb[-2][0], bb[-1][0], bb[0][0], "1", "135", "627", "; Not in martinize"]) #Close CP
    angles_BB.append([bb[-1][0], bb[0][0], bb[1][0], "1", "135", "627", "; Not in martinize"]) #Close CP
    angles_BB.append([bb[-1][0], bb[0][0], str(int(bb[0][0])+1), "2", "100", "25", "; Not in martinize"]) #Close CP

# Dihedrals
    dihedrals_BB=[]
    for i in range(len(bb)-3):
        dihedrals_BB.append([bb[i][0], bb[i+1][0], bb[i+2][0], bb[i+3][0], "2", "180", "418", "; Not in martinize"])
    dihedrals_BB.append([bb[-3][0], bb[-2][0], bb[-1][0], bb[0][0], "2", "180", "418", "; Not in martinize"]) #Close CP
    dihedrals_BB.append([bb[-2][0], bb[-1][0], bb[0][0], bb[1][0], "2", "180", "418", "; Not in martinize"]) #Close CP
    dihedrals_BB.append([bb[-1][0], bb[0][0], bb[1][0], bb[2][0], "2", "180", "418", "; Not in martinize"]) #Close CP

# Write output topology
    label=";Topology generated with python script starting from "+str(infile_name)+" for the molecule "+str(molname)+"\n"
    fout.write(label)
    fout.write("[ moleculetype ]\n")
    fout.write(";name nrexcl\n")
    fout.write("%s   3\n" % (molname))
    fout.write('\n')
# Write atoms
    fout.write("[ atoms ]\n")
    for i in range(len(atoms)):
        if atoms[i][1] == "residue":
            fout.write("%s   \n" % " ".join(atoms[i]))
        elif atoms[i][0] == ";":
            fout.write("; %s   \n" % "     ".join(atoms[i][1:]))
        else:
            fout.write("  %s   \n" % "\t".join(atoms[i]))
    fout.write('\n')
# Write exclusions & virtual_sites, if exists
    if len(exclusions) > 0:
        for i in range(len(exclusions)):
            fout.write("  %s   \n" % "\t".join(exclusions[i]))
    fout.write('\n')
    if len(virtual_sites) > 0:
        fout.write("[ virtual_sites2 ]\n")
        for i in range(len(virtual_sites)):
            fout.write("  %s   \n" % "\t".join(virtual_sites[i]))
    fout.write('\n')

# Write bonds, adding those need for closing CP
    fout.write("[ bonds ]\n")
    fout.write("; Backbone bonds \n")
    for i in range(len(bonds_BB)):
        fout.write("  %s   \n" % "\t".join(bonds_BB[i]))
    for i in range(len(bonds)):
        if bonds[i][0] == ";":
            fout.write("; %s   \n" % "     ".join(bonds[i][1:]))
        else:
            fout.write("  %s   \n" % "\t".join(bonds[i]))
    fout.write('\n')

# Write constraints, if exists
    if len(constraints) > 1:
        fout.write("[ constraints ]\n")
        for i in range(len(constraints)):
            if constraints[i][0] == ";":
                fout.write("; %s   \n" % "     ".join(constraints[i][1:]))
            else:
                fout.write("  %s   \n" % "\t".join(constraints[i]))
    fout.write('\n')

# Write angles, adding those need for closing CP
    fout.write("[ angles ]\n")
    fout.write("; Backbone angles \n")
    for i in range(len(angles_BB)-1):
        fout.write("  %s   \n" % "\t".join(angles_BB[i]))
    fout.write("; Backbone-sidechain angles\n")
    fout.write("  %s   \n" % "\t".join(angles_BB[-1]))
    for i in range(2,len(angles)):
        if angles[i][0] == ";":
            fout.write("; %s   \n" % "     ".join(angles[i][1:]))
        else:
            fout.write("  %s   \n" % "\t".join(angles[i]))
    fout.write('\n')

# Write dihedrals, adding those need for closing CP
    fout.write("[ dihedrals ]\n")
    fout.write("; Backbone dihedrals \n")
    for i in range(len(dihedrals_BB)):
        fout.write("  %s   \n" % "\t".join(dihedrals_BB[i]))
    fout.write("; Sidechain improper dihedrals\n")
    for i in range(len(dihedrals)):
        if len(dihedrals[i]) > 5:
            if dihedrals[i][0] == ";":
                fout.write("; %s   \n" % "     ".join(dihedrals[i][1:]))
            else:
                fout.write("  %s   \n" % "\t".join(dihedrals[i]))
    fout.write('\n')

if __name__=="__main__" :
        main()

