#!/usr/bin/env python

import optparse
import fileinput

desc='''Create Gaussian input for TS optimization from xyz and gaussian reference (freezing distances/angles/dihedrals)
The functional, charge, multiplicity, pseudopotentials, and bases are defined in this script. Modify if you want to use other options.
create_TS_optimization_gaussian.py -c coords.xyz -r reference.com -o CP.itp 
June 2016, by Martin Calvelo martin.calvelo@gmail.com
'''
def main() :
    parser = optparse.OptionParser(description=desc, version='%prog version 1.0')
    parser.add_option('-c', '--coordinates', help='input xyz optimizationfixed.xyz', action='store')
    parser.add_option('-r', '--reference', help='renference optF.com', action='store')
    parser.add_option('-o', '--outfile', help='output file opt.com', action='store')
    parser.set_defaults(coordinates='fixed.xyz', outfile='opt.com', reference='optF.com')
    options, arguments = parser.parse_args()
#**************************************************************************************************************************************
    infile=open(str(options.coordinates), "r")
    ref_infile=open(str(options.reference), "r")
    fout=open(str(options.outfile), "w")

#Read atoms from xyz
    atoms=[]
    for line in infile :
        if line.strip() :
            fields=line.split()
            if (len(fields)>1) and fields[0]!= "scf":
                atoms.append([str(fields[0]), str(fields[1]), str(fields[2]), str(fields[3])])
   
#Gaussian options. Change here functional (b3lyp), charge and multiplicity if applies.   
#Write it in the output
        fout.write("""# opt=(ts,calcfc,noeigentest,modredundant) freq=noraman gen b3lyp pseudo=read

____

1 1     
""")

#Write atoms in the output after gaussian options.
    for i in range(len(atoms)):
        fout.write("%s\t%s\t%s\t%s  \n" % (str(atoms[i][0]), str(atoms[i][1]), str(atoms[i][2]), str(atoms[i][3]) )) 
    fout.write("\n")  #This blank line is needed for Gaussian format.
    
#Read fixed bonds/angles/dihedrals from gaussian input of reference. They are marked with an F. In the new input, the will be marked with a B (Gaussian will monitorize them in this way)
    ref_fixed=[]
    for line in ref_infile :
         if line.strip() :
            fields=line.split()
            if fields[0] == "B" and fields[3] == "F":
                    ref_fixed.append([str(fields[0]), str(fields[1]), str(fields[2]), "B"])
            elif fields[0] == "A" and fields[4] == "F":
                    ref_fixed.append([str(fields[0]), str(fields[1]), str(fields[2]), str(fields[3]), "B"])
            elif fields[0] == "D" and fields[5] == "F":
                    ref_fixed.append([str(fields[0]), str(fields[1]), str(fields[2]), str(fields[3]), str(fields[4]), "B"])
#Write it in the output
    for i in range(len(ref_fixed)):
            if ref_fixed[i][0] == "B":
                    fout.write("%s\t%s\t%s\t%s  \n" % (str(ref_fixed[i][0]), str(ref_fixed[i][1]), str(ref_fixed[i][2]), str(ref_fixed[i][3]) ))
            if ref_fixed[i][0] == "A":
                    fout.write("%s\t%s\t%s\t%s\t%s  \n" % (str(ref_fixed[i][0]), str(ref_fixed[i][1]), str(ref_fixed[i][2]), str(ref_fixed[i][3]), str(ref_fixed[i][4]) ))
            if ref_fixed[i][0] == "D":
                    fout.write("%s\t%s\t%s\t%s\t%s\t%s  \n" % (str(ref_fixed[i][0]), str(ref_fixed[i][1]), str(ref_fixed[i][2]), str(ref_fixed[i][3]), str(ref_fixed[i][4]),  str(ref_fixed[i][5]) ))

#Write the pseudopotentials and the base
#One for each atom. If you have other types of atoms, add them here. If any of the atoms are not in your system, delete it (If not, Gaussian will give you and error
    fout.write("""
Ir 0
LANL2DZ
 ****
C 0
6-31g(d)
 ****
H 0
6-31g(d)
 ****
N  0
6-31g(d)
 ****
P  0
6-31g(d)
 ****
O 0
6-31g(d)
 ****
F 0
6-31g(d)
 **** 

Ir
LANL2DZ



""")

if __name__=="__main__" :
        main()
