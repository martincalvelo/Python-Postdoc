import optparse 
import MDAnalysis 

desc="""
script para obtener los index de los atomos de la parte de QM. Saca un pdb llamado QM.pdb. No saca link
Ej:
    python extract_indexes_v2.py -i frame770_vcon20.pdb -o QMindex.txt -s '((resid 116 117 and not name C O N H HA CA) or resid 993 to 996) and not (resname WAT or type 'Na+' 'Cl-') or resid 50934'
"""
def main():
    parser = optparse.OptionParser(description=desc, version='%prog version 1.0')
    parser.add_option('-i', '--input', help='pdb input', action='store')
    parser.add_option('-s', '--sele', help='selection', action='store')
    parser.add_option('-o', '--output', help='txt out', action='store')
    options, arguments = parser.parse_args()
########################

    u = MDAnalysis.Universe(str(options.input)) 

    selection = str(options.sele)

    atoms = u.select_atoms(selection)
    index=[]; types=[]
    for atom in atoms:
        types.append(atom.type); index.append(atom.index+1)
    atoms.write('QM.pdb') 

    out = open(str(options.output),'w')
    for x in index:
        out.write(str(x)+' ')

    name_out2 = str(options.output)[:-4]+'_byelement.txt'
    out2 = open(name_out2, 'w')
    out2.write('#'+selection+'\n\n')
    set_types=set(types)
    for item in set_types:
        out2.write('\n'+item+':\n')
        para_print=[]
        for i in range(len(types)):
            if item == types[i]:
                para_print.append(index[i])
        for x in para_print:
            out2.write(str(x)+' ')
        out2.write('\n')

if __name__=="__main__" :
        main()
