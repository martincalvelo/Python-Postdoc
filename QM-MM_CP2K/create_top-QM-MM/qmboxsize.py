import optparse 
import MDAnalysis 
########################
# Constants
########################
folg=7
########################

desc="""
script para obtener el tamanho de la cell de QM. Se selecciona la parte que queremos
python qmboxsize.py -i frame770_vcon20.pdb -s '((resid 116 117 and not name C O N H HA) or resid 993 to 996) and not (resname WAT or type 'Na+' 'Cl-') or resid 50934' 
"""
def main():
    parser = optparse.OptionParser(description=desc, version='%prog version 1.0')
    parser.add_option('-i', '--infile', help='pdb input', action='store')
    parser.add_option('-s', '--sele', help='sel con caps', action='store')
    options, arguments = parser.parse_args()
########################

    u = MDAnalysis.Universe(str(options.infile)) 
    atoms = u.select_atoms(str(options.sele))

    x_s=[]; y_s=[]; z_s=[]
    for atom in atoms:
        x_s.append(atom.position[0])
        y_s.append(atom.position[1])
        z_s.append(atom.position[2])

    data = []
    data.append(min(x_s))
    data.append(min(y_s))
    data.append(min(z_s))
    data.append(max(x_s))
    data.append(max(y_s))
    data.append(max(z_s))

    x1=float(data[0])
    x2=float(data[3])
    y1=float(data[1])
    y2=float(data[4])
    z1=float(data[2])
    z2=float(data[5])
    x=folg+x2-x1
    y=folg+y2-y1
    z=folg+z2-z1

    print('')
    print('Slack:',folg)
    print(x)
    print(y)
    print(z)
    print('')
    print('Done')
    print('')

if __name__=="__main__" :
        main()
