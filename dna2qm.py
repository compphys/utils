#!/usr/bin/env python

from pylab import *
import sys

gaussian_template = """%NProcLinda=1
%NProcShared=16
# opt B3LYP/cc-pVDZ
"""

basis_set = 'aug-cc-pVDZ'

atom_charges = {'S':'16.0','P':'15.0', 'Na':'11.0', 'O':'8.0', 'N':'7.0', 'C':'6.0', 'H':'1.0'}

replacement_groups = {
'OH':[
['O', 1.440000, 0.000000, 0.000000],
['H', 1.744214, 0.925273, 0.000000]]
}

def dist(c1, c2):
    return sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2)


def add_cap(bond_atom, replace_atom, replace_group, keep_distance = False):
    vect = array(replace_atom) - array(bond_atom)
    bond = norm(vect)
    vect = vect / norm(vect)
    rotation = [[1,0,0],[0,1,0],[0,0,1]]
    if vect[0] == -1:
        rotation = [[-1,0,0],[0,1,0],[0,0,1]]
    elif vect[0] != 1:
        middle = array([(vect[0] + 1)/2, vect[1]/2, vect[2]/2])
        middle = middle/norm(middle)
        x = middle[0]
        y = middle[1]
        z = middle[2]
        rotation = array([[2*x**2 - 1, 2*x*y, 2*x*z], [2*y*x, 2*y**2 - 1, 2*y*z], [2*z
*x, 2*z*y, 2*z**2 - 1]])
    cap_atoms = []
    replacements = replacement_groups[replace_group]
    if keep_distance:
        old_bond = replacements[0][1]
        delta = bond - old_bond
        for i in range(len(replacements)):
            replacements[i][1] = replacements[i][1] + delta
    for atom in replacement_groups[replace_group]:
        cap_atoms.append([atom[0], bond_atom + dot(rotation, array(atom[1:]))])
    return cap_atoms

def pdbreader(inname):
    infile = open(inname, 'r')
    line = infile.readline()
    while len(line) < 4 or line[:4] != 'ATOM':
        line = infile.readline()
    strands = []
    while 'ATOM' in line:
        resnum = line[22:26]
        resname = line[17:20]
        strand = []
        residue = {}
        while not 'TER   ' in line:            
            if line[17:20] in [' DT', ' DC', ' DA', ' DG']:
                if line[22:26] != resnum:
                    strand.append([resname[-1], residue])
                    residue = {}
                    resnum = line[22:26]
                residue[line[12:16].strip()] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                resname = line[17:20]
            line = infile.readline()
        if residue:
            strand.append([resname[-1], residue])
        if strand:
            strands.append(strand)
        line = infile.readline()
    if dist(strands[0][0][1]['N1'], strands[1][0][1]['N1']) > dist(strands[0][0][1]['N1'], strands[1][-1][1]['N1']):
        strands[1] = strands[1][::-1]
    return strands


def hydrogenator(nH, catom, batoms, patom=False):
    #Adds nH hydrogen atoms bonded to catom, which is bonded to the atoms in batoms.
    #The optional patoms gives an extra atom to form a plane with a single atom in
    #batoms and catom.
    catom = array(catom)
    batoms = [array(atom) for atom in batoms]
    Hdist = 1.09
    Hatoms = []
    vect = [0, 0, 0]
    for atom in batoms:
        avect = atom-catom
        avect = avect/norm(avect)
        vect = vect + avect
    vect = -vect/norm(vect)
    if nH == 1:
        H = catom + vect*Hdist
        Hatoms.append(['H', H.tolist()])
    elif nH == 2:
        pvect = array([-vect[1], vect[0], 0])
        pvect_temp = array([0, -vect[2], vect[1]])
        if norm(pvect_temp) > norm(pvect):
            pvect = pvect_temp
        if len(batoms) > 1:
            pvect = cross(batoms[0] - catom, batoms[1] - catom)
        if patom:
            patom = array(patom)
            pvect = patom-catom
            pvect = pvect - vect*dot(pvect, vect)
        pvect = pvect/norm(pvect)
        alpha = 109.4
        H1 = catom + Hdist * (cos((alpha/2)*2*pi/360) * vect + sin((alpha/2)*2*pi/360) * pvect)
        H2 = catom + Hdist * (cos((alpha/2)*2*pi/360) * vect - sin((alpha/2)*2*pi/360) * pvect)
        Hatoms.append(['H', H1.tolist()])
        Hatoms.append(['H', H2.tolist()])
    elif nH == 3:
        pvect = array([-vect[1], vect[0], 0])
        pvect_temp = array([0, -vect[2], vect[1]])
        if norm(pvect_temp) > norm(pvect):
            pvect = pvect_temp
        if patom:
            patom = array(patom)
            pvect = patom-catom
            pvect = pvect - vect*dot(pvect, vect)
        pvect = pvect/norm(pvect)
        opvect = cross(vect, pvect)
        f1 = 0.333333
        f2 = 0.942809
        f3 = 0.471405
        f4 = 0.816497
        H1 = catom + Hdist * (f1 * vect + f2 * pvect)
        H2 = catom + Hdist * (f1 * vect - f3 * pvect + f4 * opvect)
        H3 = catom + Hdist * (f1 * vect - f3 * pvect - f4 * opvect)
        Hatoms.append(['H', H1.tolist()])
        Hatoms.append(['H', H2.tolist()])
        Hatoms.append(['H', H3.tolist()])
        
        
    return Hatoms


def get_atoms(number, strand, position, no_sugar):
    atoms = []
    name = strand[number-1][0]
    na = strand[number-1][1]
    if name == 'A':
        atoms.append(['N', na['N1']])
        atoms.append(['C', na['C2']])
        atoms += hydrogenator(1, na['C2'], [na['N1'], na['N3']])
        atoms.append(['N', na['N3']])
        atoms.append(['C', na['C4']])
        atoms.append(['C', na['C5']])
        atoms.append(['C', na['C6']])
        atoms.append(['N', na['N6']])
        atoms += hydrogenator(2, na['N6'], [na['C6']], na['N1'])
        atoms.append(['N', na['N7']])
        atoms.append(['C', na['C8']])
        atoms += hydrogenator(1, na['C8'], [na['N7'], na['N9']])
        atoms.append(['N', na['N9']])
        if no_sugar:
            atoms.append(['C', na["C1'"]])
            atoms += hydrogenator(3, na["C1'"], [na['N9']], na['C4'])
    elif name == 'C':
        atoms.append(['N', na['N1']])
        if no_sugar:
            atoms.append(['C', na["C1'"]])
            atoms += hydrogenator(3, na["C1'"], [na['N1']], na['O2'])
        atoms.append(['C', na['C2']])
        atoms.append(['O', na['O2']])
        atoms.append(['N', na['N3']])
        atoms.append(['C', na['C4']])
        atoms.append(['N', na['N4']])
        atoms += hydrogenator(2, na['N4'], [na['C4']], na['N3'])
        atoms.append(['C', na['C5']])
        atoms += hydrogenator(1, na['C5'], [na['C4'], na['C6']])
        atoms.append(['C', na['C6']])
        atoms += hydrogenator(1, na['C6'], [na['C5'], na['N1']])
    elif name == 'G':
        atoms.append(['N', na['N1']])
        atoms += hydrogenator(1, na['N1'], [na['C2'], na['C6']])
        atoms.append(['C', na['C2']])
        atoms.append(['N', na['N2']])
        atoms += hydrogenator(2, na['N2'], [na['C2']], na['N3'])
        atoms.append(['N', na['N3']])
        atoms.append(['C', na['C4']])
        atoms.append(['C', na['C5']])
        atoms.append(['C', na['C6']])
        atoms.append(['O', na['O6']])
        atoms.append(['N', na['N7']])
        atoms.append(['C', na['C8']])
        atoms += hydrogenator(1, na['C8'], [na['N7'], na['N9']])
        atoms.append(['N', na['N9']])
        if no_sugar:
            atoms.append(['C', na["C1'"]])
            atoms += hydrogenator(3, na["C1'"], [na['N9']], na['C4'])
    elif name == 'T':
        atoms.append(['N', na['N1']])
        if no_sugar:
            atoms.append(['C', na["C1'"]])
            atoms += hydrogenator(3, na["C1'"], [na['N1']], na['C2'])
        atoms.append(['C', na['C2']])
        atoms.append(['O', na['O2']])
        atoms.append(['N', na['N3']])
        atoms += hydrogenator(1, na['N3'], [na['C2'], na['C4']])
        atoms.append(['C', na['C4']])
        atoms.append(['O', na['O4']])
        atoms.append(['C', na['C5']])
        atoms.append(['C', na['C6']])
        atoms += hydrogenator(1, na['C6'], [na['C5'], na['N1']])
        atoms.append(['C', na['C7']])
        atoms += hydrogenator(3, na['C7'], [na['C5']], na['O4'])

    if not no_sugar:
        atoms.append(['C', na["C1'"]])
        if name == 'A' or name == 'G':
            atoms += hydrogenator(1, na["C1'"], [na["C2'"], na["O4'"], na["N9"]])
        else:
            atoms += hydrogenator(1, na["C1'"], [na["C2'"], na["O4'"], na["N1"]])
        atoms.append(['C', na["C2'"]])
        atoms += hydrogenator(2, na["C2'"], [na["C1'"], na["C3'"]])
        atoms.append(['C', na["C3'"]])
        atoms += hydrogenator(1, na["C3'"], [na["C2'"], na["C4'"], na["O3'"]])
        atoms.append(['C', na["C4'"]])
        atoms += hydrogenator(1, na["C4'"], [na["C3'"], na["O4'"], na["C5'"]])
        atoms.append(['O', na["O4'"]])
        atoms.append(['C', na["C5'"]])
        atoms += hydrogenator(2, na["C5'"], [na["C4'"], na["O5'"]])
        if position == 'single':
            atoms += add_cap(na["C3'"], na["O3'"], 'OH', True)
            atoms += add_cap(na["C5'"], na["O5'"], 'OH', True)
        elif position == 'start':
            atoms += add_cap(na["C5'"], na["O5'"], 'OH', True)
            atoms.append(['O', na["O3'"]])
        elif position == 'end':
            atoms += add_cap(na["C3'"], na["O3'"], 'OH', True)
            atoms.append(['O', na["O5'"]])
            atoms.append(['P', na["P"]])
            atoms += add_cap(na["P"], na["OP1"], 'OH', True)
            atoms += add_cap(na["P"], na["OP2"], 'OH', True)
        else:
            atoms.append(['O', na["O5'"]])
            atoms.append(['O', na["O3'"]])
            atoms.append(['P', na["P"]])
            atoms += add_cap(na["P"], na["OP1"], 'OH', True)
            atoms += add_cap(na["P"], na["OP2"], 'OH', True)
    return atoms
            


def xyzwrite(atoms, outname):
    outfile = open(outname, 'w')
    outfile.write('  ' + str(len(atoms)) + '\n\n')
    for atom in atoms:
        newline = '{0:>5}{1:>12.6f}{2:>12.6f}{3:>12.6f}\n'.format(atom[0], atom[1][0], atom[1][1], atom[1][2])
        outfile.write(newline)


def molwrite(atoms, outname, comment = ''):
    atom_types = {}
    for atom in atoms:
        if atom[0] in atom_types:
            atom_types[atom[0]] = atom_types[atom[0]] + [atom[1]]
        else:
            atom_types[atom[0]] = [atom[1]]
    outfile = open(outname, 'w')
    outfile.write('ATOMBASIS\n')
    outfile.write(comment + '\n')
    outfile.write('-------------------------\n')
    outfile.write('Atomtypes=' + str(len(atom_types)) + ' Generators=0 Angstrom\n')
    for atom_type in atom_types:
        outfile.write('Charge=' + atom_charges[atom_type] + ' Atoms=' + str(len(atom_types[atom_type])) + ' Basis=' + basis_set + '\n')
        for atom in atom_types[atom_type]:
            newline = '{0:<}{1:>28.6f}{2:>14.6f}{3:>14.6f}\n'.format(atom_type, atom[0], atom[1], atom[2])
            outfile.write(newline)
    outfile.write('End of input\n')

def comwrite(atoms, outname, comment='Title Card Required\n'):
    outfile = open(outname, 'w')
    for line in gaussian_template:
        outfile.write(line)
    outfile.write('\n')
    outfile.write(comment)
    outfile.write('\n')
    outfile.write('0 1\n')
    for atom in atoms:
        newline = '{0:>2}{1:>28.8f}{2:>14.8f}{3:>14.8f}\n'.format(atom[0], atom[1][0], atom[1][1], atom[1][2])
        outfile.write(newline)
    outfile.write('\n')
    outfile.write('End of file\n')
    outfile.write('\n')
    

inname = sys.argv[1]
outname = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])
nosugar = False
if 'nosugar' in sys.argv:
    nosugar = True
pair = False
if 'pair' in sys.argv:
    pair = True
strands = pdbreader(inname)
if 'second' in sys.argv:
    strands = [strands[1], strands[0]]
if start < 1:
    print "Start position must be positive integer above zero"
    sys.exit()
if end > len(strands[0]):
    print "End position cannot be larger than strand length (" + len(strands[0]) + ")"
    sys.exit()
if end < start:
    print "Start position must be greater than or equal to end position"
    sys.exit()
atoms = []

for i in range(start, end+1):
    atoms += get_atoms(i, strands[0], 'single', nosugar)
    if pair:
        atoms += get_atoms(i, strands[1], 'single', nosugar)

sequence = ''
sequence2 = ''
for i in range(start-1, end):
    sequence += strands[0][i][0]
    sequence2 += strands[1][i][0]
               
if '.xyz' in outname:
    xyzwrite(atoms, outname)
elif '.com' in outname:
    comwrite(atoms, outname,
             sequence + '; ' + 
             inname + '; ' + 
             str(start) + '; ' + 
             str(end) + '; ' + 
             'pair=' +  str(pair) + '; ' + 
             'nosugar=' + str(nosugar))
else:
    molwrite(atoms, outname, 
             sequence + '; ' + 
             inname + '; ' + 
             str(start) + '; ' + 
             str(end) + '; ' + 
             'pair=' +  str(pair) + '; ' + 
             'nosugar=' + str(nosugar))

print "Wrote sequence:"
print sequence
if pair:
    print '|' * len(sequence)
    print sequence2
line = ''
if not nosugar:
    line += 'with sugar '
line +=  "to file " + outname
print line
