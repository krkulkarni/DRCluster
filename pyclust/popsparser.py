__author__ = 'kulkarnik'

""" Parses the pops.out and neli.out files created by the command:
 ~/software/pops-1.6.2/src/pops
 --pdb ../K9F313.S_2OHF_A_0256.pdb
 --popsOut K9.pops.out
 --residueOut
 --neighbourOut
 --atomOut

 with the program POPS"""

import csv
import sys

class Atom:
    def __init__(self,number,atom,residue,resnum,SASA,neighbors):
        self.number = number
        self.atom = atom,
        self.residue = residue
        self.resnum = resnum
        self.sasa = SASA
        self.neighbors = neighbors

        aacharge = {
        "GLY": "np", "ALA": "np", "VAL": "np", "LEU": "np", "MET": "np", "ILE": "np",
        "SER": "pol", "THR": "pol", "CYS": "pol", "PRO": "pol", "GLN": "pol", "ASN": "pol",
        "PHE": "ar", "TRP": "ar", "TYR": "ar",
        "LYS": "pos", "HIS": "pos", "ARG": "pos",
        "GLU": "neg", "ASP": "neg"
    }
        self.type = aacharge[self.residue]

class Surface:
    def __init__(self,num,type,area,atoms):
        self.num = num
        self.type = type
        self.area = area
        self.atoms = atoms

def parsepops(popsfile, nelifile, pdbfile):
    atomdict = {}
    with open(popsfile, 'rb') as f:
        parser = csv.reader(f, delimiter='\n')
        for x in range(4):
            parts = next(parser)
        try:
            while True:
                parts = next(parser)[0].split()
                num = int(parts[0])
                atom = Atom(int(parts[0]),parts[1],parts[2],int(parts[4]),float(parts[5]),int(parts[7]))
                atomdict[num] = atom

        except IndexError:
            print "beginning of residue list"
        except StopIteration:
            print "went too far"
            sys.exit(2)

    neighbordict = {}
    with open(nelifile, 'rb') as f:
        parser = csv.reader(f, delimiter='\n')
        for x in range(3):
            parts = next(parser)
        try:
            while True:
                parts = next(parser)[0].split()
                key = int(parts[0].split(":")[0])
                neighbors = []
                for neighb in parts[2:]:
                    neighbors.append(int(neighb.split(":")[0]))
                neighbordict[key] = neighbors

        except StopIteration:
            print "End of Iteration"
            print ""

    pdbdict = {}
    with open(pdbfile, 'rb') as f:
        parser = csv.reader(f, delimiter='\n')
        line = next(parser)[0]
        try:
            while True:
                parts = line.split()
                key = int(parts[1])
                pdbdict[key] = line
                line = next(parser)[0]
        except IndexError:
            "Termination reached"

    return atomdict, neighbordict, pdbdict

def findsurfaces(atomdict,neighbordict):
    surfacenum = 1
    surfacelist = []
    foundatoms = set()

    for key in neighbordict.keys():
        if key in foundatoms:
            continue
        foundatoms.add(key)
        atom = atomdict[key]
        surface = Surface(surfacenum,atom.type,atom.sasa,[key])
        surfacenum += 1
        surface,foundatoms = recfind(key,foundatoms,atomdict,neighbordict,surface)
        surfacelist.append(surface)

    return surfacelist

def recfind(currentkey,foundatoms,atomdict,neighbordict,surface):
    for neighbor in neighbordict[currentkey]:
        if neighbor in foundatoms:
            continue
        neighboratom = atomdict[neighbor]
        flag = checkmatchcrit(neighboratom,atomdict[currentkey],surface)
        if flag == True:
            foundatoms.add(neighbor)
            surface.atoms.append(neighbor)
            surface.area += neighboratom.sasa
            surface, foundatoms = recfind(neighbor,foundatoms,atomdict,neighbordict,surface)

    return surface, foundatoms

def checkmatchcrit(neighbor,atom,surface):
    if neighbor.type == surface.type:
        return True
    else:
        return False

if __name__ == "__main__":

    pdbname = sys.argv[1]
    parts = pdbname.split("/")
    dir = "/".join(parts[:-1])+""

    popsfile = dir + "pops.out"
    nelifile = dir + "neli.out"
    surfacepml = dir + "surface.pml"
    surfacetxt = dir + "info.surfaces.txt"

    atomdict, neighbordict,pdbdict = parsepops(popsfile,nelifile,pdbname)
    surfacelist = findsurfaces(atomdict,neighbordict)
    totalarea = 0
    totalnum = 0
    with open (surfacetxt, "w") as g:
        for surface in surfacelist:
            g.write("SURFACE TYPE: " + surface.type + '\n')
            g.write("SURFACE NUM: " + str(surface.num) + '\n')
            g.write("SURFACE AREA: " + str(surface.area) + '\n\n')
    with open (surfacepml, "w") as f:
        f.write("load " + pdbname)
        f.write("\n\n")
        for surface in surfacelist:
            name = str(surface.num) + '_' + surface.type
            f.write("select " + name + ', id ')
            for atom in surface.atoms:
                f.write(str(atom))
                f.write('+')
            f.write('\n')
            #f.write('show surface, ' + name + "\n")
            f.write('\n')
            totalarea += surface.area
            totalnum += len(surface.atoms)

        f.write('group positive, *_pos\n')
        f.write('set surface_color, red, positive\n')
        f.write('group negative, *_neg\n')
        f.write('set surface_color, blue, negative\n')
        f.write('group nonpolar, *_np\n')
        f.write('set surface_color, white, nonpolar\n')
        f.write('group polar, *_pol\n')
        f.write('set surface_color, green, polar\n')
        f.write('group aromatic, *_ar\n')
        f.write('set surface_color, yellow, aromatic\n')

    print "Total area:", totalarea
    print "Total num:", totalnum,
