from sys import argv
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import is_aa
from math import sqrt, atan2, degrees
from matplotlib import pyplot, rc
import operator
import argparse


# pythonic vector ops are extensible to higher dimensions but slower here
def add(v0, v1):
    return v0[0] + v1[0], v0[1] + v1[1], v0[2] + v1[2]


def add_pythonic(v0, v1):
    return tuple(map(operator.add, v0, v1))


def subtract(v0, v1):
    return v0[0] - v1[0], v0[1] - v1[1], v0[2] - v1[2]


def subtract_pythonic(v0, v1):
    return tuple(map(operator.sub, v0, v1))


def normalize(v):
    mag = sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    return v[0] / mag, v[1] / mag, v[2] / mag


def normalize_pythonic(v):
    mag = sqrt(sum([x ** 2 for x in v]))
    return tuple([x / mag for x in v])


def dot(v0, v1):
    return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2]


def dot_pythonic(v0, v1):
    return sum([pair[0] * pair[1] for pair in zip(v0, v1)])


def cross(v0, v1):
    return (v0[1] * v1[2] - v0[2] * v1[1], -(v0[0] * v1[2] - v0[2] * v1[0]),
            v0[0] * v1[1] - v0[1] * v1[0])


def torsion(a, b, c, d):
    # had to change to sub(a, b) from (b, a) to correct 180-degree rotation
    # (probably wrong)
    b1 = subtract(a, b)
    b2 = subtract(b, c)
    b3 = subtract(c, d)
    n1 = normalize(cross(b1, b2))
    n2 = normalize(cross(b2, b3))
    m1 = cross(n1, normalize(b2))
    x = dot(n1, n2)
    y = dot(m1, n2)
    return degrees(atan2(y, x))

# argparser = argparse.ArgumentParser()
# argparser.add_argument('filename', type=str, nargs=1, help='pdb filename')
# argparser.add_argument('chain', type=str, nargs=1, default='A', help='pdb filename')
script, filename = argv
# pull pdb id out of path string e.g. '../1axc.pdb'
pdb_id = filename.split('/')[-1].split('.')[0].upper()
parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure(pdb_id, filename)
phis = []
psis = []
"""
# plot all residues
for model in structure:
    # using list(chain) to avoid weird residue indices
    for lchain in [list(chain) for chain in model]:
        for i, res in enumerate(lchain):
            # next residue in chain for torsion calculation
            nex = lchain[i+1]
            # only calculate torsion if next residue is also an AA
            if is_aa(nex):
                # make groups of relevant atoms to calculate angle
                # phi: Ci, Ni+1, Ci+1alpha, Ci+1
                # psi: Ni, Cialpha, Ci, Ni+1
                phi_atoms = (res['C'], nex['N'], nex['CA'], nex['C'])
                psi_atoms = (res['N'], res['CA'], res['C'], nex['N'])
                # make group of coordinates of relevant atoms
                phi_coords = (a.get_coord() for a in phi_atoms)
                psi_coords = (a.get_coord() for a in psi_atoms)
                # calculate torsions from groups of coordinates and store
                phi = torsion(*phi_coords)
                phis.append(phi)
                psi = torsion(*psi_coords)
                psis.append(psi)
            # stop calculating torsions when AAs run out
            else:
                break
"""
# plot chain A residues
for model in structure:
    # using list(chain) to avoid weird residue indices
    lchain = list(model['A'])
    for i, res in enumerate(lchain):
        # next residue in chain for torsion calculation
        n = lchain[i+1]
        # only calculate torsion if next residue is also an AA
        if is_aa(n):
            # make groups of relevant atoms to calculate angle
            phi_atoms = (res['C'], n['N'], n['CA'], n['C'])
            psi_atoms = (res['N'], res['CA'], res['C'], n['N'])
            # make group of coordinates of relevant atoms
            phi_coords = (a.get_coord() for a in phi_atoms)
            psi_coords = (a.get_coord() for a in psi_atoms)
            # calculate torsions from groups of coordinates and store
            phi = torsion(*phi_coords)
            phis.append(phi)
            psi = torsion(*psi_coords)
            psis.append(psi)
        # stop calculating torsions when AAs run out
        else:
            break

rc('text', usetex=True)
pyplot.plot(phis, psis, 'bo')
pyplot.axis([-180, 180, -180, 180])
pyplot.gca().set_aspect('equal')
pyplot.xticks([-180, -90, 0, 90, 180])
pyplot.yticks([-180, -90, 0, 90, 180])
pyplot.xlabel(r'$\phi$', size=14)
pyplot.ylabel(r'$\psi$', size=14)
pyplot.title('Ramachandran plot for PDB {}'.format(pdb_id))
pyplot.grid()
pyplot.show()
