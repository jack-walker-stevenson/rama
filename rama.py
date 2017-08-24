from sys import argv
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import is_aa
from math import sqrt, atan2, degrees
from matplotlib import pyplot, rc


def add_non_pythonic(v0, v1):
    return v0[0] + v1[0], v0[1] + v1[1], v0[2] + v1[2]


def add(v0, v1):
    return tuple([pair[0] + pair[1] for pair in zip(v0, v1)])


def subtract_non_pythonic(v0, v1):
    return v0[0] + v1[0], v0[1] + v1[1], v0[2] + v1[2]


def subtract(v0, v1):
    return tuple([pair[0] - pair[1] for pair in zip(v0, v1)])


def normalize_non_pythonic(v):
    mag = sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    return v[0] / mag, v[1] / mag, v[2] / mag


def normalize(v):
    mag = sqrt(sum([x ** 2 for x in v]))
    return [x / mag for x in v]


def dot_non_pythonic(v0, v1):
    return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2]


def dot(v0, v1):
    return sum([pair[0] * pair[1] for pair in zip(v0, v1)])


def cross(v0, v1):
    return (v0[1] * v1[2] - v0[2] * v1[1], -(v0[0] * v1[2] - v0[2] * v1[0]),
            v0[0] * v1[1] - v0[1] * v1[0])


def torsion(a, b, c, d):
    b1 = subtract(b, a)
    b2 = subtract(c, b)
    b3 = subtract(d, c)
    n1 = normalize(cross(b1, b2))
    n2 = normalize(cross(b2, b3))
    m1 = cross(n1, normalize(b2))
    x = dot(n1, n2)
    y = dot(m1, n2)
    return degrees(atan2(y, x))

# phi: Ci, Ni+1, Ci+1alpha, Ci+1
# psi: Ni, Cialpha, Ci, Ni+1
script, filename = argv
parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure("test_id", filename)  # actual id later
phis = []
psis = []
"""
for model in structure:
    # using list(chain) to avoid weird residue indices
    for lchain in [list(chain) for chain in model]:
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
"""
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
pyplot.axis('equal')
pyplot.plot(phis, psis, 'o')
pyplot.axis([-180, 180, -180, 180])
pyplot.xticks([-180, -90, 0, 90, 180])
pyplot.yticks([-180, -90, 0, 90, 180])
pyplot.xlabel(r'$\phi$', size=14)
pyplot.ylabel(r'$\psi$', size=14)
pdb_name = filename.split('.')[0].upper()
pyplot.title('Ramachandran plot for PDB {}'.format(pdb_name))
pyplot.grid()
pyplot.show()
