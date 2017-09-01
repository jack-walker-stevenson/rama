import argparse
from sys import argv  # remove once argparse works

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import is_aa
from matplotlib import pyplot, rc

from vector_ops import torsion


argparser = argparse.ArgumentParser()
argparser.add_argument('filename', type=str, help='name of pdb file')
argparser.add_argument('-v', '--verbose', action='store_true',
                       help='verbose mode')
argparser.add_argument('--chain', type=str,
                       help='specific chain to plot (default is all chains)')
args = argparser.parse_args()

# pull pdb id out of path string e.g. '../1axc.pdb'
pdb_id = args.filename.split('/')[-1].split('.')[0].upper()
parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure(pdb_id, args.filename)
if args.verbose:
    print("Verbose mode enabled. Parsing PDB {}.".format(pdb_id))
phis = []
psis = []

if args.chain:  # plot residues from a single chain
    if args.verbose:
        print("Using residues from chain {}.".format(args.chain))
    for model in structure:
        # using list() to avoid weird residue indices
        lchain = list(model[args.chain])
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
                phis.append(torsion(*phi_coords))
                psis.append(torsion(*psi_coords))
            # stop calculating torsions when AAs run out
            else:
                break
else:  # plot residues from all chains
    if args.verbose:
        print("Using all residues.")
    for model in structure:
        # using list(chain) to avoid weird residue indices
        for lchain in [list(chain) for chain in model]:
            # if args.verbose: ...
            for i, res in enumerate(lchain):
                # next residue in chain for torsion calculation
                nex = lchain[i+1]
                # only calculate torsion if next residue is also an AA
                # (assume first residue is an AA)
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
                    phis.append(torsion(*phi_coords))
                    psis.append(torsion(*psi_coords))
                # stop calculating torsions when AAs run out
                else:
                    break

rc('text', usetex=True)
# throw out first phi and last psi to plot correct phi, psi pairs
pyplot.plot(phis[:-1], psis[1:], 'bo')
pyplot.axis([-180, 180, -180, 180])
pyplot.gca().set_aspect('equal')  # equal axis scale -> square plot
pyplot.xticks([-180, -90, 0, 90, 180])
pyplot.yticks([-180, -90, 0, 90, 180])
pyplot.xlabel(r'$\phi$', size=14)
pyplot.ylabel(r'$\psi$', size=14, rotation=0)
pyplot.title('Ramachandran plot for PDB {}'.format(pdb_id))
pyplot.grid()
pyplot.show()
