from ase.io import read
from ase.io.bader import attach_charges


to_read = "POSCAR"
atoms = read(to_read)


attach_charges(atoms)

for atom in atoms:
    print('Atom', atom.symbol, 'Bader charge', atom.charge)
