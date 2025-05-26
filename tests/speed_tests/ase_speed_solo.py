import sys
from ase.io import read, write

cif_file = sys.argv[1]
atoms = read(cif_file)
poscar_name = cif_file.replace('.cif', '_POSCAR')
write(poscar_name, atoms, format='vasp')
