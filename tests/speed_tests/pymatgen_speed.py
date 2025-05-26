import sys
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp import Poscar

cif_file = sys.argv[1]

parser = CifParser(cif_file, occupancy_tolerance=100)
structure = parser.get_structures()[0] 
poscar = Poscar(structure)

poscar_file = cif_file.replace('.cif', '_POSCAR')
poscar.write_file(poscar_file)
