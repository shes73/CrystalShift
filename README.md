# CrystalShift
Crystal cell operating tool: basis changer, format converter, atom coordinates editor and analyzer of molecular structure.

## Installation
You can download the repository using the command below or manually by downloading the archive:
```bash
git clone https://github.com/shes73/CrystalShift.git
```

The program is distributed with open source code, to work with CrystalShift you need to build it. To do this, you need to go to the folder *"main/crystalshift"* or *"main/crystalshift_layers"*, then you can simply compile it via Makefile:
```bash
make all
```
Or do it manually. For the main CrystalShift module use the command below:
```bash
gcc -std=c99 main.c structures.c poscar_parser.c poscar_writer.c cif_parser.c cif_writer.c xyz_parser.c xyz_writer.c atomic_coords_editor.c basis_changer.c -lm -o crystalshift
```
For CrystalShift layers module:
```bash
gcc -std=c99 crystalshift_layers.c -lm -o crystalshift_layers.exe
```

## Shortly about implementation
Let's have a close look at each module. Here are listed modules source files with their features and some points of work.

***structures.c*** is the prime source file of the program. Structures with all information about crystal lattice and atoms is located there. When you run any of the modules: converters, basis changer or coordinate editor, it will definitely write information during parsing into structures, update it, save it in a new file.

***atomic_coords_editor.c*** is module for editing atomic coordinates by adding a vector (x, y, z coordinates) that will be simply added to the coordinates of all atoms in the file.

***FORMAT_parser.c*** are source files for reading and writing all information about lattice and atoms into structures in *structure.c*. In case of POSCAR and xyz formats everything seems to be easy, but CIF files can be varied and even creative.  

## Usage
I strongly recommend not to change the basis and create a supercell at the same time. When you change the basis, Crystal Shift does not write new atomic coordinates to the structure, while when you create a supercell, the new atomic positions are written and saved.

## Examples of usage
### 1. Basis change
![basis_change_visualisation](https://github.com/shes73/CrystalShift/blob/main/images/basis_change.jpg)
```bash
> crystalshift
Choose the option you need:
1. Change basis of unit cell
2. Edit atomic coordinates
3. Converter CIF -> POSCAR
4. Converter CIF -> xyz
5. Converter POSCAR -> CIF
6. Converter POSCAR -> xyz
> 1
Change basis of unit cell
Please, choose the format of the input file: 
1. CIF
2. POSCAR
>1
Enter the path to the CIF file:
> C:\cif_editor\VEWSIC.cif
Enter the vectors of the basis transformation matrix.
a vector:
> 1
> 0
> 1
b vector:
> 0
> 1
> 0
c vector:
> 0
> 0
> 1
Enter the path to save the new CIF file: 
> C:\cif_editor\VEWSIC_edited.cif
```

### 2. Supercell generation
![basis_change_visualisation](https://github.com/shes73/CrystalShift/blob/main/images/supercell_visualisation.jpg)
```bash
> crystalshift
Choose the option you need:
1. Change basis of unit cell
2. Edit atomic coordinates
3. Converter CIF -> POSCAR
4. Converter CIF -> xyz
5. Converter POSCAR -> CIF
6. Converter POSCAR -> xyz
> 1
Change basis of unit cell
Please, choose the format of the input file: 
1. CIF
2. POSCAR
> 1
Enter the path to the CIF file:
> C:\cif_editor\VEWSIC.cif
Enter the vectors of the basis transformation matrix.
a vector:
> 2
> 0
> 0
b vector:
> 0
> 1
> 0
c vector:
> 0
> 0
> 2
Supercell detected. Do you want to add new atom coordinates in the supercell? (yes/no):
> yes
Supercell factors: [2, 1, 2]
Original atom count: 100
Updated atom count: 400
Enter the path to save the new CIF file:
> C:\cif_editor\VEWSIC_supercell.cif
```

### 3. Atomic coordinates editor
![atomic_coords_edition_visualisation](https://github.com/shes73/CrystalShift/blob/main/images/atomic_coords_edition_visualisation.jpg)
```bash
Choose the option you need:
1. Change basis of unit cell
2. Edit atomic coordinates
3. Converter CIF -> POSCAR
4. Converter CIF -> xyz
5. Converter POSCAR -> CIF
6. Converter POSCAR -> xyz
> 2
Edit atomic coordinates in the unit cell
Please, choose the format of the input file:
1. CIF
2. POSCAR
> 1
Enter the path to the CIF file:
> C:\cif_editor\VEWSIC_edited.cif  
Enter the three float numbers separated by spaces that will be added to the coordinates of all atoms in the file (e.g. 0.25 0.0 0.0):
> 0.375 0.0 0.0
Enter the path to save the new CIF file:
> C:\cif_editor\VEWSIC_layer_moved.cif
```

### 4. Converters
![CIF to POSCAR example](https://github.com/shes73/CrystalShift/blob/main/images/format_converter_example.jpg)
![CIF to POSCAR example](https://github.com/shes73/CrystalShift/blob/main/images/CIF_to_POSCAR.jpg)
```bash
> crystalshift
Choose the option you need:
1. Change basis of unit cell
2. Edit atomic coordinates
3. Converter CIF -> POSCAR
4. Converter CIF -> xyz
5. Converter POSCAR -> CIF
6. Converter POSCAR -> xyz
> 3
Converter CIF -> POSCAR
Enter the path to the CIF file:
> C:\cif_editor\VEWSIC.cif
Choose the order for writing atoms to POSCAR:
1. From lightest to heaviest;
2. From heaviest to lightest;
3. Custom order (input symbols separated by spaces).
Enter your choice:
> 1
Enter the path to the POSCAR file:
> C:\cif_editor\VEWSIC_POSCAR_1
```
For custom order of atoms:
```bash
Enter your choice:
> 3
Enter custom order:
> H C O Br
Enter the path to the POSCAR file:
> C:\cif_editor\VEWSIC_POSCAR_3
```

![CIF to xyz example](https://github.com/shes73/CrystalShift/blob/main/images/CIF_to_xyz.jpg)
```bash
> crystalshift
Choose the option you need:
1. Change basis of unit cell
2. Edit atomic coordinates
3. Converter CIF -> POSCAR
4. Converter CIF -> xyz
5. Converter POSCAR -> CIF
6. Converter POSCAR -> xyz
> 4
Converter CIF -> xyz
Enter the path to the CIF file:
> C:\cif_editor\VEWSIC.cif
Choose the format for writing xyz file:
1. Standard xyz format;
2. Extended xyz format (with lattice parameters in the comment line);
Enter your choice:
> 2
Enter the path to the xyz file:
> C:\cif_editor\VEWSIC.xyz
```

## Requirements
Main CrystalShift modules require only standard C libraries.

## License
This project is licensed under the GNU General Public License, Ver. 3 - see the LICENSE.md file for details.

## Contact information
Email: i.isupova@g.nsu.ru

Github: @shes73

Telegram: @shes73
