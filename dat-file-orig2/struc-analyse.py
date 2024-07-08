from ase.io import read
import numpy as np
import pandas as pd

to_read = "POSCAR"
lay_sep = 1.0 # Angstoms

slab = read(to_read)
atoms= read(to_read)


z_coords = atoms.positions[:, 2]
atomic_num = atoms.numbers

sorted_z = np.sort(z_coords)
print(z_coords,sorted_z)

crrctd_to_base = sorted_z-sorted_z[0] #to bring base layer to the bottom of the cell
int_rep        = np.rint(crrctd_to_base) #to represent the z values of the atoms as integers	
print(int_rep)

data = [[z_coords[i],atomic_num[i]] for i in range(len(atoms)) ]
df_z_symb = pd.DataFrame(data,columns=["Z values","Chemical Symbols"])
dfzsy_sort = df_z_symb.copy()
dfzsy_sort = dfzsy_sort.sort_values(by=["Z values"])

print(df_z_symb,dfzsy_sort)

























#for z in crrctd_to_base:

# Sort the z-coordinates in ascending order
#sorted_z_coords = np.sort(np.unique(z_coords))

# Group the atoms into layers based on their z-coordinate values
#layer_indices = []
#for z in sorted_z_coords:
#    indices = np.where(np.isclose(z_coords, z,atol=5))[0]
#    layer_indices.append(indices)

# Print the indices of atoms in each layer
#for i, indices in enumerate(layer_indices):
#    print(f'Layer {i+1}: {indices}')
