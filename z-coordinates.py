#!/usr/bin/python3
# Get z-coordinates from FHI-aims geometry
# Give geometry (.in) file as argument and a list of species in the molecule
# (used for calculating total charge of molecule and its height/corrugation)

import numpy as np
import sys

print('#############################')
print('# z-coordinate analysis')
print('#############################')

# Read geometry file (only lattice vectors and atom coordinates)
a = [] # Output data variable
f = open(sys.argv[1], 'r')
for line in f:
    if line.startswith('lattice_vector') or line.startswith('atom'):
        a.append(line)
f.close()

# Find what are the species
species = np.array([], dtype='str')
for line in a:
    if line.startswith('atom'):
        line = line.split()
        s = line[4]
        if s not in species:
            species = np.append(species, s)
N_species = np.size(species) # Number of species (for creating data matrices)


# Get z-coordinates of each species
# 2D array of z-coordinates: rows = species, columns = atom z-coordinates
z_coords = np.empty((N_species, 9999), dtype='float') # Arbitrary large size
z_coords.fill(999999) # Fill unused matrix elements with '999999'

# Write z-coords to 2D matrix
for ind, s in enumerate(species):
    s_number = 0 # initialize atom counter
    for line in a:
        if line.startswith('atom'):
            line = line.split()
            if line[4] == s: # If correct species
                z_coords[ind, s_number] = line[3] # Add to matrix 
                s_number = s_number + 1

# Convert z-coords if fractional (check from first atom, after lattice)
if a[3].split()[0] == 'atom_frac':
    box_height = float(a[2].split()[3]) # Get box z height from 3rd line
    for i, line in enumerate(z_coords):
        for j, z in enumerate(line):
            if z != 999999:
                z_coords[i, j] = z_coords[i, j] * box_height


##################### MOLECULE ANALYSIS STARTS HERE ##########################

# If arguments given for molecule species
if np.size(sys.argv) > 2:
    molspec = sys.argv[2:]
    for i in molspec:
        if i not in species:
            print('*** Error: Given molecule species ' + i + \
                    ' not found in geometry')
            sys.exit()

    # Gather indices of molecule species into an array
    molind = np.array([], dtype='int') # Indices of molecule species
    for i in molspec: # Gather indices of molecule species in array
        for ind, s in enumerate(species):
            if i == s:
                molind = np.append(molind, ind)

    # Define a new array of molecule z coordinates
    z_coord_mol = z_coords[molind, :]

    # Find lowest and highest molecule atoms
    mol_sum = 0.0 # For calculating average
    mol_min = 999999.0 # Arbitrary large
    mol_max = -999999.0 # Arbitrary small
    at_number = 0 # Initialize atom counter
    for i, line in enumerate(z_coord_mol):
        for j, z in enumerate(line):
            if z != 999999:
                c = z_coord_mol[i, j]
                mol_sum = mol_sum + c
                if c < mol_min:
                    mol_min = c
                if c > mol_max:
                    mol_max = c
                at_number = at_number + 1
    mol_avg = mol_sum / at_number # Average height

##################### MOLECULE ANALYSIS ENDS HERE ##########################


# Calculate average, min, and max height for each species
N_per_species = np.array([], dtype='int') # N atoms per species
z_avg = np.array([], dtype='float')
z_min = np.array([], dtype='float')
z_max = np.array([], dtype='float')
for i, line in enumerate(z_coords):
    s_number = 0 # Initialize atom counter
    coord_sum = 0.0 # For calculating average
    coord_min = 999999.0 # Arbitrary large
    coord_max = -999999.0 # Arbitrary small
    for j, z in enumerate(line):
        if z != 999999:
            c = z_coords[i, j]
            coord_sum = coord_sum + c
            if c < coord_min:
                coord_min = c
            if c > coord_max:
                coord_max = c
            s_number = s_number + 1
    avg = coord_sum / s_number # Average z-coordinate
    N_per_species = np.append(N_per_species, s_number)
    z_avg = np.append(z_avg, avg)
    z_min = np.append(z_min, coord_min)
    z_max = np.append(z_max, coord_max)

# Calculate deviation of z-coordinates (i.e. corrugation) for each species
z_diff = z_max - z_min

# Sort data by average layer height (might not be in order in the geometry file)
inds = z_avg.argsort() # Check order from avg height
sort_species = species[inds]
sort_N_per_species = N_per_species[inds]
sort_z_avg = z_avg[inds]
sort_z_min = z_min[inds]
sort_z_max = z_max[inds]
sort_z_diff = z_diff[inds]

# Calculate layer separations from average coordinates
#sort_z_sep = np.diff(sort_z_avg)

# Calculate molecule height from top layer (from their average heights)
# Top layer found by subtracting the number of molecule species in the array
# (counting backwards)
if np.size(sys.argv) > 2:
    N_mol_species = np.size(molind) # Number of molecule species
    top_substrate_ind = N_mol_species + 1 # How many counted back from array end
    mol_height = mol_avg - sort_z_avg[-top_substrate_ind]

# Print z-coordinates
print('# Species in the order of average z-coordinate, from top to bottom')
print('# N_at.  z_avg    z_min    z_max    z_diff')
for i in np.flip(np.arange(np.size(sort_species))):
    print(sort_species[i] + ' ' + str(sort_N_per_species[i]) + '     ' + \
        '{:.4f}'.format(sort_z_avg[i]) + '  ' + \
        '{:.4f}'.format(sort_z_min[i]) + '  ' + \
        '{:.4f}'.format(sort_z_max[i]) + '  ' + \
        '{:.4f}'.format(sort_z_diff[i]))

# Print layer separations
print('\n# Layer separations based on average z-coordinates')
#for i in np.flip(np.arange(np.size(sort_z_sep))):
#    print(sort_species[i+1] + '-' + sort_species[i] + \
#            ': ' + '{:.4f}'.format(sort_z_sep[i]))
loopsize = np.size(sort_z_avg)
for i in np.arange(loopsize):
    for j in np.arange(i+1, loopsize):
        separation = sort_z_avg[j] - sort_z_avg[i]
        print(sort_species[i] + '-' + sort_species[j] + ': ' + \
                    '{:.4f}'.format(separation))

# Print molecule data if arguments given for molecule species
if np.size(sys.argv) > 2:
    print()
    print('# Molecule defined with species: ' + str(molspec))
    print('Molecule has ' + str(at_number) + ' atoms')
    print('Molecule top height: ' + '{:.4f}'.format(mol_max))
    print('Molecule bottom height: ' + '{:.4f}'.format(mol_min))
    print('Molecule average height: ' + '{:.4f}'.format(mol_avg))
    print('Molecule corrugation (bend = top - bottom): ' \
            + '{:.4f}'.format(mol_max - mol_min))
    print('Molecule height from ' + sort_species[-top_substrate_ind] \
            + ' (calculated from avg): ' + '{:.4f}'.format(mol_height))
    print()
    print('Molecule height from all species (calculated from avg):')
    for i in np.arange(loopsize - N_mol_species):
        mol_spec_height = mol_avg - sort_z_avg[i]
        print(sort_species[i] + '-mol: ' + '{:.4f}'.format(mol_spec_height))

