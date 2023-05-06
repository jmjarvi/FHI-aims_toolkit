#!/usr/bin/python3
# Get Mulliken charges from FHI-aims output
# Give input file as argument
# This script also writes charges to output file 'mulliken_partial_charges.dat'

import numpy as np
import sys

print('#############################')
print('# Mulliken charge analysis')
print('# Positive = losing electrons')

# Read output file
a = [] # Output data variable
f = open(sys.argv[1], 'r')
for line in f:
    a.append(line)
f.close()

# Find how many species there are
for line in a:
    if '| Number of species' in line:
        line = line.split()
        N_species = int(line[5])
        break

# Find what are the species
species = np.array([], dtype='str')
i = 0
while i < N_species:
    for line in a:
        if 'Reading configuration options for species' in line:
            line = line.split()
            species = np.append(species, line[5])
            print('# Species ' + str(i+1) + ': ' + species[i])
            i = i + 1
print('#############################')

# Find which species have which atom numbers
# 2D array of atom numbers: rows = species, columns = atom numbers
a_numbers = np.empty((N_species, 9999), dtype='int') # Arbitrary large size
a_numbers.fill(9999) # Fill unused matrix elements with '9999'

# Find start and end lines of atomic structure
for ind, line in enumerate(a):
    if '| Atomic structure:' in line:
        a_start = ind + 2
    if 'Lattice parameters for 3D lattice' in line:
        a_end = ind - 2
        break

# Write atom numbers to a_numbers matrix
for ind, s in enumerate(species):
    s_number = 0 # initialize atom counter
    for i in np.arange(a_end - a_start + 1) + a_start: # Loop line numbers
        line = a[i].split()
        if line[3] == s: # If correct species
            a_numbers[ind, s_number] = int(line[1][:-1]) # Add to matrix (no ':')
            s_number = s_number + 1

# Calculate numbers of atoms per species
N_atoms = np.array([], dtype='int')
for ind, line in enumerate(a_numbers):
    count = 0
    for i in line:
        if i == 9999:
            break
        else:
            count = count + 1
    N_atoms = np.append(N_atoms, count)


# Find start and end lines of Mulliken analysis
for ind, line in enumerate(a):
    if 'Starting Mulliken Analysis' in line:
        mull_start = ind + 8
    if 'Writing Mulliken decomposition to disk' in line:
        mull_end = ind - 4
        break

# Write Mulliken charges to output file (in original atom order)
of = open('mulliken_partial_charges.dat', 'w')
of.write('# Atom number, partial charge (positive = losing electrons)\n')
for i in np.arange(mull_end - mull_start + 1) + mull_start:
    line = a[i].split()
    of.write(str(line[1]) + ' ' + str(line[3]) + '\n')
of.close()

# Write charges to 2D matrix
a_charges = np.empty((N_species, 9999), dtype='float') # Arbitrary large size
a_charges.fill(9999) # Fill unused matrix elements with '9999'
for i, line in enumerate(a_numbers): # Loop all atom numbers
    for j, n in enumerate(line):
        if n != 9999:
            # Loop mulliken analysis lines
            for k in np.arange(mull_end - mull_start + 1) + mull_start:
                line = a[k].split()
                if int(line[1]) == n: # If correct atom number
                    a_charges[i, j] = float(line[3])
                    break

# Sum charges
tot_charges = np.array([], dtype='float')
for line in a_charges:
    charge = 0.0
    for c in line:
        if c != 9999:
            charge = charge + c
    tot_charges = np.append(tot_charges, charge)

# Print charges
print()
print('# Total charges')
for ind, s in enumerate(species):
    print(s + ': ' + '{:.6f}'.format(tot_charges[ind]))

print()
print('# Average charge per atom')
for ind, s in enumerate(species):
    avg_charge = tot_charges[ind] / N_atoms[ind]
    print(s + ': ' + '{:.6f}'.format(avg_charge))

# Info about output to file
print('\n# Mulliken charges written to file "mulliken_partial_charges.dat"')

