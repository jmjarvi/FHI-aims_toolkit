#!/usr/bin/python3
# Get energy data (corrected energy) and the number of relaxation steps 
# from multiple FHI-aims output files (names must be 'output').
# Input file should be a list of directory names (must be numbers only)
# Writes output file to 'energies_from_aims_output.dat'

import numpy as np
import sys

# Read subdirectory numbers from file
dirs = np.loadtxt(sys.argv[1], dtype='int')
N = np.size(dirs)

# Initialize data arrays
steps = np.array([], dtype='int') # No of relaxation steps
E_init = np.array([]) # Energy before relaxation
E_final = np.array([]) # Energy after relaxation

for i in np.arange(N):

    cnum = dirs[i] # Calculation number

    # Read output file
    a = [] # Data variable
    f = open(str(cnum) + '/output', 'r')
    for line in f:
        a.append(line)
    f.close()

    # Check if file ends with the Happy Days
    niceday = 0
    for line in a:
        if 'Have a nice day.' in line:
            niceday = 1
    if niceday == 0:
        print('*** Error: Calculation ' + str(cnum) + \
                ' not completed! Check end of output file.')
        sys.exit()

    # Get number of relaxation steps
    for line in a:
        if 'Number of relaxation steps' in line:
            line = line.split()
            steps = np.append(steps, int(line[6]))

    # Get initial and final energies
    energies = np.array([])
    for line in a:
        if '| Total energy corrected' in line:
            line = line.split()
            energies = np.append(energies, float(line[5]))
    E_init = np.append(E_init, energies[0])
    E_final = np.append(E_final, energies[-1])

# Write output to file
of = open('energies_from_aims_output.dat', 'w')
of.write('# Structure no., Relax. steps, E_initial, E_final\n')
for i in np.arange(N):
    of.write(str(dirs[i]) + ' ' + \
            str(steps[i]) + ' ' + \
            str(E_init[i]) + ' ' + \
            str(E_final[i]) + '\n')
of.close()
