#!/usr/bin/python3
# Sum atom-projected DOS per element from FHI-aims output
# Give input directory as argument (end with /), 
# the sums are written into that directory

import numpy as np
import sys
import os
import re

indir = sys.argv[1] # Input (and output) directory
itot = 0 # Atom counter (total)

# Get element names
elements = np.array([], dtype='str')
for f in os.listdir(indir):
    if 'atom_projected_dos' in f and 'raw' not in f:
        name = f[:-8]
        name = name[19:]
        if name not in elements:
            elements = np.append(elements, name)

for e in elements:

    i = 0 # Atom counter (by element)

    # Get input file names
    infiles = np.array([])
    for f in os.listdir(indir):
        regex = re.compile('_dos_' + str(e) + r'\d{4}\.dat') # Get filename
        infile = regex.findall(f)
        if infile:
            infiles = np.append(infiles, f)

    # Find data array size from first file, initialize
    aa = np.loadtxt(indir + infiles[0])
    data = np.zeros(np.shape(aa))

    # Read and sum atom-projected dos files
    for f in infiles:
        a = np.loadtxt(indir + f)
        data[:,0] = a[:,0] # Energy values
        data[:,1:] = data[:,1:] + a[:,1:] # DOS values
        i = i + 1
    itot = itot + i

    print(str(i) + ' ' + e + ' atoms read')

    # Write to file
    of = open(indir + 'element-projected_dos_' + e + '.dat', 'w')
    of.write('# Energy (eV)  Total DOS    l=0          l=1          l=2 \
              ...\n')
    for line in data:
        for val in line:
            of.write('{:12.8f}'.format(val) + ' ')
        of.write('\n')
    of.close()

print(str(itot) + ' atoms in total')

##### Calculate total dos from output files
##### (for comparison with KS_DOS_total.dat)

j = 0 # Element counter
npoints = np.shape(aa)[0]
sumarray = np.zeros((npoints, 2))
for f in os.listdir(indir):
    if 'element-projected_dos_' in f:
        b = np.loadtxt(indir + f)
        sumarray[:,0] = b[:,0] # Energy values
        sumarray[:,1] = sumarray[:,1] + b[:,1] # DOS values
        j = j + 1

# Write to file
of = open(indir + 'element-projected_dos_total-sum.dat', 'w')
of.write('# Energy (eV)  Total DOS\n')
for line in sumarray:
    for val in line:
        of.write('{:12.8f}'.format(val) + ' ')
    of.write('\n')
of.close()

print(str(j) + ' elements summed in total dos')
