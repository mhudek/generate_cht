#!/usr/bin/env python3

"""
This scripts reads pdb file and creates plumed script for dihedral angles in chitosan
psi = C1-O4-C4'-C3'
phi = O5-C1-O4 -C4'
"""

__author__ = "Magdalena Hudek"
__version__ = "1.0"
__email__ = "magdalena.hudek@strath.ac.uk"

"""
Import basic modules
"""
import sys 
import os
import io
import timeit
from timeit import default_timer as timer

start_time = timer()

"""
Try to import numpy
"""
try:
    import numpy as np
except:
    print("Numpy not found. Exiting.")
    sys.exit(1)

if len(sys.argv) != 2:
    print("Usage: $python plumed_input_gen.py <pdb_filename>")
    sys.exit(2)

def read_pdb():
    """
    Scans pdb file for psi and phi atom ids, and return psi and phi angles in terms of their atoms ids 
    """
    try:
        infile = open(sys.argv[1],'r')
    except:
        print("Could not open file {:s}. Aborting.".format(infile))
        sys.exit(2)
    
    o4s, c1s, o5s, c4s, c3s = [], [], [], [], [] 
    lines = infile.readlines()
    for line in lines:
        list1 = line.split()
        if list1[0] != 'ATOM':
            continue
        elif list1[2] == 'O5':
            o4s.append(list1[1])
        elif list1[2] == 'C5':
            c1s.append(list1[1])
        elif list1[2] == 'O6':
            o5s.append(list1[1])
        elif list1[2] == 'C2':
            c4s.append(list1[1])
        elif list1[2] == 'C3':
            c3s.append(list1[1])
            
    #sanity check

    if len(o4s) != len(c3s):
        print("something went wrong ...")
    
    n = len(o5s)
    i = 0
    phi = np.zeros((n,4))
    psi = np.zeros((n,4))
    while i < n-1:
        
        phi[i][0] = o5s[i]
        phi[i][1] = c1s[i]
        phi[i][2] = o4s[i]
        phi[i][3] = c4s[i+1]

        psi[i][0] = c1s[i]
        psi[i][1] = o4s[i]
        psi[i][2] = c4s[i+1]
        psi[i][3] = c3s[i+1]

        i += 1

    return(phi, psi, n-1)

def write_plumed_input(phi,psi, n):
    """
    write input for PLUMED-GUI in VMD, save the angles in tabular format
    """
    try:
        out = open('phi_psi.plumed', 'w+')
    except:
        print("Could not create output file. Aborting.")
        sys.exit(2)
    
    s_psi_n = ''
    s_phi_n = ''
    for i in range(n):
        phi_n = 'phi' + str(i)
        psi_n = 'psi' + str(i)
        if i == 0:
            s_phi_n = phi_n
            s_psi_n = psi_n
        else:
            s_phi_n = s_phi_n + ',' + phi_n
            s_psi_n = s_psi_n + ',' + psi_n
        s1 = '{}: TORSION ATOMS={:d},{:d},{:d},{:d}\n'.format(phi_n,int(phi[i][0]), int(phi[i][1]), int(phi[i][2]), int(phi[i][3]))
        s2 = '{}: TORSION ATOMS={:d},{:d},{:d},{:d}\n'.format(psi_n,int(psi[i][0]), int(psi[i][1]), int(psi[i][2]), int(psi[i][3]))
        out.write(s1)
        out.write(s2)

    out.write('\n')    
    sp1 = 'PRINT ARG='+ s_phi_n  + ' STRIDE=10 FILE=COLVAR_PHI\n'
    sp2 = 'PRINT ARG='+ s_psi_n  + ' STRIDE=10 FILE=COLVAR_PSI\n'
    out.write(sp1)
    out.write(sp2)

    out.close()     
  

#Main part of this script
#Execute the above functions
phi, psi, n = read_pdb()
write_plumed_input(phi, psi, n)

end_time=timer()
runtime = end_time-start_time
hours = runtime/3600
minutes = (runtime-np.rint(hours)*3600)/60
seconds = (runtime-np.rint(hours)*3600-np.rint(minutes)*60)%60
print("## Total runtime %ih:%im:%.2fs" % (hours,minutes,seconds), file=sys.stdout)

