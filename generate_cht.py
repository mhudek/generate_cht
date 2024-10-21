#!/usr/bin/env python3

"""
This script creates a chitosan molecule chain cusing coordinates provided
by Naumov & Ignatov Citation: Naumov, V.S. & Ignatov, S.K. Modification of 56ACARBO
force field for molecular dynamic calculations of chitosan and 
its derivatives  // J Mol Model (2017) 23: 244. 

The script needs to be accompanied with .gro files in the same
directorty which contain residue coordinates with C1 atom 
positioned at (0,0,0).

Output file is 'data.pdb' by default.

"""

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
Try to import numpy; if failed, import a local version mynumpy
which needs to be provided
"""
try:
    import numpy as np
except:
    print("numpy not found. Exiting.")
    sys.exit(1)

gro_to_charm = {'CHT' : 'SDN', 'CHT0' : 'SDN', \
'CHTP' : 'SDP', 'CHTN' : 'SDN', 'ACE' : 'BGL'}

def rotate(n,v):
    '''rotate vector, v around vector n by 180 deg'''
    R = np.array([[2*n[0]**2-1, 2*n[0]*n[1], 2*n[0]*n[2]], \
                  [2*n[0]*n[1], 2*n[1]**2-1, 2*n[1]*n[2]], \
                  [2*n[0]*n[2], 2*n[1]*n[2], 2*n[2]**2-1]])
    out = np.dot(R,v)
    return out

def zrotate_ang(v, theta):
    '''rotate vector, v theta degrees around z-axis'''
    theta = np.deg2rad(theta)
    
    R = np.array([[np.cos(theta), -np.sin(theta), 0], 
    \[np.sin(theta),np.cos(theta), 0],[0, 0, 1]])
    v = np.array(v)
    out = np.dot(R,v)
    return out


def zrotate(v):
    '''rotate vector, v 180 degrees around z-axis'''
    R = np.array([[-1, 0, 0], [0,-1, 0],[0, 0, 1]])
    v = np.array(v)
    out = np.dot(R,v)
    return out

def read_monomer(res):
    '''read individual monomers from the working directory'''
    res_input = str(res) + ".gro"
    read_residue = open(res_input,'r')
    n = 0
    x, y, z = [], [], []
    atom_id, atom_name = [], []
    atom_type, charm_name = [], []
    lines = read_residue.readlines()
    for line in lines:
        list1 = line.split()
        if len(list1) != 9:
            continue
        else:
            atom_name.append(list1[2])
            atom_id.append(list1[3])
            x.append(list1[4])
            y.append(list1[5])
            z.append(list1[6])
            atom_type.append(list1[7])
            charm_name.append(list1[8])

            n += 1
    
    return(atom_name, atom_id, x, y, z, atom_type, charm_name, n)

def create_polymer(res_id, res_name, a_name,\
a_id, x, y, z, a_type, charm_n ):
    '''correct atom postions and residue ids to create a polymer''' 
    max_res = int(res_id[-1])
    n = np.zeros((max_res,3))
    #from gromos unit to charm units
    mul = 10.0
    
    for i in range( len(x)):
        if res_id[i] == 1:

            x[i], y[i], z[i] = float(x[i]), float(y[i]), float(z[i])
            x[i] = mul * x[i]
            y[i] = mul * y[i]
            z[i] = mul * z[i]

        else:
            # correct atom_id
            a_id[i] = str(i+1)
            u = np.array([ float(x[i]), float(y[i]), float(z[i])])
            
            #rotate every second monomer by 180 deg around z-axis
            if res_id[i]%2 == 0:
                u = zrotate(u)   
            # translation vector - distance O4 - O1 
            t = np.array([ 0, 0, -0.52178])
            
            #calculate new coordinates
            v = (res_id[i]-1)*t
            v = v + u
            
            x[i], y[i], z[i] = v[0], v[1], v[2]
            x[i] = mul * x[i]
            y[i] = mul * y[i]
            z[i] = mul * z[i]
    
    return(res_id, res_name, a_name, a_id,\
    x, y, z, a_type, charm_n)          
           
def create_crystal(res_id, res_name, a_name, \
a_id, x, y, z, a_type, charm_n):
    """
    Create chain by duplicating created polymer
    """
    # number of chains in x and y directions
    a = 1
    b = 0
    xs, ys, zs = [], [], []
            
  
    res_ids, res_names = [], []
    a_names, a_ids = [], []
    a_types, charm_ns = [], [] 
   
    #tilt of the chain with respect to the z-axis
    theta = 7.5
    dx = 4.85
    dy = 9.26
    i,j = 0, 0
    n = 0
    
   
    #duplicate chains to the number of 
    # chains required for the crystal
    while i < a:
        j = 0
        while j < b:
            
                               
            k = 0
            while k <len (x):
                u = np.zeros(3)
                u[0], u[1], u[2]  = float(x[k]), float(y[k]), float(z[k])
                u = zrotate_ang(u, theta)
                u[0] = u[0] + dx*i
                u[1] = u[1] + dy*j
                
                
                f, d, c = u[0], u[1], u[2] 
                
                xs.append(f)
                ys.append(d)
                zs.append(c)
                
                k+= 1        
            for r_id in res_id:
                res_ids.append(r_id+n*res_id[-1])
            for a_i in a_id:
                a_ids.append( str (int(a_i)+(n*int(a_id[-1]))))
            for a_n in a_name:
                a_names.append(a_n)
            for r_n in res_name:
                res_names.append(r_n)
            for a_t in a_type:
                a_types.append(a_t)
            for c_n in charm_n:
                charm_ns.append(c_n)
           
            n +=1            
            j +=1
        i+=1            
    
    
    return(res_ids, res_names, a_names, \
    a_ids, xs, ys,zs, a_types, charm_ns)            
        
        
        

def read_structure():
    """
    Main function for this script.
    Reads a text file with the following format:
    CHT0 CHT CHTP CHTN
    every line must begin with CHT0 and end with CHTN
    """
    try:
        #infile = open(sys.argv[1],'r')
        infile = open("input")
    except:
        print("Could not open file {:s}. Aborting.".format(infile))
        sys.exit(2)

    #Read in the monemers and see how many polymers
    nmonomeres, nchains = 0, 0
    count = 0
    strandnum = []
    restype = []
    lines = infile.readlines()
    #read file to polymer length and number of chains
    for line in lines:
        line = line.upper()
        if len(line) == 0:
            continue  
        else:
            line = line.split()
            length = len(line)
            nmonomeres += length
            nchains += 1
                  
    # rewind the sequence input file
    infile.seek(0)

    #generate the data file in GROMACS format
    try:
        #out = open("data.gro", "w+")
        out = open("data.pdb", "w+")
    except:
        print("Could not open data file for writing. Aborting.",\
                file=sys.stderr)
        sys.exit(2)

    lines = infile.readlines()
    nlines = len(lines)
         
    chainnum = 0
    
    xs, ys, zs = [], [], []
    counter = 0
    res_id, res_name  = [], []
    a_ids, els, a_names = [], [], []
    a_types, charm_ns = [], [] 
    for line in lines:
        line = line.upper()

        # skip empty lines
        if len(line) == 0:
            continue

        else:
            res = line.split()
            length = len(res)
            print("Found " + str(length) + " monomers." )
            
            #loop over residues
            for b in range(length):
                #print("b = " + str(b))
                restype.append(res[b])
                strandnum.append(chainnum)
                a_name, a_id, x, y, z, a_type, \
                charm_n, n = read_monomer(res[b])
                #print( "n = " + str(n))
                #loop over atoms in residues
                for c in range(n):
                    
                    res_id.append(b+1)
                    res_name.append(res[b])
                    xs.append(x[c])
                    ys.append(y[c])
                    zs.append(z[c])
                    a_ids.append(a_id[c])
                    a_names.append(a_name[c])
                    a_types.append(a_type[c])
                    charm_ns.append(charm_n[c])
                count += n

    size = len(xs)
    
    #create polymer
    res_id, res_name, a_name, a_id, x, y, z, a_type, charm_n = \
            create_polymer(res_id, res_name, a_names, \
            a_ids, xs, ys, zs, a_types, charm_ns)
    
    
    res_id, res_name, a_name, a_id,\
    x, y, z, a_type, charm_n = \
        create_crystal(res_id, res_name, a_name, \
        a_id, x, y, z, a_type, charm_n)
    
    #test if nested arrays
    print("a_id " + str(a_id[0]) + " type" + str(type(a_id[0]))\
    + "expect str")
    print("charm_n: " + str(charm_n[0]) + "type: " \
    + str(type(charm_n[0])) + "exp str")
    print("fake ACE")
    print("res_name: " + str (res_name[0]) + " type: " \
    + str(type(res_name[0])))
    print("res_id: " + str(res_id[0]) + " type: " \
    + str(type(res_id[0])) + " exp int")
    
    print("a_name: " + str(a_name[0]) + "t ype: " \
    + str(type(a_name[0])))
    
    print("x " + str(x[0]) + "type" + str(type(x[0]))) 
    print("y " + str(y[0]))
    print("z " + str(z[0]))
    print("a_type: " + str(a_type[0]))
             
    
    count = count*18
    #sanity check
    if size != count:
        print('Warning! Numbers do not match')
        print("Size should be " + str(size) + ", but I count " + str(count))

    """
    gromacs output format
    out.write('!comment \n')
    out.write( str(count) + '\n')
    
    for j in range(count):
        # create .gro output format
        # res_num, res_name, atom_name, atom_num, x, y, z 
        s = '{:>5d} {:>4s} {:>4s} {:>4s} {:>7.3f} {:>7.3f}\
        {:>7.3f}'.format(res_id[j], res_name[j], \
              a_name[j],a_id[j], x[j], y[j], z[j])
        out.write(s)
        out.write("\n")
    
    #dummy box size - later corrected with gromacs preprocessing
    out.write("10 10 10 \n")
    out.close()
    """

    #pdb output format
    out.write("REMARK generated using generate_cht.py script \n")
    s1 = 'ATOM   '
    
    for j in range(count):
        res_name[j] = gro_to_charm[res_name[j]]
        s2 = '{:>4s}  {:<4s}{:>3s} {:1s} {:>3d}    \
        {:>8.3f}{:>8.3f}{:>8.3f}{:>21s}{:>3s}'.\
            format(a_id[j], charm_n[j],res_name[j], 'A',res_id[j],x[j],y[j],z[j],'CHT',a_type[j])
        s = s1 + s2 + '\n'
        out.write(s)

    out.write('END\n')
    out.close()
            
# call the above main() function, which executes the program
read_structure()

end_time=timer()
runtime = end_time-start_time
hours = runtime/3600
minutes = (runtime-np.rint(hours)*3600)/60
seconds = (runtime-np.rint(hours)*3600-np.rint(minutes)*60)%60
print("## Total runtime %ih:%im:%.2fs" % \
(hours,minutes,seconds), file=sys.stdout)
