# Genrerate chitin/chitosan polymer .pdb (.gro) file 

This script creates a chitosan molecule chain compatible with
GROMOS 56Acarbo force field developed by Naumov & Ignatov
Citation: Naumov, V.S. & Ignatov, S.K. Modification of 56ACARBO
force field for molecular dynamic calculations of chitosan and 
its derivatives  // J Mol Model (2017) 23: 244. 

The script needs to be accompanied with .gro files in the same
directorty which contain residue coordinates with C1 atom 
positioned at (0,0,0). The necessary files for neutral and 
protonated chitosan monomers and chitin monomers are provided. 
This can be expanded with modified monomers if desired. 

Output file is 'data.pdb' by default.

Usage: python generate_cht.py demo_in 
