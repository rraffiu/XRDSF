"""
This program reads lattice vectors and positions
of the atoms (basis) and returns structure factor
for any general system.
"""

import sys
from ase.io import read
from numpy import *
from operator import itemgetter
from matplotlib import pyplot as plt

"""
How to run:

python xrdfs.py inputfile format

    inputfile: contains lattice vectors and basis (atomic positions) in some standard format.
    format:    Format of the input file. For example 'vasp' for CONTCAR.
"""

file_in = sys.argv[1]
format_in = sys.argv[2]

if format_in == 'vasp':
    obj = read(file_in,format=format_in)
else:
    exit('ERROR:The program only reads vasp format for now.')

# Get the reciprocal cell/vectors
rcell    = 2*pi*linalg.inv(obj.cell).T
amp_rvec = linalg.norm(rcell,axis=1)
Gcut = 10   # G_max in [A-1]
natoms = obj.get_number_of_atoms()
hklmax = [int(x) for x in Gcut/amp_rvec]
nmax = (2*hklmax[0]+1)*(2*hklmax[1]+1)*(2*hklmax[2]+1)
ints = empty([nmax,5])
indx = 0
for h in  range(-hklmax[0],hklmax[0]+1):
    for k in range(-hklmax[1],hklmax[1]+1):
        for l in range(-hklmax[2],hklmax[2]+1):
            G = h*rcell[0]+k*rcell[1]+l*rcell[2]
            strfac = 0.0
            for d in range(natoms):
                proj = dot(G,obj.positions[d])
                strfac += exp(proj*1j)
            ints[indx,:] = [h,k,l,linalg.norm(G), abs(strfac)**2]
            indx +=1

idx = argsort(ints[:,3])
ints = ints[idx,:]
gdat, ind, cnt = unique(ints[:,3].round(decimals=2),return_index=True,return_counts=True)

def pdf(x, mu, sig):
    return exp(-power((x - mu)/sig, 2.)/2)

q = linspace(0,10,200)
sf = 0.0*q

final = empty([len(gdat),5])
for i in range(len(gdat)):
    final[i,0:4] = ints[ind[i],0:4]
    psum = 0.0
    for j in range(ind[i],ind[i]+cnt[i]):
        psum += ints[j,4]
        sf += pdf(q,gdat[i],0.05)*ints[j,4]
    final[i,4] = psum

final[:,0:3] = -sort(-abs(final[:,0:3]))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(q,sf)
for i, p in enumerate(final[:,4]):
    if p > 4.1:
        lbl = '['+str(final[i,0].astype(int))+str(final[i,1].astype(int))+str(final[i,2].astype(int))+']'
        plt.annotate(lbl, (final[i,3], final[i,4]+10),rotation=90)
plt.xlabel(r"q (\AA$^{-1}$)")
plt.ylabel('Intensity (arb. units)')
plt.margins(0.05,0.1)
plt.show(block=False)
