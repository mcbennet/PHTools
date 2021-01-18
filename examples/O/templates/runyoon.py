#!/usr/bin/env python

import os
from shutil import copyfile
from numpy import array, square
from libpp import *
import socket


print('\nBEGIN RUNYOON.PY')
print(' ','CWD:',os.getcwd())

os.makedirs('workdir',exist_ok=True)

# Directory of reference yoon run
refdir='/home/mcbennett/Repositories/GitHub/mcbennet/L2/atoms/Co/ecps/ccECP_soft/'

# Print DONLP2 parameters
params=[x.split(' ')[0] for x in open('params.in').readlines()]
print('  PARAMS.IN')
print(' ','\n  '.join(params))
params=array(params,dtype=float)
print(params)

# Copy params.in to workdir
copyfile('params.in','workdir/params.in')

f = open('tmp.pp.d','r')
ftext = f.read()
f.close()
for i,v in enumerate(params):
    ftext = ftext.replace('v{:03d}'.format(i),str(v))        
#end for
f = open('workdir/minimal','w')
f.write(ftext)
f.close()

if 'ZXE' in ftext:
    replace_zxe('./workdir/minimal')
#end if

transform_yoon('workdir/minimal',1,'workdir/pp.d')

# Copy params.in to workdir
copyfile('templates/tmp.ip.d','workdir/ip.d')

# Enter workdir and run yoon
os.chdir('workdir')
if socket.gethostname() == 'hilbert':
    yooncomm = os.system("${HOME}/Programs/src/yoon/hfmesh &> hfmesh.log")
else:
    yooncomm = os.system("${HOME}/Codes/src/yoon/hfmesh &> hfmesh.log")
#end if


# EIGENVALUE RESIDUALS
eigs=get_yoon_eigs('out')
ref_eigs=get_yoon_eigs(refdir+'out')

# BOUNDEDNESS VIOLATION
unb = unboundedness('pp.d') 

# NORM CONSERVATION
norm_diffs = get_yoon_norm_diff('fort.26',refdir+'fort.26',[0.8,1.6,0.8,0.8])


os.chdir('..')
f = open('results.out','w')
f.write('{}\n'.format(sum(square(eigs-ref_eigs))))
f.write('{}\n'.format(unb))
f.write('{}\n'.format(sum(square(norm_diffs))))
f.close()

