#!/usr/bin/env python3

import os
import subprocess
from shutil import copyfile
from numpy import array, square, array2string
from libpp import *
import socket

try:
    opts, args = getopt.getopt(sys.argv[1:],"hw:r:t:d:c:v:",["weights=","rcuts=","trans=","refdir=","cc_dir=","cc_wght="])
except getopt.GetoptError:
    print('runyoon.py -w <weights> -r <rcuts> -t <transform> -d <refdir>')
    sys.exit(2)
cc_dir=None
cc_wght=None
for opt, arg in opts:
    if opt == '-h':
        print('runyoon.py -w <weights> -r <rcuts> -t <transform> -d <refdir>')
        sys.exit()
    elif opt in ("-w", "--weights"):
        weights = arg.split()
    elif opt in ("-r", "--rcuts"):
        rcuts = arg.split()
    elif opt in ("-t", "--trans"):
        trans = int(arg)
    elif opt in ("-d", "--refdir"):
        refdir = arg
    elif opt in ("-c", "--cc_dir"):
        cc_dir = arg
    elif opt in ("-v", "--cc_wght"):
        cc_wght = float(arg)

weights = array(weights,dtype=float)
rcuts = array(rcuts,dtype=float)

print('\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
print('\nBEGIN RUNYOON.PY')
print('\n ','CWD:',os.getcwd())

os.makedirs('workdir',exist_ok=True)

## Directory of reference yoon run
#if socket.gethostname() == 'hilbert':
#    refdir='/home/mcbennett/Repositories/GitHub/mcbennet/L2/atoms/Co/ecps/ccECP_soft/'
#else:
#    refdir='/home/8gb/Projects/L2_Production/atoms/Co/ecps/ccECP_soft/'
##0end if

# Print DONLP2 parameters
params=[x.split(' ')[0] for x in open('params.in').readlines()]
print('\n  PARAMS.IN')
print(' ','\n  '.join(params))
params=array(params,dtype=float)

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

transform_yoon('workdir/minimal',trans,'workdir/pp.d')

# Copy params.in to workdir
copyfile('templates/tmp.ip.d','workdir/ip.d')

# Enter workdir and run yoon
os.chdir('workdir')
if socket.gethostname() == 'hilbert':
    with open('hfmesh.log', 'w') as f, open('hfmesh.err', 'w') as ferr:
        subprocess.run(['/home/mcbennett/Programs/src/yoon/hfmesh'],stdout=f,stderr=ferr)
    #end with
elif socket.gethostname().startswith('cori') or socket.gethostname().startswith('nid'):
    with open('hfmesh.log', 'w') as f, open('hfmesh.err', 'w') as ferr:
        os.system('module swap PrgEnv-intel PrgEnv-gnu; /global/homes/m/mcbennet/Codes/src/yoon/hfmesh')
        #subprocess.run(['/global/homes/m/mcbennet/Codes/src/yoon/hfmesh'],stdout=f,stderr=ferr)
    #end with
else:
    with open('hfmesh.log', 'w') as f, open('hfmesh.err', 'w') as ferr:
        subprocess.run(['/home/8gb/Codes/src/yoon/hfmesh'],stdout=f,stderr=ferr)
    #end with
#end if


# EIGENVALUE RESIDUALS
eigs=get_yoon_eigs('out',valence=True)
ref_eigs=get_yoon_eigs(refdir+'out',valence=True)
print(ref_eigs)

# BOUNDEDNESS VIOLATION
unb = unboundedness('pp.d') 

# NORM CONSERVATION
norm_diffs = get_yoon_norm_diff('fort.26',refdir+'fort.26',rcuts)

# Coupled Cluster Energetics
if cc_dir is not None:

    if cc_wght is None:
        print('Weight for cc energies must be given via --cc_wght')
        sys.exit()
    #end if
    # Assume 1st two CCSD(T) IPs and EA are contained in file [cc_dir]/ips_and_ea.txt
    with open(cc_dir+'/ips_and_ea.txt', 'r') as fcc:
        ref_cc = fcc.read().split()
    #end with
    ref_cc = np.array(ref_cc,dtype=float)

    # Shifts
    ref_cc = ref_cc+[0.003,0.0,0.0]

    # Filename of run
    #file_base = '5z_ccecp'
    file_base = '5z_ccecp_HF'
    file_name = '{}.com'.format(file_base)

    # Truncation level of PP
    maxlval = 6

    # Filename of basis
    basis_name = 'Co.ccecp.aug-cc-pV5Z.molpro'

    # Transform minimal basis
    transform_molpro('minimal',1,zae=27,outfile='pp.molpro',maxl=maxlval)

    copyfile('../templates/{}'.format(file_name),file_name)
    copyfile('../templates/{}'.format(basis_name),basis_name)
    copyfile('../templates/states.proc','states.proc')
    quit()
   
    print('Running MOLPRO')
    with open('5z_ccecp.py.log', 'w') as f, open('5z_ccecp.err', 'w') as ferr:
        os.system('module swap PrgEnv-gnu PrgEnv-intel; module load molpro; module list; molpro -n 24 {}'.format(file_name))
        os.system('module list')
    #end with

    with open('{}.table1.txt'.format(file_base.lower())) as tbf:
        cc_text = tbf.read()
        print(cc_text)
    #end with
    scf_energies = cc_text.split()[1:]
    scf_energies = np.array(scf_energies,dtype=float)
    scf_ips = []
    for i in range(1,len(scf_energies)):
        scf_ips.append(scf_energies[i]-scf_energies[i-1])
    #end for
    scf_ips = np.array(scf_ips)
    print('IPs Residuals')
    print(scf_ips-ref_cc)
    ip_residuals_soq = sum(square((scf_ips-ref_cc)/ref_cc))


#end if


# Square residuals
eig_sosq = sum(square((eigs-ref_eigs)/ref_eigs))
norm_sosq = sum(square(norm_diffs))

os.chdir('..')
f = open('results.out','w')
f.write('{}\n'.format(weights[0]*eig_sosq))
f.write('{}\n'.format(weights[1]*unb))
f.write('{}\n'.format(weights[2]*norm_sosq))
if cc_dir is not None:
    f.write('{}\n'.format(cc_wght*ip_residuals_soq))
#end if
f.close()

print('\n  RESIDUALS')
print('  EIGENVALUES: {}'.format(     ' '.join('{:.5f}'.format(n) for n in (eigs-ref_eigs)/ref_eigs)       ))
print('  UNBOUND:     {:.5f}'.format(unb))
print('  NORM:        {}'.format(     ' '.join('{:.5f}'.format(n) for n in norm_diffs)          ))
if cc_dir is not None:
    print('  CC:          {}'.format(     ' '.join('{:.5f}'.format(n) for n in (scf_ips-ref_cc)/ref_cc)       ))
#end if

print('\n  WEIGHTED SQUARES')
print('  EIGENVALUES: {:.8f}'.format(weights[0]*eig_sosq))
print('  UNBOUND:     {:.8f}'.format(weights[1]*unb))
print('  NORM:        {:.8f}'.format(weights[2]*norm_sosq))
if cc_dir is not None:
    print('  CC:          {:.8f}'.format(cc_wght*ip_residuals_soq))
#end if

if cc_dir is None:
    print('\n  SUM of WEIGHTED SQUARES:  {}'.format(weights[0]*eig_sosq+weights[1]*unb+weights[2]*norm_sosq))
else:
    print('\n  SUM of WEIGHTED SQUARES:  {}'.format(weights[0]*eig_sosq+weights[1]*unb+weights[2]*norm_sosq+cc_wght*ip_residuals_soq))
#end if

os.system('rm -r workdir')

print()
