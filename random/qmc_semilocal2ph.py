#!/usr/bin/env python3

from pseudopotential import GaussianPP,SemilocalPP,QmcpackPP
import sys, getopt

inputfile = ''
outputfile=''
try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["ifile=","ofile="])
except getopt.GetoptError:
   print('qmc_semilocal2ph.py -i <inputppfile>')
   sys.exit(2)
for opt, arg in opts:
   if opt == '-h':
      print('qmc_semilocal2ph.py -i <inputppfile>')
      sys.exit()
   elif opt in ("-i", "--ifile"):
      inputfile = arg
   elif opt in ("-o", "--ofile"):
      outputfile = arg
#end for

print('\nInput file is: ', inputfile)

filepath=inputfile
#ppf = SemilocalPP(filepath,format='qmcpack')
ppf_ref = QmcpackPP(filepath)


ppf = ppf_ref.copy()


print(ppf)
ppf.add_L2(v=0.5*(ppf_ref.components.p-ppf_ref.components.s))
ppf.set_component('s',v=ppf_ref.components.d+ppf_ref.components.s)
ppf.local='s'
ppf.remove_nonlocal(l='p')
ppf.remove_nonlocal(l='d')
print(ppf)
ppf.write_qmcpack(outputfile)
