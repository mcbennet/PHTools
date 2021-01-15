#!/usr/bin/env python3

from pseudopotential import GaussianPP,gamessPPFile
import sys, getopt

inputfile = ''
outputfile=''
lchan=-1
try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:l:",["ifile=","ofile=","local="])
except getopt.GetoptError:
   print('gms2xml -i <inputppfile>')
   sys.exit(2)
for opt, arg in opts:
   if opt == '-h':
      print('gms2xml -i <inputppfile>')
      sys.exit()
   elif opt in ("-i", "--ifile"):
      inputfile = arg
   elif opt in ("-o", "--ofile"):
      outputfile = arg
   elif opt in ("-l", "--local"):
      lchan = int(arg)
#end for

if lchan<0:
    print('Please provide local channel via \'-l\'')
    quit()
#end if

print('\nInput file is: ', inputfile)

filepath=inputfile
ppf = gamessPPFile(filepath)
pp_ref = GaussianPP()
pp_ref.read_text(ppf.pp_text,format='gamess')


pp = pp_ref.copy()

pp.set_component('s',v=pp_ref.components.d+pp_ref.components.s)








print(pp)
quit()


pp.write_qmcpack(outputfile)
