#!/usr/bin/env python3

from pseudopotential import GaussianPP,gamessPPFile
import sys, getopt

inputfile = ''
outputfile=''
try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["ifile=","ofile="])
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
#end for

print('\nInput file is: ', inputfile)

filepath=inputfile
ppf = gamessPPFile(filepath)
pp = GaussianPP()
pp.read_text(ppf.pp_text,format='gamess')
pp.write_qmcpack(outputfile)
