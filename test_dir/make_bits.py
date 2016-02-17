#!/usr/bin/env python
# Takes a file named on the command line, and writes a Path fingerprint as a
# bitstring to stdout.

import sys
from openeye.oechem import *
from openeye.oegraphsim import *

ims = oemolistream( sys.argv[1] )
mol = OEMol()

for mol in ims.GetOEMols() :
    fp = OEFingerPrint()
    OEMakePathFP( fp , mol )
    bitstring = ''
    for b in range(0, fp.GetSize()):
        if fp.IsBitOn(b):
            bitstring += '1'
        else:
            bitstring += '0'

    print '%s %s' % ( mol.GetTitle() , bitstring )
    
