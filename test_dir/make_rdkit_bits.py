 #!/usr/bin/env python

# Produce RDKit fingerprints in binary file that could be used with Satan
# etc.

from __future__ import print_function
from sys import exit

from rdkit import Chem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit import DataStructs
from rdkit.Chem.rdmolops import RDKFingerprint

import argparse

################################################################################
def parse_arguments():

    parser = argparse.ArgumentParser(description="Generate RDKit fingerprints in ASCII format")
    parser.add_argument('-I', '--input-mol', dest='infile',
                        required=True, help='Input molecule file.')
    parser.add_argument('-O', '--output-file', dest='outfile',
                        required=True, help='Output filename')
    parser.add_argument('-F', '--fingerprint-type', dest='fptype',
                        choices=['Morgan', 'Path'],
                        default='Morgan', help='Fingerprint type (Path or Circular/Morgan). Default=Morgan')
    parser.add_argument('-N', '--num-bits', dest='numbits',
                        default=4096, type=int,
                        help='Number of bits in fingerprint. Default=4096')
    parser.add_argument('-S', '--size', dest='size',
                        default=-1, type=int,
                        help='Max. path length for Path fingerprints (default=7),'
                        ' max. radius for Morgan fingerprints (default=3)')

    args = parser.parse_args()
    if -1 == args.size:
        if 'Morgan' == args.fptype:
            args.size = 3
        elif 'Path' == args.fptype:
            args.size = 7

    return args

################################################################################
def create_mol_supplier(infile):
    """
    Take the input filename and based on extension, return the relevant RDKit
    supplier object.

    Args:
        infile (string):  name of input file

    Returns:
        RDKit mol supplier
    """

    if infile.endswith('.smi') or infile.endswith('.smi.gz'):
        return Chem.SmilesMolSupplier(infile)
    elif infile.endswith('.sdf') or infile.endswith('.sdf.gz'):
        return Chem.SDMolSupplier(infile)

    return None

################################################################################
def main():

    args = parse_arguments()
    print('Read from {}'.format(args.infile))
    print('Write to {}'.format(args.outfile))

    suppl = create_mol_supplier(args.infile)
    if not suppl:
        print('ERROR - unrecognised file format for {}'.format(infile))
        exit(1)

    with open(args.outfile, 'w') as of:
        for mol in suppl:
            if args.fptype == 'Morgan':
                fp = GetMorganFingerprintAsBitVect(mol, args.size, nBits=args.numbits)
            elif args.fptype == 'Path':
                fp = RDKFingerprint(mol, fpSize=args.numbits, maxPath=args.size)
            of.write('{} {}\n'.format(mol.GetProp('_Name'), DataStructs.BitVectToText(fp)))

################################################################################
if __name__ == "__main__":
    main()
    
