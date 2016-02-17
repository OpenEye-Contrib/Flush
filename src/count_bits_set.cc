// take a flush fp file and write out the number of bits set in each cpd.
// takes 1 command line argument, the name of the fp file.

#include <fstream>
#include <iostream>
#include <stdio.h>

#include "ByteSwapper.H"
#include "Fingerprint.H"
using namespace std;

#include "AbstractPoint.H"
int AbstractPoint::next_seq_num = 0; 

// on a bigendian machine, spells Dave.
static const int MAGIC_INT = 0x65766144;
// as it appears on a littleendian machine
static const int BUGGERED_MAGIC_INT = 0x44617665;

// ***************************************************************************
int main( int argc , char **argv ) {

  if( argc < 2 ) {
    cerr << " Error : need the name of a fingerprints file" << endl;
    exit( 1 );
  }

  FILE *infile = fopen( argv[1] , "rb" );
  if( !infile ) {
    cerr << "Failed to open " << argv[1] << " for reading" << endl;
    return false;
  }
  // take the integer off the top to get things set up correctly.
  // in a new fp file, the very first integer will be either MAGIC_INT
  // or BUGGERED_MAGIC_INT and indicates whether the machine reading and the
  // machine writing were both in the same bigendian/littleendian format.
  // if the first integer is neither of these, indicates that the file is in
  // the old format - we'll assume no byteswapping.
  bool byte_swapping;
  int i , num_chars;
  fread( &i , sizeof( int ) , 1 , infile );
  if( i == MAGIC_INT ) {
    fread( &num_chars , sizeof( int ) , 1 , infile );
    byte_swapping = false;
  } else if( i == BUGGERED_MAGIC_INT ) {
    fread( &num_chars , sizeof( int ) , 1 , infile );
    dac_byte_swapper<int>( num_chars );
    byte_swapping = true;
  } else {
    byte_swapping = false;
    num_chars = i;
  }
  fread( &i , sizeof( int ) , 1 , infile );
  
  int max_mol_name_len = 0 , len;
    
  char *mol_name = 0;
    
  unsigned char *finger_chars = 0;

  while( 1 ) {

    if( !finger_chars )
      finger_chars = new unsigned char[num_chars];

    if( !fread( &len , sizeof( int ) , 1 , infile ) )
      return false;

    if( byte_swapping )
      dac_byte_swapper<int>( len );

    if( len > max_mol_name_len ) {
      delete [] mol_name;
      mol_name = new char[len + 1];
      max_mol_name_len = len + 1;
    }
    if( !fread( mol_name , sizeof( char ) , len + 1 , infile ) )
      break;
    if( !fread( finger_chars , sizeof( unsigned char ) ,
		num_chars , infile ) )
      break;
    Fingerprint finger( mol_name , num_chars , finger_chars );
    cout << mol_name << " " << finger.get_num_bits_set() << endl;
      
  }

}
