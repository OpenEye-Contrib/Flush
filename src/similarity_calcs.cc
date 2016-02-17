//
// file similarity_calcs.cc
//
// D Cosgrove
// Zeneca Pharms
// 22 July 1997
//
// This file contains the functions that are used to calculate similarity
// between two fingerprints.  They will be used by the object FlushDataset,
// the particular one to use being established at runtime.

#include <iostream.h>
#include <iomanip.h>
#include "bit_count.h"

// *************************************************************************
// calculate the dissimilarity between the two fingerprints using
// the Tanimoto measure
float calc_tanimoto( unsigned short *pus_a , unsigned short *pus_b ,
		     int i_num_shorts , int i_num_a_bits , int i_num_b_bits ,
		     float r_alpha , float r_beta , float r_theta ) {

  double         r_div , r_ret;

  int            i_num_in_common , j;

  unsigned short          us_common;

  i_num_in_common = 0;
  for( j = 0 ; j < i_num_shorts ; j++ ) {
    us_common = pus_a[j] & pus_b[j];
    i_num_in_common += BIT_COUNT_ARRAY[us_common];
  }
  // the distance is a tanimoto dissimilarity - 0.0 is identical,
  // 1.0 is completely dissimilar
  r_div = double( i_num_a_bits ) + double( i_num_b_bits ) -
    double( i_num_in_common );
  r_ret = 1.0 - ( double( i_num_in_common ) / r_div );
  return( float( r_ret ) );

}

// **************************************************************************
// calculate the dissimilarity between the two fingerprints using
// the Tversky measure
float calc_tversky( unsigned short *pus_a , unsigned short *pus_b ,
		    int i_num_shorts , int i_num_a_bits , int i_num_b_bits ,
		    float r_alpha , float r_beta , float r_theta ) {

  double         r_div , r_ret;

  int            i_num_in_a_and_not_b , i_num_in_b_and_not_a ,
                 i_num_in_common , j;

  unsigned short i_a_and_not_b , i_b_and_not_a , i_common;

  i_num_in_a_and_not_b = i_num_in_b_and_not_a = i_num_in_common = 0;
  for( j = 0 ; j < i_num_shorts ; j++ ) {
    i_common = pus_a[j] & pus_b[j];
    i_num_in_common += BIT_COUNT_ARRAY[i_common];
    
    i_a_and_not_b = pus_a[j] &~ pus_b[j];
    i_num_in_a_and_not_b += BIT_COUNT_ARRAY[i_a_and_not_b];
    
    i_b_and_not_a = pus_b[j] &~ pus_a[j];
    i_num_in_b_and_not_a += BIT_COUNT_ARRAY[i_b_and_not_a];
  }
  // the distance is a dissimilarity - 0.0 is identical,
  // 1.0 is completely dissimilar
  r_div = double( i_num_in_a_and_not_b ) * r_alpha +
    double( i_num_in_b_and_not_a ) * r_beta +
    double( i_num_in_common ) * r_theta;
  r_ret = 1.0 - double( i_num_in_common ) / r_div;
  return( float( r_ret ) );
  
}

