//
// file histogram.cc
// David Cosgrove
// AstraZeneca
// 27th January 2014
//
// Takes 2 fingerprint files and does a histogram of distances between them.

#include <functional>
#include <numeric>
#include <vector>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "stddefs.H"
#include "FileExceptions.H"
#include "FingerprintBase.H"
#include "HashedFingerprint.H"
#include "NotHashedFingerprint.H"

using namespace boost;
using namespace std;
using namespace DAC_FINGERPRINTS;

// ****************************************************************************
int main( int argc , char **argv ) {

  string usage( "./histogram {FLUSH_FPS|BITSTIRNGS} <PROBE_FILE> <TARGET_FILE> {start_num} {finish_num}" );

  if( argc < 4 ) {
    cout << usage << endl;
    exit( 1 );
  }

  FP_FILE_FORMAT fp_format;
  bool binary_file;
  string bitstring_separator;
  decode_format_string( string( argv[1] ) , fp_format , binary_file , bitstring_separator );

  vector<FingerprintBase *> probe_fps , target_fps;
  read_fp_file( string( argv[2] ) , fp_format , bitstring_separator , probe_fps );
  read_fp_file( string( argv[3] ) , fp_format , bitstring_separator , target_fps );

  cerr << probe_fps.size() << " probe fps and " << target_fps.size() << " target fps." << endl;

  unsigned int start =  0;
  if( argc >= 5 ) {
    start = lexical_cast<unsigned int>( argv[4] );
  }
  unsigned int finish = probe_fps.size();
  if( argc >= 6 ) {
    finish = lexical_cast<unsigned int>( argv[5] );
  }
  if( finish > probe_fps.size() ) {
    finish = probe_fps.size();
  }

  vector<double> hist_fracs( 21 , 0.0 );
  for( unsigned int i = start ; i < finish ; ++i ) {
    vector<unsigned int> dist_counts( 21 , 0 );
    BOOST_FOREACH( FingerprintBase *tfp , target_fps ) {
      double dist = probe_fps[i]->calc_distance( *tfp );
      int i_dist = int( 20.0 * dist );
      ++dist_counts[i_dist];
    }
    for( int j = 0 , js = dist_counts.size() ; j < js ; ++j ) {
      hist_fracs[j] += double( dist_counts[j] ) / double( target_fps.size() );
    }
    cout << i << " : ";
    copy( hist_fracs.begin() , hist_fracs.end() , doubleOut );
    cout << endl;
  }

  double frac_tot = 0.0;
  frac_tot = accumulate( hist_fracs.begin() , hist_fracs.end() , frac_tot );
  cerr << "fraction total : " << frac_tot << endl;

}
