//
// file cad.cc
// David Cosgrove
// AstraZeneca
// 4th February 2010
//
// Takes an existing set of clusters, with fingerprints, and computes the average
// tanimoto distance in each cluster.

#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/tuple/tuple.hpp>

#include "CadSettings.H"
#include "FileExceptions.H"
#include "FingerprintBase.H"
#include "HashedFingerprint.H"

using namespace boost;
using namespace std;
using namespace DAC_FINGERPRINTS;

extern string BUILD_TIME;

// ****************************************************************************
void read_samples_file( istream &is ,
                        vector<vector<string> > &clusters ) {
  string input_line;
  getline( is , input_line ); // first line is headings
  while( 1 ) {
    getline( is , input_line );
    if( is.eof() )
      break;
    istringstream iss( input_line );
    string tmp;
    iss >> tmp >> tmp >> tmp >> tmp;
    clusters.push_back( vector<string>( istream_iterator<string>( iss ) ,
                                        istream_iterator<string>() ) );
  }

}

// ****************************************************************************
void read_csv_file( istream &is ,
                    vector<vector<string> > &clusters ) {

  string next_line;
  while( 1 ) {
    getline( is , next_line );
    if( is.eof() )
      break;
    vector<string> splits;
    boost::algorithm::split( splits , next_line , boost::algorithm::is_any_of( "," ) );
    unsigned int clus_num = boost::lexical_cast<int>( splits.front() );
    while( clus_num > clusters.size() )
      clusters.push_back( vector<string>() );
    --clus_num; // file counts from 1
    clusters[clus_num].push_back( splits[3] );
  }

}

// ****************************************************************************
void read_cluster_file( const string &cluster_input_file ,
                        CLUS_FILE_FORMAT clus_format ,
                        vector<vector<string> > &clusters ) {

  ifstream ifs( cluster_input_file.c_str() );
  if( !ifs.good() ) {
    cerr << "Error reading " << cluster_input_file << " for reading."
	 << endl;
    exit( 1 );
  }

  if( SAMPLES_FORMAT == clus_format ) {
    read_samples_file( ifs , clusters );
  } else if( CSV_FORMAT == clus_format ) {
    read_csv_file( ifs , clusters );
  }

}

// ****************************************************************************
FingerprintBase *find_fingerprint( vector<FingerprintBase *> &fps ,
                                   const string &fp_name ) {

  // equal_range needs a pointer to an actual fingerprint, and FingerprintBase is
  // abstract.
  boost::scoped_ptr<FingerprintBase> search_fp( new HashedFingerprint( fp_name ) );
  pair<vector<FingerprintBase *>::iterator,vector<FingerprintBase *>::iterator > seed_fp =
    equal_range( fps.begin() , fps.end() , search_fp ,
                 bind( greater<string>() ,
                       bind( &FingerprintBase::get_name , _1 ) ,
                       bind( &FingerprintBase::get_name , _2 ) ) );

  if( seed_fp.first == seed_fp.second ) {
    cerr << "Program cad error : fingerprint for member " << fp_name << " not found."
        << endl
        << "Program aborts with error." << endl;
    exit( 1 );
  }

  if( fp_name == (*seed_fp.first)->get_name() )
    return *seed_fp.first;
  else
    return *seed_fp.second;

}

// ****************************************************************************
void get_cluster_fps( vector<FingerprintBase *> &cluster_fps ,
                      const vector<string> &cluster ,
                      vector<FingerprintBase *> &clus_fps ) {

  clus_fps.clear();
  for( int i = 0 , is = cluster.size() ; i < is ; ++i ) {
    clus_fps.push_back( find_fingerprint( cluster_fps , cluster[i] ) );
  }

}

// ****************************************************************************
void generate_cads( vector<FingerprintBase *> &cluster_fps ,
                    const vector<vector<string> > &clusters ,
                    vector<tuple<double,double,double> > &cads ) {

  vector<FingerprintBase *> curr_clus;
  for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
    if( clusters[i].size() > 1 ) {
      get_cluster_fps( cluster_fps , clusters[i] , curr_clus );
      int num_dists = 0;
      double sum_dist = 0.0 , min_dist = 1.0 , max_dist = 0.0;
      for( int j = 0 , js = curr_clus.size() - 1 ; j < js ; ++j ) {
        for( int k = j + 1 , ks = js + 1 ; k < ks ; ++k , ++num_dists ) {
          double dist = curr_clus[j]->calc_distance( *curr_clus[k] );
          sum_dist += dist;
          if( dist > max_dist ) {
            max_dist = dist;
          }
          if( dist < min_dist ) {
            min_dist = dist;
          }
        }
      }
      cads.push_back( make_tuple( sum_dist / double( num_dists ) , min_dist , max_dist ) );
    } else {
      cads.push_back( 0.0 );
    }
  }

}

// ****************************************************************************
void write_cads( const vector<tuple<double,double,double> > &cads ,
                 const vector<vector<string> > &clusters ,
                 const string &filename ) {

  ofstream ofs( filename.c_str() );
  if( !ofs || !ofs.good() ) {
    throw DACLIB::FileReadOpenError( filename.c_str() );
  }

  for( int i = 0 , is = cads.size() ; i < is ; ++i ) {
    ofs << i << " " << clusters[i].front() << " " << cads[i].get<0>()
        << " " << cads[i].get<1>() << " " << cads[i].get<2>()
        << " " << clusters[i].size() << endl;
  }

}

// ****************************************************************************
int main( int argc , char **argv ) {

  cout << "cad - built " << BUILD_TIME << endl;

  CadSettings cs( argc , argv );
  if( !cs ) {
    cout << cs.error_message() << endl << cs.usage_text() << endl;
    cerr << cs.error_message() << endl << cs.usage_text() << endl;
  }

  vector<vector<string> > clusters;

  read_cluster_file( cs.cluster_file() , cs.clus_file_format() , clusters );

  vector<FingerprintBase *> cluster_fps;
  try {
    read_fp_file( cs.cluster_fp_file() , cs.input_format() ,
                  cs.bitstring_separator() , cluster_fps );
  } catch( DACLIB::FileReadOpenError &e ) {
    cerr << e.what() << endl;
    cout << e.what() << endl;
    exit( 1 );
  } catch( FingerprintFileError &e ) {
    cerr << e.what() << endl;
    cout << e.what() << endl;
    exit( 1 );
  }
  sort( cluster_fps.begin() , cluster_fps.end() ,
        bind( &FingerprintBase::get_name , _1 ) >
        bind( &FingerprintBase::get_name , _2 ) );

  vector<tuple<double,double,double> > cads; // mean, min, max
  generate_cads( cluster_fps , clusters , cads );

  write_cads( cads , clusters , cs.output_file() );

}
