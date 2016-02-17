//
// file amtec.cc
// David Cosgrove
// AstraZeneca
// 14th May 2007
//
// Takes an existing set of clusters, with fingerprints, and drops new
// fingerprints/molecules into those clusters, adding them to the cluster whose
// seed they are nearest. Can overwrite existing cluster file, or write a new one.

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
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/scoped_ptr.hpp>

#include "AmtecSettings.H"
#include "FileExceptions.H"
#include "FingerprintBase.H"
#include "HashedFingerprint.H"
#include "NotHashedFingerprint.H"

using namespace std;
using namespace DAC_FINGERPRINTS;

namespace po = boost::program_options;

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

  boost::scoped_ptr<FingerprintBase> search_fp( new HashedFingerprint( fp_name ) );
  pair<vector<FingerprintBase *>::iterator,vector<FingerprintBase *>::iterator > seed_fp =
      equal_range( fps.begin() , fps.end() , search_fp ,
                   bind( greater<string>() ,
                         bind( &FingerprintBase::get_name , _1 ) ,
                         bind( &FingerprintBase::get_name , _2 ) ) );

  if( seed_fp.first == seed_fp.second )
    return 0;

  if( fp_name == (*seed_fp.first)->get_name() )
    return *seed_fp.first;
  else
    return *seed_fp.second;

}

// ****************************************************************************
int find_nearest_seed( double threshold , vector<FingerprintBase *> &cluster_seeds ,
                       const FingerprintBase &fp ) {

  int nearest_seed = -1;
  double nearest_dist = threshold;
  for( int i = 0 , is = cluster_seeds.size() ; i < is ; ++i ) {
    double dist = cluster_seeds[i]->calc_distance( fp , nearest_dist );
    if( dist < nearest_dist ) {
      nearest_seed = i;
      nearest_dist = dist;
    }
  }

  return nearest_seed;

}

// ****************************************************************************
void add_fps_to_clusters( double threshold ,
                          vector<FingerprintBase *> &cluster_fps ,
                          vector<FingerprintBase *> &new_fps ,
                          vector<FingerprintBase *> &cluster_seed_fps ,
                          vector<vector<string> > &clusters ,
                          vector<int> &additions_dests ) {

  cluster_seed_fps.reserve( clusters.size() );
  for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
    if( clusters[i].empty() ) {
      continue; // Sam sometimes has the clusters out of order and non-consecutive
    }
    cluster_seed_fps.push_back( find_fingerprint( cluster_fps , clusters[i][0] ) );
    if( !cluster_seed_fps.back() ) {
      cerr << "ERROR : Cluster seed " << clusters[i][0]
           << " not found in cluster fingerprints." << endl;
      exit( 1 );
    }
  }

  for( int i = 0 , is = new_fps.size() ; i < is ; ++i ) {
    // find the nearest seed to this fp
    int nearest_seed = find_nearest_seed( threshold , cluster_seed_fps , *new_fps[i] );
    if( -1 == nearest_seed ) {
      cout << new_fps[i]->get_name() << " was beyond " << threshold
           << " from any existing cluster seed." << endl;
      additions_dests.push_back( -1 );
    } else {
      clusters[nearest_seed].push_back( new_fps[i]->get_name() );
      additions_dests.push_back( nearest_seed );
    }
  }

}

// ****************************************************************************
// clusters are in descending order of distance from seed. Any clusters that have
// grown will need this re-establishing.
void re_sort_amended_clusters( const vector<unsigned int> &orig_sizes ,
                               vector<FingerprintBase *> &cluster_seed_fps ,
                               vector<FingerprintBase *> &clus_fps ,
                               vector<FingerprintBase *> &new_fps ,
                               vector<vector<string> > &clusters ) {

  for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
    if( clusters[i].size() == orig_sizes[i] ) {
      continue;
    }
    vector<pair<string,double> > new_clus;
    new_clus.reserve( clusters[i].size() );
    for( int j = 0 , js = clusters[i].size() ; j < js ; ++j ) {
      FingerprintBase *fp = find_fingerprint( clus_fps , clusters[i][j] );
      if( !fp )
        fp = find_fingerprint( new_fps , clusters[i][j] );
      new_clus.push_back( make_pair( clusters[i][j] ,
                                     cluster_seed_fps[i]->calc_distance( *fp ) ) );
    }
    // want the sort to be on distance, with name order as a tie-breaker. The lazy
    // way to do this is:
    sort( new_clus.begin() , new_clus.end() ,
          boost::bind( less<string>() ,
                       boost::bind( &pair<string,double>::first , _1 ) ,
                       boost::bind( &pair<string,double>::first , _2 ) ) );
    stable_sort( new_clus.begin() , new_clus.end() ,
                 boost::bind( less<double>() ,
                              boost::bind( &pair<string,double>::second , _1 ) ,
                              boost::bind( &pair<string,double>::second , _2 ) ) );
  }

}

// ****************************************************************************
void apply_subset_file( const string &subset_file ,
                        vector<FingerprintBase *> &fps ) {

  vector<string> subset_names;
  ifstream ifs( subset_file.c_str() );
  if( !ifs || !ifs.good() ) {
    cerr << "Failed to open " << subset_file << " for reading." << endl;
    exit( 1 );
  }

  copy( istream_iterator<string>( ifs ) , istream_iterator<string>() ,
        back_inserter( subset_names ) );
  cout << "Read " << subset_names.size() << " subset names." << endl;
  sort( subset_names.begin() , subset_names.end() );

  for( int i = 0 , is = fps.size() ; i < is ; ++i ) {
    if( !binary_search( subset_names.begin() , subset_names.end() , fps[i]->get_name() ) ) {
      delete fps[i];
      fps[i] = 0;
    }
  }

  fps.erase( remove( fps.begin() , fps.end() , static_cast<FingerprintBase *>( 0 ) ) ,
             fps.end() );

}

#ifdef NOTYET
// ****************************************************************************
void output_samples_clusters( ostream &os ,
                              const vector<vector<string> > &clusters ,
                              const vector<unsigned int> &orig_sizes ) {

  os << "Molecule name : Cluster size : Cluster Members" << endl;

  for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
    os << clusters[i][0] << " : " << clusters[i].size() << "("
                         << orig_sizes[i] << ") : ";
    for( int j = 0 , js = clusters[i].size() ; j < js ; ++j ) {
      os << clusters[i][j] << "  ";
    }
    os << endl;
  }

}

// ****************************************************************************
void output_csv_clusters( ostream &os ,
                          const vector<vector<string> > &clusters ,
                          const vector<unsigned int> &orig_sizes ) {

  for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
    for( int j = 0 , js = clusters[i].size() ; j < js ; ++j ) {
      os << i + 1 << "," << clusters[i].size() << "," << clusters[i].front()
         << "," << clusters[i][j] << "," << orig_sizes[i] << endl;
    }
  }

}
#endif

// *******************************************************************************
void write_cluster( ostream &os , CLUS_FILE_FORMAT output_format ,
                    const vector<string> &cluster ,
                    int orig_nn_size ) {

  if( SAMPLES_FORMAT == output_format ) {
    os << cluster.front() << " : " << cluster.size() << "("
       << orig_nn_size << ") : ";
    for( int i = 0 , is = cluster.size() ; i < is ; ++i ) {
      os << cluster[i] << "  ";
    }
    os << endl;
  } else if( CSV_FORMAT == output_format ) {
    for( int j = 0 , js = cluster.size() ; j < js ; ++j ) {
      os << j + 1 << "," << cluster.size() << "," << cluster.front()
         << "," << cluster[j] << "," << orig_nn_size << endl;
    }
  }

}

// ****************************************************************************
void output_new_clusters( const string &cluster_output_file ,
                          const vector<vector<string> > &clusters ,
                          const vector<unsigned int> &orig_sizes ,
                          CLUS_FILE_FORMAT output_format ) {

  ofstream ofs( cluster_output_file.c_str() );
  if( !ofs.good() ) {
    cerr << "ERROR : can't open " << cluster_output_file << " for writing."
         << endl;
    exit( 1 );
  }

  if( SAMPLES_FORMAT == output_format ) {
    ofs << "Molecule name : Cluster size : Cluster Members" << endl;
  }
  for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
    write_cluster( ofs , output_format , clusters[i] , orig_sizes[i] );
  }

}

// ****************************************************************************
void output_additions_file( const string &add_file ,
                            const vector<FingerprintBase *> &cluster_seed_fps ,
                            const vector<FingerprintBase *> &new_fps ,
                            const vector<int> &additions_dests ) {

  ofstream ofs( add_file.c_str() );
  if( !ofs || !ofs.good() ) {
    cerr << "Couldn't open " << add_file << " for writing." << endl;
    return;
  }

  for( int i = 0 , is = additions_dests.size() ; i < is ; ++i ) {
    ofs << new_fps[i]->get_name() << " ";
    if( -1 == additions_dests[i] ) {
      ofs << "NO_CLUSTER ";
    } else {
      ofs << cluster_seed_fps[additions_dests[i]]->get_name() << " ";
    }
    ofs << additions_dests[i] << endl;
  }

}

// ****************************************************************************
int main( int argc , char **argv ) {

  cout << "amtec - built " << BUILD_TIME << endl;

  AmtecSettings as( argc , argv );
  if( !as ) {
    cout << as.error_message() << endl << as.usage_text() << endl;
    cerr << as.error_message() << endl << as.usage_text() << endl;
  }

  if( TVERSKY == as.similarity_calc() ) {
    FingerprintBase::set_tversky_alpha( as.tversky_alpha() );
    HashedFingerprint::set_similarity_calc( as.similarity_calc() );
    NotHashedFingerprint::set_similarity_calc( as.similarity_calc() );
  }

  vector<vector<string> > clusters;

  read_cluster_file( as.input_cluster_file() , as.clus_input_format() ,
                     clusters );

  vector<FingerprintBase *> cluster_fps;
  try {
    read_fp_file( as.existing_cluster_fp_file() , as.input_format() ,
                  as.bitstring_separator() , cluster_fps );
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

  vector<FingerprintBase *> new_fps;
  try {
    read_fp_file( as.incoming_cluster_fp_file() , as.input_format() ,
                  as.bitstring_separator() , new_fps );
  } catch( DACLIB::FileReadOpenError &e ) {
    cerr << e.what() << endl;
    cout << e.what() << endl;
    exit( 1 );
  } catch( FingerprintFileError &e ) {
    cerr << e.what() << endl;
    cout << e.what() << endl;
    exit( 1 );
  }

  if( !as.new_subset_file().empty() ) {
    apply_subset_file( as.new_subset_file() , new_fps );
  }
  sort( new_fps.begin() , new_fps.end() ,
        bind( &FingerprintBase::get_name , _1 ) >
        bind( &FingerprintBase::get_name , _2 ) );

  vector<unsigned int> orig_sizes;
  orig_sizes.reserve( clusters.size() );
  for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
    orig_sizes.push_back( clusters[i].size() );
  }

  cout << "Adding " << new_fps.size() << " fingerprints"
       << " to " << clusters.size() << " clusters." << endl;
  vector<FingerprintBase *> cluster_seed_fps;
  vector<int> additions_dests; // where the new fps ended up
  add_fps_to_clusters( as.threshold() , cluster_fps , new_fps ,
                       cluster_seed_fps , clusters ,
                       additions_dests );

  re_sort_amended_clusters( orig_sizes , cluster_seed_fps ,
                            cluster_fps , new_fps , clusters );

  output_new_clusters( as.output_cluster_file() , clusters , orig_sizes ,
                       as.clus_output_format() );

  if( !as.additions_file().empty() ) {
    output_additions_file( as.additions_file() , cluster_seed_fps , new_fps , additions_dests );
  }

}
