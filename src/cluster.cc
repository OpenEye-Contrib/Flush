//
// file cluster.cc
// David Cosgrove
// AstraZeneca
// 23rd February 2007
//
// Does a sphere-exclusion clustering on a fingerprint file.

#include <algorithm>
#include <functional>
#include <iterator>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <mpi.h>

#include "stddefs.H"
#include "ClusterSettings.H"
#include "HashedFingerprint.H"
#include "NotHashedFingerprint.H"
#include "FileExceptions.H"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

using namespace boost;
using namespace std;

extern string BUILD_TIME;

class SortNbsByDist :
    public binary_function<pair<int,float> , pair<int,float> , bool> {
public :
  result_type operator()( first_argument_type a , second_argument_type b ) const {
    if( a.second == b.second )
      return a.first > b.first;
    else
      return a.second < b.second;
  }
};

namespace DACLIB {
// in eponymous file
string get_cwd();

// In mpi_string_subs.cc
void mpi_send_string( const string &str , int dest_rank );
void mpi_rec_string( int source_rank , string &str );
}

using namespace DAC_FINGERPRINTS;

typedef boost::shared_ptr<FingerprintBase> pFB;

// *******************************************************************************
void read_subset_file( const string &filename , vector<string> &subset_names ) {

  ifstream ifs( filename.c_str() );
  if( !ifs || !ifs.good() ) {
    throw DACLIB::FileReadOpenError( filename.c_str() );
  }

  string next_name;
  while( 1 ) {
    ifs >> next_name;
    if( ifs.eof() || !ifs.good() ) {
      break;
    }
    subset_names.push_back( next_name );
  }

  sort( subset_names.begin() , subset_names.end() );

}

// *******************************************************************************
void apply_subset_names( const vector<string> &subset_names , vector<pFB> &fps ) {

  int j = 0;
  for( int i = 0 , is = fps.size() ; i < is ; ++i ) {
    if( !binary_search( subset_names.begin() , subset_names.end() ,
                        fps[i]->get_name() ) ) {
      fps[j++] = fps[i];
    }
  }
  fps.erase( fps.begin() + j , fps.end() );

}

// *******************************************************************************
void open_fp_file( const string &filename ,
                   DAC_FINGERPRINTS::FP_FILE_FORMAT input_format ,
                   bool &byteswapping  , gzFile &pfile) {

  byteswapping = false;
  try {
    open_fp_file_for_reading( filename , input_format , byteswapping , pfile );
  } catch( DACLIB::FileReadOpenError &e ) {
    cerr << e.what() << endl;
    cout << e.what() << endl;
    exit( 1 );
  } catch( FingerprintFileError &e ) {
    cerr << e.what() << endl;
    cout << e.what() << endl;
    exit( 1 );
  }

}

// *******************************************************************************
string fix_spaces_in_fp_name( const string &fp_name , bool verbose ) {

  string new_name = fp_name;
  for( unsigned int j = 0 , js = new_name.size() ; j < js ; ++j ) {
    if( ' ' == new_name[j] ) {
      new_name[j] = '_';
    }
  }
  if( verbose ) {
    cout << "Fingerprint name " << fp_name
         << " has space(s). Changing to " << new_name << endl;
  }
  return new_name;

}

// *******************************************************************************
// read the fingerprint file, keeping only seeds and singleton fps
void read_fp_file( ClusterSettings cs , const vector<string> &seed_names ,
                   const vector<string> &singleton_names ,
                   vector<pFB> &seed_fps , vector<pFB> &singleton_fps ) {

  gzFile fpfile;
  bool byteswapping;
  open_fp_file( cs.input_file() , cs.input_format() , byteswapping , fpfile );

  while( 1 ) {
    pFB fp( read_next_fp_from_file( fpfile , byteswapping , cs.input_format() ,
                                    cs.bitstring_separator() ) );
    if( !fp ) {
      break;
    }
    if( SAMPLES_FORMAT == cs.output_format() ) {
      // if we've got this far and there's a space in the name, we need to fix it.
      // We'd have stopped by now if that was not the case.
      fp->set_name( fix_spaces_in_fp_name( fp->get_name() , false ) );
    }

    if( binary_search( seed_names.begin() , seed_names.end() ,
                       fp->get_name() ) ) {
      seed_fps.push_back( fp );
    }
    if( binary_search( singleton_names.begin() , singleton_names.end() ,
                       fp->get_name() ) ) {
      singleton_fps.push_back( fp );
    }
  }

}

// *******************************************************************************
void fps_to_name( vector<pFB> &fps , vector<string> &fp_names ) {

  fp_names.reserve( fps.size() );
  for( int i = 0 , is = fps.size() ; i < is ; ++i ) {
    fp_names.push_back( fps[i]->get_name() );
  }
  fps.clear();

}

// *******************************************************************************
void apply_subset( ClusterSettings &cs , vector<pFB> &fps ) {

  if( !cs.subset_file().empty() ) {
    vector<string> subset_names;
    read_subset_file( cs.subset_file() , subset_names );
    apply_subset_names( subset_names , fps );
  }

}

// *******************************************************************************
void make_nnlists( bool warm_feeling , double threshold ,
                   unsigned int start_num , unsigned int stop_num ,
                   const vector<pFB> &fps ,
                   vector<vector<int> > &nns ) {

  stop_num = stop_num > fps.size() ? fps.size() : stop_num;
  if( warm_feeling ) {
    cout << "Creating neighbour lists for fps " << start_num
         << " to " << stop_num << endl;
  }

  for( unsigned int i = start_num ; i < stop_num ; ++i ) {

    vector<pair<int,float> > nbs;
    nbs.push_back( make_pair( i , 0.0F ) );
    for( unsigned int j = 0 , js = fps.size() ; j < js ; ++j ) {
      if( i == j ) {
        continue;
      }
      double dist = fps[i]->calc_distance( *fps[j] , threshold );
      if( dist < threshold ) {
        nbs.push_back( make_pair( j , dist ) );
      }
    }
    if( nbs.size() > 1 ) {
      sort( nbs.begin() + 1 , nbs.end() , SortNbsByDist() );
    }
    nns.push_back( vector<int>() );
    transform( nbs.begin() , nbs.end() , back_inserter( nns.back() ) ,
               bind( &pair<int,float>::first , _1 ) );
    if( warm_feeling && (i - start_num) && !( ( i-start_num ) % 1000 ) ) {
      cout << "Generated " << i - start_num << " near-neighbour lists."
           << endl;
    }

  }

  if( warm_feeling ) {
    cout << "Generated all " << stop_num - start_num << " near-neighbour lists."
         << endl;
  }

}

// *******************************************************************************
void check_for_spaces_in_fp_names( bool fix_spaces , vector<pFB> &fps ) {

  for( unsigned int i = 0 , is = fps.size() ; i < is ; ++i ) {
    if( string::npos != fps[i]->get_name().find( ' ' ) ) {
      if( !fix_spaces ) {
        cout << "Fingerprint " << i << " name " << fps[i]->get_name()
             << " has space(s) in its name.  Either use output format CSV"
             << " or --fix-spaces-in-names" << endl;
        exit( 1 );
      } else {
        fps[i]->set_name( fix_spaces_in_fp_name( fps[i]->get_name() , true ) );
      }
    }
  }

}

// *******************************************************************************
void check_for_spaces_in_fp_names( bool fix_spaces , vector<string> &fp_names ) {

  for( unsigned int i = 0 , is = fp_names.size() ; i < is ; ++i ) {
    if( string::npos != fp_names[i].find( ' ' ) ) {
      if( !fix_spaces ) {
        cout << "Fingerprint " << i << " name " << fp_names[i]
             << " has space(s) in its name.  Either use output format CSV"
             << " or --fix-spaces-in-names" << endl;
        exit( 1 );
      } else {
        fp_names[i] = fix_spaces_in_fp_name( fp_names[i] , true );
      }
    }
  }

}

// *******************************************************************************
void make_nnlists( ClusterSettings &cs , unsigned int start_fp ,
                   unsigned int &num_fps_to_do , vector<string> &fp_names ,
                   vector<vector<int> > &nns ) {

  gzFile gzfp;
  bool byteswapping;
  open_fp_file( cs.input_file() , cs.input_format() , byteswapping , gzfp );

  // read all the fps from the file, which we'll need even if we're only
  // doing a portion of the nnlists
  vector<FingerprintBase *> raw_fps;
  read_fps_from_file( gzfp , byteswapping , cs.input_format() , cs.bitstring_separator() ,
                      0 , numeric_limits<unsigned int>::max() , raw_fps );

  vector<pFB> fps;
  fps.reserve( raw_fps.size() );
  for( int i = 0 , is = raw_fps.size() ; i < is ; ++i ) {
    fps.push_back( pFB( raw_fps[i] ) );
  }
  apply_subset( cs , fps );

  if( SAMPLES_FORMAT == cs.output_format() ) {
    // this will stop the program if there are some and cs.fix_spaces_in_names()
    // is false
    check_for_spaces_in_fp_names( cs.fix_spaces_in_names() , fps );
  }

  nns.reserve( fps.size() );
  unsigned int stop_fp = start_fp + num_fps_to_do;
  if( stop_fp > fps.size() ) {
    num_fps_to_do = fps.size() - start_fp;
    stop_fp = fps.size();
#ifdef NOTYET
    cout << "revised num_fps_to_do to " << num_fps_to_do << endl;
#endif
  }
  make_nnlists( cs.warm_feeling() , cs.threshold() , start_fp , stop_fp , fps ,
                nns );

  // pull the names out of the fingerprints and delete
  fps_to_name( fps , fp_names );

#ifdef NOTYET
  cout << "leaving make_nnlists" << endl;
  for( int i = 0 , is = nns.size() ; i < is ; ++i ) {
    cout << fp_names[nns[i].front()] << " : ";
    for( int j = 0 , js = nns[i].size() ; j < js ; ++j ) {
      cout << nns[i][j] << " ";
    }
    cout << endl;
  }
#endif

}

// *******************************************************************************
void make_orig_nn_sizes( unsigned int num_fps , unsigned int start_fp ,
                         unsigned int num_fps_to_do ,
                         const vector<vector<int> > &nns ,
                         vector<int> &orig_nn_sizes ) {

  orig_nn_sizes = vector<int>( num_fps , 0 );
  unsigned int stop_fp = start_fp + num_fps_to_do;
  for( unsigned int i = start_fp , j = 0 ; i < stop_fp ; ++i , ++j ) {
    orig_nn_sizes[i] = nns[j].size();
  }

}

// *******************************************************************************
// find the next cluster, returning the sequence number. The next cluster is the
// one with the largest nn list, with the original neighbour list size as first
// tie-breaker, sequence number of seed as second tie-breaker - the lower being
// preferred. orig_nn_sizes is not shortened as nns is, so need to take care.
int find_next_seed( const vector<vector<int> > &nns ,
                    const vector<int> &orig_nn_sizes ) {

  int next_seed = 0;
  for( int i = 1 , is = nns.size() ; i < is ; ++i ) {
    if( nns[i].size() > nns[next_seed].size() ) {
      next_seed = i;
    } else if( nns[i].size() == nns[next_seed].size() ) {
      if( orig_nn_sizes[nns[i].front()] >= orig_nn_sizes[nns[next_seed].front()] ) {
        next_seed = i;
      }
    }
  }

  return next_seed;

}

// *******************************************************************************
void write_cluster( OUTPUT_FORMAT output_format , int clus_num ,
                    const vector<string> &fp_names , const vector<int> &clus ,
                    int orig_nn_size , ostream &output_stream ,
                    vector<string> &seed_names ,
                    vector<string> &singleton_names ) {

  switch( output_format ) {
  case SAMPLES_FORMAT :
    output_stream << fp_names[clus.front()] << " : "
                                            << clus.size() << "(" << orig_nn_size << ") : ";

    // it would be tidier and more elegant not to put a "  " at the end of the
    // line, but flush_clus used to so I have left it in in case people have
    // relied upon it. At the very least, it makes it easier to compare output
    // when debugging.
    for( int i = 0 , is = clus.size() ; i < is ; ++i ) {
      output_stream << fp_names[clus[i]] << "  ";
    }
    output_stream << endl;
    break;
  case CSV_FORMAT :
    for( int i = 0 , is = clus.size() ; i < is ; ++i ) {
      output_stream << clus_num << "," << is << ","
                    << fp_names[clus.front()] << ","
                    << fp_names[clus[i]] << "," << orig_nn_size << endl;
    }
    break;
  }

  seed_names.push_back( fp_names[clus.front()] );
  if( 1 == clus.size() ) {
    singleton_names.push_back( fp_names[clus.front()] );
  }

}

// *******************************************************************************
void remove_cluster_from_nns( unsigned int num_fps , const vector<int> &cluster ,
                              vector<vector<int> > &nns ) {

  static vector<char> in_cluster;
  if( in_cluster.empty() ) {
    in_cluster = vector<char>( num_fps , 0 );
  }

  for(int i = 0 , is = cluster.size() ; i < is ; ++i ) {
    in_cluster[cluster[i]] = 1;
  }

  for( int i = 0 , is = nns.size() ; i < is ; ++i ) {
    if( nns[i].empty() ) {
      continue;
    }
    // if the seed is in the cluster, take out the entire cluster
    if( in_cluster[nns[i][0]] ) {
      nns[i].clear();
    }
    // set neighbours in the cluster to -1 and tidy them out at the end
    for( int j = 0 , js = nns[i].size() ; j < js ; ++j ) {
      if( in_cluster[nns[i][j]] ) {
        nns[i][j] = -1;
      }
    }
    nns[i].erase( remove( nns[i].begin() , nns[i].end() , -1 ) , nns[i].end() );
  }

  // take out any empty clusters
  nns.erase( remove_if( nns.begin() , nns.end() ,
                        bind( &vector<int>::empty , _1 ) ) , nns.end() );

  // reset in_cluster for next time round
  for(int i = 0 , is = cluster.size() ; i < is ; ++i ) {
    in_cluster[cluster[i]] = 0;
  }

}

// *******************************************************************************
// do the clustering and output as we go
void output_clusters( bool warm_feeling , vector<string> &fp_names ,
                      vector<vector<int> > &nns , OUTPUT_FORMAT output_format ,
                      ostream &output_stream , vector<string> &seed_names ,
                      vector<string> &singleton_names ) {

  if( SAMPLES_FORMAT == output_format ) {
    output_stream << "Molecule name : Cluster size : Cluster Members" << endl;
  }

  vector<int> orig_nn_sizes( nns.size() , -1 );
  for( int i = 0 , is = nns.size() ; i < is ; ++i ) {
    orig_nn_sizes[i] = nns[i].size();
  }

  int num_written = 0 , tot = 0;
  while( !nns.empty() ) {

    int next_seed_num = find_next_seed( nns , orig_nn_sizes );

    write_cluster( output_format , num_written + 1 , fp_names ,
                   nns[next_seed_num] ,
                   orig_nn_sizes[nns[next_seed_num].front()] , output_stream ,
        seed_names , singleton_names );
    ++num_written;
    tot += nns[next_seed_num].size();
    vector<int> cluster( nns[next_seed_num] );
    remove_cluster_from_nns( fp_names.size() , cluster , nns );

    if( warm_feeling && !( num_written % 100 ) ) {
      cout << "Written " << num_written << " clusters, average size "
           << tot / num_written << "." << endl;
    }
  }

  string fp_out = " fingerprint";
  if( fp_names.size() > 1 ) {
    fp_out += "s";
  }
  string clus_out = " cluster";
  if( num_written > 1 ) {
    clus_out += "s";
  }
  cout << "Clustered " << fp_names.size() << fp_out << " into "
       << num_written << clus_out << "." << endl;

}

// *******************************************************************************
void serial_run( ClusterSettings &cs , vector<string> &seed_names ,
                 vector<string> &singleton_names ) {

  // open the output stream right away, in case we can't. It's best to find
  // out before we've done a potentially long job.
  ofstream output_stream( cs.output_file().c_str() );
  if( !output_stream.good() ) {
    cerr << "Couldn't open " << cs.output_file() << " for writing." << endl;
    exit( 1 );
  }

  vector<vector<int> > nns;
  vector<string> fp_names;

  unsigned num_fps_to_do = numeric_limits<unsigned int>::max();
  make_nnlists( cs , 0 , num_fps_to_do , fp_names , nns );
  output_clusters( cs.warm_feeling() , fp_names , nns , cs.output_format() ,
                   output_stream , seed_names , singleton_names );

}

// *******************************************************************************
void send_cluster_to_master( const vector<vector<int> > &nnlists ) {


  int clus_num;
  MPI_Recv( &clus_num , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );

  unsigned int i = nnlists[clus_num].size();
  MPI_Send( &i , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD );
  MPI_Send( const_cast<int *>( &(nnlists[clus_num])[0] ) , i , MPI_INT , 0 , 0 , MPI_COMM_WORLD );

}

// *******************************************************************************
void send_best_cluster_details_to_master( const vector<vector<int> > &nnlists ,
                                          const vector<int> &orig_nn_sizes ) {

#ifdef NOTYET
  int world_rank;
  MPI_Comm_rank( MPI_COMM_WORLD , &world_rank );
  cout << "send_best_cluster_details_to_master for " << world_rank << endl;
#endif

  if( nnlists.empty() ) {
    int i = -1;
    MPI_Send( &i , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD );
    return;
  }

  int best_clus = find_next_seed( nnlists , orig_nn_sizes );
  MPI_Send( &best_clus , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD );

  // send the size, the original size and the first member in one go
  unsigned int best_clus_details[3];
  best_clus_details[0] = nnlists[best_clus].size();
  best_clus_details[1] = orig_nn_sizes[nnlists[best_clus].front()];
  best_clus_details[2] = nnlists[best_clus].front();
  MPI_Send( &best_clus_details , 3 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD );

}

// *******************************************************************************
void send_best_cluster_to_master( const vector<vector<int> > &nnlists ,
                                  const vector<int> &orig_nn_sizes ) {

  if( nnlists.empty() ) {
    int i = -1;
    MPI_Send( &i , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD );
    return;
  }

  int best_clus = find_next_seed( nnlists , orig_nn_sizes );
  int i = nnlists[best_clus].size();
  MPI_Send( &i , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD );
  MPI_Send( const_cast<int *>( &nnlists[best_clus][0] ) , i , MPI_INT , 0 , 0 , MPI_COMM_WORLD );
  i = orig_nn_sizes[nnlists[best_clus][0]];
  MPI_Send( &i , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD );

#ifdef NOTYET
  cout << "sending best_clus : " << best_clus << " size " << nnlists[best_clus].size()
       << " orig nn size : " << i << endl;
#endif

}

// *******************************************************************************
void cross_off_cluster( unsigned int num_fps ,
                        vector<vector<int> > &nnlists ) {

  int clus_size = 0;
  MPI_Recv( &clus_size , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  vector<int> cluster( clus_size , -1 );
  MPI_Recv( &cluster[0] , clus_size , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );

  remove_cluster_from_nns( num_fps , cluster , nnlists );

}

// *******************************************************************************
void receive_best_clusters_from_slaves( int world_size ,
                                        vector<int> &cluster , int &orig_nn_size ) {

  for( int i = 1 ; i < world_size ; ++i ) {
    DACLIB::mpi_send_string( string( "Send_Best_Cluster_Details" ) , i );
  }

  int best_clus = -1 , best_clus_slave = 0;
  unsigned int best_clus_size = 0 , best_orig_nn_size = 0 , best_clus_front = 0;

  int num_to_report = 0;
  int num_slaves = world_size - 1; // process 0 is the master

  while( num_to_report < num_slaves ) {

    MPI_Status status;
    MPI_Probe( MPI_ANY_SOURCE , 0 , MPI_COMM_WORLD , &status );
    int this_best_clus;
    MPI_Recv( &this_best_clus , 1 , MPI_INT , status.MPI_SOURCE , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
    ++num_to_report;
#ifdef NOTYET
    cout << "slave " << status.MPI_SOURCE << " reports this_best_clus = " << this_best_clus << endl;
#endif
    if( -1 == this_best_clus ) {
      continue; // this slave has no more nnlists
    }
    unsigned int this_best_clus_details[3];
    MPI_Recv( &this_best_clus_details , 3 , MPI_UNSIGNED , status.MPI_SOURCE , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
#ifdef NOTYET
    cout << "details : " << this_best_clus_details[0]
            << " : " << this_best_clus_details[1] << " : " << this_best_clus_details[2] << endl;
#endif
    // this_best_clus_details contains, in order, the size of the best cluster
    // on slave i, the orginal nnlists size of that cluster, and the first
    // member of the cluster.
    if( this_best_clus_details[0] > best_clus_size ) {
      best_clus_size = this_best_clus_details[0];
      best_orig_nn_size = this_best_clus_details[1];
      best_clus_front = this_best_clus_details[2];
      best_clus = this_best_clus;
      best_clus_slave = status.MPI_SOURCE;
    } else if( this_best_clus_details[0] == best_clus_size ) {
      if( this_best_clus_details[1] > best_orig_nn_size ) {
        best_orig_nn_size = this_best_clus_details[1];
        best_clus_front = this_best_clus_details[2];
        best_clus = this_best_clus;
        best_clus_slave = status.MPI_SOURCE;
      } else if( this_best_clus_details[1] == best_orig_nn_size ) {
        if( this_best_clus_details[2] > best_clus_front ) {
          best_clus_front = this_best_clus_details[2];
          best_clus = this_best_clus;
          best_clus_slave = status.MPI_SOURCE;
        }
      }
    }
  }

#ifdef NOTYET
  cout << "best_clus = " << best_clus << endl;
#endif

  if( -1 == best_clus ) {
    // we're done, flag it as such
    cluster.clear();
    orig_nn_size = 0;
  } else {
    // need to get the best cluster off the appropriate slave
    DACLIB::mpi_send_string( string( "Send_Cluster" ) , best_clus_slave );
    MPI_Send( &best_clus , 1 , MPI_INT , best_clus_slave , 0 , MPI_COMM_WORLD );

    cluster.resize( best_clus_size );
    unsigned int clus_size;
    MPI_Recv( &clus_size , 1 , MPI_UNSIGNED , best_clus_slave , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
    if( clus_size != best_clus_size ) {
      cerr << "AWOOGA - slave " << best_clus_slave
           << " has sent the wrong cluster. It's clearly messed up, so we're done."
           << endl;
      MPI_Finalize();
      exit( 1 );
    }
    MPI_Recv( &cluster[0] , clus_size , MPI_INT , best_clus_slave , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
    orig_nn_size = best_orig_nn_size;
  }

}

// *******************************************************************************
void tell_slaves_to_cross_off_cluster( int world_size ,
                                       const vector<int> &cluster ) {

  for( int i = 1 ; i < world_size ; ++i ) {
    DACLIB::mpi_send_string( string( "Cross_Off_Cluster" ) , i );

    int j = cluster.size();
    MPI_Send( &j , 1 , MPI_INT , i , 0 , MPI_COMM_WORLD );
    MPI_Send( const_cast<int *>( &cluster[0] ) , j , MPI_INT , i , 0 , MPI_COMM_WORLD );
  }

}

// *******************************************************************************
void send_search_details( ClusterSettings &cs , unsigned int num_fps ,
                          int world_size , unsigned int &slave_does ) {

  // each slave needs to do something like num_fps / num_procs fps
  // each against all targets
  slave_does = num_fps / ( world_size - 1 );
  while( slave_does * ( world_size - 1 ) < num_fps ) {
    ++slave_does;
  }
  if( cs.warm_feeling() ) {
    cout << "Each slave does " << slave_does << " fps " << endl;
  }

  for( int i = 1 ; i < world_size ; ++i ) {

    DACLIB::mpi_send_string( string( "Search_Details" ) , i );
    cs.send_contents_via_mpi( i );
    // send the number of fps each slave must do, and the slave number,
    // so it knows where to start
    MPI_Send( &slave_does , 1 , MPI_UNSIGNED , i , 0 , MPI_COMM_WORLD );
    MPI_Send( &i , 1 , MPI_INT , i , 0 , MPI_COMM_WORLD );

  }

}

// *******************************************************************************
void receive_search_details( ClusterSettings &cs , unsigned int &num_fps_to_do ,
                             unsigned int &start_fp ) {

  cs.receive_contents_via_mpi();
  MPI_Recv( &num_fps_to_do , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  int slave_rank;
  MPI_Recv( &slave_rank , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  --slave_rank; // first slave has rank 1.
  start_fp = slave_rank * num_fps_to_do;

  if( cs.warm_feeling() ) {
    cout << "This slave to do " << num_fps_to_do << " nnlists, starting at fp "
         << start_fp << endl;
  }

}

// ********************************************************************
void send_cwd_to_slaves( int world_size ) {

  string cwd = DACLIB::get_cwd();
  if( cwd.length() ) {
    // process 0 is the master
    for( int i = 1 ; i < world_size ; ++i ) {
      DACLIB::mpi_send_string( string( "New_CWD" ) , i );
      DACLIB::mpi_send_string( cwd , i );
    }
  }

}

// *******************************************************************************
void tell_master_slave_has_done_nnlists() {

  DACLIB::mpi_send_string( string( "NNLists_Done" ) , 0 );

}

// *******************************************************************************
void parallel_run( ClusterSettings &cs , int world_size ,
                   vector<string> &seed_names ,
                   vector<string> &singleton_names ) {

  // open the output stream right away, in case we can't. It's best to find
  // out before we've done a potentially long job.
  ofstream output_stream( cs.output_file().c_str() );
  if( !output_stream.good() ) {
    cerr << "Couldn't open " << cs.output_file() << " for writing." << endl;
    exit( 1 );
  }
  if( SAMPLES_FORMAT == cs.output_format() ) {
    output_stream << "Molecule name : Cluster size : Cluster Members" << endl;
  }

  unsigned int num_fps = 0;
  try {
    num_fps = count_fps_in_file( cs.input_file() , cs.input_format() ,
                                 cs.bitstring_separator() );
  } catch( DACLIB::FileReadOpenError &e ) {
    cerr << e.what() << endl;
    cout << e.what() << endl;
    exit( 1 );
  } catch( FingerprintFileError &e ) {
    cerr << e.what() << endl;
    cout << e.what() << endl;
    exit( 1 );
  }

  if( num_fps ) {
    send_cwd_to_slaves( world_size );
    // send_search_details also fires off the jobs on the slaves
    unsigned int chunk_size;
    send_search_details( cs , num_fps , world_size , chunk_size );
    if( cs.warm_feeling() ) {
      cout << "NN list requirements all sent. Each slave will produce "
           << chunk_size << " NN lists." << endl;
    }
    // whilst the slaves are making the nnlists, there's time to re-read
    // the file and get the fingerprint names out
    vector<string> fp_names;
    get_fp_names( cs.input_file() , cs.input_format() ,
                  cs.bitstring_separator() , fp_names );
    check_for_spaces_in_fp_names( cs.fix_spaces_in_names() , fp_names );

    // wait for all slaves to announce they're done. This is because we don't
    // want the master sending messages to a slave until it's finished the
    // nnlists, because that interferes with messages if the master dies
    // and we want the slave to get those so it stops early if necessary.
    int num_finished = 0;
    int num_slaves = world_size - 1; // process 0 is the master
    while( num_finished < num_slaves ) {
      MPI_Status status;
      MPI_Probe( MPI_ANY_SOURCE , 0 , MPI_COMM_WORLD , &status );

      string msg;
      DACLIB::mpi_rec_string( status.MPI_SOURCE , msg );
      if( cs.warm_feeling() ) {
        cout << "Slave " << status.MPI_SOURCE << " has finished nnlists." << endl;
      }
      ++num_finished;
    }

    int num_written = 0 , tot = 0;
    while( 1 ) {

      vector<int> cluster;
      int orig_nn_size = 0;
      receive_best_clusters_from_slaves( world_size , cluster , orig_nn_size );
      if( cluster.empty() ) {
        break;
      }
      write_cluster( cs.output_format() , num_written + 1 , fp_names , cluster ,
                     orig_nn_size , output_stream , seed_names ,
                     singleton_names );
      ++num_written;
      tot += cluster.size();
      if( cs.warm_feeling() && !( num_written % 100 ) ) {
        cout << "Written " << num_written << " clusters, average size "
             << tot / num_written << "." << endl;
      }

      tell_slaves_to_cross_off_cluster( world_size , cluster );

    }

    for( int i = 1 ; i < world_size ; ++i ) {
      DACLIB::mpi_send_string( string( "Finished" ) , i );
    }

    string fp_out = " fingerprint";
    if( tot > 1 ) {
      fp_out += "s";
    }
    string clus_out = " cluster";
    if( num_written > 1 ) {
      clus_out += "s";
    }
    cout << "Clustered " << tot << fp_out << " into "
         << num_written << clus_out << "." << endl;
  }

}

// ********************************************************************
void receive_new_cwd() {

  string new_cwd;
  DACLIB::mpi_rec_string( 0 , new_cwd );
  if(chdir( new_cwd.c_str() ) < 0){
    cerr << "ERROR : couldn't change to directory " << new_cwd << endl;
    exit(1);
  }

}

// *******************************************************************************
void slave_event_loop() {

  string msg;

  ClusterSettings cs;
  unsigned int num_fps_to_do , start_fp;
  vector<vector<int> > nns;
  vector<int> orig_nn_sizes;
  vector<string> fp_names;

#ifdef NOTYET
  int wr;
  MPI_Comm_rank( MPI_COMM_WORLD , &wr );
  cout << "Slave event loop for : " << wr << endl;
#endif

  while( 1 ) {

    DACLIB::mpi_rec_string( 0 , msg );
#if DEBUG == 1
    cout << world_rank << " : received message : " << msg << endl;
#endif
    if( string( "Finished" ) == msg ) {
      break;
    } else if( string( "Search_Details" ) == msg ) {
      receive_search_details( cs , num_fps_to_do , start_fp );
      make_nnlists( cs , start_fp , num_fps_to_do , fp_names , nns );
      // orig_nn_sizes needs to be indexed for the original fp set.
      make_orig_nn_sizes( fp_names.size() , start_fp , num_fps_to_do ,
                          nns , orig_nn_sizes );
      tell_master_slave_has_done_nnlists();
    } else if( string( "Send_Best_Cluster" ) == msg ) {
      send_best_cluster_to_master( nns , orig_nn_sizes );
    } else if( string( "Send_Best_Cluster_Details" ) == msg ) {
      send_best_cluster_details_to_master( nns , orig_nn_sizes );
    } else if( string( "Send_Cluster" ) == msg ) {
      // reads the cluster number off the pvm buffer, sends that cluster to
      // the master
      send_cluster_to_master( nns );
    } else if( string( "Cross_Off_Cluster" ) == msg ) {
      cross_off_cluster( fp_names.size() , nns );
    } else if( string( "New_CWD" ) == msg ) {
      receive_new_cwd();
    } else {
      int world_rank;
      MPI_Comm_rank( MPI_COMM_WORLD , &world_rank );
      cout << world_rank << " received suspect message " << msg << endl;
    }

  }

}

// *******************************************************************************
void read_next_samples_cluster( istream &is ,
                                vector<string> &cluster , int &orig_nn_size ) {

  string next_line;
  getline( is , next_line );
  if( is.eof() ) {
    return;
  }
  istringstream iss( next_line );
  string tmp , nns;
  iss >> tmp >> tmp >> nns >> tmp;
  cluster = ( vector<string>( istream_iterator<string>( iss ) ,
                               istream_iterator<string>() ) );
  vector<string> splits;
  boost::algorithm::split( splits , nns , boost::algorithm::is_any_of( "()" ) );
  if( splits.size() > 1 ) {
    orig_nn_size = lexical_cast<int>( splits[1] );
  } else {
    orig_nn_size = -1; // but there's a problem.
  }

}

// *******************************************************************************
void read_next_csv_cluster( istream &is ,
                            vector<string> &cluster , int &orig_nn_size ) {

  string next_line;
  getline( is , next_line );
  if( is.eof() ) {
    return;
  }
  vector<string> splits;
  boost::algorithm::split( splits , next_line , boost::algorithm::is_any_of( "," ) );
  if( splits.size() < 5 ) {
    return;
  }
  cluster.push_back( splits[3] );
  orig_nn_size = lexical_cast<int>( splits[4] );
  unsigned int clus_size = boost::lexical_cast<int>( splits[1] );
  // we've already read the first line, so read clus_size - 1 more.
  for( unsigned int i = 0 ; i < clus_size - 1 ; ++i ) {
    getline( is , next_line );
    if( is.eof() ) {
      return;
    }
    vector<string> splits1;
    boost::algorithm::split( splits1 , next_line , boost::algorithm::is_any_of( "," ) );
    if( splits1.size() < 5 ) {
      return;
    }
    cluster.push_back( splits1[3] );
  }

}

// *******************************************************************************
void read_next_cluster( istream &is , ClusterSettings &cs ,
                        vector<string> &cluster , int &orig_nn_size ) {

  static bool first_call( true );
  if( SAMPLES_FORMAT == cs.output_format() ) {
    if( first_call ) {
      // 1st line is headings
      string next_line;
      getline( is , next_line );
      first_call = false;
    }
    read_next_samples_cluster( is , cluster , orig_nn_size );
  } else if( CSV_FORMAT == cs.output_format() ) {
    read_next_csv_cluster( is , cluster , orig_nn_size );
  }

}

// *******************************************************************************
void write_cluster( ostream &os , OUTPUT_FORMAT output_format ,
                    const vector<string> &cluster ,
                    int orig_nn_size , int clus_num ) {

  if( SAMPLES_FORMAT == output_format ) {
    os << cluster.front() << " : " << cluster.size() << "("
       << orig_nn_size << ") : ";
    for( int i = 0 , is = cluster.size() ; i < is ; ++i ) {
      os << cluster[i] << "  ";
    }
    os << endl;
  } else if( CSV_FORMAT == output_format ) {
    for( int j = 0 , js = cluster.size() ; j < js ; ++j ) {
      os << clus_num << "," << cluster.size() << "," << cluster.front()
         << "," << cluster[j] << "," << orig_nn_size << endl;
    }
  }

}

// *******************************************************************************
void update_clusters_file( ClusterSettings &cs ,
                          const vector<vector<pair<int,float> > > &seed_nbs ,
                          const map<string,int> &seed_fps_map ,
                          const vector<pFB> &seed_fps ,
                          const vector<pFB> &singleton_fps ) {

  ifstream ifs( cs.output_file().c_str() );
  filesystem::path temp = boost::filesystem::unique_path();
  const string tmp_clus_file = temp.native();  // optional
  ofstream ofs( tmp_clus_file.c_str() );
  if( SAMPLES_FORMAT == cs.output_format() ) {
    ofs << "Molecule name : Cluster size : Cluster Members" << endl;
  }

  int orig_nn_size , clus_num = 1;
  while( 1 ) {
    vector<string> cluster;
    read_next_cluster( ifs , cs , cluster , orig_nn_size );
    if( cluster.empty() ) {
      break;
    }
#ifdef NOTYET
    cout << "Next cluster : ";
    copy( cluster.begin() , cluster.end() , stringOut );
    cout << " : " << orig_nn_size << endl;
#endif
    map<string,int>::const_iterator p = seed_fps_map.find( cluster.front() );
    // if seed_fps[p->second] is empty, this is a singleton that's been put
    // in a different cluster, so don't write it
    if( p == seed_fps_map.end() ) {
      cout << "ERROR : failed to find seed " << cluster.front()
           << " in seeds map." << endl;
      cout << "seeds size : " << seed_fps_map.size() << endl;
      map<string,int>::const_iterator pp , pps;
      for( pp = seed_fps_map.begin() , pps = seed_fps_map.end() ; pp != pps ; ++pp ) {
        cout << pp->first << endl;
      }
      exit( 1 );
    }
    if( seed_fps[p->second] ) {
      if( p != seed_fps_map.end() && !seed_nbs[p->second].empty() ) {
        for( int i = 0 , is = seed_nbs[p->second].size() ; i < is ; ++i ) {
          cluster.push_back( singleton_fps[seed_nbs[p->second][i].first]->get_name() );
        }
      }
      write_cluster( ofs , cs.output_format() , cluster , orig_nn_size ,
                     clus_num );
      ++clus_num;
    }
  }

  ifs.close();
  filesystem::remove( filesystem::path( cs.output_file() ) );
  try {
    filesystem::rename( temp , filesystem::path( cs.output_file() ) );
  } catch( filesystem::filesystem_error &e ) {
    filesystem::copy_file( temp , filesystem::path( cs.output_file() ) );
    filesystem::remove( filesystem::path( temp ) );
  }

}

// *******************************************************************************
void collapse_singletons( ClusterSettings &cs ,
                          vector<string> &seed_names ,
                          vector<string> &singleton_names ) {

  if( cs.warm_feeling() ) {
    cout << "Collapse_singletons at " << cs.singletons_threshold() << "." << endl;
    if( 1 == singleton_names.size() ) {
      cout << "There is 1 singleton";
    } else {
      cout << "There are " << singleton_names.size() << " singletons";
    }
    cout << " to slot into " << seed_names.size() << " cluster";
    if( seed_names.size() > 1 ) {
      cout << "s";
    }
    cout << "." << endl;
  }

  sort( seed_names.begin() , seed_names.end() );
  sort( singleton_names.begin() , singleton_names.end() );

  // re-read the fp file for the seeds and singletons
  vector<pFB> seed_fps , singleton_fps;
  read_fp_file( cs , seed_names , singleton_names , seed_fps , singleton_fps );

  vector<vector<pair<int,float> > >seed_nbs( seed_names.size() );
  map<string,int> seed_fps_map , singleton_fps_map;
  for( int i = 0 , is = seed_fps.size() ; i < is ; ++i ) {
    seed_fps_map.insert( make_pair( seed_fps[i]->get_name() , i ) );
  }
  for( int i = 0 , is = singleton_fps.size() ; i < is ; ++i ) {
    singleton_fps_map.insert( make_pair( singleton_fps[i]->get_name() , i ) );
  }

  for( int i = 0 , is = singleton_fps.size() ; i < is ; ++i ) {
    if( !singleton_fps[i] ) {
      continue; // it might have been promoted by now
    }
    double nearest_dist = cs.singletons_threshold();
    int nearest_seed = -1;
    for( int j = 0 , js = seed_fps.size() ; j < js ; ++j ) {
      if( !seed_fps[j] || singleton_fps[i]->get_name() == seed_fps[j]->get_name() ) {
        continue;
      }
      double dist = seed_fps[j]->calc_distance( *singleton_fps[i] , cs.singletons_threshold() );
      if( dist < nearest_dist ) {
        nearest_dist = dist;
        nearest_seed = j;
      }
    }
    if( nearest_seed != -1 ) {
      if( cs.warm_feeling() ) {
        cout << "Singleton " << singleton_fps[i]->get_name() << " goes into cluster of "
                << seed_fps[nearest_seed]->get_name() << " at distance " << nearest_dist << endl;
      }
      seed_nbs[nearest_seed].push_back( make_pair( i , nearest_dist ) );
      // remove this singleton from seeds
      map<string,int>::const_iterator p = seed_fps_map.find( singleton_fps[i]->get_name() );
      seed_fps[p->second].reset( static_cast<FingerprintBase *>( 0 ) );
      // if the seed was a singleton it isn't any more, so take it out
      p = singleton_fps_map.find( seed_fps[nearest_seed]->get_name() );
      if( p != singleton_fps_map.end() ) {
        singleton_fps[p->second].reset( static_cast<FingerprintBase *>( 0 ) );
      }
    }
  }

#ifdef NOTYET
  for( int i = 0 , is = seed_nbs.size() ; i < is ; ++i ) {
    sort( seed_nbs[i].begin() , seed_nbs[i].end() , SortNbsByDist() );
    if( seed_fps[i] ) {
      cout << seed_fps[i]->get_name() << " : ";
      for( int j = 0 , js = seed_nbs[i].size() ; j < js ; ++j ) {
        cout << seed_nbs[i][j].first << " ";
      }
      cout << endl;
    } else {
      cout << i << " is empty" << endl;
    }
  }
#endif

  update_clusters_file( cs , seed_nbs , seed_fps_map , seed_fps , singleton_fps );

}

// *******************************************************************************
int main( int argc , char **argv ) {

  cout << "cluster - built " << BUILD_TIME << endl;

  // sort out the  MPI environment
  MPI_Init( NULL , NULL );
  int world_rank , world_size;
  MPI_Comm_rank( MPI_COMM_WORLD , &world_rank );
  MPI_Comm_size( MPI_COMM_WORLD , &world_size );

#ifdef NOTYET
  cout << "world_rank : " << world_rank << endl
       << "world size : " << world_size << endl;
#endif

  // if we're not world_rank 0, we're a slave.
  if( world_rank && world_size > 1 ) {
    slave_event_loop();
    MPI_Finalize();
    exit( 0 );
  }

  try {
    ClusterSettings cs( argc , argv );
    if( !cs ) {
      cout << cs.error_message() << endl;
      cout << cs.usage_text() << endl;
      cerr << cs.error_message() << endl;
      cerr << cs.usage_text() << endl;
      MPI_Finalize();
      exit( 1 );
    }

    vector<string> seed_names , singleton_names;
    if( 1 == world_size ) {
      serial_run( cs , seed_names , singleton_names );
    } else {
      parallel_run( cs , world_size , seed_names , singleton_names );
    }

    if( cs.singletons_threshold() > cs.threshold() ) {
      collapse_singletons( cs , seed_names , singleton_names );
    }
  } catch( ClusterInputFormatError &e ) {
    cerr << e.what() << endl;
    MPI_Finalize();
    exit( 1 );
  } catch( ClusterOutputFormatError &e ) {
    cerr << e.what() << endl;
    MPI_Finalize();
    exit( 1 );
  } catch( DACLIB::FileReadOpenError &e ) {
    cerr << e.what() << endl;
    MPI_Finalize();
    exit( 1 );
  }

  MPI_Finalize();
  exit( 0 );

}
