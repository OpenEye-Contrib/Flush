//
// file satan.cc - So Are There Any Neighbours?
// David Cosgrove
// AstraZeneca
// 16th February 2006
//
// Takes 2 fingerprint files and checks to see if there are any fingerprints
// in the first that have at least a given number of fingerprints in the second
// within a threshold tanimoto distance.

#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric> // for the accumulate algorithm
#include <sstream>
#include <vector>

#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/scoped_ptr.hpp>

#include "FileExceptions.H"
#include "FingerprintBase.H"
#include "HashedFingerprint.H"
#include "NotHashedFingerprint.H"
#include "SatanSettings.H"
#include "chrono.h"

#include <mpi.h>

using namespace boost;
using namespace std;
using namespace DAC_FINGERPRINTS;

class SortNbsByDist :
    public binary_function<pair<string,double> , pair<string,double> , bool> {
public :
  result_type operator()( first_argument_type a , second_argument_type b ) const {
    if( a.second == b.second )
      return a.first < b.first;
    else
      return a.second < b.second;
  }
};

// in eponymous file
namespace DACLIB {
string get_cwd();

// In mpi_string_subs.cc
void mpi_send_string( const string &str , int dest_rank );
void mpi_rec_string( int source_rank , string &str );
}

extern string BUILD_TIME; // in build_time.cc

// static const int FP_CHUNK_SIZE = 500000;

// ****************************************************************************
void output_neighbours_satan( unsigned int min_count ,
                              ostream &output_stream ,
                              vector<pair<string,vector<pair<string,double> > > > &nbs ) {

  for( int i = 0 , is = nbs.size() ; i < is ; ++i ) {
    if( min_count && (nbs[i].second).size() >= min_count ) {
      for( unsigned int j = 0 ; j < min_count ; ++j ) {
        output_stream << nbs[i].first << " "
                      << (nbs[i].second)[j].first << " "
                      << (nbs[i].second)[j].second << endl;
      }
    }
  }

}

// ****************************************************************************
void output_neighbours_satan( ostream &output_stream ,
                              vector<pair<string,vector<pair<string,double> > > > &nbs ) {

  for( int i = 0 , is = nbs.size() ; i < is ; ++i ) {
    for( unsigned int j = 0 , js = (nbs[i].second).size() ; j < js ; ++j ) {
      output_stream << nbs[i].first << " "
                    << (nbs[i].second)[j].first << " "
                    << (nbs[i].second)[j].second << endl;
    }
  }

}

// ****************************************************************************
unsigned int max_probe_name_len( const vector<pair<string,vector<pair<string,double> > > > &nbs ) {

  unsigned int max_len = 0;
  for( int i = 0 , is = nbs.size() ; i < is ; ++i ) {
    if( (nbs[i].first).length() > max_len ) {
      max_len = (nbs[i].first).length();
    }
  }

  return max_len;

}

// ****************************************************************************
unsigned int max_target_name_len( const vector<pair<string,vector<pair<string,double> > > > &nbs ) {

  unsigned int max_len = 0;
  for( int i = 0 , is = nbs.size() ; i < is ; ++i ) {
    for( int k = 0 , ks = (nbs[i].second).size() ; k < ks ; ++k ) {
      if( (nbs[i].second)[k].first.length() > max_len ) {
        max_len = (nbs[i].second)[k].first.length();
      }
    }
  }

  return max_len;

}

// ****************************************************************************
unsigned int max_count_name_len( const vector<pair<string,vector<unsigned int> > > &counts ) {

  unsigned int max_len = 0;
  typedef pair<string,vector<unsigned int> > SVUI;
  BOOST_FOREACH( SVUI cnt , counts ) {
    if( cnt.first.length() > max_len ) {
      max_len = cnt.first.length();
    }
  }

  return max_len;
}

// ****************************************************************************
void pad_spaces( unsigned int str_len , unsigned int max_len ,
                 bool with_colon , ostream &output_stream ) {

  // this is the way that snailflush does it - I can't remember why, but I
  // know there ought to be a better way.
  for( unsigned int i = str_len ; i <= max_len ; ++i ) {
    output_stream << " ";
  }
  if( with_colon ) {
    output_stream << "  : ";
  }

}

// ****************************************************************************
void output_neighbours_nnlists( unsigned int min_count ,
                                ostream &output_stream ,
                                vector<pair<string,vector<pair<string,double> > > > &nbs ) {

  // snailflush pads the nnlists with spaces so all records are the same length.
  // Do the same here for consistency, but bearing in mind that it might
  // not be completely kosher as the probes are being done in batches.
  unsigned int max_probe_len = max_probe_name_len( nbs );
  unsigned int max_target_len = max_target_name_len( nbs );

  for( int i = 0 , is = nbs.size() ; i < is ; ++i ) {
    if( min_count && (nbs[i].second).size() >= min_count ) {
      output_stream << nbs[i].first;
      pad_spaces( nbs[i].first.length() , max_probe_len , true , output_stream );
      output_stream << endl;

      for( unsigned int j = 0 ; j < min_count ; ++j ) {
        // for nnlists, by tradition we don't output the neighbour if it appears
        // to be the same molecule
        if( (nbs[i].second)[j].second == 0.0F &&
            (nbs[i].second)[j].first == nbs[i].first ) {
          continue;
        }
        output_stream << "      " << (nbs[i].second)[j].first;
        pad_spaces( (nbs[i].second)[j].first.length() , max_target_len ,
                    true , output_stream );
        output_stream << setw( 6 ) << setprecision( 4 )
                      << setiosflags( ios::showpoint )
                      << (nbs[i].second)[j].second << endl;
      }
      output_stream << endl;
    }
  }

}

// ****************************************************************************
void output_neighbours_nnlists( ostream &output_stream ,
                                vector<pair<string,vector<pair<string,double> > > > &nbs ) {

  // snailflush pads the nnlists with spaces so all records are the same length.
  // Do the same here for consistency, but bearing in mind that it might
  // not be completely kosher as the probes are being done in batches.
  unsigned int max_probe_len = max_probe_name_len( nbs );
  unsigned int max_target_len = max_target_name_len( nbs );

  for( int i = 0 , is = nbs.size() ; i < is ; ++i ) {
    output_stream << nbs[i].first;
    pad_spaces( nbs[i].first.length() , max_probe_len , true , output_stream );
    output_stream << endl;

    for( unsigned int j = 0 , js = (nbs[i].second).size() ; j < js ; ++j ) {
      // for nnlists, by tradition we don't output the neighbour if it appears
      // to be the same molecule
      if( (nbs[i].second)[j].second == 0.0F &&
          (nbs[i].second)[j].first == nbs[i].first ) {
        continue;
      }
      output_stream << "      " << (nbs[i].second)[j].first;
      pad_spaces( (nbs[i].second)[j].first.length() , max_target_len ,
                  true , output_stream );
      output_stream << setw( 6 ) << setprecision( 4 )
                    << setiosflags( ios::showpoint )
                    << (nbs[i].second)[j].second << endl;
    }
    output_stream << endl;
  }

}

// ****************************************************************************
void output_neighbours( unsigned int min_count , const string &output_format ,
                        ostream &output_stream ,
                        vector<pair<string,vector<pair<string,double> > > > &nbs ) {

  if( min_count ) {
    if( string( "SATAN" ) == output_format ) {
      output_neighbours_satan( min_count , output_stream , nbs );
    } else {
      output_neighbours_nnlists( min_count , output_stream , nbs );
    }
  } else {
    if( string( "SATAN" ) == output_format ) {
      output_neighbours_satan( output_stream , nbs );
    } else {
      output_neighbours_nnlists( output_stream , nbs );
    }
  }

}

// ****************************************************************************
void output_counts( ostream &output_stream ,
                    vector<pair<string,vector<unsigned int> > > &counts ) {

  // snailflush pads the nnlists with spaces so all records are the same length.
  // Do the same here for consistency, but bearing in mind that it might
  // not be completely kosher as the probes are being done in batches.
  unsigned int max_name_len = max_count_name_len( counts );

  for( int i = 0 , is = counts.size() ; i < is ; ++i ) {
    output_stream << counts[i].first;
    pad_spaces( counts[i].first.length() , max_name_len + 3 ,
                false , output_stream );
    int sum_count = 0;
    for( int j = 0 ; j < 10 ; ++j ) {
      sum_count += counts[i].second[j];
      output_stream << setw( 6 ) << sum_count << "  ";
    }
    output_stream << endl;
  }

}

// ****************************************************************************
void dump_fps( vector<FingerprintBase *> &fps ) {

  for( int i = 0 , is = fps.size() ; i < is ; ++i ) {
    delete fps[i];
  }
  fps.clear();

}

// ****************************************************************************
void open_fp_file( const string &filename ,
                   DAC_FINGERPRINTS::FP_FILE_FORMAT input_format ,
                   bool &byteswapping  , gzFile &pfile) {

  byteswapping = false;
  try {
    open_fp_file_for_reading( filename , input_format ,
                              byteswapping , pfile );
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

// ****************************************************************************
void target_against_probes( const FingerprintBase *target_fp ,
                            const vector<FingerprintBase *> &probe_fps ,
                            double threshold , unsigned int min_count ,
                            vector<pair<string,vector<pair<string,double> > > > &nbs ) {

  for( int i = 0 , is = probe_fps.size() ; i < is ; ++i ) {
    if( !min_count || nbs[i].second.size() < min_count ) {
      double dist = target_fp->calc_distance( *(probe_fps[i] ) , threshold );
      if( dist <= threshold ) {
        nbs[i].second.push_back( make_pair( target_fp->get_name() , dist ) );
      }
    }
  }

}

// ****************************************************************************
// the counts version.  If dist is 0.44, then counts[4] will be incremented
void target_against_probes( const FingerprintBase *target_fp ,
                            const vector<FingerprintBase *> &probe_fps ,
                            vector<pair<string,vector<unsigned int> > > &counts ) {

  for( int i = 0 , is = probe_fps.size() ; i < is ; ++i ) {
    // traditionally, we don't report the compound with itself, even though
    // the test is going to slow things down badly.
    if( probe_fps[i]->get_name() != target_fp->get_name() ) {
      double dist = 10.0 * target_fp->calc_distance( *(probe_fps[i]) );
      int cbin = int( dist );
      cbin = 10 == cbin ? 9 : cbin;
      // bins go to <= dist, so on the border is in the previous bin
      if( cbin && fabs( double( cbin ) - dist ) < 1.0e-16 ) {
        --cbin;
      }
      ++(counts[i].second[cbin]);
    }
  }

}

// ****************************************************************************
void process_fingerprints( const SatanSettings &ss ,
                           unsigned int num_probe_fps , int chunk_num ,
                           vector<pair<string,vector<pair<string,double> > > > &nbs ,
                           vector<pair<string,vector<unsigned int> > > &counts ) {

  gzFile pfile , tfile;
  bool target_byteswapping , probe_byteswapping;

  // read next lot of probe fps
  open_fp_file( ss.probe_file() , ss.input_format() , probe_byteswapping , pfile );
  vector<FingerprintBase *> probe_fps;
  unsigned int start_probe_fp = num_probe_fps * chunk_num;
  read_fps_from_file( pfile , probe_byteswapping , ss.input_format() ,
                      ss.bitstring_separator() , start_probe_fp ,
                      num_probe_fps , probe_fps );

  if( probe_fps.empty() ) {
    cerr << "Error : premature end of file " << ss.probe_file() << endl;
    gzclose( pfile );
    exit( 1 );
  }
  if( ss.warm_feeling() ) {
    cout << "Read " << probe_fps.size() << " probes." << endl;
  }

  if( string( "COUNTS" ) == ss.output_format() ) {
    counts.reserve( probe_fps.size() );
    BOOST_FOREACH( FingerprintBase *pfp , probe_fps ) {
      counts.push_back( make_pair( pfp->get_name() , vector<unsigned int>( 10 , 0 ) ) );
    }
  } else {
    nbs.reserve( probe_fps.size() );
    BOOST_FOREACH( FingerprintBase *pfp , probe_fps ) {
      nbs.push_back( make_pair( pfp->get_name() , vector<pair<string,double> >() ) );
    }
  }

  open_fp_file( ss.target_file() , ss.input_format() , target_byteswapping , tfile );

  int num_targets = 0;
  bool counts_output = string( "COUNTS" ) == ss.output_format() ? true : false;

  while( 1 ) {
    FingerprintBase *target_fp = read_next_fp_from_file( tfile , target_byteswapping ,
                                                         ss.input_format() ,
                                                         ss.bitstring_separator() );
    if( !target_fp ) {
      break;
    }
    ++num_targets;

    if( counts_output ) {
      target_against_probes( target_fp , probe_fps , counts );
    } else {
      target_against_probes( target_fp , probe_fps , ss.threshold() ,
                             ss.min_count() , nbs );
    }
    delete target_fp;
  }

  gzclose( pfile );
  gzclose( tfile );

  // we're done with probes
  dump_fps( probe_fps );

  // sort the neighbour lists ready for output
  for( int i = 0 , is = nbs.size() ; i < is ; ++i ) {
    sort( nbs[i].second.begin() , nbs[i].second.end() , SortNbsByDist() );
  }

}

// ****************************************************************************
void serial_run( const SatanSettings &ss ) {

  // open the output stream right away, in case we can't. It's best to find
  // out before we've done a potentially long job.
  ofstream output_stream( ss.output_file().c_str() ) ;
  if( !output_stream.good() ) {
    cerr << "Couldn't open " << ss.output_file() << " for writing." << endl;
    exit( 1 );
  }

    vector<pair<string,vector<pair<string,double> > > > nbs;
    vector<pair<string,vector<unsigned int> > > counts;
    process_fingerprints( ss , numeric_limits<unsigned int>::max() , 0 ,
                          nbs , counts );
    if( !nbs.empty() ) {
      output_neighbours( ss.min_count() , ss.output_format() , output_stream , nbs );
    }
    if( !counts.empty( ) ){
      output_counts( output_stream , counts );
    }

}

// ****************************************************************************
void tell_master_slave_has_done_nnlists() {

  // tell the master we're finished.  It'll then pick the results up in
  // the correct order.
  DACLIB::mpi_send_string( string( "Slave_Finished" ) , 0 );

}

// ****************************************************************************
void send_results_to_master( int chunk_num ,
                             vector<pair<string,vector<pair<string,double> > > > &nbs ) {

  MPI_Send( &chunk_num , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD );
  unsigned int num_to_send = nbs.size();
  MPI_Send( &num_to_send , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD );

  for( unsigned int i = 0 , is = nbs.size() ; i < is ; ++i ) {
    DACLIB::mpi_send_string( nbs[i].first , 0 );
    unsigned int j = (nbs[i].second).size();
    MPI_Send( &j , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD );
    for( int k = 0 , ks = (nbs[i].second).size() ; k < ks ; ++k ) {
      DACLIB::mpi_send_string( (nbs[i].second)[k].first , 0 );
      MPI_Send( &(nbs[i].second)[k].second , 1 , MPI_DOUBLE , 0 , 0 , MPI_COMM_WORLD );
    }
  }

}

// ****************************************************************************
void send_results_to_master( int chunk_num ,
                             vector<pair<string,vector<unsigned int> > > &counts ) {

  MPI_Send( &chunk_num , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD );
  unsigned int num_to_send = counts.size();
  MPI_Send( &num_to_send , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD );

  for( unsigned int i = 0 , is = counts.size() ; i < is ; ++i ) {
    DACLIB::mpi_send_string( counts[i].first , 0 );
    MPI_Send( &counts[i].second[0] , 10 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD );
  }

}

// ****************************************************************************
void wait_till_all_slaves_done( bool warm_feeling , int world_size ) {

  int slaves_running = world_size - 1;
  while( slaves_running > 0 ) {
    MPI_Status status;
    MPI_Probe( MPI_ANY_SOURCE , 0 , MPI_COMM_WORLD , &status );

    string msg;
    DACLIB::mpi_rec_string( status.MPI_SOURCE , msg );
    if( string( "Slave_Finished" ) != msg ) {
      cerr << "Error, expected message Slave_Finished from slave, but got "
           << msg << ". Can't go on." << endl;
      MPI_Finalize();
      exit( 1 );
    }
    --slaves_running;
    if( warm_feeling ) {
      cout << "Slave " << status.MPI_SOURCE << " is finished.  " << slaves_running
           << " still running." << endl;
    }
  }

}

// ****************************************************************************
void receive_slave_counts_results( int world_size ,
                                   ostream &output_stream ) {

  vector<pair<string,vector<unsigned int> > > counts;
  // get the results from the slaves in order, because that's important
  // to keep the output in probe input order.  The chunks have been sent off
  // in slave_tids order.
  for( int i = 1 ; i < world_size ; ++i ) {
    DACLIB::mpi_send_string( string( "Send_Results" ) , i );

    int chunk_num;
    MPI_Recv( &chunk_num , 1 , MPI_INT , i , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
    unsigned int num_to_rec;
    MPI_Recv( &num_to_rec , 1 , MPI_UNSIGNED , i , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
    counts.resize( num_to_rec );
    for( unsigned int k = 0 ; k < num_to_rec ; ++k ) {
      DACLIB::mpi_rec_string( i , counts[k].first );
      counts[k].second = vector<unsigned int>( 10 , 0 );
      MPI_Recv( &counts[k].second[0] , 10 , MPI_UNSIGNED , i , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
    }
    output_counts( output_stream , counts );
    counts.clear();
  }

}

// ****************************************************************************
void receive_slave_nnlists_results( int world_size , int min_count ,
                                    const string &output_format ,
                                    ostream &output_stream ) {

  vector<pair<string,vector<pair<string,double> > > > nbs;
  // get the results from the slaves in order, because that's important
  // to keep the output in probe input order.  The chunks have been sent off
  // in slave_tids order.
  for( int i = 1 ; i < world_size ; ++i ) {
    DACLIB::mpi_send_string( string( "Send_Results" ) , i );

    int chunk_num;
    unsigned int num_to_rec , num_nbs;
    string target_name;
    double target_dist;
    MPI_Recv( &chunk_num , 1 , MPI_INT , i , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
    MPI_Recv( &num_to_rec , 1 , MPI_UNSIGNED , i , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
    nbs.resize( num_to_rec );
    for( unsigned int k = 0 ; k < num_to_rec ; ++k ) {
      DACLIB::mpi_rec_string( i , nbs[k].first );
      MPI_Recv( &num_nbs , 1 , MPI_UNSIGNED , i , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
      nbs[i].second.reserve( num_nbs );
      for( unsigned int j = 0 ; j < num_nbs ; ++j ) {
        DACLIB::mpi_rec_string( i , target_name );
        MPI_Recv( &target_dist , 1 , MPI_DOUBLE , i , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
        nbs[k].second.push_back( make_pair( target_name , target_dist ) );
      }
    }
    output_neighbours( min_count , output_format , output_stream , nbs );
    nbs.clear();
  }

}

// ****************************************************************************
void receive_slave_results( bool warm_feeling , int world_size ,
                            int min_count , const string &output_format ,
                            ostream &output_stream ) {

  // first, wait for all the slaves to report they're done. Doing it this way
  // will allow us to spot dead slaves early and take appropriate action.
  wait_till_all_slaves_done( warm_feeling , world_size );

  if( string( "COUNTS" ) == output_format ) {
    receive_slave_counts_results( world_size , output_stream );
  } else {
    receive_slave_nnlists_results( world_size , min_count ,
                                   output_format , output_stream );
  }

}

// ****************************************************************************
void send_search_details( SatanSettings &ss , unsigned int num_probe_fps ,
                          int world_size ) {

  // each slave needs to do something like num_probe_fps / num_procs probe fps
  // each against all targets
  unsigned int slave_does = num_probe_fps / ( world_size - 1 );
  while( slave_does * ( world_size - 1 ) < num_probe_fps ) {
    ++slave_does;
  }
  cout << "Each slave does " << slave_does << " probe fps " << endl;

  for( int i = 1 ; i < world_size ; ++i ) {
    DACLIB::mpi_send_string( string( "Search_Details" ) , i );

    ss.send_contents_via_mpi( i );

    // send the number of fps each slave must do, and the slave number,
    // so it knows where to start
    MPI_Send( &slave_does , 1 , MPI_INT , i , 0 , MPI_COMM_WORLD );
    MPI_Send( &i , 1 , MPI_INT , i , 0 , MPI_COMM_WORLD );
  }

}

// ****************************************************************************
void receive_search_details( SatanSettings &ss ,
                             unsigned int &num_probe_fps_to_do ,
                             int &chunk_num ) {

  ss.receive_contents_via_mpi();
  if( TVERSKY == ss.similarity_calc() ) {
    FingerprintBase::set_tversky_alpha( ss.tversky_alpha() );
    HashedFingerprint::set_similarity_calc( ss.similarity_calc() );
    NotHashedFingerprint::set_similarity_calc( ss.similarity_calc() );
  }

  MPI_Recv( &num_probe_fps_to_do , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &chunk_num , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  --chunk_num; // chunk_num is based on the MPI rank of the slave, which starts from 1

}

// ****************************************************************************
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

// ****************************************************************************
void receive_new_cwd() {

  string new_cwd;
  DACLIB::mpi_rec_string( 0 , new_cwd );
  if(chdir( new_cwd.c_str() ) < 0){
    cerr << "ERROR : couldn't change to directory " << new_cwd << endl;
    exit(1);
  }

}

// ****************************************************************************
void parallel_run( SatanSettings &ss , int world_size ) {

  // open the output stream right away, in case we can't. It's best to find
  // out before we've done a potentially long job.
  ofstream output_stream( ss.output_file().c_str() ) ;
  if( !output_stream.good() ) {
    cerr << "Couldn't open " << ss.output_file() << " for writing." << endl;
    exit( 1 );
  }

  unsigned int num_probe_fps = 0;
  try {
    num_probe_fps = count_fps_in_file( ss.probe_file() , ss.input_format() ,
                                       ss.bitstring_separator() );
  } catch( DACLIB::FileReadOpenError &e ) {
    cerr << e.what() << endl;
    cout << e.what() << endl;
    MPI_Finalize();
    exit( 1 );
  } catch( FingerprintFileError &e ) {
    cerr << e.what() << endl;
    cout << e.what() << endl;
    MPI_Finalize();
    exit( 1 );
  }

  if( num_probe_fps ) {
    send_cwd_to_slaves( world_size );
    // send_search_details also fires off the jobs on the slaves
    send_search_details( ss , num_probe_fps , world_size );
    // get the results and write directly to file. This way, we don't ever
    // have to hold the whole, potentially enormous, neighbour list in
    // memory
    receive_slave_results( ss.warm_feeling() , world_size , ss.min_count() ,
                           ss.output_format() , output_stream );
  }

  for( int i = 1 ; i < world_size ; ++i ) {
    DACLIB::mpi_send_string( string( "Finished" ) , i );
  }

}

// ****************************************************************************
void slave_event_loop() {

#if DEBUG == 1
  int wr;
  MPI_Comm_rank( MPI_COMM_WORLD , &wr );
  cout << "top of slave_event_loop for " << wr << endl;
#endif

  SatanSettings ss;
  unsigned int num_probe_fps_to_do = 0;
  int chunk_num = 0;
  vector<pair<string,vector<pair<string,double> > > > nbs;
  vector<pair<string,vector<unsigned int> > > counts;

  while( 1 ) {
    
    string msg;
    DACLIB::mpi_rec_string( 0 , msg );
#if DEBUG == 1
    cout << "received message : " << msg << endl;
#endif
    if( string( "Finished" ) == msg ) {
      break;
    } else if( string( "Search_Details" ) == msg ) {
      receive_search_details( ss , num_probe_fps_to_do , chunk_num );
      process_fingerprints( ss , num_probe_fps_to_do , chunk_num ,  nbs , counts );
      tell_master_slave_has_done_nnlists();
    } else if( string( "Send_Results" ) == msg ) {
      if( string( "COUNTS" ) == ss.output_format() ) {
        send_results_to_master( chunk_num , counts );
      } else {
        send_results_to_master( chunk_num , nbs );
      }
    } else if( string( "New_CWD" ) == msg  ) {
      receive_new_cwd();
    } else {
      int world_rank;
      MPI_Comm_rank( MPI_COMM_WORLD , &world_rank );
      cout << world_rank << " received suspect message " << msg << endl;
    }

  }

}

// ****************************************************************************
int main( int argc , char **argv ) {

  cout << "satan - built " << BUILD_TIME << endl;

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

  SatanSettings ss( argc , argv );
  if( !ss ) {
    cout << "ERROR : " << ss.error_message() << endl << ss.usage_text() << endl;
    cerr << "ERROR : " << ss.error_message() << endl << ss.usage_text() << endl;
  }

  if( TVERSKY == ss.similarity_calc() ) {
    FingerprintBase::set_tversky_alpha( ss.tversky_alpha() );
    HashedFingerprint::set_similarity_calc( ss.similarity_calc() );
    NotHashedFingerprint::set_similarity_calc( ss.similarity_calc() );
  }

  if( 1 == world_size ) {
    serial_run( ss );
  } else {
    parallel_run( ss , world_size );
  }

  MPI_Finalize();
  exit( 0 );

}
