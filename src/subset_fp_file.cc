//
// file subset_fp_file.cc
// David Cosgrove
// AstraZeneca
// 28th February 2007
//
// Does what it says on the tin.

#include <iostream>
#include <fstream>
#include <vector>

#include <boost/regex.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "FileExceptions.H"
#include "FingerprintBase.H"
#include "HashedFingerprint.H"
#include "NotHashedFingerprint.H"

using namespace std;
using namespace DAC_FINGERPRINTS;
using namespace boost;
namespace po = boost::program_options;

extern string BUILD_TIME;

// ***********************************************************************
void build_program_options( po::options_description &desc ,
                            string &input_fp_file ,
                            string &subset_names_file , string &output_file ,
                            string &format_string , string &bitstring_separator ,
                            bool &warm_feeling ) {

  desc.add_options()
      ( "help" , "Produce help text." )
      ( "output-file,O" , po::value<string>( &output_file ) , "Output filename" )
      ( "input-fp-file,I" , po::value<string>( &input_fp_file ) ,
        "Input filename" )
      ( "input-format,F" , po::value<string>( &format_string ) ,
        "Input format : FLUSH_FPS|BITSTRINGS|BIN_FRAG_NUMS|FRAG_NUMS (default FLUSH_FPS)" )
      ( "subset-names-file,S" , po::value<string>( &subset_names_file ) ,
        "Name of file containing names for subset." )
      ( "verbose,V" , po::value<bool>( &warm_feeling )->zero_tokens() ,
        "Verbose mode" )
      ( "warm-feeling,W" , po::value<bool>( &warm_feeling )->zero_tokens() ,
        "Verbose mode" )
      ( "bitstring-separator" , po::value<string>( &bitstring_separator ) ,
        "For bitstrings input, the separator between bits (defaults to no separator)." )
      ( "frag-num-separator" , po::value<string>( &bitstring_separator ) ,
        "For fragment numbers input, the separator between numbers (defaults to space)." );

}

// *************************************************************************
void verify_program_options( po::options_description &desc ,
                             po::variables_map &vm , int argc ,
                             bool &verbose ) {

  if( 1 == argc || vm.count( "help" ) ) {
    cout << desc << endl;
    exit( 1 );
  }

  if( !vm.count( "input-fp-file" ) ) {
    cerr << "Need an input fingerprint file." << endl << desc << endl;
    exit( 1 );
  }

  if( !vm.count( "output-file" ) ) {
    cerr << "Need an output_file." << endl << desc << endl;
    exit( 1 );
  }

  if( !vm.count( "subset-names-file" ) ) {
    cerr << "Need a subset names file." << endl << desc << endl;
    exit( 1 );
  }

  if( vm.count( "verbose" ) || vm.count( "warm-feeling" ) )
    verbose = true;

}

// *******************************************************************************
void read_subset_names( const string &subset_names_file , bool warm_feeling ,
                        vector<string> &subset_names ) {

  ifstream ifs( subset_names_file.c_str() );
  if( !ifs || !ifs.good() ) {
    cerr << "Error opening " << subset_names_file << " for reading." << endl;
  }

  string next_name;
  while( 1 ) {
    ifs >> next_name;
    if( ifs.eof() || ifs.fail() )
      break;
    subset_names.push_back( next_name );
  }

  if( warm_feeling )
    cout << "Read " << subset_names.size() << " from file " << subset_names_file
         << endl;

  bool needs_sorting = false;
  for( int i = 0 , j = 1 , is = subset_names.size() - 1 ; i < is ; ++i , ++j )
    if( subset_names[i] > subset_names[j] ) {
      needs_sorting = true;
      break;
    }

  if( needs_sorting ) {
    if( warm_feeling )
      cout << "Sorting subset names." << endl;
    sort( subset_names.begin() , subset_names.end() );
  }

}

// *************************************************************************
void open_output_file( const string &output_file , FP_FILE_FORMAT fp_file_format ,
                       gzFile &gzfp , FILE *&ucfp ) {

  boost::regex gzip( ".*\\.gz" );
  bool compr = boost::regex_match( output_file , gzip );

  if( compr ) {
    open_fp_file_for_writing( output_file , HashedFingerprint::num_ints() * sizeof( unsigned int ) ,
                              fp_file_format , gzfp );
    ucfp = 0;
  } else {
    open_fp_file_for_writing( output_file , HashedFingerprint::num_ints() * sizeof( unsigned int ) ,
                              fp_file_format , ucfp );
    gzfp = 0;
  }

}

// *************************************************************************
void write_fp_to_file( const FingerprintBase &fp , gzFile &gzfp ,
                       FILE *ucfp , FP_FILE_FORMAT fp_file_format ,
                       const string &bitstring_separator ) {

  switch( fp_file_format ) {
  case FLUSH_FPS : case BIN_FRAG_NUMS :
    if( gzfp ) {
      fp.binary_write( gzfp );
    } else {
      fp.binary_write( ucfp );
    }
    break;
  case BITSTRINGS : case FRAG_NUMS :
    if( gzfp ) {
      fp.ascii_write( gzfp , bitstring_separator );
    } else {
      fp.ascii_write( ucfp , bitstring_separator );
    }
    break;
  }

}

// *******************************************************************************
int main( int argc , char **argv ) {

  string input_fp_file , subset_names_file , output_file;
  string format_string , bitstring_separator;
  bool warm_feeling( false ) , binary_file( false );
  po::options_description desc( "Allowed Options" );
  build_program_options( desc , input_fp_file , subset_names_file ,
                         output_file , format_string , bitstring_separator ,
                         warm_feeling );

  po::variables_map vm;
  po::store( po::parse_command_line( argc , argv , desc ) , vm );
  po::notify( vm );

  verify_program_options( desc , vm , argc , warm_feeling );

  FP_FILE_FORMAT fp_file_format( FLUSH_FPS );
  if( format_string.empty() ) {
    format_string = "FLUSH_FPS";
  }

  decode_format_string( format_string , fp_file_format , binary_file ,
                        bitstring_separator );

  bool byteswapping = false;
  gzFile gzfp = 0;
  try {
    if( binary_file ) {
      open_fp_file_for_reading( input_fp_file , fp_file_format ,
                                byteswapping , gzfp );
    } else {
      open_fp_file_for_reading( input_fp_file , gzfp );
    }
  } catch( DACLIB::FileReadOpenError &e ) {
    cerr << e.what() << endl;
    exit( 1 );
  } catch( FingerprintFileError &e ) {
    cerr << e.what() << endl;
    exit( 1 );
  }

  vector<FingerprintBase *> fps;
  read_fp_file( gzfp , byteswapping , fp_file_format ,
                bitstring_separator , fps );
  if( warm_feeling ) {
    cout << "Read " << fps.size() << " fingerprints" << endl;
  }

  gzclose( gzfp );

  vector<string> subset_names;
  read_subset_names( subset_names_file , warm_feeling , subset_names );

  gzfp = 0; // for binary formats
  FILE *ucfp = 0;
  open_output_file( output_file , fp_file_format , gzfp , ucfp );

  for( int i = 0 , is = fps.size() ; i < is ; ++i ) {
    if( binary_search( subset_names.begin() , subset_names.end() ,
                       fps[i]->get_name() ) ) {
      write_fp_to_file( *fps[i] , gzfp , ucfp , fp_file_format ,
                        bitstring_separator );
    }
  }

  if( gzfp ) {
    gzclose( gzfp );
  }
  if( ucfp ) {
    fclose( ucfp );
  }

}
