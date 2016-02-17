//
// file merge_fp_files.cc
// David Cosgrove
// AstraZeneca
// 16th May 2007
//
// Combines 2 or more fp files into a new one.

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
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
                            vector<string> &input_files , string &output_file ,
                            string &input_format_string , string &output_format_string , 
			    string &bitstring_separator , bool &warm_feeling ) {

  desc.add_options()
      ( "help" , "Produce help text." )
      ( "output-file,O" , po::value<string>( &output_file ) , "Output filename" )
      ( "input-file,I" , po::value<vector<string> >( &input_files ) ,
        "Input filename" )
      ( "input-format" , po::value<string>( &input_format_string ) ,
        "Input format : FLUSH_FPS|BITSTRINGS|BIN_FRAG_NUMS|FRAG_NUMS (default FLUSH_FPS)" )
      ( "output-format" , po::value<string>( &output_format_string ) ,
        "Output format : FLUSH_FPS|BITSTRINGS|BIN_FRAG_NUMS|FRAG_NUMS (default FLUSH_FPS)" )
      ( "verbose,V" , po::value<bool>( &warm_feeling )->zero_tokens() ,
        "Verbose mode" )
      ( "warm-feeling" , po::value<bool>( &warm_feeling )->zero_tokens() ,
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

  if( !vm.count( "input-file" ) ) {
    cerr << "Need an input_file." << endl << desc << endl;
    exit( 1 );
  }

  if( !vm.count( "output-file" ) ) {
    cerr << "Need an output_file." << endl << desc << endl;
    exit( 1 );
  }

  if( vm.count( "verbose" ) || vm.count( "warm-feeling" ) )
    verbose = true;

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
void write_fps_to_file( gzFile &gzfp , FILE *ucfp , FP_FILE_FORMAT fp_file_format ,
                        const string &bitstring_separator ,
                        const vector<FingerprintBase *> &fps ) {

  for( int i = 0 , is = fps.size() ; i < is ; ++i ) {

    switch( fp_file_format ) {
    case FLUSH_FPS : case BIN_FRAG_NUMS :
      if( gzfp ) {
        fps[i]->binary_write( gzfp );
      } else {
        fps[i]->binary_write( ucfp );
      }
      break;
    case BITSTRINGS : case FRAG_NUMS :
      if( gzfp ) {
        fps[i]->ascii_write( gzfp , bitstring_separator );
      } else {
        fps[i]->ascii_write( ucfp , bitstring_separator );
      }
      break;
    }
  }

}

// *************************************************************************
int main( int argc , char **argv ) {

  cout << "merge_fp_files - built " << BUILD_TIME << endl;

  vector<string> input_files;
  string input_format_string( "FLUSH_FPS" );
  string output_file , output_format_string( "FLUSH_FPS" );
  string bitstring_separator;
  bool   warm_feeling( false ) , binary_file( false );

  po::options_description desc( "Allowed Options" );
  build_program_options( desc , input_files , output_file ,
			 input_format_string , output_format_string ,
                         bitstring_separator , warm_feeling );

  po::variables_map vm;
  po::store( po::parse_command_line( argc , argv , desc ) , vm );
  po::notify( vm );

  verify_program_options( desc , vm , argc , warm_feeling );

  FP_FILE_FORMAT in_fp_file_format( FLUSH_FPS );
  decode_format_string( input_format_string , in_fp_file_format ,
			binary_file , bitstring_separator );
  FP_FILE_FORMAT out_fp_file_format( FLUSH_FPS );
  decode_format_string( output_format_string , out_fp_file_format ,
			binary_file , bitstring_separator );

  gzFile gzfp = 0; // for binary formats
  FILE *ucfp = 0;

  int num_fps_read = 0;
  for( int i = 0 , is = input_files.size() ; i < is ; ++i ) {
    vector<FingerprintBase *> next_fps;
    if( warm_feeling ) {
      cout << "Reading fingerprint file " << input_files[i] << endl;
    }
    try {
      read_fp_file( input_files[i] , in_fp_file_format , bitstring_separator ,
                    next_fps );
    } catch( DACLIB::FileReadOpenError &e ) {
      cerr << e.what() << endl;
      cout << e.what() << endl;
      exit( 1 );
    } catch( FingerprintFileError &e ) {
      cerr << e.what() << endl;
      cout << e.what() << endl;
      exit( 1 );
    } catch( DAC_FINGERPRINTS::HashedFingerprintLengthError &e ) {
      cerr << e.what() << endl;
      cout << e.what() << endl;
      exit( 1 );
    }

    num_fps_read += next_fps.size();
    if( warm_feeling ) {
      cout << "Read " << next_fps.size() << " fingerprint";
      if( next_fps.size() > 1 )
        cout << "s";
      cout << " from file number " << i + 1 << " : " << input_files[i] << endl;
    }
    // delay opening file until we've read a fingerprint and know how
    // many chars there are in it.
    if( !gzfp && !ucfp ) {
      open_output_file( output_file , out_fp_file_format , gzfp , ucfp );
    }

    write_fps_to_file( gzfp , ucfp , out_fp_file_format , bitstring_separator ,
                       next_fps );
    for( int j = 0 , js = next_fps.size() ; j < js ; ++j ) {
      delete next_fps[j];
    }
  }

  if( warm_feeling ) {
    cout << "Written " << num_fps_read << " fingerprint";
    if( num_fps_read > 1 ) cout << "s";
    cout << " from " << input_files.size() << " file";
    if( input_files.size() > 1 ) cout << "s";
    cout << "." << endl;
  }

  if( gzfp ) {
    gzclose( gzfp );
  }
  if( ucfp ) {
    fclose( ucfp );
  }

}
