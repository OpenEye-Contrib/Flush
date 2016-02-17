//
// file CadSettings.cc
// David Cosgrove
// 4th February 2010
//
// This class parses the command-line arguments for program cad and holds the
// corresponding settings.

#include <iostream>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "CadSettings.H"

using namespace std;
using namespace DAC_FINGERPRINTS;
namespace po = boost::program_options;

// ***************************************************************************
CadSettings::CadSettings( int argc , char **argv ) :
  binary_file_( false ) , input_format_( FLUSH_FPS ) ,
  clus_file_format_( SAMPLES_FORMAT ) , fp_format_string_( "FLUSH_FPS" ) ,
  clus_format_string_( "SAMPLES_FORMAT" ) {

  po::options_description desc( "Allowed Options" );
  build_program_options( desc );

  po::variables_map vm;
  try {
    po::store( po::parse_command_line( argc , argv , desc ) , vm );
  } catch( po::error &e ) {
    cerr << "Error parsing command line : " << e.what() << endl
        << "cad aborts." << endl;
    exit( 1 );
  }
  po::notify( vm );

  if( argc < 2 || vm.count( "help" ) ) {
    cout << desc << endl;
    exit( 1 );
  }

  decode_formats();

  ostringstream oss;
  oss << desc;
  usage_text_ = oss.str();

}

// ***************************************************************************
bool CadSettings::operator!() const {

  if( clus_file_.empty() ) {
    error_msg_ = "No cluster file specified.";
    return true;
  } else if( fp_file_.empty() ) {
    error_msg_ = "No fingerprint file specified.";
    return true;
  } else if( output_file_.empty() ) {
    error_msg_ = "No output file specified.";
    return true;
  }

  return false;

}

// ****************************************************************************
void CadSettings::build_program_options( po::options_description &desc ) {

  desc.add_options()
    ( "help" , "Produce this help text." )
    ( "cluster-file,C" , po::value<string>( &clus_file_ ) ,
      "Name of cluster file." )
    ( "output-file,O" , po::value<string>( &output_file_ ) ,
      "Name of output file file." )
    ( "cluster-fp-file,F" , po::value<string>( &fp_file_ ) ,
      "Name of the fingerprint file for the input clusters." )
    ( "fingerprint-format" , po::value<string>( &fp_format_string_ ) ,
      "Fingerprint file format : FLUSH_FPS|BITSTRINGS|BIN_FRAG_NUMS|FRAG_NUMS (default FLUSH_FPS)" )
    ( "cluster-file-format" , po::value<string>( &clus_format_string_ ) ,
      "Clusters input format : CSV_FORMAT|SAMPLES_FORMAT (default SAMPLES_FORMAT)" )
    ( "bitstring-separator" , po::value<string>( &bitstring_separator_ ) ,
      "For bitstrings input, the separator between bits (defaults to no separator)." )
    ( "frag-num-separator" , po::value<string>( &bitstring_separator_ ) ,
      "For fragment numbers input, the separator between numbers (defaults to space)." );
  
}

// ***************************************************************************
void CadSettings::decode_formats() {

  if( fp_format_string_ == "FLUSH_FPS" ) {
    input_format_ = FLUSH_FPS;
    binary_file_ = true;
  } else if( fp_format_string_ == "BITSTRINGS" ) {
    input_format_ = BITSTRINGS;
  } else if( fp_format_string_ == "BIN_FRAG_NUMS" ) {
    input_format_ = BIN_FRAG_NUMS;
    binary_file_ = true;
  } else if( fp_format_string_ == "FRAG_NUMS" ) {
    input_format_ = FRAG_NUMS;
    if( bitstring_separator_.empty() ) {
      bitstring_separator_ = " ";
    }
  } else {
    throw FingerprintInputFormatError( fp_format_string_ );
  }

  if( clus_format_string_ == "SAMPLES_FORMAT" ) {
    clus_file_format_ = SAMPLES_FORMAT;
  } else if( clus_format_string_ == "CSV_FORMAT" ) {
    clus_file_format_ = CSV_FORMAT;
  } else {
    throw ClusterFileFormatError( clus_format_string_ );
  }

}
