//
// file AmtecSettings.cc
// David Cosgrove
// 3rd March 2009
//
// This class parses the command-line arguments for program amtec and holds the
// corresponding settings.

#include <iostream>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "AmtecSettings.H"

using namespace std;
using namespace DAC_FINGERPRINTS;
namespace po = boost::program_options;

// ***************************************************************************
AmtecSettings::AmtecSettings( int argc , char **argv ) :
  threshold_( 0.3F ) , tversky_alpha_( 0.5F ) , binary_file_( false ) ,
  input_format_( FLUSH_FPS ) , sim_calc_( TANIMOTO ) ,
  clus_input_format_( SAMPLES_FORMAT ) , clus_output_format_( SAMPLES_FORMAT ) ,
  input_format_string_( "FLUSH_FPS" ) , sim_calc_string_( "TANIMOTO" ) ,
  clus_input_format_string_( "SAMPLES_FORMAT" ) , clus_output_format_string_( "SAMPLES_FORMAT" ) {

  po::options_description desc( "Allowed Options" );
  build_program_options( desc );

  po::variables_map vm;
  try {
    po::store( po::parse_command_line( argc , argv , desc ) , vm );
  } catch( po::error &e ) {
    cerr << "Error parsing command line : " << e.what() << endl
        << "amtec aborts." << endl;
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
bool AmtecSettings::operator!() const {

  if( input_clus_file_.empty() ) {
    error_msg_ = "No existing cluster file specified.";
    return true;
  } else if( output_clus_file_.empty() ) {
    error_msg_ = "No output cluster file specified.";
    return true;
  } else if( clus_fp_file_.empty() ) {
    error_msg_ = "No fingerprint file for existing clusters specified.";
    return true;
  } else if( new_fp_file_.empty() ) {
    error_msg_ = "No file for incoming fingerprints specified.";
    return true;
  } else if( threshold_ < 0.0F || threshold_ > 1.0F ) {
    error_msg_ = string( "Invalid distance threshold " ) +
      boost::lexical_cast<string>( threshold_ ) + string( "." );
    return true;
  } else if( tversky_alpha_ < 0.0F || tversky_alpha_ > 1.0F ) {
    error_msg_ = string( "Invalid tversky_alpha " ) +
      boost::lexical_cast<string>( threshold_ ) + string( "." );
    return true;
  }

  return false;

}

// ****************************************************************************
void AmtecSettings::build_program_options( po::options_description &desc ) {

  desc.add_options()
    ( "help" , "Produce this help text." )
    ( "input-cluster-file,I" , po::value<string>( &input_clus_file_ ) ,
      "Name of existing cluster file." )
    ( "output-cluster-file,O" , po::value<string>( &output_clus_file_ ) ,
      "Name of output cluster file file." )
    ( "existing-cluster-fp-file,E" , po::value<string>( &clus_fp_file_ ) ,
      "Name of the fingerprint file for the input clusters." )
    ( "new-fingerprint-file,N" , po::value<string>( &new_fp_file_ ) ,
      "Name of new fingerprints file." )
    ( "new-fingerprint-subset" , po::value<string>( &new_subset_file_ ) ,
      "Name of file of subset of new fingerprints file.")
    ( "additions-file" , po::value<string>( &additions_file_ ) ,
      "Output file showing which clusters the new fingerprints ended up in.")
    ( "threshold,T" , po::value<double>( &threshold_ ) ,
      "Clustering threshold (default 0.3)" )
    ( "input-format,F" , po::value<string>( &input_format_string_ ) ,
      "Input format : FLUSH_FPS|BITSTRINGS|BIN_FRAG_NUMS|FRAG_NUMS (default FLUSH_FPS)" )
    ( "clus-output-format" , po::value<string>( &clus_output_format_string_ ) ,
      "Clusters output format : CSV_FORMAT|SAMPLES_FORMAT (default SAMPLES_FORMAT)" )
    ( "clus-input-format" , po::value<string>( &clus_input_format_string_ ) ,
      "Clusters input format : CSV_FORMAT|SAMPLES_FORMAT (default SAMPLES_FORMAT)" )
    ( "distance-calculation" , po::value<string>( &sim_calc_string_ ) ,
      "Distance calculation : TANIMOTO|TVERSKY (default TANIMOTO)" )
    ( "tversky-alpha" , po::value<double>( &tversky_alpha_ ) ,
      "Tversky alpha paramter (0.0-1.0, default 0.5" )
    ( "bitstring-separator" , po::value<string>( &bitstring_separator_ ) ,
      "For bitstrings input, the separator between bits (defaults to no separator)." )
    ( "frag-num-separator" , po::value<string>( &bitstring_separator_ ) ,
      "For fragment numbers input, the separator between numbers (defaults to space)." );
  
}

// ***************************************************************************
void AmtecSettings::decode_formats() {

  if( input_format_string_ == "FLUSH_FPS" ) {
    input_format_ = FLUSH_FPS;
    binary_file_ = true;
  } else if( input_format_string_ == "BITSTRINGS" ) {
    input_format_ = BITSTRINGS;
  } else if( input_format_string_ == "BIN_FRAG_NUMS" ) {
    input_format_ = BIN_FRAG_NUMS;
    binary_file_ = true;
  } else if( input_format_string_ == "FRAG_NUMS" ) {
    input_format_ = FRAG_NUMS;
    if( bitstring_separator_.empty() ) {
      bitstring_separator_ = " ";
    }
  } else {
    throw FingerprintInputFormatError( input_format_string_ );
  }

  if( clus_input_format_string_ == "SAMPLES_FORMAT" ) {
    clus_input_format_ = SAMPLES_FORMAT;
  } else if( clus_input_format_string_ == "CSV_FORMAT" ) {
    clus_input_format_ = CSV_FORMAT;
  } else {
    throw ClusterFileFormatError( clus_input_format_string_ );
  }

  if( clus_output_format_string_ == "SAMPLES_FORMAT" ) {
    clus_output_format_ = SAMPLES_FORMAT;
  } else if( clus_output_format_string_ == "CSV_FORMAT" ) {
    clus_output_format_ = CSV_FORMAT;
  } else {
    throw ClusterFileFormatError( clus_output_format_string_ );
  }

  if( sim_calc_string_ == "TANIMOTO" ) {
    sim_calc_ = DAC_FINGERPRINTS::TANIMOTO;
  } else if( sim_calc_string_ == "TVERSKY" ) {
    sim_calc_ = DAC_FINGERPRINTS::TVERSKY;
  } else {
    throw FingerprintDistCalcError( sim_calc_string_ );
  }

}
