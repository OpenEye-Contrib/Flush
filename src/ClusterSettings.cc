//
// file ClusterSettings.cc
// David Cosgrove
// AstraZeneca
// 23rd February 2009
//

#include <iostream>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <mpi.h>

#include "ClusterSettings.H"

using namespace std;
using namespace DAC_FINGERPRINTS;
namespace po = boost::program_options;

namespace DACLIB {
// In mpi_string_subs.cc
void mpi_send_string( const string &str , int dest_rank );
void mpi_rec_string( int source_rank , string &str );
}

// ****************************************************************************
ClusterSettings::ClusterSettings( int argc , char **argv ) :
  threshold_( 0.3 ) , singletons_threshold_( -1.0 ) , warm_feeling_( false ) ,
  output_format_string_( "SAMPLES_FORMAT" ) ,
  input_format_string_( "FLUSH_FPS" ) ,
  output_format_( SAMPLES_FORMAT ) , input_format_( FLUSH_FPS ) ,
  binary_file_( false ) , fix_spaces_in_names_( false ) {

  po::options_description desc( "Allowed Options" );
  build_program_options( desc );

  po::variables_map vm;
  try {
    po::store( po::parse_command_line( argc , argv , desc ) , vm );
  }  catch( po::error &e ) {
    cerr << "Error parsing command line : " << e.what() << endl
        << "cluster aborts." << endl;
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
bool ClusterSettings::operator!() const {

  if( input_file_.empty() ) {
    error_msg_ = "No input file specified.";
    return true;
  } else if( output_file_.empty() ) {
    error_msg_ = "No output file specified.";
    return true;
  } else if( threshold_ < 0.0 || threshold_ > 1.0 ) {
    error_msg_ = string( "Invalid distance threshold " ) +
      boost::lexical_cast<string>( threshold_ ) + string( "." );
    return true;
  }

  return false;

}

// ****************************************************************************
void ClusterSettings::send_contents_via_mpi( int dest_slave ) {

  using namespace DACLIB;
  mpi_send_string( input_file_ , dest_slave );
  mpi_send_string( output_file_ , dest_slave );
  mpi_send_string( subset_file_ , dest_slave );

  MPI_Send( &threshold_ , 1 , MPI_DOUBLE , dest_slave , 0 , MPI_COMM_WORLD );
  int i( warm_feeling_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );

  mpi_send_string( input_format_string_ , dest_slave );
  i = int( output_format_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = int( input_format_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = int( binary_file_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  mpi_send_string( bitstring_separator_ , dest_slave );
  i = int( fix_spaces_in_names_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );

}

// ****************************************************************************
void ClusterSettings::receive_contents_via_mpi() {

  using namespace DACLIB;

  mpi_rec_string( 0 , input_file_ );
  mpi_rec_string( 0 , output_file_ );
  mpi_rec_string( 0 , subset_file_ );

  MPI_Recv( &threshold_ , 1 , MPI_DOUBLE , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  int i = 0;
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  warm_feeling_ = static_cast<bool>( i );

  mpi_rec_string( 0 , input_format_string_ );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  output_format_ = static_cast<OUTPUT_FORMAT>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  input_format_ = static_cast<DAC_FINGERPRINTS::FP_FILE_FORMAT>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  binary_file_ = static_cast<bool>( i );

  mpi_rec_string( 0 , bitstring_separator_ );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  fix_spaces_in_names_ = static_cast<bool>( i );

}

// ****************************************************************************
void ClusterSettings::build_program_options( po::options_description &desc ) {

  desc.add_options()
    ( "help" , "Produce this help text." )
    ( "input-file,I" , po::value<string>( &input_file_ ) ,
      "Name of input fingerprint file." )
    ( "output-file,O" , po::value<string>( &output_file_ ) ,
      "Name of the output clusters file." )
    ( "subset-file,S" , po::value<string>( &subset_file_ ) ,
      "File containing names of fingerprints giving subset to be used in clustering.")
    ( "threshold,T" , po::value<double>( &threshold_ ) ,
      "Clustering threshold (default 0.3)" )
      ( "singletons-threshold" , po::value<double>( &singletons_threshold_ ) ,
        "Threshold for collapsing singletons. Defaults to -1.0, no collapse." )
    ( "warm-feeling,W" , po::value<bool>( &warm_feeling_ )->zero_tokens() ,
      "Verbose" )
    ( "verbose,V" , po::value<bool>( &warm_feeling_ )->zero_tokens() ,
      "Verbose" )
    ( "input-format,F" , po::value<string>( &input_format_string_ ) ,
      "Input format : FLUSH_FPS|BITSTRINGS|BIN_FRAG_NUMS|FRAG_NUMS (default FLUSH_FPS)" )
    ( "output-format" , po::value<string>( &output_format_string_ ) ,
      "Output format : CSV_FORMAT|SAMPLES_FORMAT (default SAMPLES_FORMAT)" )
    ( "bitstring-separator" , po::value<string>( &bitstring_separator_ ) ,
      "For bitstrings input, the separator between bits (defaults to no separator)." )
    ( "frag-num-separator" , po::value<string>( &bitstring_separator_ ) ,
      "For fragment numbers input, the separator between numbers (defaults to space)." )
      ( "fix-spaces-in-names" , po::value<bool>( &fix_spaces_in_names_ )->zero_tokens()->default_value( false ) ,
        "Changes spaces in fingerprint names to \'_\' so as not to mess up SAMPLES format file.");
  
}

// ***************************************************************************
void ClusterSettings::decode_formats() {

  decode_format_string( input_format_string_ , input_format_ ,
                        binary_file_ , bitstring_separator_ );

  if( output_format_string_ == "CSV_FORMAT" ) {
    output_format_ = CSV_FORMAT;
  } else if(  output_format_string_ == "SAMPLES_FORMAT" ) {
    output_format_ = SAMPLES_FORMAT;
  } else {
    throw ClusterOutputFormatError( output_format_string_ );
  }

  if( FLUSH_FPS == input_format_ || BIN_FRAG_NUMS == input_format_ ) {
    binary_file_ = true;
  }

}
