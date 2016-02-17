//
// file SatanSettings.cc
// David Cosgrove
// 3rd March 2009
//
// This class parses the command-line arguments for program satan and holds the
// corresponding settings.

#include <iostream>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <mpi.h>
#include "SatanSettings.H"

using namespace std;
using namespace DAC_FINGERPRINTS;
namespace po = boost::program_options;

namespace DACLIB {
// In mpi_string_subs.cc
void mpi_send_string( const string &str , int dest_rank );
void mpi_rec_string( int source_rank , string &str );
}

// ***************************************************************************
SatanSettings::SatanSettings( int argc , char **argv ) :
  threshold_( 0.3 ) , min_count_( 0 ) , probe_chunk_size_( -1 ) ,
  tversky_alpha_( 0.5F ) ,
  warm_feeling_( false ) , binary_file_( false ) ,
  input_format_( FLUSH_FPS ) , sim_calc_( TANIMOTO ) ,
  input_format_string_( "FLUSH_FPS" ) , output_format_string_( "SATAN" ) ,
  sim_calc_string_( "TANIMOTO" ) {

  po::options_description desc( "Allowed Options" );
  build_program_options( desc );

  po::variables_map vm;
  try {
    po::store( po::parse_command_line( argc , argv , desc ) , vm );
  } catch( po::error &e ) {
    cerr << "Error parsing command line : " << e.what() << endl
         << "satan aborts." << endl;
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
bool SatanSettings::operator!() const {

  if( probe_file_.empty() ) {
    error_msg_ = "No probe file specified.";
    return true;
  } else if( target_file_.empty() ) {
    error_msg_ = "No target file specified.";
    return true;
  } else if( output_file_.empty() ) {
    error_msg_ = "No output file specified.";
    return true;
  } else if( threshold_ < 0.0 || threshold_ > 1.0 ) {
    error_msg_ = string( "Invalid distance threshold " ) +
        boost::lexical_cast<string>( threshold_ ) + string( "." );
    return true;
  } else if( tversky_alpha_ < 0.0F || tversky_alpha_ > 1.0F ) {
    error_msg_ = string( "Invalid tversky_alpha " ) +
        boost::lexical_cast<string>( threshold_ ) + string( "." );
    return true;
  }

  if( string( "SATAN" ) != output_format_string_ &&
      string( "NNLISTS" ) != output_format_string_ &&
      string( "COUNTS" ) != output_format_string_ ) {
    error_msg_ = string( "Invalid output format string : " ) + output_format_string_ +
        string( "\nMust be one of SATAN or NNLISTS or COUNTS.\n" );
    return true;
  }

  return false;

}

// ***************************************************************************
void SatanSettings::send_contents_via_mpi( int dest_rank ) {

  DACLIB::mpi_send_string( probe_file_ , dest_rank );
  DACLIB::mpi_send_string( target_file_ , dest_rank );

  MPI_Send( &threshold_ , 1 , MPI_DOUBLE , dest_rank , 0 , MPI_COMM_WORLD );
  MPI_Send( &min_count_ , 1 , MPI_INT , dest_rank , 0 , MPI_COMM_WORLD );
  MPI_Send( &probe_chunk_size_ , 1 , MPI_INT , dest_rank , 0 , MPI_COMM_WORLD );
  MPI_Send( &tversky_alpha_ , 1 , MPI_FLOAT , dest_rank , 0 , MPI_COMM_WORLD );
  int i = int( binary_file_ );
  MPI_Send( &i , 1 , MPI_INT , dest_rank , 0 , MPI_COMM_WORLD );
  i = int( input_format_ );
  MPI_Send( &i , 1 , MPI_INT , dest_rank , 0 , MPI_COMM_WORLD );
  i = int( sim_calc_ );
  MPI_Send( &i , 1 , MPI_INT , dest_rank , 0 , MPI_COMM_WORLD );

  DACLIB::mpi_send_string( input_format_string_ , dest_rank );
  DACLIB::mpi_send_string( bitstring_separator_ , dest_rank );
  DACLIB::mpi_send_string( sim_calc_string_ , dest_rank );
  DACLIB::mpi_send_string( output_format_string_ , dest_rank );

}

// ***************************************************************************
void SatanSettings::receive_contents_via_mpi() {

  DACLIB::mpi_rec_string( 0 , probe_file_ );
  DACLIB::mpi_rec_string( 0 , target_file_ );

  MPI_Recv( &threshold_ , 1 , MPI_DOUBLE , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &min_count_ , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &probe_chunk_size_ , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &tversky_alpha_ , 1 , MPI_FLOAT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  int i;
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  binary_file_ = static_cast<bool>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  input_format_ = static_cast<DAC_FINGERPRINTS::FP_FILE_FORMAT>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  sim_calc_ = static_cast<DAC_FINGERPRINTS::SIMILARITY_CALC>( i );

  DACLIB::mpi_rec_string( 0 , input_format_string_ );
  DACLIB::mpi_rec_string( 0 , bitstring_separator_ );
  DACLIB::mpi_rec_string( 0 , sim_calc_string_ );
  DACLIB::mpi_rec_string( 0 , output_format_string_ );

}

// ****************************************************************************
void SatanSettings::build_program_options( po::options_description &desc ) {

  desc.add_options()
      ( "help" , "Produce this help text." )
      ( "probe-file,P" , po::value<string>( &probe_file_ ) ,
        "Name of probe fingerprint file." )
      ( "target-file,T" , po::value<string>( &target_file_ ) ,
        "Name of target fingerprint file." )
      ( "output-file,O" , po::value<string>( &output_file_ ) ,
        "Name of the output clusters file." )
      ( "threshold" , po::value<double>( &threshold_ ) ,
        "Neighbour list distance threshold (default 0.3)" )
      ( "min-count,M" , po::value<int>( &min_count_ ) ,
        "Minimum neighbour count, defaults to 0 (report all neighbours)" )
      ( "probe-chunk-size" , po::value<int>( &probe_chunk_size_ ) ,
        "Controls the size of the pieces in which the probe is dealt with. Needs to be relatively small for large jobs." )
      ( "warm-feeling,W" , po::value<bool>( &warm_feeling_ )->zero_tokens() ,
        "Verbose" )
      ( "verbose,V" , po::value<bool>( &warm_feeling_ )->zero_tokens() ,
        "Verbose" )
      ( "input-format,F" , po::value<string>( &input_format_string_ ) ,
        "Input format : FLUSH_FPS|BITSTRINGS|BIN_FRAG_NUMS|FRAG_NUMS (default FLUSH_FPS)" )
      ( "output-format" , po::value<string>( &output_format_string_ ) ,
        "Output format : SATAN|NNLISTS|COUNTS (default SATAN)" )
      ( "distance-calculation" , po::value<string>( &sim_calc_string_ ) ,
        "Distance calculation : TANIMOTO|TVERSKY (default TANIMOTO)" )
      ( "tversky-alpha" , po::value<float>( &tversky_alpha_ ) ,
        "Tversky alpha parameter (0.0-1.0, default 0.5" )
      ( "bitstring-separator" , po::value<string>( &bitstring_separator_ ) ,
        "For bitstrings input, the separator between bits (defaults to no separator)." )
      ( "frag-num-separator" , po::value<string>( &bitstring_separator_ ) ,
        "For fragment numbers input, the separator between numbers (defaults to space)." );
  
}

// ***************************************************************************
void SatanSettings::decode_formats() {

  decode_format_string( input_format_string_ , input_format_ ,
                        binary_file_ , bitstring_separator_ );

  if( sim_calc_string_ == "TANIMOTO" ) {
    sim_calc_ = DAC_FINGERPRINTS::TANIMOTO;
  } else if( sim_calc_string_ == "TVERSKY" ) {
    sim_calc_ = DAC_FINGERPRINTS::TVERSKY;
  } else {
    throw FingerprintDistCalcError( sim_calc_string_ );
  }

}
