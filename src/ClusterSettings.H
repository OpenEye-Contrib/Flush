//
// file ClusterSettings.H
// David Cosgrove
// AstraZeneca
// 23rd February 2009
//
// This class parse and holds the command-line settings for program cluster.

#ifndef DAC_CLUSTER_SETTINGS
#define DAC_CLUSTER_SETTINGS

#include <iosfwd>
#include <string>
#include <boost/program_options/options_description.hpp>

#include "FingerprintBase.H"

// **********************************************************************

typedef enum { CSV_FORMAT , SAMPLES_FORMAT } OUTPUT_FORMAT;

class ClusterSettings {

public :

  ClusterSettings() {} // will be filled by unpack_contents_from_pvm_buffer, probably
  ClusterSettings( int argc , char **argv );
  ~ClusterSettings() {}

  bool operator!() const;

  void print_usage( std::ostream &os ) const;
  std::string input_file() const { return input_file_; }
  std::string output_file() const { return output_file_; }
  std::string subset_file() const { return subset_file_; }
  double threshold() const { return threshold_; }
  double singletons_threshold() const { return singletons_threshold_; }

  bool warm_feeling() const { return warm_feeling_; }
  OUTPUT_FORMAT output_format() const { return output_format_; }
  DAC_FINGERPRINTS::FP_FILE_FORMAT input_format() const { return input_format_; }
  bool binary_file() const { return binary_file_; }
  std::string bitstring_separator() const { return bitstring_separator_; }
  bool fix_spaces_in_names() const { return fix_spaces_in_names_; }
  std::string usage_text() const { return usage_text_; }
  std::string error_message() const { return error_msg_; }

  void send_contents_via_mpi( int dest_slave );
  void receive_contents_via_mpi();

private :

  std::string input_file_;
  std::string output_file_;
  std::string subset_file_;
  double threshold_;
  double singletons_threshold_; // for collapse singletons

  bool warm_feeling_;
  std::string output_format_string_;
  std::string input_format_string_;
  OUTPUT_FORMAT output_format_;
  DAC_FINGERPRINTS::FP_FILE_FORMAT input_format_;
  bool binary_file_;
  std::string bitstring_separator_;
  bool fix_spaces_in_names_;
  std::string usage_text_;
  mutable std::string error_msg_;

  void build_program_options( boost::program_options::options_description &desc );

  void decode_formats();

};

class ClusterInputFormatError {

public :

  explicit ClusterInputFormatError( const std::string &bad_format_string ) :
    msg_( std::string( "Bad input format specifier string : " ) +
	  bad_format_string ) {}

  const char *what() {
    return msg_.c_str();
  }

private :

  std::string msg_;

};

class ClusterOutputFormatError {

public :

  explicit ClusterOutputFormatError( const std::string &bad_format_string ) :
    msg_( std::string( "Bad output format specifier string : " ) +
	  bad_format_string ) {}

  const char *what() {
    return msg_.c_str();
  }

private :

  std::string msg_;

};

#endif
