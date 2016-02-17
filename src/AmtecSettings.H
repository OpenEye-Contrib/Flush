//
// file AmtecSettings.H
// David Cosgrove
// 3rd March 2009
//
// This class parses the command-line arguments for program amtec and holds the
// corresponding settings.

#ifndef DAC_AMTEC_SETTINGS
#define DAC_AMTEC_SETTINGS

#include <iosfwd>
#include <string>
#include <boost/program_options/options_description.hpp>

#include "FingerprintBase.H"

// *******************************************************************

typedef enum { CSV_FORMAT , SAMPLES_FORMAT } CLUS_FILE_FORMAT;

class AmtecSettings {

public :

  AmtecSettings( int argc , char **argv );
  ~AmtecSettings() {}

  bool operator!() const;

  std::string input_cluster_file() const { return input_clus_file_; }
  std::string output_cluster_file() const { return output_clus_file_; }
  std::string existing_cluster_fp_file() const { return clus_fp_file_; }
  std::string incoming_cluster_fp_file() const { return new_fp_file_; }
  std::string additions_file() const { return additions_file_; }
  std::string new_subset_file() const { return new_subset_file_; }

  double threshold() const { return threshold_; }
  double tversky_alpha() const { return tversky_alpha_; }
  DAC_FINGERPRINTS::FP_FILE_FORMAT input_format() const { return input_format_; }
  DAC_FINGERPRINTS::SIMILARITY_CALC similarity_calc() const { return sim_calc_; }
  CLUS_FILE_FORMAT clus_output_format() const { return clus_output_format_; }
  CLUS_FILE_FORMAT clus_input_format() const { return clus_input_format_; }
  std::string bitstring_separator() const { return bitstring_separator_; }
  bool binary_file() const { return binary_file_; }
  std::string usage_text() const { return usage_text_; }
  std::string error_message() const { return error_msg_; }

private :

  std::string input_clus_file_;
  std::string output_clus_file_;
  std::string clus_fp_file_;
  std::string new_fp_file_;
  std::string new_subset_file_;
  std::string additions_file_;
  double threshold_;
  double tversky_alpha_;
  bool binary_file_;
  DAC_FINGERPRINTS::FP_FILE_FORMAT input_format_;
  DAC_FINGERPRINTS::SIMILARITY_CALC sim_calc_;
  CLUS_FILE_FORMAT clus_input_format_ , clus_output_format_;
  std::string input_format_string_;
  std::string bitstring_separator_;
  std::string sim_calc_string_;
  std::string clus_input_format_string_ , clus_output_format_string_;
  std::string usage_text_;
  mutable std::string error_msg_;

  void build_program_options( boost::program_options::options_description &desc );

  void decode_formats();

};

// **********************************************************************
class ClusterFileFormatError {
public :
    explicit ClusterFileFormatError( const std::string &bad_format_string ) :
    msg_( std::string( "Bad cluster file format specifier string : " ) +
          bad_format_string ) {}
const char *what() {
  return msg_.c_str();
}
private :
  std::string msg_;
};

#endif
