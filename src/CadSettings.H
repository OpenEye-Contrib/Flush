//
// file CadSettings.H
// David Cosgrove
// 4th February 2010
//
// This class parses the command-line arguments for program cad and holds the
// corresponding settings.

#ifndef DAC_CAD_SETTINGS
#define DAC_CAD_SETTINGS

#include <iosfwd>
#include <string>
#include <boost/program_options/options_description.hpp>

#include "FingerprintBase.H"

// *******************************************************************

typedef enum { CSV_FORMAT , SAMPLES_FORMAT } CLUS_FILE_FORMAT;

class CadSettings {

public :

  CadSettings( int argc , char **argv );
  ~CadSettings() {}

  bool operator!() const;

  std::string cluster_file() const { return clus_file_; }
  std::string cluster_fp_file() const { return fp_file_; }
  std::string output_file() const { return output_file_; }

  DAC_FINGERPRINTS::FP_FILE_FORMAT input_format() const { return input_format_; }
  CLUS_FILE_FORMAT clus_file_format() const { return clus_file_format_; }
  std::string bitstring_separator() const { return bitstring_separator_; }
  bool binary_file() const { return binary_file_; }
  std::string usage_text() const { return usage_text_; }
  std::string error_message() const { return error_msg_; }

private :

  std::string clus_file_;
  std::string output_file_;
  std::string fp_file_;
  bool binary_file_;
  DAC_FINGERPRINTS::FP_FILE_FORMAT input_format_;
  CLUS_FILE_FORMAT clus_file_format_;
  std::string fp_format_string_;
  std::string bitstring_separator_;
  std::string clus_format_string_;
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
