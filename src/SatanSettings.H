//
// file SatanSettings.H
// David Cosgrove
// 3rd March 2009
//
// This class parses the command-line arguments for program satan and holds the
// corresponding settings.

#ifndef DAC_SATAN_SETTINGS
#define DAC_SATAN_SETTINGS

#include <iosfwd>
#include <string>
#include <boost/program_options/options_description.hpp>

#include "FingerprintBase.H"

// *******************************************************************

class SatanSettings {

public :

  SatanSettings() {} // will be filled by unpacked from PVM, probably
  SatanSettings( int argc , char **argv );
  ~SatanSettings() {}

  bool operator!() const;

  std::string probe_file() const { return probe_file_; }
  std::string target_file() const { return target_file_; }
  std::string output_file() const { return output_file_; }
  double threshold() const { return threshold_; }
  int min_count() const { return min_count_; }
  int probe_chunk_size() const { return probe_chunk_size_; }
  float tversky_alpha() const { return tversky_alpha_; }
  DAC_FINGERPRINTS::FP_FILE_FORMAT input_format() const { return input_format_; }
  std::string output_format() const { return output_format_string_; }
  DAC_FINGERPRINTS::SIMILARITY_CALC similarity_calc() const { return sim_calc_; }
  std::string bitstring_separator() const { return bitstring_separator_; }
  bool warm_feeling() const { return warm_feeling_; }
  bool binary_file() const { return binary_file_; }
  std::string usage_text() const { return usage_text_; }
  std::string error_message() const { return error_msg_; }

  std::string slave_name() const { return slave_name_; }
  int num_slave_procs() const { return num_slave_procs_; }
  std:: string pvm_hosts_file() const { return pvm_hosts_file_; }

  void send_contents_via_mpi( int dest_rank );
  void receive_contents_via_mpi();

private :

  std::string probe_file_;
  std::string target_file_;
  std::string output_file_;
  double threshold_;
  int min_count_;
  int probe_chunk_size_; /* how the probe should be divided up - needs to be
			    small for large jobs, defaults to FP_CHUNK_SIZE */
  float tversky_alpha_;
  bool warm_feeling_;
  bool binary_file_;
  DAC_FINGERPRINTS::FP_FILE_FORMAT input_format_;
  DAC_FINGERPRINTS::SIMILARITY_CALC sim_calc_;
  std::string input_format_string_;
  std::string output_format_string_;
  std::string bitstring_separator_;
  std::string sim_calc_string_;
  std::string usage_text_;
  mutable std::string error_msg_;

  // for parallel runs
  std::string slave_name_;
  int num_slave_procs_;
  std:: string pvm_hosts_file_;

  void build_program_options( boost::program_options::options_description &desc );

  void decode_formats();

};

#endif
