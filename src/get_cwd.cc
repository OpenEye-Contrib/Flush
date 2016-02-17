//
// file get_cwd.cc
// David Cosgrove
// AstraZeneca
// 9th December 2014
//
// This file contains a function to return the current working directory
// as a string.  It's a spell I found on t'Interweb.

#include <sys/param.h>
#include <unistd.h>

#include <string>

namespace DACLIB {

// ****************************************************************************
std::string get_cwd(){

  char *buffer = new char[MAXPATHLEN];
  getcwd(buffer,MAXPATHLEN);
  if(buffer != NULL){
    std::string ret(buffer);
    delete[] buffer;
    return ret;
  } else {
    return std::string();
  }

}

}
