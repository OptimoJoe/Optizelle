#pragma once

// File IO utilities 

// std::string
#include <string>

// std::runtime_error
#include <stdexcept>

// strerror
#include <cstring>

// std::function
#include <functional>

// Location macro
#define S1(x) #x
#define S2(x) S1(x)
#define __LOC__ std::string("File \"" __FILE__ "\", line " S2(__LINE__))

// Checks that we opened a file.  We need this as a macro to get the file
// location correct.
#define CHECK_FILE(s,fname) \
    if(!((s).is_open())) \
        throw std::runtime_error( \
            __LOC__  + ", unable to open the file " + (fname) + ": " \
                + strerror(errno));

// Checks for errors on stream objects.  We need this as a macro to get the
// file location correct.
#define CHECK_STREAM(s) \
    if((s).bad()) \
        throw std::runtime_error( \
            __LOC__ + ", error with the stream object: " \
            + strerror(errno));

// Prints out nested exceptions 
void print_exception(std::exception const & e);

// Iterate over a stream, line by line, and apply the function f to each line 
template <typename Stream>
void iter_over_stream(
    std::function <void(std::string const &)> const & f,
    Stream & in
);

// Iterate over a file, line by line, and apply the function f to each line 
void iter_over_file(
    std::function <void(std::string const &)> const & f,
    std::string const & fname
);

// Produces a function that writes to a stream line by line
template <typename Stream>
std::function <void(std::string const &)> iter_to_stream(Stream & out);

// Our templated functions
#include "file.tpp"
