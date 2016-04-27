// Our header
#include "file.h"

// std::cerr
#include <iostream>

// std::exception, std::throw_with_nested, std::rethrow_if_nested
#include <exception>

// Prints out nested exceptions 
void print_exception(std::exception const & e) {
    std::cerr << e.what() << std::endl; 
    try {
        std::rethrow_if_nested(e);
    } catch(std::exception const & nested) {
        print_exception(nested);
    } catch(...) {
        std::cerr << "Unknown exception" << std::endl;
    }
}

// Iterate over a file, line by line, and apply the function f to each line 
void iter_over_file(
    std::function <void(std::string const &)> const & f,
    std::string const & fname
) try {
    // Open the file into a stream
    std::ifstream fin(fname);
    CHECK_FILE(fin,fname);

    // Now, iterate the function over the stream
    iter_over_stream(f,fin);

} catch(...) {
    std::throw_with_nested(
        std::runtime_error(__LOC__ + 
            ", problem iterating a function over the file, " + fname));
}
