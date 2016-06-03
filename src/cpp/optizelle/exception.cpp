#include "exception.h"

// std::cerr
#include <iostream>

namespace Exception {
    // Prints out nested exceptions 
    void print_exception(std::exception const & e) {
        // Print out the top layer of the exception
        std::cerr << e.what() << std::endl; 

        // Recursively break things down if we're nested
        try {
            std::rethrow_if_nested(e);

        // If we're still a std::exception, continue unrolling 
        } catch(std::exception const & nested) {
            print_exception(nested);

        // Bail once we can't unroll any further
        } catch(...) {
            std::cerr << "Unknown exception" << std::endl;
        }
    }

    // Converts a nested expression to a string
    std::string exception_to_string(std::exception const & e) {
        // Grab the top layer of the exception
        auto s = std::string(e.what()) + '\n';

        // Recursively break things down if we're nested
        try {
            std::rethrow_if_nested(e);
            return s; 

        // If we're still a std::exception, continue unrolling 
        } catch(std::exception const & nested) {
            return s + exception_to_string(nested);

        // Bail once we can't unroll any further
        } catch(...) {
            return s + "Unknown exception" + '\n';
        }
    }
}
