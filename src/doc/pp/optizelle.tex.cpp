// Stream::macro
#include "optizelle/stream.h"

// Exception::print_exception
#include "optizelle/exception.h"

// std::ifstream
#include <fstream>

// std::cout
#include <iostream>

int main(int argc,char * argv[]) try {
    // Make sure that we have our data
    if(argc!=2)
        throw std::runtime_error(std::string("Usage: ") + argv[0] + " <data>");

    // Set our file name for the data 
    auto const data_name = argv[1];

    // Create a stream out of stdin
    auto in = Stream::of_std <Stream::cin,std::string>(
        new Stream::cin(),
        {'\n'});

    // Macro process stdin
    Stream::macro(
        {
            {R"(%-flags)",
            [&data_name](auto const & ws) {
                // Open our data 
                auto file = std::unique_ptr <std::ifstream> (
                    new std::ifstream(data_name));
                CHECK_FILE(*file,data_name);
                auto data = Stream::of_std <std::ifstream,std::string>(
                    file.release(),
                    {'\n'});

                // Output the formatted data
                Stream::iter <std::string> (
                    [&](auto const & s) {
                        std::cout
                            << ws 
                            << s
                            << R"( \)"
                            << std::endl;
                    },
                    data);
            }},
            {"default",
            [](auto const & s) {
                // Just print the string
                std::cout << s << std::endl;
            }}
        },
        in);
} catch(std::exception const & e){
    Exception::print_exception(e);
    throw;
}
