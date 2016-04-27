// File IO 
#include "file.h"

// Produces a function that indents the content and adds newlines 
template <typename Stream>
std::function <void(std::string const &)> iter_to_config(Stream & out) {
    // Produce a function that writes out the strings line by line
    return [&] (auto const & s) {
        try {
            out << "    " << s << " \\" << std::endl;
            CHECK_STREAM(out);
        } catch(...) {
            std::throw_with_nested(
                std::runtime_error(__LOC__ + ", problem writing the flags"));
        }
    };
}

// Converts CMake flags into a build script 
void config_to_build(std::string const & iname,std::string const & oname) try {

    // Open a file for output
    std::ofstream fout(oname);
    CHECK_FILE(fout,oname);

    // Write the start of the cmake command 
    fout << "cmake \\" << std::endl;

    // Iterate over the given file and push the information to the above
    // file 
    iter_over_file(iter_to_config(fout),iname); 

    // Give the final directory 
    fout << ".." << std::endl;

} catch(...) {
    std::throw_with_nested(
        std::runtime_error(__LOC__ + ", problem writing the build file")); 
}

int main(int argc,char * argv[]) try {
    if(argc!=3)
        throw std::runtime_error(__LOC__ + ", process_config <config> <out>");
    config_to_build(argv[1],argv[2]);
    return EXIT_SUCCESS;
} catch(std::exception const & e){
    print_exception(e);
    throw;
}
