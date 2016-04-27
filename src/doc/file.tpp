// std::endl
#include <iostream>

// std::ifstream
#include <fstream>

// Iterate over a stream, line by line, and apply the function f to each line 
template <typename Stream>
void iter_over_stream(
    std::function <void(std::string const &)> const & f,
    Stream & in
) {
    // Grab each line, one at a time, and apply our function to it 
    auto line=std::string("");
    while(std::getline(in,line))
        f(line);
    CHECK_STREAM(in);
}

// Produces a function that writes to a stream line by line
template <typename Stream>
std::function <void(std::string const &)> iter_to_stream(Stream & out) {
    // Produce a function that writes out the strings line by line
    return [&] (auto const & s) {
        try {
            out << "OUT: " << s << std::endl;
            CHECK_STREAM(out);
        } catch(...) {
            std::throw_with_nested(
                std::runtime_error(__LOC__ + ", problem writing to a stream"));
        }
    };
}
