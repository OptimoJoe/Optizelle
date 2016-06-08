#include "stream.h"

// std::find_if
#include <algorithm>

namespace Optizelle { namespace Stream {
    // Implements a macro processor that iterates over a stream and then
    // does macro replacement for all the strings found in a map
    void macro(
        std::map<
            std::string,
            std::function <void(std::string const &)>
        > const macros,
        t <std::string> & stream
    ) try {
        // Make sure that the macros contain an default element
        if(macros.count("default")==0)
            throw Exception::t(__LOC__
                + ", macro processor requires a default function (default)");

        // Otherwise, iterate over all of the elements and apply the appropriate
        // macro function
        iter <std::string> (
            [&macros](std::string const & s_) {
                // Divide the string into the initial white space, if any,
                // and the rest of the string
                auto divider = std::find_if(s_.begin(),s_.end(),
                    std::not1(std::ptr_fun<int,int>(isspace) ));
                auto ws = std::string(s_.begin(),divider);
                auto s = std::string(divider,s_.end()); 

                // Loop over all of the elements
                for(auto const & macro: macros) {
                    // Don't match the default element
                    if(macro.first=="default")
                        continue;

                    // See if we find our token at the beginning of the string
                    auto i = s.find(macro.first);

                    // If our token is the first thing we find, we call our
                    // corresponding routine for macro replacement with the
                    // amount of whitespace that we've seen and skip the rest
                    // of the macros
                    if(i==0) {
                        macro.second(ws);
                        return;
                    }
                }

                // If we've not matched anything, call the default element 
                macros.at("default")(s_);
            },
            stream);
    } catch(...) {
        std::throw_with_nested(
            Exception::t(__LOC__ + ", problem macro processing over the stream"));
    }

    // Wrap cin in something that we can instantiate
    cin::cin() {  
        // Set this buffer to be the same as std::cin
        rdbuf(std::cin.rdbuf());
    }
} }
