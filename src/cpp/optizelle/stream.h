// std::exception
#include <exception>

// std::function
#include <functional>

// strerror
#include <cstring>

// std::map
#include <map>

// __LOC__
#include "exception.h"

// std::istream
#include <iostream>

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

// Stream functions
namespace Stream{ 
    // A : Set -> t 
    template <typename A>
    struct t;

    // Set -> A : Set -> t A
    template <typename Stream,typename A>
    struct of_std; 

    // A : Set -> (A -> unit) -> t A -> unit
    template <typename A>
    void iter(std::function <void(A const &)> const & f,t <A> & stream);

    // A : Set -> (A -> unit) -> (A -> unit) -> t A -> unit 
    template <typename A>
    void iter_with_last(
        std::function <void(A const &)> const & f,
        std::function <void(A const &)> const & g,
        t <A> & stream
    );

    // A : Set -> B : Set -> (A -> B) -> t A -> t B 
    template <typename A,typename B>
    struct map; 

    // A : Set -> (A -> bool) -> t A -> t A  
    template <typename A>
    struct filter;

    // A : Set -> B : Set -> (A->B->B) -> t A -> B -> B  
    template <typename A,typename B>
    B fold(
        std::function<B(A const &,B const &)> const & f,
        t <A> & stream,
        B acc
    );

    // map string (string->unit) -> t string -> unit
    void macro(
        std::map<
            std::string,
            std::function <void(std::string const &)>
        > const macros,
        t <std::string> & stream
    );

    // istream
    struct cin : public std::istream {
        cin();
    };
}

#include "stream.tpp"
