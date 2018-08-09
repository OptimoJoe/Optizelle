#pragma once

// std::exception
#include <exception>

// std::string
#include <string>

#include <stdexcept>

// Location macro
#define S1(x) #x
#define S2(x) S1(x)
#define __LOC__ std::string("File \"" __FILE__ "\", line " S2(__LINE__))

//---Optizelle0---
namespace Optizelle {namespace Exception {
//---Optizelle1---
    // exception -> unit
    void to_stderr(std::exception const & e);

    // exception -> std::string
    std::string to_string(std::exception const & e);

    // t
    //---Exception0---
    struct t : public std::runtime_error {
        using std::runtime_error::runtime_error;
    };
    //---Exception1---
//---Optizelle2---
}}
//---Optizelle3---
