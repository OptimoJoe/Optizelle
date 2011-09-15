#ifndef PEOPT_H
#define PEOPT_H

#include<list>
#include<cmath>
#include<sstream>
#include<string>
#include<iomanip>
#include<limits>

/** \file peopt.h
    \brief Routines for constructing optimization algorithms for parameter
    estimation
**/

/// Parameter Estimation Optimization toolkit
namespace peopt {
  #include "core.h"
  #include "diff_operator.h"
}

#endif  // PEOPT_H
