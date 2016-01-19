#pragma once

// Some easy to use spaces for our examples
template <typename Real>
using XX = Optizelle::Rm <Real>;
typedef double Real;
typedef XX <Real> X;
typedef typename X::Vector Vector;
typedef Unit::Matrix <Real>::t Matrix;
