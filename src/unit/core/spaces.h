#pragma once

// Some easy to use spaces for our examples
typedef double Real;

template <typename Real>
using XX = Optizelle::Rm <Real>;
typedef XX <Real> X;
typedef typename X::Vector X_Vector;

template <typename Real>
using YY = Optizelle::Rm <Real>;
typedef YY <Real> X;
typedef typename X::Vector X_Vector;
