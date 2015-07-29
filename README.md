[![Optizelle](http://www.optimojoe.com/img/optizelle-gh.jpg "Optizelle")](http://www.optimojoe.com/products/optizelle)

<div align="center">
<p><strong>Brought to you by</strong></p>
<a href="http://www.optimojoe.com"><img src="http://www.optimojoe.com/img/optimojoe-gh.jpg" alt="OptimoJoe" ></a>
</div>


# Optizelle 

**Optizelle** [op-tuh-zel] is an open source software library designed to solve general purpose nonlinear optimization problems of the form

| min f(x) | min f(x) st g(x)=0 |
|:-----------|:--------------------------------------------------|
| **min f(x) st h(x)&ge;0** | **min f(x) st g(x)=0, h(x)&ge;0**  |


# It features
* **State of the art algorithms**
    * Unconstrained -- steepest descent, preconditioned nonlinear-CG (Fletcher-Reeves, Polak-Ribiere, Hestenes-Stiefel), BFGS, Newton-CG, SR1, trust-region Newton (truncated-CG and truncated-MINRES), Barzilai-Borwein two-point approximation
    * Equality constrained -- inexact composite-step SQP
    * Inequality constrained -- primal-dual interior point method for cone constraints (linear, second-order cone, and semidefinite), log-barrier method for cone constraints
    * Constrained -- any combination of the above
* **Open source**
    * Released under the 2-Clause BSD License
    * Free and ready to use with both open and closed sourced commercial codes
* **Multilanguage support**
    * Interfaces to C++, MATLAB/Octave, and Python
* **Robust computations and repeatability**
    * Can stop, archive, and restart the computation from any optimization iteration
    * Combined with the multilanguage support, the optimization can be started in one language and migrated to another.  For example, archived optimization runs that started in Python can be migrated and completed in C++.
* **User-defined parallelism**
    * Fully compatible with OpenMP, MPI, or GPUs
* **Extensible linear algebra**
    * Supports user-defined vector algebra and preconditioners
    * Enables sparse, dense, and matrix-free computations
    * Ability to define custom inner products and compatibility with preconditioners such as algebraic multigrid make Optizelle well-suited for PDE constrained optimization
* **Sophisticated Control of the Optimization Algorithms**
    * Allows the user to insert arbitrary code into the optimization algorithm, which enables custom heuristics to be embedded without modifying the source.  For example, in signal processing applications, the optimization iterates could be run through a band-pass filter at the end of each optimization iteration.

# Download

For the most current version of Optizelle, please download a [zipped archive](https://github.com/OptimoJoe/Optizelle/archive/master.zip) of our code.  For power users, we provide public access to our git repository on our [Github page](https://github.com/OptimoJoe/Optizelle).  In order to clone the Optizelle repository, use the command

```
git clone https://github.com/OptimoJoe/Optizelle.git
```

# Documentation

We provide a full set of instructions for building, installing, and using Optizelle in our manual ([letter](http://www.optimojoe.com/uploads/reports/Optizelle-1.1.2-letter.pdf),[a4](http://www.optimojoe.com/uploads/reports/Optizelle-1.1.2-a4.pdf).)

# Installation

For a full set of installation instructions, please see our manual ([letter](http://www.optimojoe.com/uploads/reports/Optizelle-1.1.2-letter.pdf),[a4](http://www.optimojoe.com/uploads/reports/Optizelle-1.1.2-a4.pdf).)  For power users, installation may be as simple as executing the following commands from the base Optizelle directory:

1. `mkdir build`
1. `cd build`
1. `ccmake ..`
1. Configure the build
1. `make install`

# Support

For general questions, please visit our [community forum](http://forum.optimojoe.com).  In addition, we provide paid support and consulting for Optizelle. If you are interested, please [contact](http://www.optimojoe.com/contact/) us.

# Contributing

We appreciate community contributions to Optizelle!  If you notice a bug, please file a report on our [issues page](https://github.com/OptimoJoe/Optizelle/issues).  Alternatively, for more in-depth contributions, clone our repository and send a pull request via our [Github page](https://github.com/OptimoJoe/Optizelle).  Finally, we appreciate help in answering general user questions on our [community forum](http://forum.optimojoe.com).
