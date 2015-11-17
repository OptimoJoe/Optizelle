/*
Copyright 2013-2014 OptimoJoe.

For the full copyright notice, see LICENSE.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Author: Joseph Young (joe@optimojoe.com)
*/

#include <Utility.h>

namespace Optizelle {
    namespace StoppingCondition { 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & opt_stop) {
            // Do the conversion
            switch(opt_stop){
            case NotConverged:
                return Python::enumToPyObject(
                    "StoppingCondition","NotConverged");
            case RelativeGradientSmall:
                return Python::enumToPyObject(
                    "StoppingCondition","RelativeGradientSmall");
            case RelativeStepSmall:
                return Python::enumToPyObject(
                    "StoppingCondition","RelativeStepSmall");
            case MaxItersExceeded:
                return Python::enumToPyObject(
                    "StoppingCondition","MaxItersExceeded");
            case InteriorPointInstability:
                return Python::enumToPyObject(
                    "StoppingCondition","InteriorPointInstability");
            case UserDefined:
                return Python::enumToPyObject(
                    "StoppingCondition","UserDefined");
            default:
                throw;
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member) {
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(member);

            if(m==Python::enumToNatural("StoppingCondition","NotConverged"))
                return NotConverged;
            else if(m==Python::enumToNatural(
                "StoppingCondition","RelativeGradientSmall")
            )
                return RelativeGradientSmall;
            else if(m==Python::enumToNatural(
                "StoppingCondition","RelativeStepSmall")
            )
                return RelativeStepSmall;
            else if(m==Python::enumToNatural(
                "StoppingCondition","MaxItersExceeded")
            )
                return MaxItersExceeded;
            else if(m==Python::enumToNatural(
                "StoppingCondition","InteriorPointInstability")
            )
                return InteriorPointInstability;
            else if(m==Python::enumToNatural("StoppingCondition","UserDefined"))
                return UserDefined;
            else
                throw;
        }
    }
    
    namespace KrylovStop { 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & krylov_stop) {
            // Do the conversion
            switch(krylov_stop){
            case NotConverged:
                return Python::enumToPyObject("KrylovStop","NotConverged");
            case NegativeCurvature:
                return Python::enumToPyObject("KrylovStop","NegativeCurvature");
            case RelativeErrorSmall:
                return Python::enumToPyObject(
                    "KrylovStop","RelativeErrorSmall");
            case MaxItersExceeded:
                return Python::enumToPyObject("KrylovStop","MaxItersExceeded");
            case TrustRegionViolated:
                return Python::enumToPyObject(
                    "KrylovStop","TrustRegionViolated");
            case NanDetected:
                return Python::enumToPyObject("KrylovStop","NanDetected");
            case LossOfOrthogonality:
                return Python::enumToPyObject(
                    "KrylovStop","LossOfOrthogonality");
            case InvalidTrustRegionOffset:
                return Python::enumToPyObject(
                    "KrylovStop","InvalidTrustRegionOffset");
            case TooManyFailedSafeguard:
                return Python::enumToPyObject(
                    "KrylovStop","TooManyFailedSafeguard");
            default:
                throw;
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member) {
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(member);

            if(m==Python::enumToNatural("KrylovStop","NotConverged"))
                return NotConverged;
            else if(m==Python::enumToNatural("KrylovStop","NegativeCurvature"))
                return NegativeCurvature;
            else if(m==Python::enumToNatural("KrylovStop","RelativeErrorSmall"))
                return RelativeErrorSmall;
            else if(m==Python::enumToNatural("KrylovStop","MaxItersExceeded"))
                return MaxItersExceeded;
            else if(m==Python::enumToNatural(
                "KrylovStop","TrustRegionViolated")
            )
                return TrustRegionViolated;
            else if(m==Python::enumToNatural("KrylovStop","NanDetected"))
                return NanDetected;
            else if(m==Python::enumToNatural(
                "KrylovStop","LossOfOrthogonality")
            )
                return LossOfOrthogonality;
            else if(m==Python::enumToNatural(
                "KrylovStop","InvalidTrustRegionOffset")
            )
                return InvalidTrustRegionOffset;
            else if(m==Python::enumToNatural(
                "KrylovStop","TooManyFailedSafeguard")
            )
                return TooManyFailedSafeguard;
            else
                throw;
        }
    }
    
    namespace KrylovSolverTruncated { 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & truncated_krylov) {
            // Do the conversion
            switch(truncated_krylov){
            case ConjugateDirection:
                return Python::enumToPyObject(
                    "KrylovSolverTruncated","ConjugateDirection");
            case MINRES:
                return Python::enumToPyObject(
                    "KrylovSolverTruncated","MINRES");
            default:
                throw;
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member) {
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(member);

            if(m==Python::enumToNatural(
                "KrylovSolverTruncated","ConjugateDirection")
            )
                return ConjugateDirection;
            else if(m==Python::enumToNatural("KrylovSolverTruncated","MINRES"))
                return MINRES;
            else
                throw;
        }
    }

    namespace AlgorithmClass { 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & algorithm_class) {
            // Do the conversion
            switch(algorithm_class){
            case TrustRegion:
                return Python::enumToPyObject("AlgorithmClass","TrustRegion");
            case LineSearch:
                return Python::enumToPyObject("AlgorithmClass","LineSearch");
            case UserDefined:
                return Python::enumToPyObject("AlgorithmClass","UserDefined");
            default:
                throw;
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member) {
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(member);

            if(m==Python::enumToNatural("AlgorithmClass","TrustRegion"))
                return TrustRegion;
            else if(m==Python::enumToNatural("AlgorithmClass","LineSearch"))
                return LineSearch;
            else if(m==Python::enumToNatural("AlgorithmClass","UserDefined"))
                return UserDefined;
            else
                throw;
        }
    }

    namespace Operators { 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & op) {
            // Do the conversion
            switch(op){
            case Identity:
                return Python::enumToPyObject("Operators","Identity");
            case ScaledIdentity:
                return Python::enumToPyObject("Operators","ScaledIdentity");
            case BFGS:
                return Python::enumToPyObject("Operators","BFGS");
            case InvBFGS:
                return Python::enumToPyObject("Operators","InvBFGS");
            case SR1:
                return Python::enumToPyObject("Operators","SR1");
            case InvSR1:
                return Python::enumToPyObject("Operators","InvSR1");
            case UserDefined:
                return Python::enumToPyObject("Operators","UserDefined");
            default:
                throw;
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member) {
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(member);

            if(m==Python::enumToNatural("Operators","Identity"))
                return Identity;
            else if(m==Python::enumToNatural("Operators","ScaledIdentity"))
                return ScaledIdentity;
            else if(m==Python::enumToNatural("Operators","BFGS"))
                return BFGS;
            else if(m==Python::enumToNatural("Operators","InvBFGS"))
                return InvBFGS;
            else if(m==Python::enumToNatural("Operators","SR1"))
                return SR1;
            else if(m==Python::enumToNatural("Operators","InvSR1"))
                return InvSR1;
            else if(m==Python::enumToNatural("Operators","UserDefined"))
                return UserDefined;
            else
                throw;
        }
    }

    namespace LineSearchDirection {
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & dir) {
            // Do the conversion
            switch(dir){
            case SteepestDescent:
                return Python::enumToPyObject("LineSearchDirection",
                    "SteepestDescent");
            case FletcherReeves:
                return Python::enumToPyObject("LineSearchDirection",
                    "FletcherReeves");
            case PolakRibiere:
                return Python::enumToPyObject("LineSearchDirection",
                    "PolakRibiere");
            case HestenesStiefel:
                return Python::enumToPyObject("LineSearchDirection",
                    "HestenesStiefel");
            case BFGS:
                return Python::enumToPyObject("LineSearchDirection","BFGS");
            case NewtonCG:
                return Python::enumToPyObject("LineSearchDirection","NewtonCG");
            default:
                throw;
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member) {
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(member);

            if(m==Python::enumToNatural("LineSearchDirection",
                "SteepestDescent")
            )
                return SteepestDescent;
            else if(m==Python::enumToNatural("LineSearchDirection",
                "FletcherReeves")
            )
                return FletcherReeves;
            else if(m==Python::enumToNatural("LineSearchDirection",
                "PolakRibiere")
            )
                return PolakRibiere;
            else if(m==Python::enumToNatural("LineSearchDirection",
                "HestenesStiefel")
            )
                return HestenesStiefel;
            else if(m==Python::enumToNatural("LineSearchDirection","BFGS"))
                return BFGS;
            else if(m==Python::enumToNatural("LineSearchDirection","NewtonCG"))
                return NewtonCG;
            else
                throw;
        }
    }

    namespace LineSearchKind { 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & kind) {
            // Do the conversion
            switch(kind){
            case Brents:
                return Python::enumToPyObject("LineSearchKind","Brents");
            case GoldenSection:
                return Python::enumToPyObject("LineSearchKind","GoldenSection");
            case BackTracking:
                return Python::enumToPyObject("LineSearchKind","BackTracking");
            case TwoPointA:
                return Python::enumToPyObject("LineSearchKind","TwoPointA");
            case TwoPointB:
                return Python::enumToPyObject("LineSearchKind","TwoPointB");
            default:
                throw;
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member) {
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(member);

            if(m==Python::enumToNatural("LineSearchKind","Brents"))
                return Brents;
            else if(m==Python::enumToNatural("LineSearchKind","GoldenSection"))
                return GoldenSection;
            else if(m==Python::enumToNatural("LineSearchKind","BackTracking"))
                return BackTracking;
            else if(m==Python::enumToNatural("LineSearchKind","TwoPointA"))
                return TwoPointA;
            else if(m==Python::enumToNatural("LineSearchKind","TwoPointB"))
                return TwoPointB;
            else
                throw;
        }
    }

    namespace OptimizationLocation { 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & loc) {
            // Do the conversion
            switch(loc){
            case BeginningOfOptimization:
                return Python::enumToPyObject(
                    "OptimizationLocation","BeginningOfOptimization");
            case BeforeInitialFuncAndGrad:
                return Python::enumToPyObject(
                    "OptimizationLocation","BeforeInitialFuncAndGrad");
            case AfterInitialFuncAndGrad:
                return Python::enumToPyObject(
                    "OptimizationLocation","AfterInitialFuncAndGrad");
            case BeforeOptimizationLoop:
                return Python::enumToPyObject(
                    "OptimizationLocation","BeforeOptimizationLoop");
            case BeginningOfOptimizationLoop:
                return Python::enumToPyObject(
                    "OptimizationLocation","BeginningOfOptimizationLoop");
            case BeforeSaveOld:
                return Python::enumToPyObject(
                    "OptimizationLocation","BeforeSaveOld");
            case BeforeStep:
                return Python::enumToPyObject(
                    "OptimizationLocation","BeforeStep");
            case BeforeGetStep:
                return Python::enumToPyObject(
                    "OptimizationLocation","BeforeGetStep");
            case GetStep:
                return Python::enumToPyObject("OptimizationLocation","GetStep");
            case AfterStepBeforeGradient:
                return Python::enumToPyObject(
                    "OptimizationLocation","AfterStepBeforeGradient");
            case AfterGradient:
                return Python::enumToPyObject(
                    "OptimizationLocation","AfterGradient");
            case BeforeQuasi:
                return Python::enumToPyObject(
                    "OptimizationLocation","BeforeQuasi");
            case AfterQuasi:
                return Python::enumToPyObject(
                    "OptimizationLocation","AfterQuasi");
            case AfterCheckStop:
                return Python::enumToPyObject(
                    "OptimizationLocation","AfterCheckStop");
            case EndOfOptimizationIteration:
                return Python::enumToPyObject(
                    "OptimizationLocation","EndOfOptimizationIteration");
            case BeforeLineSearch:
                return Python::enumToPyObject(
                    "OptimizationLocation","BeforeLineSearch");
            case AfterRejectedTrustRegion:
                return Python::enumToPyObject(
                    "OptimizationLocation","AfterRejectedTrustRegion");
            case AfterRejectedLineSearch:
                return Python::enumToPyObject(
                    "OptimizationLocation","AfterRejectedLineSearch");
            case BeforeActualVersusPredicted:
                return Python::enumToPyObject(
                    "OptimizationLocation","BeforeActualVersusPredicted");
            case EndOfKrylovIteration:
                return Python::enumToPyObject(
                    "OptimizationLocation","EndOfKrylovIteration");
            case EndOfOptimization:
                return Python::enumToPyObject(
                    "OptimizationLocation","EndOfOptimization");
            default:
                throw;
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member) {
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(member);

            if(m==Python::enumToNatural(
                "OptimizationLocation","BeginningOfOptimization"))
                return BeginningOfOptimization;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","BeforeInitialFuncAndGrad"))
                return BeforeInitialFuncAndGrad;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","AfterInitialFuncAndGrad"))
                return AfterInitialFuncAndGrad;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","BeforeOptimizationLoop"))
                return BeforeOptimizationLoop;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","BeginningOfOptimizationLoop"))
                return BeginningOfOptimizationLoop;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","BeforeSaveOld"))
                return BeforeSaveOld;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","BeforeStep"))
                return BeforeStep;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","BeforeGetStep"))
                return BeforeGetStep;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","GetStep"))
                return GetStep;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","AfterStepBeforeGradient"))
                return AfterStepBeforeGradient;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","AfterGradient"))
                return AfterGradient;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","BeforeQuasi"))
                return BeforeQuasi;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","AfterQuasi"))
                return AfterQuasi;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","AfterCheckStop"))
                return AfterCheckStop;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","EndOfOptimizationIteration"))
                return EndOfOptimizationIteration;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","BeforeLineSearch"))
                return BeforeLineSearch;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","AfterRejectedTrustRegion"))
                return AfterRejectedTrustRegion;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","AfterRejectedLineSearch"))
                return AfterRejectedLineSearch;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","BeforeActualVersusPredicted"))
                return BeforeActualVersusPredicted;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","EndOfKrylovIteration"))
                return EndOfKrylovIteration;
            else if(m==Python::enumToNatural(
                "OptimizationLocation","EndOfOptimization"))
                return EndOfOptimization;
            else
                throw;
        }
    }

    namespace FunctionDiagnostics { 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & diag) {
            // Do the conversion
            switch(diag){
            case NoDiagnostics:
                return Python::enumToPyObject("FunctionDiagnostics",
                    "NoDiagnostics");
            case FirstOrder:
                return Python::enumToPyObject("FunctionDiagnostics",
                    "FirstOrder");
            case SecondOrder:
                return Python::enumToPyObject("FunctionDiagnostics",
                    "SecondOrder");
            default:
                throw;
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member) {
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(member);

            if(m==Python::enumToNatural("FunctionDiagnostics","NoDiagnostics"))
                return NoDiagnostics;
            else if(m==Python::enumToNatural("FunctionDiagnostics",
                "FirstOrder")
            )
                return FirstOrder;
            else if(m==Python::enumToNatural("FunctionDiagnostics",
                "SecondOrder")
            )
                return SecondOrder;
            else
                throw;
        }
    }

    namespace VectorSpaceDiagnostics { 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & diag) {
            // Do the conversion
            switch(diag){
            case NoDiagnostics:
                return Python::enumToPyObject("VectorSpaceDiagnostics",
                    "NoDiagnostics");
            case Basic:
                return Python::enumToPyObject("VectorSpaceDiagnostics",
                    "Basic");
            case EuclideanJordan:
                return Python::enumToPyObject("VectorSpaceDiagnostics",
                    "EuclideanJordan");
            default:
                throw;
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member) {
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(member);

            if(m==Python::enumToNatural("VectorSpaceDiagnostics",
                "NoDiagnostics")
            )
                return NoDiagnostics;
            else if(m==Python::enumToNatural("VectorSpaceDiagnostics",
                "Basic")
            )
                return Basic;
            else if(m==Python::enumToNatural("VectorSpaceDiagnostics",
                "EuclideanJordan")
            )
                return EuclideanJordan;
            else
                throw;
        }
    }

    namespace DiagnosticScheme { 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & dscheme) {
            // Do the conversion
            switch(dscheme){
            case Never:
                return Python::enumToPyObject("DiagnosticScheme","Never");
            case DiagnosticsOnly:
                return Python::enumToPyObject("DiagnosticScheme",
                    "DiagnosticsOnly");
            case EveryIteration:
                return Python::enumToPyObject("DiagnosticScheme",
                    "EveryIteration");
            default:
                throw;
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member) {
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(member);

            if(m==Python::enumToNatural("DiagnosticScheme","Never"))
                return Never;
            else if(m==Python::enumToNatural("DiagnosticScheme",
                "DiagnosticsOnly")
            )
                return DiagnosticsOnly;
            else if(m==Python::enumToNatural("DiagnosticScheme",
                "EveryIteration")
            )
                return EveryIteration;
            else
                throw;
        }
    }

    namespace json {
        // Serialization utility for the Rm vector space
        template <>
        struct Serialization <double,Python::PythonVS> {
            static std::string serialize (
                Python::Vector const & x,
                std::string const & name_,
                Natural const & iter_
            ) {
                // Grab the serialization module 
                Python::PyObjectPtr module(PyImport_ImportModule(
                    "Optizelle.json.Serialization")); 

                // Now, get the serialize routine
                Python::PyObjectPtr serialize(PyObject_GetAttrString(
                    module.get(),"serialize"));
                
                // Make a Python object of the name and iteration 
                Python::PyObjectPtr name(PyString_FromString(name_.c_str()));
                Python::PyObjectPtr iter(PyInt_FromSize_t(iter_));

                // Call the serialize routine on the vector
                Python::PyObjectPtr x_json(
                    Python::PyObject_CallObject3(
                        serialize.get(),
                        const_cast <Python::Vector &> (x).get(),
                        name.get(),
                        iter.get()));
            
                // Check errors
                if(x_json.get()==nullptr) {
                    std::string msg(
                        "Evaluation of the serialize function failed.\n");
                    PySys_WriteStderr("%s",msg.c_str());
                    Python::PyErr_SetString_Optizelle(msg);
                    throw Optizelle::Python::Exception();
                }

                // Convert the serialized vector to a string and return it 
                return std::string(PyString_AsString(x_json.get()));
            }

            static Python::Vector deserialize (
                Python::Vector const & x_,
                std::string const & x_json_
            ) {
                // Grab the serialization module 
                Python::PyObjectPtr module(PyImport_ImportModule(
                    "Optizelle.json.Serialization")); 

                // Now, get the deserialize routine
                Python::PyObjectPtr deserialize(PyObject_GetAttrString(
                    module.get(),"deserialize"));

                // Convert the inputed string into Python
                Python::PyObjectPtr x_json(
                    PyString_FromString(x_json_.c_str()));

                // Allocate memory for a new Python vector
                Python::Vector x(const_cast <Python::Vector &> (x_).init());
                
                // Call the deserialize routine on the reference vector and the
                // json vector
                Python::PyObjectPtr x_raw(Python::PyObject_CallObject2(
                    deserialize.get(),
                    x.get(),
                    x_json.get()));
            
                // Check errors
                if(x_raw.get()==nullptr) {
                    std::string msg(
                        "Evaluation of the deserialize function failed.\n");
                    PySys_WriteStderr("%s",msg.c_str());
                    Python::PyErr_SetString_Optizelle(msg);
                    throw Optizelle::Python::Exception();
                }

                // Move the raw information into the Python vector
                x.reset(x_raw.release());

                // Move out the new vector
                return std::move(x);
            }
        };
    }

    namespace Python {
        // Converts Py_ssize_t to Natural
        Natural Py_ssize_t_to_Natural(Py_ssize_t const & x) {
            return x < 0 ? 0 : x;
        }

        // A function to alter the behavior of PyTuple_SetItem so that we don't
        // have to hand increment the reference to the object since SetItem
        // takes control of its arguments.
        void MyPyTuple_SetItem(PyObject * p,Natural const & pos,PyObject * o) {
            Py_INCREF(o);
            PyTuple_SetItem(p,pos,o);
        }

        // Calls a Python function with one argument 
        PyObject * PyObject_CallObject1(
            PyObject * const fn,
            PyObject * const arg1
        ) {
            PyObjectPtr args(PyTuple_New(1)); 
            MyPyTuple_SetItem(args.get(),0,arg1);
            return PyObject_CallObject(fn,args.get()); 
        }
        
        // Calls a Python function with two arguments
        PyObject * PyObject_CallObject2(
            PyObject * const fn,
            PyObject * const arg1,
            PyObject * const arg2
        ) {
            PyObjectPtr args(PyTuple_New(2)); 
            MyPyTuple_SetItem(args.get(),0,arg1); 
            MyPyTuple_SetItem(args.get(),1,arg2); 
            return PyObject_CallObject(fn,args.get()); 
        }
        
        // Calls a Python function with three arguments
        PyObject * PyObject_CallObject3(
            PyObject * const fn,
            PyObject * const arg1,
            PyObject * const arg2,
            PyObject * const arg3
        ) {
            PyObjectPtr args(PyTuple_New(3)); 
            MyPyTuple_SetItem(args.get(),0,arg1); 
            MyPyTuple_SetItem(args.get(),1,arg2); 
            MyPyTuple_SetItem(args.get(),2,arg3); 
            return PyObject_CallObject(fn,args.get()); 
        }
        
        // Calls a Python function with four arguments
        PyObject * PyObject_CallObject4(
            PyObject * const fn,
            PyObject * const arg1,
            PyObject * const arg2,
            PyObject * const arg3,
            PyObject * const arg4
        ) {
            PyObjectPtr args(PyTuple_New(4)); 
            MyPyTuple_SetItem(args.get(),0,arg1); 
            MyPyTuple_SetItem(args.get(),1,arg2); 
            MyPyTuple_SetItem(args.get(),2,arg3); 
            MyPyTuple_SetItem(args.get(),3,arg4); 
            return PyObject_CallObject(fn,args.get()); 
        }

        // Used to catch Python exceptions
        Exception::Exception() {}

        // Deep copy of a Python object and return the result
        PyObject * deepcopy(PyObject * const in) {
            // Grab the deepcopy function from the copy module 
            PyObjectPtr module(PyImport_ImportModule("copy")); 
            PyObjectPtr deepcopy(PyObject_GetAttrString(module.get(),
                "deepcopy")); 

            // Call deepcopy on vec and return the result
            PyObjectPtr args(PyTuple_New(1)); 
            MyPyTuple_SetItem(args.get(),0,in); 
            return PyObject_CallObject(
                deepcopy.get(),
                args.get()); 
        }

        // On construction, initialize the pointer and figure out if
        // we're capturing the pointer or attaching to it
        PyObjectPtr::PyObjectPtr(
            PyObject * const ptr_,
            PyObjectPtrMode::t const mode 
        ) : ptr(ptr_) {
            switch(mode) {
            case PyObjectPtrMode::Capture:
                break;
            case PyObjectPtrMode::Attach:
                Py_XINCREF(ptr);
                break;
            }
        }
            
        // Move constructor
        PyObjectPtr::PyObjectPtr(PyObjectPtr&& ptr_) noexcept
            : ptr(ptr_.release()) {}
        
        // Move assignment operator
        PyObjectPtr const & PyObjectPtr::operator=(PyObjectPtr&& ptr_)noexcept {
            ptr=ptr_.release();
            return *this;
        }

        // For a reset, we decrement the pointer and then assign a new
        // value.
        void PyObjectPtr::reset(PyObject * const ptr_) {
            Py_XDECREF(ptr);
            ptr=ptr_;
        }

        // For an attach, we decrement the pointer, assign a new value,
        // and then increment the reference count.
        void PyObjectPtr::attach(PyObject * const ptr_) {
            Py_XDECREF(ptr);
            ptr=ptr_;
            Py_XINCREF(ptr);
        }

        // On a get, we simply return the pointer.
        PyObject * PyObjectPtr::get() {
            return ptr;
        }
    
        // On a release, we return the underlying pointer and then clear
        // the vector.  This will prevent a decrement later.
        PyObject * PyObjectPtr::release() {
            PyObject * ptr_=ptr;
            ptr=nullptr;
            return ptr_;
        }

        // On destruction, decrement the Python reference counter and do
        // not delete the pointer.
        PyObjectPtr::~PyObjectPtr() {
            Py_XDECREF(ptr);
            ptr=nullptr;
        }
            
        // On construction, we just grab the pointer to the messaging object
        Messaging::Messaging(
            PyObject * const ptr_,
            PyObjectPtrMode::t const mode
        ) : PyObjectPtr(ptr_,mode) {}
            
        // Move constructor
        Messaging::Messaging(Messaging && msg) noexcept
            : PyObjectPtr(msg.release()) {}

        // Move assignment operator
        Messaging const & Messaging::operator = (Messaging && msg) noexcept {
            ptr = msg.release();
            return *this;
        }
            
        // Prints a message
        void Messaging::print(std::string const & msg_) const {
            // Call the print function on msg
            PyObjectPtr print(PyObject_GetAttrString(ptr,"print"));
            PyObjectPtr msg(PyString_FromString(msg_.c_str()));
            PyObjectPtr ret(PyObject_CallObject1(
                print.get(),
                msg.get()));

            // Check errors
            if(ret.get()==nullptr)
                error("Evaluation of the print function in the Messaging "
                    "object failed.");
        }

        // Prints an error
        void Messaging::error(std::string const & msg_) const {
            // Call the error function on msg
            PyObjectPtr error(PyObject_GetAttrString(ptr,"error"));
            PyObjectPtr msg(PyString_FromString(msg_.c_str()));
            PyObjectPtr ret(PyObject_CallObject1(
                error.get(),
                msg.get()));

            // Check errors
            if(ret.get()==nullptr) {
                std::string msg2="Evaluation of the error function in the "
                    "Messaging object failed.\n";
                PySys_WriteStderr("%s",msg2.c_str());
                PyErr_SetString_Optizelle(msg2);
                throw Exception();
            }

            // Raise a Python exception
            PyErr_SetString_Optizelle(msg_);
            throw Exception();
        }

        // Create a vector with the appropriate messaging and vector space 
        Vector::Vector(
            PyObject * const msg_,
            PyObject * const vs_,
            PyObject * const vec,
            PyObjectPtrMode::t mode
        ) : 
            PyObjectPtr(vec,mode),
            msg(msg_,PyObjectPtrMode::Attach),
            vs(vs_,PyObjectPtrMode::Attach)
        {}
            
        // Create a move constructor so we can interact with stl objects
        Vector::Vector(Vector && vec) noexcept :
            PyObjectPtr(std::move(vec)),
            msg(std::move(vec.msg)),
            vs(std::move(vec.vs))
        { }
            
        // Move assignment operator
        Vector const & Vector::operator = (Vector && vec) noexcept {
            ptr = vec.release(); 
            msg = std::move(vec.msg);
            vs = std::move(vec.vs);
            return *this;
        }

        // Memory allocation and size setting 
        Vector Vector::init() { 
            // Call the init function on the internal and store in y 
            PyObjectPtr init(PyObject_GetAttrString(vs.get(),"init"));
            PyObjectPtr y(PyObject_CallObject1(
                init.get(),
                get()));

            // Check errors
            if(y.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function init failed.");

            // Create and return a new vector based on y
            return std::move(Vector(msg.get(),vs.get(),y.release()));
        } 
        
        // y <- x (Shallow.  No memory allocation.)  Internal is y.
        void Vector::copy(Vector & x) { 
            // Call the copy function on x and the internal 
            PyObjectPtr copy(PyObject_GetAttrString(vs.get(),"copy"));
            PyObjectPtr ret(PyObject_CallObject2(
                copy.get(),
                x.get(),
                get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function copy failed.");
        } 

        // x <- alpha * x.  Internal is x.
        void Vector::scal(double const & alpha_) { 
            // Call the scal function on alpha and the internal storage 
            PyObjectPtr scal(PyObject_GetAttrString(vs.get(),"scal"));
            PyObjectPtr alpha(PyFloat_FromDouble(alpha_));
            PyObjectPtr ret(PyObject_CallObject2(
                scal.get(),
                alpha.get(),
                get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function scal failed.");
        } 

        // x <- 0.  Internal is x. 
        void Vector::zero() { 
            // Call the zero function on this vector.
            PyObjectPtr zero(PyObject_GetAttrString(vs.get(),"zero"));
            PyObjectPtr ret(PyObject_CallObject1(
                zero.get(),
                get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function zero failed.");
        } 

        // y <- alpha * x + y.   Internal is y.
        void Vector::axpy(double const & alpha_,Vector & x) { 
            // Call the axpy function on alpha, x, and the internal storage.
            PyObjectPtr axpy(PyObject_GetAttrString(vs.get(),"axpy"));
            PyObjectPtr alpha(PyFloat_FromDouble(alpha_));
            PyObjectPtr ret(PyObject_CallObject3(
                axpy.get(),
                alpha.get(),
                x.get(),
                get()));
           
            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function axpy failed.");
        } 

        // innr <- <x,y>.  Internal is y.
        double Vector::innr(Vector & x) { 
            // Call the innr function on x and the internal.  Store in z. 
            PyObjectPtr innr(PyObject_GetAttrString(vs.get(),"innr"));
            PyObjectPtr z(PyObject_CallObject2(
                innr.get(),
                x.get(),
                get()));

            // Check errors
            if(z.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function innr failed.");

            // Return the result 
            return PyFloat_AsDouble(z.get()); 
        } 

        // x <- random.  Internal is x. 
        void Vector::rand() { 
            // Call the rand function on this vector.
            PyObjectPtr rand(PyObject_GetAttrString(vs.get(),"rand"));
            PyObjectPtr ret(PyObject_CallObject1(
                rand.get(),
                get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function rand failed.");
        } 

        // Jordan product, z <- x o y.  Internal is z.
        void Vector::prod(Vector & x,Vector & y) { 
            // Call the prod function on x, y, and the internal 
            PyObjectPtr prod(PyObject_GetAttrString(vs.get(),"prod"));
            PyObjectPtr ret(PyObject_CallObject3(
                prod.get(),
                x.get(),
                y.get(),
                get()));

            // Check errors
            if(get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function prod failed.");
        } 

        // Identity element, x <- e such that x o e = x .  Internal is x.
        void Vector::id() { 
            // Call the id function on the internal.
            PyObjectPtr id(PyObject_GetAttrString(vs.get(),"id"));
            PyObjectPtr ret(PyObject_CallObject1(
                id.get(),
                get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function id failed.");
        } 

        // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y.
        // Internal is z.
        void Vector::linv(Vector& x, Vector& y) { 
            // Call the linv function on x, y, and the internal
            PyObjectPtr linv(PyObject_GetAttrString(vs.get(),"linv"));
            PyObjectPtr ret(PyObject_CallObject3(
                linv.get(),
                x.get(),
                y.get(),
                get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function linv failed.");
        } 

        // Barrier function, barr <- barr(x) where x o grad barr(x) = e.
        // Internal is x.
        double Vector::barr() { 
            // Call the barr function on the internal.  Store in z.
            PyObjectPtr barr(PyObject_GetAttrString(vs.get(),"barr"));
            PyObjectPtr z(PyObject_CallObject1(
                barr.get(),
                get()));

            // Check errors
            if(z.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function barr failed.");

            // Return the result 
            return PyFloat_AsDouble(z.get());
        } 

        // Line search, srch <- argmax {alpha in Real >= 0 : alpha x + y >= 0} 
        // where y > 0.  Internal is y.
        double Vector::srch(Vector& x) {  
            // Call the srch function on x and the internal.  Store in z.
            PyObjectPtr srch(PyObject_GetAttrString(vs.get(),"srch"));
            PyObjectPtr z(PyObject_CallObject2(
                srch.get(),
                x.get(),
                get()));

            // Check errors
            if(z.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function srch failed.");

            // Return the result 
            return PyFloat_AsDouble(z.get());
        } 

        // Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
        // operator.  Internal is x.
        void Vector::symm() { 
            // Call the symm function on the internal.
            PyObjectPtr symm(PyObject_GetAttrString(vs.get(),"symm"));
            PyObjectPtr ret(PyObject_CallObject1(
                symm.get(),
                get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function symm failed.");

        } 
        
        // Converts (copies) a value into Python.  This assumes memory
        // has been allocated both in the vector as well as Python.
        void Vector::toPython(PyObject * const ptr) {
            // Call the copy function on the internal and x
            PyObjectPtr copy(PyObject_GetAttrString(vs.get(),"copy"));
            PyObjectPtr ret(PyObject_CallObject2(
                copy.get(),
                get(),
                ptr));

            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function copy failed.");
        } 
        
        // Converts (copies) a value from Python.  This assumes memory
        // has been allocated both in the vector as well as Python.
        void Vector::fromPython(PyObject * const ptr) {
            // Call the copy function on ptr and the internal 
            PyObjectPtr copy(PyObject_GetAttrString(vs.get(),"copy"));
            PyObjectPtr ret(PyObject_CallObject2(
                copy.get(),
                ptr,
                get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function copy failed.");
        } 
            
        // Convert a C++ state to a Python state 
        template <>
        void State <PyUnconstrained>::toPython(
            typename PyUnconstrained::State::t const & state
        ) {
            Unconstrained::State::toPython(state,ptr);
        }
        template <>
        void State <PyEqualityConstrained>::toPython(
            typename PyEqualityConstrained::State::t const & state
        ) {
            EqualityConstrained::State::toPython(state,ptr);
        }
        template <>
        void State <PyInequalityConstrained>::toPython(
            typename PyInequalityConstrained::State::t const & state
        ) {
            InequalityConstrained::State::toPython(state,ptr);
        }
        template <>
        void State <PyConstrained>::toPython(
            typename PyConstrained::State::t const & state
        ) {
            Constrained::State::toPython(state,ptr);
        }

        // Convert a Python state to C++ 
        template <>
        void State <PyUnconstrained>::fromPython(
            typename PyUnconstrained::State::t & state
        ) {
            Unconstrained::State::fromPython(ptr,state);
        }
        template <>
        void State <PyEqualityConstrained>::fromPython(
            typename PyEqualityConstrained::State::t & state
        ) {
            EqualityConstrained::State::fromPython(ptr,state);
        }
        template <>
        void State <PyInequalityConstrained>::fromPython(
            typename PyInequalityConstrained::State::t & state
        ) {
            InequalityConstrained::State::fromPython(ptr,state);
        }
        template <>
        void State <PyConstrained>::fromPython(
            typename PyConstrained::State::t & state
        ) {
            Constrained::State::fromPython(ptr,state);
        }
        
        // Convert a Python bundle to C++ 
        template <>
        void Functions <PyUnconstrained>::fromPython(
            typename PyUnconstrained::Functions::t & fns 
        ) {
            Unconstrained::Functions::fromPython(
                msg.get(),ptr,pystate.get(),state,fns);
        }
        template <>
        void Functions <PyEqualityConstrained>::fromPython(
            typename PyEqualityConstrained::Functions::t & fns 
        ) {
            EqualityConstrained::Functions::fromPython(
                msg.get(),ptr,pystate.get(),state,fns);
        }
        template <>
        void Functions <PyInequalityConstrained>::fromPython(
            typename PyInequalityConstrained::Functions::t & fns 
        ) {
            InequalityConstrained::Functions::fromPython(
                msg.get(),ptr,pystate.get(),state,fns);
        }
        template <>
        void Functions <PyConstrained>::fromPython(
            typename PyConstrained::Functions::t & fns 
        ) {
            Constrained::Functions::fromPython(
                msg.get(),ptr,pystate.get(),state,fns);
        }

        // Create a function 
        ScalarValuedFunction::ScalarValuedFunction(
            PyObject * const msg_,
            PyObject * const f,
            PyObjectPtrMode::t mode
        ) :
            PyObjectPtr(f,mode),
            msg(msg_,PyObjectPtrMode::Attach)
        { }

        // <- f(x) 
        double ScalarValuedFunction::eval(Vector const & x) const { 
            // Call the objective function on x.  Store in z.
            PyObjectPtr eval(PyObject_GetAttrString(ptr,"eval"));
            PyObjectPtr z(PyObject_CallObject1(
                eval.get(),
                const_cast <Vector &> (x).get()));

            // Check errors
            if(z.get()==nullptr)
                msg.error("Evaluation of the objective f failed.");

            // Return the result
            return PyFloat_AsDouble(z.get());
        }

        // grad = grad f(x) 
        void ScalarValuedFunction::grad(
            Vector const & x,
            Vector & grad
        ) const { 
            // Call the gradient function on x and grad. 
            PyObjectPtr pygrad(PyObject_GetAttrString(ptr,"grad"));
            PyObjectPtr ret(PyObject_CallObject2(
                pygrad.get(),
                const_cast <Vector &>(x).get(),
                grad.get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error("Evaluation of the gradient of f failed.");
        }

        // H_dx = hess f(x) dx 
        void ScalarValuedFunction::hessvec(
            Vector const & x,
            Vector const & dx,
            Vector & H_dx
        ) const {
            // Call the hessvec function on x, dx, and H_dx.
            PyObjectPtr hessvec(PyObject_GetAttrString(ptr,"hessvec"));
            PyObjectPtr ret(PyObject_CallObject3(
                hessvec.get(),
                const_cast <Vector &> (x).get(),
                const_cast <Vector &> (dx).get(),
                H_dx.get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error("Evaluation of the Hessian-vector product"
                    " of f failed.");
        }

        // Create a function 
        VectorValuedFunction::VectorValuedFunction(
            std::string const & name_,
            PyObject * const msg_,
            PyObject * const f,
            PyObjectPtrMode::t mode
        ) :
            PyObjectPtr(f,mode),
            msg(msg_,PyObjectPtrMode::Attach),
            name(name_)
        {}

        // y=f(x)
        void VectorValuedFunction::eval(
            Vector const & x,
            VectorValuedFunction::Y_Vector& y
        ) const {
            // Call the evaluate function on x and y.
            PyObjectPtr eval(PyObject_GetAttrString(ptr,"eval"));
            PyObjectPtr ret(PyObject_CallObject2(
                eval.get(),
                const_cast <Vector &> (x).get(),
                y.get()));

            // Check errors
            if(ret.get()==nullptr) {
                std::stringstream ss;
                ss << "Evaluation of the constraint " << name << " failed.";
                msg.error(ss.str());
            }
        }

        // y=f'(x)dx 
        void VectorValuedFunction::p(
            Vector const & x,
            Vector const & dx,
            VectorValuedFunction::Y_Vector& y
        ) const {
            // Call the prime function on x, dx, and y
            PyObjectPtr p(PyObject_GetAttrString(ptr,"p"));
            PyObjectPtr ret(PyObject_CallObject3(
                p.get(),
                const_cast <Vector &> (x).get(),
                const_cast <Vector &> (dx).get(),
                y.get()));
           
            // Check errors
            if(ret.get()==nullptr) {
                std::stringstream ss;
                ss << "Evaluation of the derivative of the constraint "
                    << name << " failed.";
                msg.error(ss.str());
            }
        }

        // z=f'(x)*dy
        void VectorValuedFunction::ps(
            Vector const & x,
            Vector const & dy,
            VectorValuedFunction::X_Vector& z
        ) const {
            // Call the prime-adjoint function on x, dy, and z
            PyObjectPtr ps(PyObject_GetAttrString(ptr,"ps"));
            PyObjectPtr ret(PyObject_CallObject3(
                ps.get(),
                const_cast <Vector &> (x).get(),
                const_cast <Vector &> (dy).get(),
                z.get()));

            // Check errors
            if(ret.get()==nullptr) {
                std::stringstream ss;
                ss << "Evaluation of the derivative-adjoint of the constraint "
                    << name << " failed.";
                msg.error(ss.str());
            }
        }
             
        // z=(f''(x)dx)*dy
        void VectorValuedFunction::pps(
            Vector const & x,
            Vector const & dx,
            Vector const & dy,
            X_Vector& z
        ) const { 
            // Call the prime-adjoint function on x, dx, dy, and z
            PyObjectPtr pps(PyObject_GetAttrString(ptr,"pps"));
            PyObjectPtr ret(PyObject_CallObject4(
                pps.get(),
                const_cast <Vector &> (x).get(),
                const_cast <Vector &> (dx).get(),
                const_cast <Vector &> (dy).get(),
                z.get()));

            // Check errors
            if(ret.get()==nullptr) {
                std::stringstream ss;
                ss << "Evaluation of the second derivative-adjoint of the "
                    "constraint " << name << " failed.";
                msg.error(ss.str());
            }
        }

        // Calls the Optizelle exception with a string
        void PyErr_SetString_Optizelle(std::string const& msg) {
            Python::PyObjectPtr module(PyImport_ImportModule("Optizelle")); 
            Python::PyObjectPtr exception(PyObject_GetAttrString(module.get(),
                "Exception"));
            PyErr_SetString(exception.get(),msg.c_str());
        }

        // Converts an Optizelle enumerated type to a PyObject * 
        PyObject * enumToPyObject(
            std::string const & type,
            std::string const & member 
        ) {
            // Grab the enumerated type object from the Optizelle module.
            // We just use simple classes in Python to represent the
            // enumerated type
            Python::PyObjectPtr module(PyImport_ImportModule("Optizelle")); 
            Python::PyObjectPtr pyclass(PyObject_GetAttrString(module.get(),
                type.c_str()));

            // Grab and return the member
            return PyObject_GetAttrString(pyclass.get(),member.c_str());
        }
       
        // Converts an Optizelle enumerated type to a Natural
        Natural enumToNatural(
            std::string const & type,
            std::string const & member 
        ) {
            // Grab the PyObject * for the type and member requested
            PyObjectPtr obj(enumToPyObject(type,member));

            // Convert and return the member
            return PyInt_AsSsize_t(obj.get());
        }
        
        // Converts elements from C++ to Python 
        namespace toPython {
                    
            // Sets a real in a Python state 
            void Real(
                std::string const & name,
                double const & value,
                PyObject * const obj 
            ) {
                PyObjectPtr item(PyFloat_FromDouble(value));
                PyObject_SetAttrString(obj,name.c_str(),item.get());
            }
        
            // Sets a natural in a Python state 
            void Natural(
                std::string const & name,
                Optizelle::Natural const & value,
                PyObject * const obj 
            ) {
                PyObjectPtr item(PyInt_FromSize_t(value));
                PyObject_SetAttrString(obj,name.c_str(),item.get());
            }
        
            // Sets a vector in a Python state 
            void Vector(
                std::string const & name,
                Python::Vector const & value,
                PyObject * const obj 
            ) {
                PyObjectPtr item(PyObject_GetAttrString(obj,name.c_str()));
                const_cast <Python::Vector &> (value).toPython(item.get());
            }
        
            // Sets a list of vectors in a Python state 
            void VectorList(
                std::string const & name,
                std::list <Python::Vector> const & values,
                PyObject * const obj 
            ) {
                // Create a new Python list that we insert elements into
                PyObjectPtr items(PyList_New(0));

                // Loop over all of the items inside values and then insert 
                // them into items 
                for(std::list <Python::Vector>::const_iterator value
                        = values.cbegin();
                    value!=values.cend();
                    value++
                ) {
                    // Allocate memory for a new vector
                    Python::Vector item(
                        const_cast <Python::Vector &> (*value).init());

                    // Copy the information from the current iterator into this
                    // new vector
                    item.copy(const_cast <Python::Vector &> (*value));

                    // Release the pointer into the Python list
                    PyList_Append(items.get(),item.release());
                }
                
                // Insert the items into obj
                PyObject_SetAttrString(obj,name.c_str(),items.get());
            }

            // Sets restart vectors in Python 
            void Vectors(
                Python::Vectors const & values,
                PyObject * const pyvalues 
            ) {
            
                // Loop over all of the items inside values and then insert 
                // them into pyvalues 
                for(typename Python::Vectors::const_iterator value
                        = values.cbegin();
                    value!=values.cend();
                    value++
                ) {
                    // Allocate memory for a new vector
                    Python::Vector pyvalue(
                        const_cast <Python::Vector &>(value->second).init());

                    // Copy the information from the current iterator into this
                    // new vector
                    pyvalue.copy(const_cast <Python::Vector &> (value->second));

                    // Release the pointer into the Python list
                    PyList_Append(pyvalues,PyTuple_Pack(2,
                        PyString_FromString(value->first.c_str()),
                        pyvalue.release()));
                }
            }
        
            // Sets restart reals in Python 
            void Reals(
                Python::Reals const & values,
                PyObject * const pyvalues 
            ) {
                // Loop over all of the items inside values and then insert 
                // them into pyvalues 
                for(typename Python::Reals::const_iterator value
                        = values.cbegin();
                    value!=values.cend();
                    value++
                ) {
                    // Insert the double into the Python list 
                    PyList_Append(pyvalues,PyTuple_Pack(2,
                        PyString_FromString(value->first.c_str()),
                        PyFloat_FromDouble(value->second)));
                }
            }
        
            // Converts a list of naturals to a Python list 
            void Naturals(
                Python::Naturals const & values,
                PyObject * const pyvalues 
            ) {
                // Loop over all of the items inside values and then insert 
                // them into pyvalues 
                for(typename Python::Naturals::const_iterator value
                        = values.cbegin();
                    value!=values.cend();
                    value++
                ) {
                    // Insert the double into the Python list 
                    PyList_Append(pyvalues,PyTuple_Pack(2,
                        PyString_FromString(value->first.c_str()),
                        PyInt_FromSize_t(value->second)));
                }
            }
        
            // Sets restart parameters in Python 
            void Params(
                Python::Params const & values,
                PyObject * const pyvalues 
            ) {
                // Loop over all of the items inside values and then insert 
                // them into pyvalues 
                for(typename Python::Params::const_iterator value
                        = values.cbegin();
                    value!=values.cend();
                    value++
                ) {
                    // Insert the double into the Python list 
                    PyList_Append(pyvalues,PyTuple_Pack(2,
                        PyString_FromString(value->first.c_str()),
                        PyString_FromString(value->second.c_str())));
                }
            }
        }
        
        // Converts elements from Python to C++ 
        namespace fromPython {
        
            // Sets a real in a C++ state 
            void Real(
                std::string const & name,
                PyObject * const obj,
                double & value
            ) {
                PyObjectPtr item(PyObject_GetAttrString(obj,name.c_str()));
                value=PyFloat_AsDouble(item.get());
            }
            
            // Sets a natural in a C++ state 
            void Natural(
                std::string const & name,
                PyObject * const obj,
                Optizelle::Natural & value
            ) {
                PyObjectPtr item(PyObject_GetAttrString(obj,name.c_str()));
                value=PyInt_AsSsize_t(item.get());
            }
            
            // Sets a list of vectors in a C++ state 
            void VectorList(
                std::string const & name,
                PyObject * const obj,
                Python::Vector const & vec,
                std::list <Python::Vector> & values
            ) {
                // Grab the list of items
                PyObjectPtr items(PyObject_GetAttrString(obj,name.c_str()));

                // Loop over all the elements in items and insert them one
                // at a time into values
                values.clear();
                for(Optizelle::Natural i=0;
                    i<Py_ssize_t_to_Natural(PyList_Size(items.get()));
                    i++
                ) {
                    // Grab the current item from Python
                    PyObject * item(PyList_GetItem(items.get(),i));

                    // Create a new vector in values 
                    values.emplace_back(std::move(
                        const_cast<Python::Vector &>(vec).init()));

                    // Copy the Python item into the new value
                    values.back().fromPython(item);
                }
            }
        
            // Sets a scalar-valued function in a C++ function bundle 
            void ScalarValuedFunction(
                std::string const & name,
                PyObject * const msg,
                PyObject * const obj,
                std::unique_ptr <PyScalarValuedFunction> & value
            ) {
                value.reset(new Python::ScalarValuedFunction(msg,
                    PyObject_GetAttrString(obj,name.c_str())));
            }
            
            // Sets a vector-valued function in a C++ function bundle 
            void VectorValuedFunction(
                std::string const & name,
                PyObject * const msg,
                PyObject * const obj,
                std::unique_ptr <PyVectorValuedFunction> & value
            ) {
                value.reset(new Python::VectorValuedFunction(name,msg,
                    PyObject_GetAttrString(obj,name.c_str())));
            }
            
            // Sets a vector in a C++ state 
            void Vector(
                std::string const & name,
                PyObject * const obj,
                Python::Vector & value
            ) {
                PyObjectPtr item(PyObject_GetAttrString(obj,name.c_str()));
                value.fromPython(item.get());
            }
        
            // Sets restart vectors in C++ 
            void Vectors(
                Python::Vector const & vec,
                PyObject * const pyvalues,
                Python::Vectors & values
            ) {
                // Loop over all the elements in pyvalues and insert them one
                // at a time into values
                values.clear();
                for(Optizelle::Natural i=0;
                    i<Py_ssize_t_to_Natural(PyList_Size(pyvalues));
                    i++
                ) {
                    // Grab the current item from Python
                    PyObject * pyvalue(PyList_GetItem(pyvalues,i));

                    // Create the elements in values 
                    values.emplace_back(
                        PyString_AsString(PyTuple_GetItem(pyvalue,0)),
                        std::move(const_cast<Python::Vector &>(vec).init()));

                    // Copy the Python value into the C++ value
                    values.back().second.fromPython(PyTuple_GetItem(pyvalue,1));
                }
            }
            
            // Sets restart reals in C++ 
            void Reals(
                PyObject * const pyvalues,
                Python::Reals & values
            ) {
                // Loop over all the elements in pyvalues and insert them one
                // at a time into values
                values.clear();
                for(Optizelle::Natural i=0;
                    i<Py_ssize_t_to_Natural(PyList_Size(pyvalues));
                    i++
                ) {
                    // Grab the current item from Python
                    PyObject * pyvalue(PyList_GetItem(pyvalues,i));
                    
                    // Create the elements in values 
                    values.emplace_back(
                        PyString_AsString(PyTuple_GetItem(pyvalue,0)),
                        PyFloat_AsDouble(PyTuple_GetItem(pyvalue,1)));
                }
            }
            
            // Sets restart naturals in C++ 
            void Naturals(
                PyObject * const pyvalues,
                Python::Naturals & values
            ) {
                // Loop over all the elements in pyvalues and insert them one
                // at a time into values
                values.clear();
                for(Optizelle::Natural i=0;
                    i<Py_ssize_t_to_Natural(PyList_Size(pyvalues));
                    i++
                ) {
                    // Grab the current item from Python
                    PyObject * pyvalue(PyList_GetItem(pyvalues,i));
                    
                    // Create the elements in values 
                    values.emplace_back(
                        PyString_AsString(PyTuple_GetItem(pyvalue,0)),
                        PyInt_AsSsize_t(PyTuple_GetItem(pyvalue,1)));
                }
            }
            
            // Sets restart parameters in C++ 
            void Params(
                PyObject * const pyvalues,
                Python::Params & values
            ) {
                // Loop over all the elements in pyvalues and insert them one
                // at a time into values
                values.clear();
                for(Optizelle::Natural i=0;
                    i<Py_ssize_t_to_Natural(PyList_Size(pyvalues));
                    i++
                ) {
                    // Grab the current item from Python
                    PyObject * pyvalue(PyList_GetItem(pyvalues,i));
                    
                    // Create the elements in values 
                    values.emplace_back(
                        PyString_AsString(PyTuple_GetItem(pyvalue,0)),
                        PyString_AsString(PyTuple_GetItem(pyvalue,1)));
                }
            }
        }
       
        // Routines that manipulate and support problems of the form
        // 
        // min_{x \in X} f(x)
        //
        // where f : X -> R.
        namespace Unconstrained {

            // Routines that manipulate the internal state of the optimization 
            // algorithm
            namespace State {
                // Convert a C++ state to a Python state 
                void toPython_(
                    typename PyUnconstrained::State::t const & state,
                    PyObject * const pystate
                ){
                    // Set each of the required items in the Python state
                    toPython::Real("eps_grad",state.eps_grad,pystate);
                    toPython::Real("eps_dx",state.eps_dx,pystate);
                    toPython::Natural("stored_history",
                        state.stored_history,pystate);
                    toPython::Natural("history_reset",
                        state.history_reset,pystate);
                    toPython::Natural("iter",state.iter,pystate);
                    toPython::Natural("iter_max",state.iter_max,pystate);
                    toPython::Param <StoppingCondition::t> (
                        "opt_stop",
                        StoppingCondition::toPython,
                        state.opt_stop,
                        pystate);
                    toPython::Natural("krylov_iter",state.krylov_iter,pystate);
                    toPython::Natural("krylov_iter_max",
                        state.krylov_iter_max,pystate);
                    toPython::Natural("krylov_iter_total",
                        state.krylov_iter_total,pystate);
                    toPython::Natural("krylov_orthog_max",
                        state.krylov_orthog_max,pystate);
                    toPython::Param <KrylovStop::t> (
                        "krylov_stop",
                        KrylovStop::toPython,
                        state.krylov_stop,
                        pystate);
                    toPython::Real("krylov_rel_err",
                        state.krylov_rel_err,pystate);
                    toPython::Real("eps_krylov",state.eps_krylov,pystate);
                    toPython::Param <KrylovSolverTruncated::t> (
                        "krylov_solver",
                        KrylovSolverTruncated::toPython,
                        state.krylov_solver,
                        pystate);
                    toPython::Param <AlgorithmClass::t> (
                        "algorithm_class",
                        AlgorithmClass::toPython,
                        state.algorithm_class,
                        pystate);
                    toPython::Param <Operators::t> (
                        "PH_type",
                        Operators::toPython,
                        state.PH_type,
                        pystate);
                    toPython::Param <Operators::t> (
                        "H_type",
                        Operators::toPython,
                        state.H_type,
                        pystate);
                    toPython::Real("norm_gradtyp",state.norm_gradtyp,pystate);
                    toPython::Real("norm_dxtyp",state.norm_dxtyp,pystate);
                    toPython::Vector("x",state.x,pystate);
                    toPython::Vector("grad",state.grad,pystate);
                    toPython::Vector("dx",state.dx,pystate);
                    toPython::Vector("x_old",state.x_old,pystate);
                    toPython::Vector("grad_old",state.grad_old,pystate);
                    toPython::Vector("dx_old",state.dx_old,pystate);
                    toPython::VectorList("oldY",state.oldY,pystate);
                    toPython::VectorList("oldS",state.oldS,pystate);
                    toPython::Real("f_x",state.f_x,pystate);
                    toPython::Real("f_xpdx",state.f_xpdx,pystate);
                    toPython::Natural("msg_level",state.msg_level,pystate);
                    toPython::Natural("failed_safeguard_max",
                        state.failed_safeguard_max,pystate);
                    toPython::Natural("failed_safeguard",
                        state.failed_safeguard,pystate);
                    toPython::Natural("failed_safeguard_total",
                        state.failed_safeguard_total,pystate);
                    toPython::Real("alpha_x",state.alpha_x,pystate);
                    toPython::Real("alpha_x_qn",state.alpha_x_qn,pystate);
                    toPython::Real("delta",state.delta,pystate);
                    toPython::Real("eta1",state.eta1,pystate);
                    toPython::Real("eta2",state.eta2,pystate);
                    toPython::Real("ared",state.ared,pystate);
                    toPython::Real("pred",state.pred,pystate);
                    toPython::Natural("rejected_trustregion",
                        state.rejected_trustregion,pystate);
                    toPython::Real("alpha0",state.alpha0,pystate);
                    toPython::Real("alpha",state.alpha,pystate);
                    toPython::Real("c1",state.c1,pystate);
                    toPython::Natural("linesearch_iter",
                        state.linesearch_iter,pystate);
                    toPython::Natural("linesearch_iter_max",
                        state.linesearch_iter_max,pystate);
                    toPython::Natural("linesearch_iter_total",
                        state.linesearch_iter_total,pystate);
                    toPython::Real("eps_ls",state.eps_ls,pystate);
                    toPython::Param <LineSearchDirection::t> (
                        "dir",
                        LineSearchDirection::toPython,
                        state.dir,
                        pystate);
                    toPython::Param <LineSearchKind::t> (
                        "kind",
                        LineSearchKind::toPython,
                        state.kind,
                        pystate);
                    toPython::Param <FunctionDiagnostics::t> (
                        "f_diag",
                        FunctionDiagnostics::toPython,
                        state.f_diag,
                        pystate);
                    toPython::Param <FunctionDiagnostics::t> (
                        "L_diag",
                        FunctionDiagnostics::toPython,
                        state.L_diag,
                        pystate);
                    toPython::Param <VectorSpaceDiagnostics::t> (
                        "x_diag",
                        VectorSpaceDiagnostics::toPython,
                        state.x_diag,
                        pystate);
                    toPython::Param <DiagnosticScheme::t> (
                        "dscheme",
                        DiagnosticScheme::toPython,
                        state.dscheme,
                        pystate);
                }
                void toPython(
                    typename PyUnconstrained::State::t const & state,
                    PyObject * const pystate
                ){
                    Unconstrained::State::toPython_(state,pystate);
                }
                
                // Convert a Python state to C++ 
                void fromPython_(
                    PyObject * const pystate,
                    typename PyUnconstrained::State::t & state
                ){
                    // Set each of the required items in the Python state
                    fromPython::Real("eps_grad",pystate,state.eps_grad);
                    fromPython::Real("eps_dx",pystate,state.eps_dx);
                    fromPython::Natural("stored_history",
                        pystate,state.stored_history);
                    fromPython::Natural("history_reset",
                        pystate,state.history_reset);
                    fromPython::Natural("iter",pystate,state.iter);
                    fromPython::Natural("iter_max",pystate,state.iter_max);
                    fromPython::Param <StoppingCondition::t> (
                        "opt_stop",
                        StoppingCondition::fromPython,
                        pystate,
                        state.opt_stop);
                    fromPython::Natural("krylov_iter",
                        pystate,state.krylov_iter);
                    fromPython::Natural("krylov_iter_max",
                        pystate,state.krylov_iter_max);
                    fromPython::Natural("krylov_iter_total",
                        pystate,state.krylov_iter_total);
                    fromPython::Natural("krylov_orthog_max",
                        pystate,state.krylov_orthog_max);
                    fromPython::Param <KrylovStop::t> (
                        "krylov_stop",
                        KrylovStop::fromPython,
                        pystate,
                        state.krylov_stop);
                    fromPython::Real("krylov_rel_err",
                        pystate,state.krylov_rel_err);
                    fromPython::Real("eps_krylov",pystate,state.eps_krylov);
                    fromPython::Param <KrylovSolverTruncated::t> (
                        "krylov_solver",
                        KrylovSolverTruncated::fromPython,
                        pystate,
                        state.krylov_solver);
                    fromPython::Param <AlgorithmClass::t> (
                        "algorithm_class",
                        AlgorithmClass::fromPython,
                        pystate,
                        state.algorithm_class);
                    fromPython::Param <Operators::t> (
                        "PH_type",
                        Operators::fromPython,
                        pystate,
                        state.PH_type);
                    fromPython::Param <Operators::t> (
                        "H_type",
                        Operators::fromPython,
                        pystate,
                        state.H_type);
                    fromPython::Real("norm_gradtyp",
                        pystate,state.norm_gradtyp);
                    fromPython::Real("norm_dxtyp",pystate,state.norm_dxtyp);
                    fromPython::Vector("x",pystate,state.x);
                    fromPython::Vector("grad",pystate,state.grad);
                    fromPython::Vector("dx",pystate,state.dx);
                    fromPython::Vector("x_old",pystate,state.x_old);
                    fromPython::Vector("grad_old",pystate,state.grad_old);
                    fromPython::Vector("dx_old",pystate,state.dx_old);
                    fromPython::VectorList("oldY",pystate,state.x,state.oldY);
                    fromPython::VectorList("oldS",pystate,state.x,state.oldS);
                    fromPython::Real("f_x",pystate,state.f_x);
                    fromPython::Real("f_xpdx",pystate,state.f_xpdx);
                    fromPython::Natural("msg_level",pystate,state.msg_level);
                    fromPython::Natural("failed_safeguard_max",
                        pystate,state.failed_safeguard_max);
                    fromPython::Natural("failed_safeguard",
                        pystate,state.failed_safeguard);
                    fromPython::Natural("failed_safeguard_total",
                        pystate,state.failed_safeguard_total);
                    fromPython::Real("alpha_x",pystate,state.alpha_x);
                    fromPython::Real("alpha_x_qn",pystate,state.alpha_x_qn);
                    fromPython::Real("delta",pystate,state.delta);
                    fromPython::Real("eta1",pystate,state.eta1);
                    fromPython::Real("eta2",pystate,state.eta2);
                    fromPython::Real("ared",pystate,state.ared);
                    fromPython::Real("pred",pystate,state.pred);
                    fromPython::Natural("rejected_trustregion",
                        pystate,state.rejected_trustregion);
                    fromPython::Real("alpha0",pystate,state.alpha0);
                    fromPython::Real("alpha",pystate,state.alpha);
                    fromPython::Real("c1",pystate,state.c1);
                    fromPython::Natural("linesearch_iter",
                        pystate,state.linesearch_iter);
                    fromPython::Natural("linesearch_iter_max",
                        pystate,state.linesearch_iter_max);
                    fromPython::Natural("linesearch_iter_total",pystate,
                        state.linesearch_iter_total);
                    fromPython::Real("eps_ls",pystate,state.eps_ls);
                    fromPython::Param <LineSearchDirection::t> (
                        "dir",
                        LineSearchDirection::fromPython,
                        pystate,
                        state.dir);
                    fromPython::Param <LineSearchKind::t> (
                        "kind",
                        LineSearchKind::fromPython,
                        pystate,
                        state.kind);
                    fromPython::Param <FunctionDiagnostics::t> (
                        "f_diag",
                        FunctionDiagnostics::fromPython,
                        pystate,
                        state.f_diag);
                    fromPython::Param <FunctionDiagnostics::t> (
                        "L_diag",
                        FunctionDiagnostics::fromPython,
                        pystate,
                        state.L_diag);
                    fromPython::Param <VectorSpaceDiagnostics::t> (
                        "x_diag",
                        VectorSpaceDiagnostics::fromPython,
                        pystate,
                        state.x_diag);
                    fromPython::Param <DiagnosticScheme::t> (
                        "dscheme",
                        DiagnosticScheme::fromPython,
                        pystate,
                        state.dscheme);
                }
                void fromPython(
                    PyObject * const pystate,
                    typename PyUnconstrained::State::t & state
                ){
                    Unconstrained::State::fromPython_(pystate,state);
                }

                // Creates a state and inserts the elements into pystate 
                PyObject * create(
                    PyObject * self,
                    PyObject * args
                ){
                    // Calling convention should be (pystate,X,msg,x) 
                    PyObject *pystate_,*X,*msg,*x_;
                    if(!PyArg_ParseTuple(args,"OOOO",&pystate_,&X,&msg,&x_))
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a vector from the user input
                        Vector x(msg,X,x_,PyObjectPtrMode::Attach);

                        // Create a Python state 
                        Python::State <PyUnconstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);

                        // Create a new C++ state
                        typename PyUnconstrained::State::t state(x);

                        // Convert the state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
        
                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be (X,msg,fname,pystate) 
                    PyObject *X,*msg_,*fname_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOO",
                        &X,&msg_,&fname_,&pystate_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Grab the file name
                        std::string fname(PyString_AsString(fname_));

                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);

                        // Create a Python state 
                        Python::State <PyUnconstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                    
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));

                        // Create a new C++ state
                        typename PyUnconstrained::State::t state(x);
                        
                        // Convert the Python state to a C++ state
                        pystate.fromPython(state);

                        // Read the JSON file into the C++ state
                        PyJsonUnconstrained::read(msg,fname,state);

                        // Convert the C++ state to a Python state
                        pystate.toPython(state);
                                
                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Python bundle to C++ 
                void fromPython(
                    PyObject * const msg,
                    PyObject * const pyfns,
                    PyObject * const pystate,
                    typename PyUnconstrained::State::t const & state,
                    typename PyUnconstrained::Functions::t & fns 
                ) {
                    Unconstrained::Functions::fromPython_
                        <PyUnconstrained> (msg,pyfns,pystate,state,fns);
                }
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem 
                PyObject * getMin(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be
                    // (X,msg,fns,state,smanip)
                    PyObject *X,*msg_,*pyfns_,*pystate_,*smanip_;
                    if(!PyArg_ParseTuple(args,"OOOOO",
                        &X,&msg_,&pyfns_,&pystate_,&smanip_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);
                            
                        // Create a Python state 
                        Python::State <PyUnconstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));

                        // Create a C++ state
                        typename PyUnconstrained::State::t state(x);
                        
                        // Convert the Python state to a C++ state
                        pystate.fromPython(state);

                        // Create a Python bundle of functions
                        Python::Functions <PyUnconstrained> pyfns(
                            msg.get(),
                            pystate.get(),
                            state,
                            pyfns_,
                            PyObjectPtrMode::Attach);

                        // Create a C++ bundle of functions
                        typename PyUnconstrained::Functions::t fns;
                        
                        // Convert the Python bundle of functions to C++ 
                        pyfns.fromPython(fns);
                        
                        // Create a state manipulator 
                        Python::StateManipulator <PyUnconstrained> smanip(
                            msg.get(),
                            pystate.get(),
                            pyfns.get(),
                            smanip_,
                            PyObjectPtrMode::Attach);
                       
                        // Minimize
                        PyUnconstrained::Algorithms::getMin(
                            msg,fns,state,smanip);
                        
                        // Convert the C++ state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
            }
            
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                PyObject * release(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be
                    // (X,msg,state,xs,reals,nats,params)
                    PyObject *X,*msg,*pystate_,*pyxs,*pyreals,*pynats,*pyparams;
                    if(!PyArg_ParseTuple(args,"OOOOOOO",
                        &X,&msg,&pystate_,&pyxs,&pyreals,&pynats,&pyparams)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a Python state 
                        Python::State <PyUnconstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg,X,
                            PyObject_GetAttrString(pystate.get(),"x"));

                        // Create a C++ state
                        typename PyUnconstrained::State::t state(x);
                        
                        // Convert the Python state to a C++ state
                        pystate.fromPython(state);

                        // Do a release 
                        PyUnconstrained::Restart::X_Vectors xs;
                        PyUnconstrained::Restart::Reals reals;
                        PyUnconstrained::Restart::Naturals nats;
                        PyUnconstrained::Restart::Params params;
                        PyUnconstrained::Restart
                            ::release(state,xs,reals,nats,params);

                        // Convert the restart information to Python 
                        toPython::Vectors(xs,pyxs);
                        toPython::Reals(reals,pyreals);
                        toPython::Naturals(nats,pynats);
                        toPython::Params(params,pyparams);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }

                // Capture data from structures controlled by the user.  
                PyObject * capture(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be
                    // (X,msg,state,xs,reals,nats,params)
                    PyObject *X,*msg_,*pystate_,*pyxs,*pyreals,*pynats,
                        *pyparams;
                    if(!PyArg_ParseTuple(args,"OOOOOOO",
                        &X,&msg_,&pystate_,&pyxs,&pyreals,&pynats,&pyparams)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);

                        // Create a Python state 
                        Python::State <PyUnconstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));

                        // Create a C++ state
                        typename PyUnconstrained::State::t state(x);
                       
                        // Allocate memory for the released vectors
                        PyUnconstrained::Restart::X_Vectors xs;
                        PyUnconstrained::Restart::Reals reals;
                        PyUnconstrained::Restart::Naturals nats;
                        PyUnconstrained::Restart::Params params;
                        
                        // Convert the restart information from Python 
                        fromPython::Vectors(x,pyxs,xs);
                        fromPython::Reals(pyreals,reals);
                        fromPython::Naturals(pynats,nats);
                        fromPython::Params(pyparams,params);

                        // Do a capture 
                        PyUnconstrained::Restart
                            ::capture(msg,state,xs,reals,nats,params);

                        // Convert the C++ state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
                
                // Writes a json restart file
                PyObject * write_restart(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be (X,msg,fname,state)
                    PyObject *X,*msg_,*fname_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOO",
                        &X,&msg_,&fname_,&pystate_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the file name
                        std::string fname(PyString_AsString(fname_));

                        // Create a Python state 
                        Python::State <PyUnconstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        
                        // Create a C++ state
                        typename PyUnconstrained::State::t state(x);
                        
                        // Convert Python state to C++ 
                        pystate.fromPython(state);

                        // Write the restart file
                        PyJsonUnconstrained::write_restart(msg,fname,state);
                        
                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
                
                // Reads a json restart file
                PyObject * read_restart(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be (X,msg,fname,x,state)
                    PyObject *X,*msg_,*fname_,*x_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOOO",
                        &X,&msg_,&fname_,&x_,&pystate_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the file name
                        std::string fname(PyString_AsString(fname_));

                        // Create a Python state 
                        Python::State <PyUnconstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the reference vector 
                        Vector x(msg_,X,x_,PyObjectPtrMode::Attach);
                        
                        // Create a C++ state
                        typename PyUnconstrained::State::t state(x);

                        // Read the restart file into the C++ state 
                        PyJsonUnconstrained::read_restart(msg,fname,x,state);
                        
                        // Convert the C++ state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
            }
        }

        // Routines that manipulate and support problems of the form
        // 
        // min_{x \in X} f(x) st g(x) = 0
        //
        // where f : X -> R and g : X -> Y
        namespace EqualityConstrained {
            // Routines that manipulate the internal state of the optimization 
            // algorithm
            namespace State {
                // Convert a C++ state to a Python state 
                void toPython_(
                    typename PyEqualityConstrained::State::t const & state,
                    PyObject * const pystate
                ){
                    toPython::Vector("y",state.y,pystate);
                    toPython::Vector("dy",state.dy,pystate);
                    toPython::Real("zeta",state.zeta,pystate);
                    toPython::Real("eta0",state.eta0,pystate);
                    toPython::Real("rho",state.rho,pystate);
                    toPython::Real("rho_old",state.rho_old,pystate);
                    toPython::Real("rho_bar",state.rho_bar,pystate);
                    toPython::Real("eps_constr",state.eps_constr,pystate);
                    toPython::Real("xi_qn",state.xi_qn,pystate);
                    toPython::Real("xi_pg",state.xi_pg,pystate);
                    toPython::Real("xi_proj",state.xi_proj,pystate);
                    toPython::Real("xi_tang",state.xi_tang,pystate);
                    toPython::Real("xi_lmh",state.xi_lmh,pystate);
                    toPython::Real("xi_lmg",state.xi_lmg,pystate);
                    toPython::Real("xi_4",state.xi_4,pystate);
                    toPython::Real("rpred",state.rpred,pystate);
                    toPython::Param <Operators::t> (
                        "PSchur_left_type",
                        Operators::toPython,
                        state.PSchur_left_type,
                        pystate);
                    toPython::Param <Operators::t> (
                        "PSchur_right_type",
                        Operators::toPython,
                        state.PSchur_right_type,
                        pystate);
                    toPython::Natural("augsys_iter_max",
                        state.augsys_iter_max,pystate);
                    toPython::Natural("augsys_rst_freq",
                        state.augsys_rst_freq,pystate);
                    toPython::Natural("augsys_qn_iter",
                        state.augsys_qn_iter,pystate);
                    toPython::Natural("augsys_pg_iter",
                        state.augsys_pg_iter,pystate);
                    toPython::Natural("augsys_proj_iter",
                        state.augsys_proj_iter,pystate);
                    toPython::Natural("augsys_tang_iter",
                        state.augsys_tang_iter,pystate);
                    toPython::Natural("augsys_lmh_iter",
                        state.augsys_lmh_iter,pystate);
                    toPython::Natural("augsys_qn_iter_total",
                        state.augsys_qn_iter_total,pystate);
                    toPython::Natural("augsys_pg_iter_total",
                        state.augsys_pg_iter_total,pystate);
                    toPython::Natural("augsys_proj_iter_total",
                        state.augsys_proj_iter_total,pystate);
                    toPython::Natural("augsys_tang_iter_total",
                        state.augsys_tang_iter_total,pystate);
                    toPython::Natural("augsys_lmh_iter_total",
                        state.augsys_lmh_iter_total,pystate);
                    toPython::Real("augsys_qn_err",
                        state.augsys_qn_err,pystate);
                    toPython::Real("augsys_pg_err",
                        state.augsys_pg_err,pystate);
                    toPython::Real("augsys_proj_err",
                        state.augsys_proj_err,pystate);
                    toPython::Real("augsys_tang_err",
                        state.augsys_tang_err,pystate);
                    toPython::Real("augsys_lmh_err",
                        state.augsys_lmh_err,pystate);
                    toPython::Real("augsys_qn_err_target",
                        state.augsys_qn_err_target,pystate);
                    toPython::Real("augsys_pg_err_target",
                        state.augsys_pg_err_target,pystate);
                    toPython::Real("augsys_proj_err_target",
                        state.augsys_proj_err_target,pystate);
                    toPython::Real("augsys_tang_err_target",
                        state.augsys_tang_err_target,pystate);
                    toPython::Real("augsys_lmh_err_target",
                        state.augsys_lmh_err_target,pystate);
                    toPython::Natural("augsys_iter_total",
                        state.augsys_iter_total,pystate);
                    toPython::Vector("g_x",state.g_x,pystate);
                    toPython::Real("norm_gxtyp",state.norm_gxtyp,pystate);
                    toPython::Vector("gpxdxn_p_gx",state.gpxdxn_p_gx,pystate);
                    toPython::Vector("gpxdxt",state.gpxdxt,pystate);
                    toPython::Real("norm_gpxdxnpgx",
                        state.norm_gpxdxnpgx,pystate);
                    toPython::Vector("dx_n",state.dx_n,pystate);
                    toPython::Vector("dx_ncp",state.dx_ncp,pystate);
                    toPython::Vector("dx_t",state.dx_t,pystate);
                    toPython::Vector("dx_t_uncorrected",
                        state.dx_t_uncorrected,pystate);
                    toPython::Vector("dx_tcp_uncorrected",
                        state.dx_tcp_uncorrected,pystate);
                    toPython::Vector("H_dxn",state.H_dxn,pystate);
                    toPython::Vector("W_gradpHdxn",state.W_gradpHdxn,pystate);
                    toPython::Vector("H_dxtuncorrected",
                        state.H_dxtuncorrected,pystate);
                    toPython::Param <FunctionDiagnostics::t> (
                        "g_diag",
                        FunctionDiagnostics::toPython,
                        state.g_diag,
                        pystate);
                    toPython::Param <VectorSpaceDiagnostics::t> (
                        "y_diag",
                        VectorSpaceDiagnostics::toPython,
                        state.y_diag,
                        pystate);
                }
                void toPython(
                    typename PyEqualityConstrained::State::t const & state,
                    PyObject * const pystate
                ){
                    Unconstrained::State::toPython_(state,pystate);
                    EqualityConstrained::State::toPython_(state,pystate);
                }
                
                // Convert a Python state to C++ 
                void fromPython_(
                    PyObject * const pystate,
                    typename PyEqualityConstrained::State::t & state
                ){
                    fromPython::Vector("y",pystate,state.y);
                    fromPython::Vector("dy",pystate,state.dy);
                    fromPython::Real("zeta",pystate,state.zeta);
                    fromPython::Real("eta0",pystate,state.eta0);
                    fromPython::Real("rho",pystate,state.rho);
                    fromPython::Real("rho_old",pystate,state.rho_old);
                    fromPython::Real("rho_bar",pystate,state.rho_bar);
                    fromPython::Real("eps_constr",pystate,state.eps_constr);
                    fromPython::Real("xi_qn",pystate,state.xi_qn);
                    fromPython::Real("xi_pg",pystate,state.xi_pg);
                    fromPython::Real("xi_proj",pystate,state.xi_proj);
                    fromPython::Real("xi_tang",pystate,state.xi_tang);
                    fromPython::Real("xi_lmh",pystate,state.xi_lmh);
                    fromPython::Real("xi_lmg",pystate,state.xi_lmg);
                    fromPython::Real("xi_4",pystate,state.xi_4);
                    fromPython::Real("rpred",pystate,state.rpred);
                    fromPython::Param <Operators::t> (
                        "PSchur_left_type",
                        Operators::fromPython,
                        pystate,
                        state.PSchur_left_type);
                    fromPython::Param <Operators::t> (
                        "PSchur_right_type",
                        Operators::fromPython,
                        pystate,
                        state.PSchur_right_type);
                    fromPython::Natural("augsys_iter_max",
                        pystate,state.augsys_iter_max);
                    fromPython::Natural("augsys_rst_freq",
                        pystate,state.augsys_rst_freq);
                    fromPython::Natural("augsys_qn_iter",
                        pystate,state.augsys_qn_iter);
                    fromPython::Natural("augsys_pg_iter",
                        pystate,state.augsys_pg_iter);
                    fromPython::Natural("augsys_proj_iter",
                        pystate,state.augsys_proj_iter);
                    fromPython::Natural("augsys_tang_iter",
                        pystate,state.augsys_tang_iter);
                    fromPython::Natural("augsys_lmh_iter",
                        pystate,state.augsys_lmh_iter);
                    fromPython::Natural("augsys_qn_iter_total",
                        pystate,state.augsys_qn_iter_total);
                    fromPython::Natural("augsys_pg_iter_total",
                        pystate,state.augsys_pg_iter_total);
                    fromPython::Natural("augsys_proj_iter_total",
                        pystate,state.augsys_proj_iter_total);
                    fromPython::Natural("augsys_tang_iter_total",
                        pystate,state.augsys_tang_iter_total);
                    fromPython::Natural("augsys_lmh_iter_total",
                        pystate,state.augsys_lmh_iter_total);
                    fromPython::Real("augsys_qn_err",
                        pystate,state.augsys_qn_err);
                    fromPython::Real("augsys_pg_err",
                        pystate,state.augsys_pg_err);
                    fromPython::Real("augsys_proj_err",
                        pystate,state.augsys_proj_err);
                    fromPython::Real("augsys_tang_err",
                        pystate,state.augsys_tang_err);
                    fromPython::Real("augsys_lmh_err",
                        pystate,state.augsys_lmh_err);
                    fromPython::Real("augsys_qn_err_target",
                        pystate,state.augsys_qn_err_target);
                    fromPython::Real("augsys_pg_err_target",
                        pystate,state.augsys_pg_err_target);
                    fromPython::Real("augsys_proj_err_target",
                        pystate,state.augsys_proj_err_target);
                    fromPython::Real("augsys_tang_err_target",
                        pystate,state.augsys_tang_err_target);
                    fromPython::Real("augsys_lmh_err_target",
                        pystate,state.augsys_lmh_err_target);
                    fromPython::Natural("augsys_iter_total",
                        pystate,state.augsys_iter_total);
                    fromPython::Vector("g_x",pystate,state.g_x);
                    fromPython::Real("norm_gxtyp",pystate,state.norm_gxtyp);
                    fromPython::Vector("gpxdxn_p_gx",pystate,state.gpxdxn_p_gx);
                    fromPython::Vector("gpxdxt",pystate,state.gpxdxt);
                    fromPython::Real("norm_gpxdxnpgx",
                        pystate,state.norm_gpxdxnpgx);
                    fromPython::Vector("dx_n",pystate,state.dx_n);
                    fromPython::Vector("dx_ncp",pystate,state.dx_ncp);
                    fromPython::Vector("dx_t",pystate,state.dx_t);
                    fromPython::Vector("dx_t_uncorrected",
                        pystate,state.dx_t_uncorrected);
                    fromPython::Vector("dx_tcp_uncorrected",
                        pystate,state.dx_tcp_uncorrected);
                    fromPython::Vector("H_dxn",pystate,state.H_dxn);
                    fromPython::Vector("W_gradpHdxn",pystate,state.W_gradpHdxn);
                    fromPython::Vector("H_dxtuncorrected",
                        pystate,state.H_dxtuncorrected);
                    fromPython::Param <FunctionDiagnostics::t> (
                        "g_diag",
                        FunctionDiagnostics::fromPython,
                        pystate,
                        state.g_diag);
                    fromPython::Param <VectorSpaceDiagnostics::t> (
                        "y_diag",
                        VectorSpaceDiagnostics::fromPython,
                        pystate,
                        state.y_diag);
                }
                void fromPython(
                    PyObject * const pystate,
                    typename PyEqualityConstrained::State::t & state
                ){
                    Unconstrained::State::fromPython_(pystate,state);
                    EqualityConstrained::State::fromPython_(pystate,state);
                }
        
                // Creates a state and inserts the elements into pystate 
                PyObject * create(
                    PyObject * self,
                    PyObject * args
                ){
                    // Calling convention should be (pystate,X,Y,msg,x,y) 
                    PyObject *pystate_,*X,*Y,*msg,*x_,*y_;
                    if(!PyArg_ParseTuple(
                        args,"OOOOOO",&pystate_,&X,&Y,&msg,&x_,&y_)
                    )
                        return nullptr; 


                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create vectors from the user input
                        Vector x(msg,X,x_,PyObjectPtrMode::Attach);
                        Vector y(msg,Y,y_,PyObjectPtrMode::Attach);

                        // Create a Python state 
                        Python::State <PyEqualityConstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);

                        // Create a new C++ state
                        typename PyEqualityConstrained::State::t state(x,y);

                        // Convert the state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
        
                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be (X,Y,msg,fname,pystate) 
                    PyObject *X,*Y,*msg_,*fname_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOOO",
                        &X,&Y,&msg_,&fname_,&pystate_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Grab the file name
                        std::string fname(PyString_AsString(fname_));

                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);

                        // Create a Python state 
                        Python::State <PyEqualityConstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                    
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        Vector y(msg_,Y,
                            PyObject_GetAttrString(pystate.get(),"y"));

                        // Create a new C++ state
                        typename PyEqualityConstrained::State::t state(x,y);
                        
                        // Convert the Python state to a C++ state
                        pystate.fromPython(state);

                        // Read the JSON file into the C++ state
                        PyJsonEqualityConstrained::read(msg,fname,state);

                        // Convert the C++ state to a Python state
                        pystate.toPython(state);
                                
                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Python bundle to C++ 
                void fromPython(
                    PyObject * const msg,
                    PyObject * const pyfns,
                    PyObject * const pystate,
                    typename PyEqualityConstrained::State::t const & state,
                    typename PyEqualityConstrained::Functions::t & fns 
                ) {
                    Unconstrained::Functions::fromPython_
                        <PyEqualityConstrained> (msg,pyfns,pystate,state,fns);
                    EqualityConstrained::Functions::fromPython_
                        <PyEqualityConstrained> (msg,pyfns,pystate,state,fns);
                }
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem 
                PyObject * getMin(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be
                    // (X,Y,msg,fns,state,smanip)
                    PyObject *X,*Y,*msg_,*pyfns_,*pystate_,*smanip_;
                    if(!PyArg_ParseTuple(args,"OOOOOO",
                        &X,&Y,&msg_,&pyfns_,&pystate_,&smanip_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);
                            
                        // Create a Python state 
                        Python::State <PyEqualityConstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        Vector y(msg_,Y,
                            PyObject_GetAttrString(pystate.get(),"y"));

                        // Create a C++ state
                        typename PyEqualityConstrained::State::t state(x,y);
                        
                        // Convert the Python state to a C++ state
                        pystate.fromPython(state);

                        // Create a Python bundle of functions
                        Python::Functions <PyEqualityConstrained> pyfns(
                            msg.get(),
                            pystate.get(),
                            state,
                            pyfns_,
                            PyObjectPtrMode::Attach);

                        // Create a C++ bundle of functions
                        typename PyEqualityConstrained::Functions::t fns;
                        
                        // Convert the Python bundle of functions to C++ 
                        pyfns.fromPython(fns);
                        
                        // Create a state manipulator 
                        Python::StateManipulator <PyEqualityConstrained> smanip(
                            msg.get(),
                            pystate.get(),
                            pyfns.get(),
                            smanip_,
                            PyObjectPtrMode::Attach);
                       
                        // Minimize
                        PyEqualityConstrained::Algorithms::getMin(
                            msg,fns,state,smanip);
                        
                        // Convert the C++ state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
            }
            
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                PyObject * release(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be
                    // (X,Y,msg,state,xs,ys,reals,nats,params)
                    PyObject *X,*Y,*msg,*pystate_,*pyxs,*pyys,*pyreals,
                        *pynats,*pyparams;
                    if(!PyArg_ParseTuple(args,"OOOOOOOOO",
                        &X,&Y,&msg,&pystate_,&pyxs,&pyys,&pyreals,
                        &pynats,&pyparams)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a Python state 
                        Python::State <PyEqualityConstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        Vector y(msg,Y,
                            PyObject_GetAttrString(pystate.get(),"y"));

                        // Create a C++ state
                        typename PyEqualityConstrained::State::t state(x,y);
                        
                        // Convert the Python state to a C++ state
                        pystate.fromPython(state);

                        // Do a release 
                        PyEqualityConstrained::Restart::X_Vectors xs;
                        PyEqualityConstrained::Restart::Y_Vectors ys;
                        PyEqualityConstrained::Restart::Reals reals;
                        PyEqualityConstrained::Restart::Naturals nats;
                        PyEqualityConstrained::Restart::Params params;
                        PyEqualityConstrained::Restart
                            ::release(state,xs,ys,reals,nats,params);

                        // Convert the restart information to Python 
                        toPython::Vectors(xs,pyxs);
                        toPython::Vectors(ys,pyys);
                        toPython::Reals(reals,pyreals);
                        toPython::Naturals(nats,pynats);
                        toPython::Params(params,pyparams);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }

                // Capture data from structures controlled by the user.  
                PyObject * capture(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be
                    // (X,Y,msg,state,xs,ys,reals,nats,params)
                    PyObject *X,*Y,*msg_,*pystate_,*pyxs,*pyys,
                        *pyreals,*pynats,*pyparams;
                    if(!PyArg_ParseTuple(args,"OOOOOOOOO",
                        &X,&Y,&msg_,&pystate_,&pyxs,&pyys,
                        &pyreals,&pynats,&pyparams)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);

                        // Create a Python state 
                        Python::State <PyEqualityConstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        Vector y(msg_,Y,
                            PyObject_GetAttrString(pystate.get(),"y"));

                        // Create a C++ state
                        typename PyEqualityConstrained::State::t state(x,y);
                       
                        // Allocate memory for the released vectors
                        PyEqualityConstrained::Restart::X_Vectors xs;
                        PyEqualityConstrained::Restart::Y_Vectors ys;
                        PyEqualityConstrained::Restart::Reals reals;
                        PyEqualityConstrained::Restart::Naturals nats;
                        PyEqualityConstrained::Restart::Params params;
                        
                        // Convert the restart information from Python 
                        fromPython::Vectors(x,pyxs,xs);
                        fromPython::Vectors(y,pyys,ys);
                        fromPython::Reals(pyreals,reals);
                        fromPython::Naturals(pynats,nats);
                        fromPython::Params(pyparams,params);

                        // Do a capture 
                        PyEqualityConstrained::Restart
                            ::capture(msg,state,xs,ys,reals,nats,params);

                        // Convert the C++ state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
                
                // Writes a json restart file
                PyObject * write_restart(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be (X,Y,msg,fname,state)
                    PyObject *X,*Y,*msg_,*fname_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOOO",
                        &X,&Y,&msg_,&fname_,&pystate_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the file name
                        std::string fname(PyString_AsString(fname_));

                        // Create a Python state 
                        Python::State <PyEqualityConstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        Vector y(msg_,Y,
                            PyObject_GetAttrString(pystate.get(),"y"));
                        
                        // Create a C++ state
                        typename PyEqualityConstrained::State::t state(x,y);
                        
                        // Convert Python state to C++ 
                        pystate.fromPython(state);

                        // Write the restart file
                        PyJsonEqualityConstrained::write_restart(
                            msg,fname,state);
                        
                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
                
                // Reads a json restart file
                PyObject * read_restart(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be (X,Y,msg,fname,x,y,state)
                    PyObject *X,*Y,*msg_,*fname_,*x_,*y_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOOOOO",
                        &X,&Y,&msg_,&fname_,&x_,&y_,&pystate_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the file name
                        std::string fname(PyString_AsString(fname_));

                        // Create a Python state 
                        Python::State <PyEqualityConstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the reference vector 
                        Vector x(msg_,X,x_,PyObjectPtrMode::Attach);
                        Vector y(msg_,Y,y_,PyObjectPtrMode::Attach);
                        
                        // Create a C++ state
                        typename PyEqualityConstrained::State::t state(x,y);

                        // Read the restart file into the C++ state 
                        PyJsonEqualityConstrained::read_restart(
                            msg,fname,x,y,state);
                        
                        // Convert the C++ state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
            }
        }
        
        // Routines that manipulate and support problems of the form
        // 
        // min_{x \in X} f(x) st h(x) >=_K 0
        //
        // where f : X -> R and h : X -> Z
        namespace InequalityConstrained {
            // Routines that manipulate the internal state of the optimization 
            // algorithm
            namespace State {
                // Convert a C++ state to a Python state 
                void toPython_(
                    typename PyInequalityConstrained::State::t const & state,
                    PyObject * const pystate
                ){
                    toPython::Vector("z",state.z,pystate);
                    toPython::Vector("dz",state.dz,pystate);
                    toPython::Vector("h_x",state.h_x,pystate);
                    toPython::Real("mu",state.mu,pystate);
                    toPython::Real("mu_est",state.mu_est,pystate);
                    toPython::Real("mu_typ",state.mu_typ,pystate);
                    toPython::Real("eps_mu",state.eps_mu,pystate);
                    toPython::Real("sigma",state.sigma,pystate);
                    toPython::Real("gamma",state.gamma,pystate);
                    toPython::Real("alpha_z",state.alpha_z,pystate);
                    toPython::Param <FunctionDiagnostics::t> (
                        "h_diag",
                        FunctionDiagnostics::toPython,
                        state.h_diag,
                        pystate);
                    toPython::Param <VectorSpaceDiagnostics::t> (
                        "z_diag",
                        VectorSpaceDiagnostics::toPython,
                        state.z_diag,
                        pystate);
                }
                void toPython(
                    typename PyInequalityConstrained::State::t const & state,
                    PyObject * const pystate
                ){
                    Unconstrained::State::toPython_(state,pystate);
                    InequalityConstrained::State::toPython_(state,pystate);
                }
                
                // Convert a Python state to C++ 
                void fromPython_(
                    PyObject * const pystate,
                    typename PyInequalityConstrained::State::t & state
                ){
                    fromPython::Vector("z",pystate,state.z);
                    fromPython::Vector("dz",pystate,state.dz);
                    fromPython::Vector("h_x",pystate,state.h_x);
                    fromPython::Real("mu",pystate,state.mu);
                    fromPython::Real("mu_est",pystate,state.mu_est);
                    fromPython::Real("mu_typ",pystate,state.mu_typ);
                    fromPython::Real("eps_mu",pystate,state.eps_mu);
                    fromPython::Real("sigma",pystate,state.sigma);
                    fromPython::Real("gamma",pystate,state.gamma);
                    fromPython::Real("alpha_z",pystate,state.alpha_z);
                    fromPython::Param <FunctionDiagnostics::t> (
                        "h_diag",
                        FunctionDiagnostics::fromPython,
                        pystate,
                        state.h_diag);
                    fromPython::Param <VectorSpaceDiagnostics::t> (
                        "z_diag",
                        VectorSpaceDiagnostics::fromPython,
                        pystate,
                        state.z_diag);
                }
                void fromPython(
                    PyObject * const pystate,
                    typename PyInequalityConstrained::State::t & state
                ){
                    Unconstrained::State::fromPython_(pystate,state);
                    InequalityConstrained::State::fromPython_(pystate,state);
                }

                // Creates a state and inserts the elements into pystate 
                PyObject * create( 
                    PyObject * self,
                    PyObject * args
                ){
                    // Calling convention should be (pystate,X,Z,msg,x,z) 
                    PyObject *pystate_,*X,*Z,*msg,*x_,*z_;
                    if(!PyArg_ParseTuple(
                        args,"OOOOOO",&pystate_,&X,&Z,&msg,&x_,&z_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create vectors from the user input
                        Vector x(msg,X,x_,PyObjectPtrMode::Attach);
                        Vector z(msg,Z,z_,PyObjectPtrMode::Attach);

                        // Create a Python state 
                        Python::State<PyInequalityConstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);

                        // Create a new C++ state
                        typename PyInequalityConstrained::State::t state(x,z);

                        // Convert the state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
        
                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be (X,Z,msg,fname,pystate) 
                    PyObject *X,*Z,*msg_,*fname_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOOO",
                        &X,&Z,&msg_,&fname_,&pystate_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Grab the file name
                        std::string fname(PyString_AsString(fname_));

                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);

                        // Create a Python state 
                        Python::State<PyInequalityConstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                    
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        Vector z(msg_,Z,
                            PyObject_GetAttrString(pystate.get(),"z"));

                        // Create a new C++ state
                        typename PyInequalityConstrained::State::t state(x,z);
                        
                        // Convert the Python state to a C++ state
                        pystate.fromPython(state);

                        // Read the JSON file into the C++ state
                        PyJsonInequalityConstrained::read(msg,fname,state);

                        // Convert the C++ state to a Python state
                        pystate.toPython(state);
                                
                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
            }
            
            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Python bundle to C++ 
                void fromPython(
                    PyObject * const msg,
                    PyObject * const pyfns,
                    PyObject * const pystate,
                    typename PyInequalityConstrained::State::t const & state,
                    typename PyInequalityConstrained::Functions::t & fns 
                ) {
                    Unconstrained::Functions::fromPython_
                        <PyInequalityConstrained> (msg,pyfns,pystate,state,fns);
                    InequalityConstrained::Functions::fromPython_
                        <PyInequalityConstrained> (msg,pyfns,pystate,state,fns);
                }
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem 
                PyObject * getMin(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be
                    // (X,Z,msg,fns,state,smanip)
                    PyObject *X,*Z,*msg_,*pyfns_,*pystate_,*smanip_;
                    if(!PyArg_ParseTuple(args,"OOOOOO",
                        &X,&Z,&msg_,&pyfns_,&pystate_,&smanip_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);
                            
                        // Create a Python state 
                        Python::State<PyInequalityConstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        Vector z(msg_,Z,
                            PyObject_GetAttrString(pystate.get(),"z"));

                        // Create a C++ state
                        typename PyInequalityConstrained::State::t state(x,z);
                        
                        // Convert the Python state to a C++ state
                        pystate.fromPython(state);

                        // Create a Python bundle of functions
                        Python::Functions <PyInequalityConstrained> pyfns(
                            msg.get(),
                            pystate.get(),
                            state,
                            pyfns_,
                            PyObjectPtrMode::Attach);

                        // Create a C++ bundle of functions
                        typename PyInequalityConstrained::Functions::t fns;
                        
                        // Convert the Python bundle of functions to C++ 
                        pyfns.fromPython(fns);
                        
                        // Create a state manipulator 
                        Python::StateManipulator<PyInequalityConstrained>smanip(
                            msg.get(),
                            pystate.get(),
                            pyfns.get(),
                            smanip_,
                            PyObjectPtrMode::Attach);
                       
                        // Minimize
                        PyInequalityConstrained::Algorithms::getMin(
                            msg,fns,state,smanip);
                        
                        // Convert the C++ state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
            }
            
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                PyObject * release(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be
                    // (X,Z,msg,state,xs,zs,reals,nats,params)
                    PyObject *X,*Z,*msg,*pystate_,*pyxs,*pyzs,*pyreals,
                        *pynats,*pyparams;
                    if(!PyArg_ParseTuple(args,"OOOOOOOOO",
                        &X,&Z,&msg,&pystate_,&pyxs,&pyzs,&pyreals,
                        &pynats,&pyparams)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a Python state 
                        Python::State <PyInequalityConstrained>pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        Vector z(msg,Z,
                            PyObject_GetAttrString(pystate.get(),"z"));

                        // Create a C++ state
                        typename PyInequalityConstrained::State::t state(x,z);
                        
                        // Convert the Python state to a C++ state
                        pystate.fromPython(state);

                        // Do a release 
                        PyInequalityConstrained::Restart::X_Vectors xs;
                        PyInequalityConstrained::Restart::Z_Vectors zs;
                        PyInequalityConstrained::Restart::Reals reals;
                        PyInequalityConstrained::Restart::Naturals nats;
                        PyInequalityConstrained::Restart::Params params;
                        PyInequalityConstrained::Restart
                            ::release(state,xs,zs,reals,nats,params);

                        // Convert the restart information to Python 
                        toPython::Vectors(xs,pyxs);
                        toPython::Vectors(zs,pyzs);
                        toPython::Reals(reals,pyreals);
                        toPython::Naturals(nats,pynats);
                        toPython::Params(params,pyparams);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }

                // Capture data from structures controlled by the user.  
                PyObject * capture(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be
                    // (X,Y,msg,state,xs,ys,reals,nats,params)
                    PyObject *X,*Z,*msg_,*pystate_,*pyxs,*pyzs,
                        *pyreals,*pynats,*pyparams;
                    if(!PyArg_ParseTuple(args,"OOOOOOOOO",
                        &X,&Z,&msg_,&pystate_,&pyxs,&pyzs,
                        &pyreals,&pynats,&pyparams)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);

                        // Create a Python state 
                        Python::State <PyInequalityConstrained>pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        Vector z(msg_,Z,
                            PyObject_GetAttrString(pystate.get(),"z"));

                        // Create a C++ state
                        typename PyInequalityConstrained::State::t state(x,z);
                       
                        // Allocate memory for the released vectors
                        PyInequalityConstrained::Restart::X_Vectors xs;
                        PyInequalityConstrained::Restart::Z_Vectors zs;
                        PyInequalityConstrained::Restart::Reals reals;
                        PyInequalityConstrained::Restart::Naturals nats;
                        PyInequalityConstrained::Restart::Params params;
                        
                        // Convert the restart information from Python 
                        fromPython::Vectors(x,pyxs,xs);
                        fromPython::Vectors(z,pyzs,zs);
                        fromPython::Reals(pyreals,reals);
                        fromPython::Naturals(pynats,nats);
                        fromPython::Params(pyparams,params);

                        // Do a capture 
                        PyInequalityConstrained::Restart
                            ::capture(msg,state,xs,zs,reals,nats,params);

                        // Convert the C++ state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
                
                // Writes a json restart file
                PyObject * write_restart(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be (X,Z,msg,fname,state)
                    PyObject *X,*Z,*msg_,*fname_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOOO",
                        &X,&Z,&msg_,&fname_,&pystate_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the file name
                        std::string fname(PyString_AsString(fname_));

                        // Create a Python state 
                        Python::State <PyInequalityConstrained>pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        Vector z(msg_,Z,
                            PyObject_GetAttrString(pystate.get(),"z"));
                        
                        // Create a C++ state
                        typename PyInequalityConstrained::State::t state(x,z);
                        
                        // Convert Python state to C++ 
                        pystate.fromPython(state);

                        // Write the restart file
                        PyJsonInequalityConstrained::write_restart(
                            msg,fname,state);
                        
                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
                
                // Reads a json restart file
                PyObject * read_restart(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be (X,Z,msg,fname,x,z,state)
                    PyObject *X,*Z,*msg_,*fname_,*x_,*z_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOOOOO",
                        &X,&Z,&msg_,&fname_,&x_,&z_,&pystate_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the file name
                        std::string fname(PyString_AsString(fname_));

                        // Create a Python state 
                        Python::State <PyInequalityConstrained>pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the reference vector 
                        Vector x(msg_,X,x_,PyObjectPtrMode::Attach);
                        Vector z(msg_,Z,z_,PyObjectPtrMode::Attach);
                        
                        // Create a C++ state
                        typename PyInequalityConstrained::State::t state(x,z);

                        // Read the restart file into the C++ state 
                        PyJsonInequalityConstrained::read_restart(
                            msg,fname,x,z,state);
                        
                        // Convert the C++ state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
            }
        }
        
        // Routines that manipulate and support problems of the form
        // problem of the form
        // 
        // min_{x \in X} f(x) st g(x) = 0, h(x) >=_K 0
        //
        // where f : X -> R, g : X -> Y, and h : X -> Z
        namespace Constrained {
            // Routines that manipulate the internal state of the optimization 
            // algorithm
            namespace State {
                // Convert a C++ state to a Python state 
                void toPython(
                    typename PyConstrained::State::t const & state,
                    PyObject * const pystate
                ){
                    Unconstrained::State::toPython_(state,pystate);
                    EqualityConstrained::State::toPython_(state,pystate);
                    InequalityConstrained::State::toPython_(state,pystate);
                }
                
                // Convert a Python state to C++ 
                void fromPython(
                    PyObject * const pystate,
                    typename PyConstrained::State::t & state
                ){
                    Unconstrained::State::fromPython_(pystate,state);
                    EqualityConstrained::State::fromPython_(pystate,state);
                    InequalityConstrained::State::fromPython_(pystate,state);
                }

                // Creates a state and inserts the elements into pystate 
                PyObject * create(
                    PyObject * self,
                    PyObject * args
                ){
                    // Calling convention should be (pystate,X,Y,Z,msg,x,y,z) 
                    PyObject *pystate_,*X,*Y,*Z,*msg,*x_,*y_,*z_;
                    if(!PyArg_ParseTuple(args,"OOOOOOOO",
                        &pystate_,&X,&Y,&Z,&msg,&x_,&y_,&z_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create vectors from the user input
                        Vector x(msg,X,x_,PyObjectPtrMode::Attach);
                        Vector y(msg,Y,y_,PyObjectPtrMode::Attach);
                        Vector z(msg,Z,z_,PyObjectPtrMode::Attach);

                        // Create a Python state 
                        Python::State <PyConstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);

                        // Create a new C++ state
                        typename PyConstrained::State::t state(x,y,z);

                        // Convert the state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }

                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be (X,Y,Z,msg,fname,pystate) 
                    PyObject *X,*Y,*Z,*msg_,*fname_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOOOO",
                        &X,&Y,&Z,&msg_,&fname_,&pystate_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Grab the file name
                        std::string fname(PyString_AsString(fname_));

                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);

                        // Create a Python state 
                        Python::State <PyConstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                    
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        Vector y(msg_,Y,
                            PyObject_GetAttrString(pystate.get(),"y"));
                        Vector z(msg_,Z,
                            PyObject_GetAttrString(pystate.get(),"z"));

                        // Create a new C++ state
                        typename PyConstrained::State::t state(x,y,z);
                        
                        // Convert the Python state to a C++ state
                        pystate.fromPython(state);

                        // Read the JSON file into the C++ state
                        PyJsonConstrained::read(msg,fname,state);

                        // Convert the C++ state to a Python state
                        pystate.toPython(state);
                                
                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
            }
            
            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Python bundle to C++ 
                void fromPython(
                    PyObject * const msg,
                    PyObject * const pyfns,
                    PyObject * const pystate,
                    typename PyConstrained::State::t const & state,
                    typename PyConstrained::Functions::t & fns 
                ) {
                    Unconstrained::Functions::fromPython_
                        <PyConstrained> (msg,pyfns,pystate,state,fns);
                    EqualityConstrained::Functions::fromPython_
                        <PyConstrained> (msg,pyfns,pystate,state,fns);
                    InequalityConstrained::Functions::fromPython_
                        <PyConstrained> (msg,pyfns,pystate,state,fns);
                }
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem 
                PyObject * getMin(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be
                    // (X,Y,Z,msg,fns,state,smanip)
                    PyObject *X,*Y,*Z,*msg_,*pyfns_,*pystate_,*smanip_;
                    if(!PyArg_ParseTuple(args,"OOOOOOO",
                        &X,&Y,&Z,&msg_,&pyfns_,&pystate_,&smanip_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);
                            
                        // Create a Python state 
                        Python::State<PyConstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        Vector y(msg_,Y,
                            PyObject_GetAttrString(pystate.get(),"y"));
                        Vector z(msg_,Z,
                            PyObject_GetAttrString(pystate.get(),"z"));

                        // Create a C++ state
                        typename PyConstrained::State::t state(x,y,z);
                        
                        // Convert the Python state to a C++ state
                        pystate.fromPython(state);

                        // Create a Python bundle of functions
                        Python::Functions <PyConstrained> pyfns(
                            msg.get(),
                            pystate.get(),
                            state,
                            pyfns_,
                            PyObjectPtrMode::Attach);

                        // Create a C++ bundle of functions
                        typename PyConstrained::Functions::t fns;
                        
                        // Convert the Python bundle of functions to C++ 
                        pyfns.fromPython(fns);
                        
                        // Create a state manipulator 
                        Python::StateManipulator<PyConstrained>smanip(
                            msg.get(),
                            pystate.get(),
                            pyfns.get(),
                            smanip_,
                            PyObjectPtrMode::Attach);
                       
                        // Minimize
                        PyConstrained::Algorithms::getMin(
                            msg,fns,state,smanip);
                        
                        // Convert the C++ state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
            }

            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                PyObject * release(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be
                    // (X,Y,Z,msg,state,xs,ys,zs,reals,nats,params)
                    PyObject *X,*Y,*Z,*msg,*pystate_,*pyxs,*pyys,*pyzs,
                        *pyreals,*pynats,*pyparams;
                    if(!PyArg_ParseTuple(args,"OOOOOOOOOOO",
                        &X,&Y,&Z,&msg,&pystate_,&pyxs,&pyys,&pyzs,
                        &pyreals,&pynats,&pyparams)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a Python state 
                        Python::State <PyConstrained> pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        Vector y(msg,Y,
                            PyObject_GetAttrString(pystate.get(),"y"));
                        Vector z(msg,Z,
                            PyObject_GetAttrString(pystate.get(),"z"));

                        // Create a C++ state
                        typename PyConstrained::State::t state(x,y,z);
                        
                        // Convert the Python state to a C++ state
                        pystate.fromPython(state);

                        // Do a release 
                        PyConstrained::Restart::X_Vectors xs;
                        PyConstrained::Restart::Y_Vectors ys;
                        PyConstrained::Restart::Z_Vectors zs;
                        PyConstrained::Restart::Reals reals;
                        PyConstrained::Restart::Naturals nats;
                        PyConstrained::Restart::Params params;
                        PyConstrained::Restart
                            ::release(state,xs,ys,zs,reals,nats,params);

                        // Convert the restart information to Python 
                        toPython::Vectors(xs,pyxs);
                        toPython::Vectors(ys,pyys);
                        toPython::Vectors(zs,pyzs);
                        toPython::Reals(reals,pyreals);
                        toPython::Naturals(nats,pynats);
                        toPython::Params(params,pyparams);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }

                // Capture data from structures controlled by the user.  
                PyObject * capture(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be
                    // (X,Y,msg,state,xs,ys,reals,nats,params)
                    PyObject *X,*Y,*Z,*msg_,*pystate_,*pyxs,*pyys,*pyzs,
                        *pyreals,*pynats,*pyparams;
                    if(!PyArg_ParseTuple(args,"OOOOOOOOOOO",
                        &X,&Y,&Z,&msg_,&pystate_,&pyxs,&pyys,&pyzs,
                        &pyreals,&pynats,&pyparams)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);

                        // Create a Python state 
                        Python::State <PyConstrained>pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        Vector y(msg_,Y,
                            PyObject_GetAttrString(pystate.get(),"y"));
                        Vector z(msg_,Z,
                            PyObject_GetAttrString(pystate.get(),"z"));

                        // Create a C++ state
                        typename PyConstrained::State::t state(x,y,z);
                       
                        // Allocate memory for the released vectors
                        PyConstrained::Restart::X_Vectors xs;
                        PyConstrained::Restart::Y_Vectors ys;
                        PyConstrained::Restart::Z_Vectors zs;
                        PyConstrained::Restart::Reals reals;
                        PyConstrained::Restart::Naturals nats;
                        PyConstrained::Restart::Params params;
                        
                        // Convert the restart information from Python 
                        fromPython::Vectors(x,pyxs,xs);
                        fromPython::Vectors(y,pyys,ys);
                        fromPython::Vectors(z,pyzs,zs);
                        fromPython::Reals(pyreals,reals);
                        fromPython::Naturals(pynats,nats);
                        fromPython::Params(pyparams,params);

                        // Do a capture 
                        PyConstrained::Restart
                            ::capture(msg,state,xs,ys,zs,reals,nats,params);

                        // Convert the C++ state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
                
                // Writes a json restart file
                PyObject * write_restart(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be (X,Y,Z,msg,fname,state)
                    PyObject *X,*Y,*Z,*msg_,*fname_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOOOO",
                        &X,&Y,&Z,&msg_,&fname_,&pystate_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the file name
                        std::string fname(PyString_AsString(fname_));

                        // Create a Python state 
                        Python::State <PyConstrained>pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the base vectors from the Python state
                        Vector x(msg_,X,
                            PyObject_GetAttrString(pystate.get(),"x"));
                        Vector y(msg_,Y,
                            PyObject_GetAttrString(pystate.get(),"y"));
                        Vector z(msg_,Z,
                            PyObject_GetAttrString(pystate.get(),"z"));
                        
                        // Create a C++ state
                        typename PyConstrained::State::t state(x,y,z);
                        
                        // Convert Python state to C++ 
                        pystate.fromPython(state);

                        // Write the restart file
                        PyJsonConstrained::write_restart(msg,fname,state);
                        
                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
                
                // Reads a json restart file
                PyObject * read_restart(
                    PyObject * self,
                    PyObject * args
                ) {
                    // Calling convention should be(X,Y,Z,msg,fname,x,y,z,state)
                    PyObject *X,*Y,*Z,*msg_,*fname_,*x_,*y_,*z_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOOOOOOO",
                        &X,&Y,&Z,&msg_,&fname_,&x_,&y_,&z_,&pystate_)
                    )
                        return nullptr; 

                    // Make sure we bail if we detect a Python exception
                    try {
                        // Create a messaging object
                        Optizelle::Python::Messaging msg(msg_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the file name
                        std::string fname(PyString_AsString(fname_));

                        // Create a Python state 
                        Python::State <PyConstrained>pystate(pystate_,
                            PyObjectPtrMode::Attach);
                        
                        // Grab the reference vector 
                        Vector x(msg_,X,x_,PyObjectPtrMode::Attach);
                        Vector y(msg_,Y,y_,PyObjectPtrMode::Attach);
                        Vector z(msg_,Z,z_,PyObjectPtrMode::Attach);
                        
                        // Create a C++ state
                        typename PyConstrained::State::t state(x,y,z);

                        // Read the restart file into the C++ state 
                        PyJsonConstrained::read_restart(
                            msg,fname,x,y,z,state);
                        
                        // Convert the C++ state to a Python state
                        pystate.toPython(state);

                        // Return nothing 
                        return Py_None; 

                    // In theory, we should have set the appropriate error
                    } catch (Exception& exc){
                        return nullptr;
                    }
                }
            }
        }
    }
}

// Collects all the methods 
PyMethodDef methods[] = {
    { const_cast <char*> ("UnconstrainedStateCreate"),
        (PyCFunction)Optizelle::Python::Unconstrained::State::create,
        METH_VARARGS,
        const_cast <char*> ("Creates an unconstrained state")},

    { const_cast <char*> ("UnconstrainedStateReadJson"),
        (PyCFunction)Optizelle::Python::Unconstrained::State::readJson,
        METH_VARARGS,
        const_cast <char*> ("Reads unconstrained state parameters from file")},

    { const_cast <char*> ("UnconstrainedAlgorithmsGetMin"),
        (PyCFunction)Optizelle::Python::Unconstrained::Algorithms::getMin,
        METH_VARARGS,
        const_cast <char*> ("Solves an unconstrained optimization problem")},

    { const_cast <char*> ("UnconstrainedRestartRelease"),
        (PyCFunction)Optizelle::Python::Unconstrained::Restart::release,
        METH_VARARGS,
        const_cast <char*> (
            "Release the state in an unconstrained optimization problem")},

    { const_cast <char*> ("UnconstrainedRestartCapture"),
        (PyCFunction)Optizelle::Python::Unconstrained::Restart::capture,
        METH_VARARGS,
        const_cast <char*> (
            "Capture the state in an unconstrained optimization problem")},

    { const_cast <char*> ("UnconstrainedRestartWriteRestart"),
        (PyCFunction)Optizelle::Python::Unconstrained::Restart::write_restart,
        METH_VARARGS,
        const_cast <char*> ("Writes a json restart file")},

    { const_cast <char*> ("UnconstrainedRestartReadRestart"),
        (PyCFunction)Optizelle::Python::Unconstrained::Restart::read_restart,
        METH_VARARGS,
        const_cast <char*> ("Reads a json restart file")},
        
    { const_cast <char*> ("EqualityConstrainedStateCreate"),
        (PyCFunction)Optizelle::Python::EqualityConstrained::State::create,
        METH_VARARGS,
        const_cast <char*> ("Creates an equality constrained state")},

    { const_cast <char*> ("EqualityConstrainedStateReadJson"),
        (PyCFunction)Optizelle::Python::EqualityConstrained::State::readJson,
        METH_VARARGS,
        const_cast <char*> (
            "Reads equality constrained state parameters from file")},

    { const_cast <char*> ("EqualityConstrainedAlgorithmsGetMin"),
        (PyCFunction)Optizelle::Python::EqualityConstrained::Algorithms::getMin,
        METH_VARARGS,
        const_cast <char*> (
            "Solves an equality constrained optimization problem")},
    
    { const_cast <char*> ("EqualityConstrainedRestartRelease"),
        (PyCFunction)Optizelle::Python::EqualityConstrained::Restart::release,
        METH_VARARGS,
        const_cast <char*> (
            "Release the state in an equality constrained optimization "
            "problem")},

    { const_cast <char*> ("EqualityConstrainedRestartCapture"),
        (PyCFunction)Optizelle::Python::EqualityConstrained::Restart::capture,
        METH_VARARGS,
        const_cast <char*> (
            "Capture the state in an equality constrained optimization "
            "problem")},
    
    { const_cast <char*> ("EqualityConstrainedRestartWriteRestart"),
        (PyCFunction)Optizelle::Python::EqualityConstrained::Restart
            ::write_restart,
        METH_VARARGS,
        const_cast <char*> ("Writes a json restart file")},

    { const_cast <char*> ("EqualityConstrainedRestartReadRestart"),
        (PyCFunction)Optizelle::Python::EqualityConstrained::Restart
            ::read_restart,
        METH_VARARGS,
        const_cast <char*> ("Reads a json restart file")},

    { const_cast <char*> ("InequalityConstrainedStateCreate"),
        (PyCFunction)Optizelle::Python::InequalityConstrained::State::create,
        METH_VARARGS,
        const_cast <char*> ("Creates an inequality constrained state")},

    { const_cast <char*> ("InequalityConstrainedStateReadJson"),
        (PyCFunction)Optizelle::Python::InequalityConstrained::State::readJson,
        METH_VARARGS,
        const_cast <char*> (
            "Reads inequality constrained state parameters from file")},

    { const_cast <char*> ("InequalityConstrainedAlgorithmsGetMin"),
        (PyCFunction)Optizelle::Python::InequalityConstrained::Algorithms
            ::getMin,
        METH_VARARGS,
        const_cast <char*> (
            "Solves an inequality constrained optimization problem")},
    
    { const_cast <char*> ("InequalityConstrainedRestartRelease"),
        (PyCFunction)Optizelle::Python::InequalityConstrained::Restart::release,
        METH_VARARGS,
        const_cast <char*> (
            "Release the state in an inequality constrained optimization "
            "problem")},

    { const_cast <char*> ("InequalityConstrainedRestartCapture"),
        (PyCFunction)Optizelle::Python::InequalityConstrained::Restart::capture,
        METH_VARARGS,
        const_cast <char*> (
            "Capture the state in an inequality constrained optimization "
            "problem")},
    
    { const_cast <char*> ("InequalityConstrainedRestartWriteRestart"),
        (PyCFunction)Optizelle::Python::InequalityConstrained::Restart
            ::write_restart,
        METH_VARARGS,
        const_cast <char*> ("Writes a json restart file")},

    { const_cast <char*> ("InequalityConstrainedRestartReadRestart"),
        (PyCFunction)Optizelle::Python::InequalityConstrained::Restart
            ::read_restart,
        METH_VARARGS,
        const_cast <char*> ("Reads a json restart file")},

    { const_cast <char*> ("ConstrainedStateCreate"),
        (PyCFunction)Optizelle::Python::Constrained::State::create,
        METH_VARARGS,
        const_cast <char*> ("Creates a constrained state")},

    { const_cast <char*> ("ConstrainedStateReadJson"),
        (PyCFunction)Optizelle::Python::Constrained::State::readJson,
        METH_VARARGS,
        const_cast <char*> (
            "Reads constrained state parameters from file")},

    { const_cast <char*> ("ConstrainedAlgorithmsGetMin"),
        (PyCFunction)Optizelle::Python::Constrained::Algorithms::getMin,
        METH_VARARGS,
        const_cast <char*> ("Solves a constrained optimization problem")},
    
    { const_cast <char*> ("ConstrainedRestartRelease"),
        (PyCFunction)Optizelle::Python::Constrained::Restart::release,
        METH_VARARGS,
        const_cast <char*> (
            "Release the state in a constrained optimization problem")},

    { const_cast <char*> ("ConstrainedRestartCapture"),
        (PyCFunction)Optizelle::Python::Constrained::Restart::capture,
        METH_VARARGS,
        const_cast <char*> (
            "Capture the state in a constrained optimization problem")},
    
    { const_cast <char*> ("ConstrainedRestartWriteRestart"),
        (PyCFunction)Optizelle::Python::Constrained::Restart::write_restart,
        METH_VARARGS,
        const_cast <char*> ("Writes a json restart file")},

    { const_cast <char*> ("ConstrainedRestartReadRestart"),
        (PyCFunction)Optizelle::Python::Constrained::Restart::read_restart,
        METH_VARARGS,
        const_cast <char*> ("Reads a json restart file")},

    {nullptr}  // Sentinel
};

PyMODINIT_FUNC initUtility() {
    PyObject * m;

    // Initilize the module
    m = Py_InitModule3(
        "Utility",
        methods,
        "Internal utility functions for Optizelle");

    if (m == nullptr)
      return;
}
