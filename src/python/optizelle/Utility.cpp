/*
Copyright 2013 OptimoJoe.

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
        PyObject* toPython::operator () (t const & opt_stop) const {
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
        t fromPython::operator () (PyObject const * const member) const{
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(const_cast <PyObject*> (member));

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
        }
    }
    
    namespace KrylovStop { 
        // Converts t to a Python enumerated type
        PyObject* toPython::operator () (t const & krylov_stop) const {
            // Do the conversion
            switch(krylov_stop){
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
            case Instability:
                return Python::enumToPyObject("KrylovStop","Instability");
            case InvalidTrustRegionCenter:
                return Python::enumToPyObject(
                    "KrylovStop","InvalidTrustRegionCenter");
            default:
                throw;
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython::operator () (PyObject const * const member) const{
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(const_cast <PyObject*> (member));

            if(m==Python::enumToNatural("KrylovStop","NegativeCurvature"))
                return NegativeCurvature;
            else if(m==Python::enumToNatural("KrylovStop","RelativeErrorSmall"))
                return RelativeErrorSmall;
            else if(m==Python::enumToNatural("KrylovStop","MaxItersExceeded"))
                return MaxItersExceeded;
            else if(m==Python::enumToNatural(
                "KrylovStop","TrustRegionViolated")
            )
                return TrustRegionViolated;
            else if(m==Python::enumToNatural("KrylovStop","Instability"))
                return Instability;
            else if(m==Python::enumToNatural(
                "KrylovStop","InvalidTrustRegionCenter")
            )
                return InvalidTrustRegionCenter;
        }
    }
    
    namespace KrylovSolverTruncated { 
        // Converts t to a Python enumerated type
        PyObject* toPython::operator () (t const & truncated_krylov) const {
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
        t fromPython::operator () (PyObject const * const member) const{
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(const_cast <PyObject*> (member));

            if(m==Python::enumToNatural(
                "KrylovSolverTruncated","ConjugateDirection")
            )
                return ConjugateDirection;
            else if(m==Python::enumToNatural("KrylovSolverTruncated","MINRES"))
                return MINRES;
        }
    }

    namespace AlgorithmClass { 
        // Converts t to a Python enumerated type
        PyObject* toPython::operator () (t const & algorithm_class) const {
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
        t fromPython::operator () (PyObject const * const member) const{
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(const_cast <PyObject*> (member));

            if(m==Python::enumToNatural("AlgorithmClass","TrustRegion"))
                return TrustRegion;
            else if(m==Python::enumToNatural("AlgorithmClass","LineSearch"))
                return LineSearch;
            else if(m==Python::enumToNatural("AlgorithmClass","UserDefined"))
                return UserDefined;
        }
    }

    namespace Operators { 
        // Converts t to a Python enumerated type
        PyObject* toPython::operator () (t const & op) const {
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
        t fromPython::operator () (PyObject const * const member) const{
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(const_cast <PyObject*> (member));

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
        }
    }

    namespace LineSearchDirection {
        // Converts t to a Python enumerated type
        PyObject* toPython::operator () (t const & dir) const {
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
        t fromPython::operator () (PyObject const * const member) const{
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(const_cast <PyObject*> (member));

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
        }
    }

    namespace LineSearchKind { 
        // Converts t to a Python enumerated type
        PyObject* toPython::operator () (t const & kind) const {
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
        t fromPython::operator () (PyObject const * const member) const{
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(const_cast <PyObject*> (member));

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
        }
    }

    namespace OptimizationLocation { 
        // Converts t to a Python enumerated type
        PyObject* toPython::operator () (t const & loc) const {
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
        t fromPython::operator () (PyObject const * const member) const{
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(const_cast <PyObject*> (member));

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
        }
    }

    namespace InteriorPointMethod { 
        // Converts t to a Python enumerated type
        PyObject* toPython::operator () (t const & ipm) const {
            // Do the conversion
            switch(ipm){
            case PrimalDual:
                return Python::enumToPyObject("InteriorPointMethod",
                    "PrimalDual");
            case PrimalDualLinked:
                return Python::enumToPyObject("InteriorPointMethod",
                    "PrimalDualLinked");
            case LogBarrier:
                return Python::enumToPyObject("InteriorPointMethod",
                    "LogBarrier");
            default:
                throw;
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython::operator () (PyObject const * const member) const{
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(const_cast <PyObject*> (member));

            if(m==Python::enumToNatural("InteriorPointMethod","PrimalDual"))
                return PrimalDual;
            else if(m==Python::enumToNatural("InteriorPointMethod",
                "PrimalDualLinked")
            )
                return PrimalDualLinked;
            else if(m==Python::enumToNatural("InteriorPointMethod",
                "LogBarrier")
            )
                return LogBarrier;
        }
    }

    namespace CentralityStrategy { 
        // Converts t to a Python enumerated type
        PyObject* toPython::operator () (t const & cstrat) const {
            // Do the conversion
            switch(cstrat){
            case Constant:
                return Python::enumToPyObject("CentralityStrategy","Constant");
            case StairStep:
                return Python::enumToPyObject("CentralityStrategy","StairStep");
            case PredictorCorrector:
                return Python::enumToPyObject("CentralityStrategy",
                    "PredictorCorrector");
            default:
                throw;
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython::operator () (PyObject const * const member) const{
            // Convert the member to a Natural 
            Natural m=PyInt_AsSsize_t(const_cast <PyObject*> (member));

            if(m==Python::enumToNatural("CentralityStrategy","Constant"))
                return Constant;
            else if(m==Python::enumToNatural("CentralityStrategy","StairStep"))
                return StairStep;
            else if(m==Python::enumToNatural("CentralityStrategy",
                "PredictorCorrector")
            )
                return PredictorCorrector;
        }
    }

    namespace Python {
        // Like PyObject_GetAttrString, but returns obj.name1.name2 
        PyObject* PyObject_GetAttrString2(
            PyObject const * const obj,
            std::string const & name1,
            std::string const & name2
        ) {
            // Grab name1 from obj and store in value 1 
            PyObjectPtr value1(PyObject_GetAttrString(
                const_cast <PyObject*> (obj),name1.c_str()));

            // Check if we grabbed name1 from obj.  If we did, evaluate
            // name 2.  Otherwise, return null
            if(value1.get())
                return PyObject_GetAttrString(value1.get(),name2.c_str());
            else
                return nullptr;
        }

        // Calls a Python function with one argument 
        PyObject* PyObject_CallObject1(
            PyObject const * const fn,
            PyObject const * const arg1
        ) {
            PyObjectPtr args(PyTuple_New(1)); 
            MyPyTuple_SetItem(args.get(),0,const_cast <PyObject*> (arg1));
            return PyObject_CallObject(const_cast <PyObject*> (fn),args.get()); 
        }
        
        // Calls a Python function with two arguments
        PyObject* PyObject_CallObject2(
            PyObject const * const fn,
            PyObject const * const arg1,
            PyObject const * const arg2
        ) {
            PyObjectPtr args(PyTuple_New(2)); 
            MyPyTuple_SetItem(args.get(),0,const_cast <PyObject*> (arg1)); 
            MyPyTuple_SetItem(args.get(),1,const_cast <PyObject*> (arg2)); 
            return PyObject_CallObject(const_cast <PyObject*> (fn),args.get()); 
        }
        
        // Calls a Python function with three arguments
        PyObject* PyObject_CallObject3(
            PyObject const * const fn,
            PyObject const * const arg1,
            PyObject const * const arg2,
            PyObject const * const arg3
        ) {
            PyObjectPtr args(PyTuple_New(3)); 
            MyPyTuple_SetItem(args.get(),0,const_cast <PyObject*> (arg1)); 
            MyPyTuple_SetItem(args.get(),1,const_cast <PyObject*> (arg2)); 
            MyPyTuple_SetItem(args.get(),2,const_cast <PyObject*> (arg3)); 
            return PyObject_CallObject(const_cast <PyObject*> (fn),args.get()); 
        }
        
        // Calls a Python function with four arguments
        PyObject* PyObject_CallObject4(
            PyObject const * const fn,
            PyObject const * const arg1,
            PyObject const * const arg2,
            PyObject const * const arg3,
            PyObject const * const arg4
        ) {
            PyObjectPtr args(PyTuple_New(4)); 
            MyPyTuple_SetItem(args.get(),0,const_cast <PyObject*> (arg1)); 
            MyPyTuple_SetItem(args.get(),1,const_cast <PyObject*> (arg2)); 
            MyPyTuple_SetItem(args.get(),2,const_cast <PyObject*> (arg3)); 
            MyPyTuple_SetItem(args.get(),3,const_cast <PyObject*> (arg4)); 
            return PyObject_CallObject(const_cast <PyObject*> (fn),args.get()); 
        }

        // Deep copy of a Python object and return the result
        PyObject* deepcopy(PyObject const * const in) {
            // Grab the deepcopy function from the copy module 
            PyObjectPtr module(PyImport_ImportModule("copy")); 
            PyObjectPtr deepcopy(PyObject_GetAttrString(module.get(),
                "deepcopy")); 

            // Call deepcopy on vec and return the result
            PyObjectPtr args(PyTuple_New(1)); 
            MyPyTuple_SetItem(args.get(),0,const_cast <PyObject*> (in)); 
            return PyObject_CallObject(deepcopy.get(),args.get()); 
        }

        // Used to catch Python exceptions
        Exception::Exception() {}

        // On construction, initialize the pointer and figure out if
        // we're capturing the pointer or attaching to it
        PyObjectPtr::PyObjectPtr(
            PyObject* ptr_,
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

        // For a reset, we decrement the pointer and then assign a new
        // value.
        void PyObjectPtr::reset(PyObject* ptr_) {
            Py_XDECREF(ptr);
            ptr=ptr_;
        }

        // For an attach, we decrement the pointer, assign a new value,
        // and then increment the reference count.
        void PyObjectPtr::attach(PyObject* ptr_) {
            Py_XDECREF(ptr);
            ptr=ptr_;
            Py_XINCREF(ptr);
        }

        // On a get, we simply return the pointer.
        PyObject const * const PyObjectPtr::get() const {
            return ptr;
        }
        PyObject* PyObjectPtr::get() {
            return ptr;
        }
    
        // On a release, we return the underlying pointer and then clear
        // the vector.  This will prevent a decrement later.
        PyObject* PyObjectPtr::release() {
            PyObject* ptr_=ptr;
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
            PyObject* ptr_,
            PyObjectPtrMode::t const mode
        ) : PyObjectPtr(ptr_,mode) { }
            
        // Move constructor
        Messaging::Messaging(Messaging&& ptr_) noexcept
            : PyObjectPtr(ptr_.release()) {}
            
        // Prints a message
        void Messaging::print(std::string const & msg_) const {
            // Call the print function on msg
            PyObjectPtr print_(PyObject_GetAttrString(ptr,"print"));
            PyObjectPtr msg(PyString_FromString(msg_.c_str()));
            PyObjectPtr ret(PyObject_CallObject1(print_.get(),msg.get()));

            // Check errors
            if(ret.get()==nullptr)
                error("Evaluation of the print function in the Messaging "
                    "object failed.");
        }

        // Prints an error
        void Messaging::error(std::string const & msg_) const {
            // Call the error function on msg
            PyObjectPtr error_(PyObject_GetAttrString(ptr,"error"));
            PyObjectPtr msg(PyString_FromString(msg_.c_str()));
            PyObjectPtr ret(PyObject_CallObject1(error_.get(),msg.get()));

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

        // Create an empty vector.  Don't use these objects until
        // initialized.
        Vector::Vector() :
            PyObjectPtr(nullptr),
            msg(nullptr),
            vs(nullptr)
        { }

        // Create a vector with the appropriate messaging and vector space 
        Vector::Vector(PyObject* msg_,PyObject* vs_,PyObject* vec,
            PyObjectPtrMode::t mode
        ) : 
            PyObjectPtr(vec,mode),
            msg(msg_,PyObjectPtrMode::Attach),
            vs(vs_,PyObjectPtrMode::Attach)
        { }
            
        // Create a move constructor so we can interact with stl objects
        Vector::Vector(Vector&& vec) :
            PyObjectPtr(std::move(vec)),
            msg(std::move(vec.msg)),
            vs(std::move(vec.vs))
        { }

        // Memory allocation and size setting 
        void Vector::init(Vector const & x) { 
            // Attach to the messaging object and vector space
            msg.attach(const_cast <Vector&> (x).msg.get());
            vs.attach(const_cast <Vector&> (x).vs.get());

            // Call the init function on x and store internally
            PyObjectPtr init_(PyObject_GetAttrString2(
                vs.get(),"init","__func__"));
            reset(PyObject_CallObject1(init_.get(),x.get()));

            // Check errors
            if(get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function init failed.");
        } 
        
        // y <- x (Shallow.  No memory allocation.) 
        void Vector::copy(Vector const & x) { 
            // Call the copy function on x and the internal 
            PyObjectPtr copy_(PyObject_GetAttrString2(
                vs.get(),"copy","__func__"));
            PyObjectPtr ret(PyObject_CallObject2(copy_.get(),x.get(),get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function copy failed.");
        } 

        // x <- alpha * x
        void Vector::scal(double const & alpha_) { 
            // Call the scal function on alpha and the internal storage 
            PyObjectPtr scal_(PyObject_GetAttrString2(
                vs.get(),"scal","__func__"));
            PyObjectPtr alpha(PyFloat_FromDouble(alpha_));
            PyObjectPtr ret(
                PyObject_CallObject2(scal_.get(),alpha.get(),get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function scal failed.");
        } 

        // x <- 0 
        void Vector::zero() { 
            // Call the zero function on this vector.
            PyObjectPtr zero_(PyObject_GetAttrString2(
                vs.get(),"zero","__func__"));
            PyObjectPtr ret(PyObject_CallObject1(zero_.get(),get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function zero failed.");
        } 

        // y <- alpha * x + y 
        void Vector::axpy(double const & alpha_,Vector const & x) { 
            // Call the axpy function on alpha, x, and the internal storage.
            PyObjectPtr axpy_(PyObject_GetAttrString2(
                vs.get(),"axpy","__func__"));
            PyObjectPtr alpha(PyFloat_FromDouble(alpha_));
            PyObjectPtr ret(
                PyObject_CallObject3(axpy_.get(),alpha.get(),x.get(),get()));
           
            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function axpy failed.");
        } 

        // innr <- <x,y> 
        double Vector::innr(Vector const & x) const { 
            // Call the innr function on x and the internal.  Store in z. 
            PyObjectPtr innr_(PyObject_GetAttrString2(
                vs.get(),"innr","__func__"));
            PyObjectPtr z(PyObject_CallObject2(innr_.get(),x.get(),get()));

            // Check errors
            if(z.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function innr failed.");

            // Return the result 
            return PyFloat_AsDouble(z.get()); 
        } 

        // Jordan product, z <- x o y
        void Vector::prod(Vector const & x,Vector const & y) { 
            // Call the prod function on x, y, and the internal 
            PyObjectPtr prod_(PyObject_GetAttrString2(
                vs.get(),"prod","__func__"));
            PyObjectPtr ret(
                PyObject_CallObject3(prod_.get(),x.get(),y.get(),get()));

            // Check errors
            if(get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function prod failed.");
        } 

        // Identity element, x <- e such that x o e = x 
        void Vector::id() { 
            // Call the id function on the internal.
            PyObjectPtr id_(PyObject_GetAttrString2(
                vs.get(),"id","__func__"));
            PyObjectPtr ret(PyObject_CallObject1(id_.get(),get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function id failed.");
        } 

        // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y 
        void Vector::linv(const Vector& x, const Vector& y) { 
            // Call the linv function on x, y, and the internal
            PyObjectPtr linv_(PyObject_GetAttrString2(
                vs.get(),"linv","__func__"));
            PyObjectPtr ret(
                PyObject_CallObject3(linv_.get(),x.get(),y.get(),get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function linv failed.");
        } 

        // Barrier function, barr <- barr(x) where x o grad barr(x) = e 
        double Vector::barr() const { 
            // Call the barr function on the internal.  Store in z.
            PyObjectPtr barr_(PyObject_GetAttrString2(
                vs.get(),"barr","__func__"));
            PyObjectPtr z(PyObject_CallObject1(barr_.get(),get()));

            // Check errors
            if(z.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function barr failed.");

            // Return the result 
            return PyFloat_AsDouble(z.get());
        } 

        // Line search, srch <- argmax {alpha in Real >= 0 : alpha x + y >= 0} 
        // where y > 0. 
        double Vector::srch(const Vector& x) const {  
            // Call the srch function on x and the internal.  Store in z.
            PyObjectPtr srch_(PyObject_GetAttrString2(
                vs.get(),"srch","__func__"));
            PyObjectPtr z(PyObject_CallObject2(srch_.get(),x.get(),get()));

            // Check errors
            if(z.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function srch failed.");

            // Return the result 
            return PyFloat_AsDouble(z.get());
        } 

        // Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
        // operator.
        void Vector::symm() { 
            // Call the symm function on the internal.
            PyObjectPtr symm_(PyObject_GetAttrString2(
                vs.get(),"symm","__func__"));
            PyObjectPtr ret(PyObject_CallObject1(symm_.get(),get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error(
                    "Evaluation of the vector space function symm failed.");

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
            PyObject* msg_,
            PyObject* f,
            PyObjectPtrMode::t mode
        ) :
            PyObjectPtr(f,mode),
            msg(msg_,PyObjectPtrMode::Attach)
        { }

        // <- f(x) 
        double ScalarValuedFunction::operator () (
            const ScalarValuedFunction::Vector& x
        ) const { 
            // Call the objective function on x.  Store in z.
            PyObjectPtr eval_(PyObject_GetAttrString(ptr,"eval"));
            PyObjectPtr z(PyObject_CallObject1(eval_.get(),x.get()));

            // Check errors
            if(z.get()==nullptr)
                msg.error("Evaluation of the objective f failed.");

            // Return the result
            return PyFloat_AsDouble(z.get());
        }

        // grad = grad f(x) 
        void ScalarValuedFunction::grad(
            const ScalarValuedFunction::Vector& x,
            ScalarValuedFunction::Vector& grad
        ) const { 
            // Call the gradient function on x and grad. 
            PyObjectPtr grad_(PyObject_GetAttrString(ptr,"grad"));
            PyObjectPtr ret(PyObject_CallObject2(
                grad_.get(),x.get(),grad.get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error("Evaluation of the gradient of f failed.");
        }

        // H_dx = hess f(x) dx 
        void ScalarValuedFunction::hessvec(
            const ScalarValuedFunction::Vector& x,
            const ScalarValuedFunction::Vector& dx,
            ScalarValuedFunction::Vector& H_dx
        ) const {
            // Call the hessvec function on x, dx, and H_dx.
            PyObjectPtr hessvec_(PyObject_GetAttrString(ptr,"hessvec"));
            PyObjectPtr ret(PyObject_CallObject3(
                hessvec_.get(),x.get(),dx.get(),H_dx.get()));

            // Check errors
            if(ret.get()==nullptr)
                msg.error("Evaluation of the Hessian-vector product"
                    " of f failed.");
        }

        // Create a function 
        VectorValuedFunction::VectorValuedFunction(
            std::string const & name_,
            PyObject* msg_,
            PyObject* f,
            PyObjectPtrMode::t mode
        ) :
            PyObjectPtr(f,mode),
            msg(msg_,PyObjectPtrMode::Attach),
            name(name_)
        { }

        // y=f(x)
        void VectorValuedFunction::operator () (
            const VectorValuedFunction::X_Vector& x,
            VectorValuedFunction::Y_Vector& y
        ) const {
            // Call the evaluate function on x and y.
            PyObjectPtr eval_(PyObject_GetAttrString(ptr,"eval"));
            PyObjectPtr ret(PyObject_CallObject2(eval_.get(),x.get(),y.get()));

            // Check errors
            if(ret.get()==nullptr) {
                std::stringstream ss;
                ss << "Evaluation of the constraint " << name << "failed.";
                msg.error(ss.str());
            }
        }

        // y=f'(x)dx 
        void VectorValuedFunction::p(
            const VectorValuedFunction::X_Vector& x,
            const VectorValuedFunction::X_Vector& dx,
            VectorValuedFunction::Y_Vector& y
        ) const {
            // Call the prime function on x, dx, and y
            PyObjectPtr p_(PyObject_GetAttrString(ptr,"p"));
            PyObjectPtr ret(PyObject_CallObject3(
                p_.get(),x.get(),dx.get(),y.get()));
           
            // Check errors
            if(ret.get()==nullptr) {
                std::stringstream ss;
                ss << "Evaluation of the derivative of the constraint "
                    << name << "failed.";
                msg.error(ss.str());
            }
        }

        // z=f'(x)*dy
        void VectorValuedFunction::ps(
            const VectorValuedFunction::X_Vector& x,
            const VectorValuedFunction::Y_Vector& dy,
            VectorValuedFunction::X_Vector& z
        ) const {
            // Call the prime-adjoint function on x, dy, and z
            PyObjectPtr ps_(PyObject_GetAttrString(ptr,"ps"));
            PyObjectPtr ret(PyObject_CallObject3(
                ps_.get(),x.get(),dy.get(),z.get()));

            // Check errors
            if(ret.get()==nullptr) {
                std::stringstream ss;
                ss << "Evaluation of the derivative-adjoint of the constraint "
                    << name << "failed.";
                msg.error(ss.str());
            }
        }
             
        // z=(f''(x)dx)*dy
        void VectorValuedFunction::pps(
            const VectorValuedFunction::X_Vector& x,
            const VectorValuedFunction::X_Vector& dx,
            const VectorValuedFunction::Y_Vector& dy,
            X_Vector& z
        ) const { 
            // Call the prime-adjoint function on x, dx, dy, and z
            PyObjectPtr pps_(PyObject_GetAttrString(ptr,"pps"));
            PyObjectPtr ret(PyObject_CallObject4(
                pps_.get(),x.get(),dx.get(),dy.get(),z.get()));

            // Check errors
            if(ret.get()==nullptr) {
                std::stringstream ss;
                ss << "Evaluation of the second derivative-adjoint of the "
                    "constraint " << name << "failed.";
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

        // Converts an Optizelle enumerated type to a PyObject* 
        PyObject* enumToPyObject(
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
            // Grab the PyObject* for the type and member requested
            PyObjectPtr obj(enumToPyObject(type,member));

            // Convert and return the member
            return PyInt_AsSsize_t(obj.get());
        }
                    
        // Sets a floating point in a Python class 
        void setFloat(
            std::string const & name,
            double const & value,
            PyObject* obj 
        ) {
            PyObjectPtr item(PyFloat_FromDouble(value));
            PyObject_SetAttrString(obj,name.c_str(),item.get());
        }
        
        // Sets a floating point in a C++ class 
        void setFloat(
            std::string const & name,
            PyObject const * const obj,
            double & value
        ) {
            PyObjectPtr item(PyObject_GetAttrString(
                const_cast <PyObject*> (obj),name.c_str()));
            value=PyFloat_AsDouble(item.get());
        }
        
        // Sets an integer in a Python class 
        void setNatural(
            std::string const & name,
            Natural const & value,
            PyObject* obj 
        ) {
            PyObjectPtr item(PyInt_FromSsize_t(value));
            PyObject_SetAttrString(obj,name.c_str(),item.get());
        }
        
        // Sets an integer in a C++ class 
        void setNatural(
            std::string const & name,
            PyObject const * const obj,
            Natural & value
        ) {
            PyObjectPtr item(PyObject_GetAttrString(
                const_cast <PyObject*>(obj),name.c_str()));
            value=PyInt_AsSsize_t(item.get());
        }
        
        // Sets a vector in a Python class 
        void setVector(
            std::string const & name,
            Vector const & value,
            PyObject* obj 
        ) {
            PyObjectPtr item(deepcopy(value.get()));
            PyObject_SetAttrString(obj,name.c_str(),item.get());
        }
        
        // Sets a vector in a C++ class 
        void setVector(
            std::string const & name,
            PyObject const * const obj,
            Vector & value
        ) {
            PyObjectPtr item(PyObject_GetAttrString(
                const_cast <PyObject*> (obj),name.c_str()));
            value.reset(deepcopy(item.get()));
        }
        
        // Sets a list of vectors in a Python class 
        void setVectors(
            std::string const & name,
            std::list <Vector> const & values,
            PyObject* obj 
        ) {
            // Create a new Python list that we insert elements into
            PyObjectPtr items(PyList_New(0));

            // Loop over all of the items inside values and then insert 
            // them into items 
            for(std::list <Vector>::const_iterator value=values.cbegin();
                value!=values.cend();
                value++
            )
                PyList_Append(items.get(),deepcopy(value->get()));
            
            // Insert the items into obj
            PyObject_SetAttrString(obj,name.c_str(),items.get());
        }
        
        // Sets a list of vectors in a C++ class 
        void setVectors(
            std::string const & name,
            PyObject const * const obj,
            std::list <Vector> & values
        ) {
            // Grab the list of items
            PyObjectPtr items(PyObject_GetAttrString(
                const_cast <PyObject*> (obj),name.c_str()));

            // Loop over all the elements in items and insert them one
            // at a time into values
            values.clear();
            for(Natural i=0;i<PyList_Size(items.get());i++) {
                // Grab the current item
                PyObjectPtr item(PyList_GetItem(items.get(),i));

                // Create a new vector in values 
                values.emplace_back(Vector());

                // Insert this item into values
                values.back().reset(deepcopy(item.get()));
            }
        }
        
        // Sets a scalar-valued function in a Python class 
        void setScalarValuedFunction(
            std::string const & name,
            PyObject * const msg,
            PyObject * const obj,
            std::unique_ptr <PyScalarValuedFunction> & value
        ) {
            value.reset(new ScalarValuedFunction(msg,
                PyObject_GetAttrString(obj,name.c_str())));
        }
        
        // Sets a vector-valued function in a Python class 
        void setVectorValuedFunction(
            std::string const & name,
            PyObject * const msg,
            PyObject * const obj,
            std::unique_ptr <PyVectorValuedFunction> & value
        ) {
            value.reset(new VectorValuedFunction(name,msg,
                PyObject_GetAttrString(obj,name.c_str())));
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
                    PyObject* pystate
                ){
                    // Set each of the required items in the Python state
                    setFloat("eps_grad",state.eps_grad,pystate);
                    setFloat("eps_dx",state.eps_dx,pystate);
                    setNatural("stored_history",state.stored_history,pystate);
                    setNatural("history_reset",state.history_reset,pystate);
                    setNatural("iter",state.iter,pystate);
                    setNatural("iter_max",state.iter_max,pystate);
                    setEnum <StoppingCondition::t> (
                        "opt_stop",
                        StoppingCondition::toPython(),
                        state.opt_stop,
                        pystate);
                    setNatural("krylov_iter",state.krylov_iter,pystate);
                    setNatural("krylov_iter_max",state.krylov_iter_max,pystate);
                    setNatural("krylov_iter_total",
                        state.krylov_iter_total,pystate);
                    setNatural("krylov_orthog_max",
                        state.krylov_orthog_max,pystate);
                    setEnum <KrylovStop::t> (
                        "krylov_stop",
                        KrylovStop::toPython(),
                        state.krylov_stop,
                        pystate);
                    setFloat("krylov_rel_err",state.krylov_rel_err,pystate);
                    setFloat("eps_krylov",state.eps_krylov,pystate);
                    setEnum <KrylovSolverTruncated::t> (
                        "krylov_solver",
                        KrylovSolverTruncated::toPython(),
                        state.krylov_solver,
                        pystate);
                    setEnum <AlgorithmClass::t> (
                        "algorithm_class",
                        AlgorithmClass::toPython(),
                        state.algorithm_class,
                        pystate);
                    setEnum <Operators::t> (
                        "PH_type",
                        Operators::toPython(),
                        state.PH_type,
                        pystate);
                    setEnum <Operators::t> (
                        "H_type",
                        Operators::toPython(),
                        state.H_type,
                        pystate);
                    setFloat("norm_gradtyp",state.norm_gradtyp,pystate);
                    setFloat("norm_dxtyp",state.norm_dxtyp,pystate);
                    setVector("x",state.x.front(),pystate);
                    setVector("grad",state.grad.front(),pystate);
                    setVector("dx",state.dx.front(),pystate);
                    setVector("x_old",state.x_old.front(),pystate);
                    setVector("grad_old",state.grad_old.front(),pystate);
                    setVector("dx_old",state.dx_old.front(),pystate);
                    setVectors("oldY",state.oldY,pystate);
                    setVectors("oldS",state.oldS,pystate);
                    setFloat("f_x",state.f_x,pystate);
                    setFloat("f_xpdx",state.f_xpdx,pystate);
                    setNatural("msg_level",state.msg_level,pystate);
                    setFloat("delta",state.delta,pystate);
                    setFloat("eta1",state.eta1,pystate);
                    setFloat("eta2",state.eta2,pystate);
                    setFloat("ared",state.ared,pystate);
                    setFloat("pred",state.pred,pystate);
                    setNatural("rejected_trustregion",
                        state.rejected_trustregion,pystate);
                    setFloat("alpha0",state.alpha0,pystate);
                    setFloat("alpha",state.alpha,pystate);
                    setFloat("c1",state.c1,pystate);
                    setNatural("linesearch_iter",state.linesearch_iter,pystate);
                    setNatural("linesearch_iter_max",
                        state.linesearch_iter_max,pystate);
                    setNatural("linesearch_iter_total",
                        state.linesearch_iter_total,pystate);
                    setFloat("eps_ls",state.eps_ls,pystate);
                    setEnum <LineSearchDirection::t> (
                        "dir",
                        LineSearchDirection::toPython(),
                        state.dir,
                        pystate);
                    setEnum <LineSearchKind::t> (
                        "kind",
                        LineSearchKind::toPython(),
                        state.kind,
                        pystate);
                }
                void toPython(
                    typename PyUnconstrained::State::t const & state,
                    PyObject* pystate
                ){
                    Unconstrained::State::toPython_(state,pystate);
                }
                
                // Convert a Python state to C++ 
                void fromPython_(
                    PyObject const * const pystate,
                    typename PyUnconstrained::State::t & state
                ){
                    // Set each of the required items in the Python state
                    setFloat("eps_grad",pystate,state.eps_grad);
                    setFloat("eps_dx",pystate,state.eps_dx);
                    setNatural("stored_history",pystate,state.stored_history);
                    setNatural("history_reset",pystate,state.history_reset);
                    setNatural("iter",pystate,state.iter);
                    setNatural("iter_max",pystate,state.iter_max);
                    setEnum <StoppingCondition::t> (
                        "opt_stop",
                        StoppingCondition::fromPython(),
                        pystate,
                        state.opt_stop);
                    setNatural("krylov_iter",pystate,state.krylov_iter);
                    setNatural("krylov_iter_max",pystate,state.krylov_iter_max);
                    setNatural("krylov_iter_total",
                        pystate,state.krylov_iter_total);
                    setNatural("krylov_orthog_max",
                        pystate,state.krylov_orthog_max);
                    setEnum <KrylovStop::t> (
                        "krylov_stop",
                        KrylovStop::fromPython(),
                        pystate,
                        state.krylov_stop);
                    setFloat("krylov_rel_err",pystate,state.krylov_rel_err);
                    setFloat("eps_krylov",pystate,state.eps_krylov);
                    setEnum <KrylovSolverTruncated::t> (
                        "krylov_solver",
                        KrylovSolverTruncated::fromPython(),
                        pystate,
                        state.krylov_solver);
                    setEnum <AlgorithmClass::t> (
                        "algorithm_class",
                        AlgorithmClass::fromPython(),
                        pystate,
                        state.algorithm_class);
                    setEnum <Operators::t> (
                        "PH_type",
                        Operators::fromPython(),
                        pystate,
                        state.PH_type);
                    setEnum <Operators::t> (
                        "H_type",
                        Operators::fromPython(),
                        pystate,
                        state.H_type);
                    setFloat("norm_gradtyp",pystate,state.norm_gradtyp);
                    setFloat("norm_dxtyp",pystate,state.norm_dxtyp);
                    setVector("x",pystate,state.x.front());
                    setVector("grad",pystate,state.grad.front());
                    setVector("dx",pystate,state.dx.front());
                    setVector("x_old",pystate,state.x_old.front());
                    setVector("grad_old",pystate,state.grad_old.front());
                    setVector("dx_old",pystate,state.dx_old.front());
                    setVectors("oldY",pystate,state.oldY);
                    setVectors("oldS",pystate,state.oldS);
                    setFloat("f_x",pystate,state.f_x);
                    setFloat("f_xpdx",pystate,state.f_xpdx);
                    setNatural("msg_level",pystate,state.msg_level);
                    setFloat("delta",pystate,state.delta);
                    setFloat("eta1",pystate,state.eta1);
                    setFloat("eta2",pystate,state.eta2);
                    setFloat("ared",pystate,state.ared);
                    setFloat("pred",pystate,state.pred);
                    setNatural("rejected_trustregion",
                        pystate,state.rejected_trustregion);
                    setFloat("alpha0",pystate,state.alpha0);
                    setFloat("alpha",pystate,state.alpha);
                    setFloat("c1",pystate,state.c1);
                    setNatural("linesearch_iter",pystate,state.linesearch_iter);
                    setNatural("linesearch_iter_max",
                        pystate,state.linesearch_iter_max);
                    setNatural("linesearch_iter_total",pystate,
                        state.linesearch_iter_total);
                    setFloat("eps_ls",pystate,state.eps_ls);
                    setEnum <LineSearchDirection::t> (
                        "dir",
                        LineSearchDirection::fromPython(),
                        pystate,
                        state.dir);
                    setEnum <LineSearchKind::t> (
                        "kind",
                        LineSearchKind::fromPython(),
                        pystate,
                        state.kind);
                }
                void fromPython(
                    PyObject const * const pystate,
                    typename PyUnconstrained::State::t & state
                ){
                    Unconstrained::State::fromPython_(pystate,state);
                }

                // Creates a state and inserts the elements into pystate 
                PyObject* create(
                    PyObject* self,
                    PyObject* args
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
                PyObject* readJson(
                    PyObject* self,
                    PyObject* args
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
                PyObject* getMin(
                    PyObject* self,
                    PyObject* args
                ) {
                    // Calling convention should be
                    // (X,msg,smanip,fns,state)
                    PyObject *X,*msg_,*smanip_,*pyfns_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOOO",
                        &X,&msg_,&smanip_,&pyfns_,&pystate_)
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
                            msg,smanip,fns,state);
                        
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
                    PyObject* pystate
                ){
                    setVector("y",state.y.front(),pystate);
                    setVector("dy",state.dy.front(),pystate);
                    setFloat("zeta",state.zeta,pystate);
                    setFloat("eta0",state.eta0,pystate);
                    setFloat("rho",state.rho,pystate);
                    setFloat("rho_old",state.rho_old,pystate);
                    setFloat("rho_bar",state.rho_bar,pystate);
                    setFloat("eps_constr",state.eps_constr,pystate);
                    setFloat("xi_qn",state.xi_qn,pystate);
                    setFloat("xi_pg",state.xi_pg,pystate);
                    setFloat("xi_proj",state.xi_proj,pystate);
                    setFloat("xi_tang",state.xi_tang,pystate);
                    setFloat("xi_lmh",state.xi_lmh,pystate);
                    setFloat("xi_lmg",state.xi_lmg,pystate);
                    setFloat("xi_4",state.xi_4,pystate);
                    setFloat("rpred",state.rpred,pystate);
                    setEnum <Operators::t> (
                        "PSchur_left_type",
                        Operators::toPython(),
                        state.PSchur_left_type,
                        pystate);
                    setEnum <Operators::t> (
                        "PSchur_right_type",
                        Operators::toPython(),
                        state.PSchur_right_type,
                        pystate);
                    setNatural("augsys_iter_max",state.augsys_iter_max,pystate);
                    setNatural("augsys_rst_freq",state.augsys_rst_freq,pystate);
                    setVector("g_x",state.g_x.front(),pystate);
                    setFloat("norm_gxtyp",state.norm_gxtyp,pystate);
                    setVector("gpxdxn_p_gx",state.gpxdxn_p_gx.front(),pystate);
                    setVector("gpxdxt",state.gpxdxt.front(),pystate);
                    setFloat("norm_gpxdxnpgx",state.norm_gpxdxnpgx,pystate);
                    setVector("dx_n",state.dx_n.front(),pystate);
                    setVector("dx_ncp",state.dx_ncp.front(),pystate);
                    setVector("dx_t",state.dx_t.front(),pystate);
                    setVector("dx_t_uncorrected",
                        state.dx_t_uncorrected.front(),pystate);
                    setVector("dx_tcp_uncorrected",
                        state.dx_tcp_uncorrected.front(),pystate);
                    setVector("H_dxn",state.H_dxn.front(),pystate);
                    setVector("W_gradpHdxn",state.W_gradpHdxn.front(),pystate);
                    setVector("H_dxtuncorrected",
                        state.H_dxtuncorrected.front(),pystate);
                }
                void toPython(
                    typename PyEqualityConstrained::State::t const & state,
                    PyObject* pystate
                ){
                    Unconstrained::State::toPython_(state,pystate);
                    EqualityConstrained::State::toPython_(state,pystate);
                }
                
                // Convert a Python state to C++ 
                void fromPython_(
                    PyObject const * const pystate,
                    typename PyEqualityConstrained::State::t & state
                ){
                    setVector("y",pystate,state.y.front());
                    setVector("dy",pystate,state.dy.front());
                    setFloat("zeta",pystate,state.zeta);
                    setFloat("eta0",pystate,state.eta0);
                    setFloat("rho",pystate,state.rho);
                    setFloat("rho_old",pystate,state.rho_old);
                    setFloat("rho_bar",pystate,state.rho_bar);
                    setFloat("eps_constr",pystate,state.eps_constr);
                    setFloat("xi_qn",pystate,state.xi_qn);
                    setFloat("xi_pg",pystate,state.xi_pg);
                    setFloat("xi_proj",pystate,state.xi_proj);
                    setFloat("xi_tang",pystate,state.xi_tang);
                    setFloat("xi_lmh",pystate,state.xi_lmh);
                    setFloat("xi_lmg",pystate,state.xi_lmg);
                    setFloat("xi_4",pystate,state.xi_4);
                    setFloat("rpred",pystate,state.rpred);
                    setEnum <Operators::t> (
                        "PSchur_left_type",
                        Operators::fromPython(),
                        pystate,
                        state.PSchur_left_type);
                    setEnum <Operators::t> (
                        "PSchur_right_type",
                        Operators::fromPython(),
                        pystate,
                        state.PSchur_right_type);
                    setNatural("augsys_iter_max",pystate,state.augsys_iter_max);
                    setNatural("augsys_rst_freq",pystate,state.augsys_rst_freq);
                    setVector("g_x",pystate,state.g_x.front());
                    setFloat("norm_gxtyp",pystate,state.norm_gxtyp);
                    setVector("gpxdxn_p_gx",pystate,state.gpxdxn_p_gx.front());
                    setVector("gpxdxt",pystate,state.gpxdxt.front());
                    setFloat("norm_gpxdxnpgx",pystate,state.norm_gpxdxnpgx);
                    setVector("dx_n",pystate,state.dx_n.front());
                    setVector("dx_ncp",pystate,state.dx_ncp.front());
                    setVector("dx_t",pystate,state.dx_t.front());
                    setVector("dx_t_uncorrected",
                        pystate,state.dx_t_uncorrected.front());
                    setVector("dx_tcp_uncorrected",
                        pystate,state.dx_tcp_uncorrected.front());
                    setVector("H_dxn",pystate,state.H_dxn.front());
                    setVector("W_gradpHdxn",pystate,state.W_gradpHdxn.front());
                    setVector("H_dxtuncorrected",
                        pystate,state.H_dxtuncorrected.front());
                }
                void fromPython(
                    PyObject const * const pystate,
                    typename PyEqualityConstrained::State::t & state
                ){
                    Unconstrained::State::fromPython_(pystate,state);
                    EqualityConstrained::State::fromPython_(pystate,state);
                }
        
                // Creates a state and inserts the elements into pystate 
                PyObject* create(
                    PyObject* self,
                    PyObject* args
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
                PyObject* readJson(
                    PyObject* self,
                    PyObject* args
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
                PyObject* getMin(
                    PyObject* self,
                    PyObject* args
                ) {
                    // Calling convention should be
                    // (X,Y,msg,smanip,fns,state)
                    PyObject *X,*Y,*msg_,*smanip_,*pyfns_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOOOO",
                        &X,&Y,&msg_,&smanip_,&pyfns_,&pystate_)
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
                            msg,smanip,fns,state);
                        
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
                    PyObject* pystate
                ){
                    setVector("z",state.z.front(),pystate);
                    setVector("dz",state.dz.front(),pystate);
                    setVector("h_x",state.h_x.front(),pystate);
                    setFloat("mu",state.mu,pystate);
                    setFloat("mu_est",state.mu_est,pystate);
                    setFloat("mu_typ",state.mu_typ,pystate);
                    setFloat("eps_mu",state.eps_mu,pystate);
                    setFloat("sigma",state.sigma,pystate);
                    setFloat("gamma",state.gamma,pystate);
                    setEnum <InteriorPointMethod::t> (
                        "ipm",
                        InteriorPointMethod::toPython(),
                        state.ipm,
                        pystate);
                    setEnum <CentralityStrategy::t> (
                        "cstrat",
                        CentralityStrategy::toPython(),
                        state.cstrat,
                        pystate);
                }
                void toPython(
                    typename PyInequalityConstrained::State::t const & state,
                    PyObject* pystate
                ){
                    Unconstrained::State::toPython_(state,pystate);
                    InequalityConstrained::State::toPython_(state,pystate);
                }
                
                // Convert a Python state to C++ 
                void fromPython_(
                    PyObject const * const pystate,
                    typename PyInequalityConstrained::State::t & state
                ){
                    setVector("z",pystate,state.z.front());
                    setVector("dz",pystate,state.dz.front());
                    setVector("h_x",pystate,state.h_x.front());
                    setFloat("mu",pystate,state.mu);
                    setFloat("mu_est",pystate,state.mu_est);
                    setFloat("mu_typ",pystate,state.mu_typ);
                    setFloat("eps_mu",pystate,state.eps_mu);
                    setFloat("sigma",pystate,state.sigma);
                    setFloat("gamma",pystate,state.gamma);
                    setEnum <InteriorPointMethod::t> (
                        "ipm",
                        InteriorPointMethod::fromPython(),
                        pystate,
                        state.ipm);
                    setEnum <CentralityStrategy::t> (
                        "cstrat",
                        CentralityStrategy::fromPython(),
                        pystate,
                        state.cstrat);
                }
                void fromPython(
                    PyObject const * const pystate,
                    typename PyInequalityConstrained::State::t & state
                ){
                    Unconstrained::State::fromPython_(pystate,state);
                    InequalityConstrained::State::fromPython_(pystate,state);
                }

                // Creates a state and inserts the elements into pystate 
                PyObject* create( 
                    PyObject* self,
                    PyObject* args
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
                        Vector x(msg,X,deepcopy(x_));
                        Vector z(msg,Z,deepcopy(z_));

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
                PyObject* readJson(
                    PyObject* self,
                    PyObject* args
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
                PyObject* getMin(
                    PyObject* self,
                    PyObject* args
                ) {
                    // Calling convention should be
                    // (X,Z,msg,smanip,fns,state)
                    PyObject *X,*Z,*msg_,*smanip_,*pyfns_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOOOO",
                        &X,&Z,&msg_,&smanip_,&pyfns_,&pystate_)
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
                            msg,smanip,fns,state);
                        
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
                    PyObject* pystate
                ){
                    Unconstrained::State::toPython_(state,pystate);
                    EqualityConstrained::State::toPython_(state,pystate);
                    InequalityConstrained::State::toPython_(state,pystate);
                }
                
                // Convert a Python state to C++ 
                void fromPython(
                    PyObject const * const pystate,
                    typename PyConstrained::State::t & state
                ){
                    Unconstrained::State::fromPython_(pystate,state);
                    EqualityConstrained::State::fromPython_(pystate,state);
                    InequalityConstrained::State::fromPython_(pystate,state);
                }

                // Creates a state and inserts the elements into pystate 
                PyObject* create(
                    PyObject* self,
                    PyObject* args
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
                        Vector x(msg,X,deepcopy(x_));
                        Vector y(msg,Y,deepcopy(y_));
                        Vector z(msg,Z,deepcopy(z_));

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
                PyObject* readJson(
                    PyObject* self,
                    PyObject* args
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
                PyObject* getMin(
                    PyObject* self,
                    PyObject* args
                ) {
                    // Calling convention should be
                    // (X,Y,Z,msg,smanip,fns,state)
                    PyObject *X,*Y,*Z,*msg_,*smanip_,*pyfns_,*pystate_;
                    if(!PyArg_ParseTuple(args,"OOOOOOO",
                        &X,&Y,&Z,&msg_,&smanip_,&pyfns_,&pystate_)
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
                            msg,smanip,fns,state);
                        
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

    {nullptr}  // Sentinel
};

PyMODINIT_FUNC initUtility() {
    PyObject* m;

    // Initilize the module
    m = Py_InitModule3(
        "Utility",
        methods,
        "Internal utility functions for Optizelle");

    if (m == nullptr)
      return;
}
