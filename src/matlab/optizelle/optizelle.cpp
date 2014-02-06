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

#include "optizelle.h"

namespace Optizelle {
    namespace StoppingCondition { 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & opt_stop) {
            // Do the conversion
            switch(opt_stop){
            case NotConverged:
                return Matlab::enumToMxArray(
                    "StoppingCondition","NotConverged");
            case RelativeGradientSmall:
                return Matlab::enumToMxArray(
                    "StoppingCondition","RelativeGradientSmall");
            case RelativeStepSmall:
                return Matlab::enumToMxArray(
                    "StoppingCondition","RelativeStepSmall");
            case MaxItersExceeded:
                return Matlab::enumToMxArray(
                    "StoppingCondition","MaxItersExceeded");
            case InteriorPointInstability:
                return Matlab::enumToMxArray(
                    "StoppingCondition","InteriorPointInstability");
            case UserDefined:
                return Matlab::enumToMxArray(
                    "StoppingCondition","UserDefined");
            default:
                throw;
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member) {
            // Convert the member to a Natural 
            Natural m(*mxGetPr(member));

            if(m==Matlab::enumToNatural("StoppingCondition","NotConverged"))
                return NotConverged;
            else if(m==Matlab::enumToNatural(
                "StoppingCondition","RelativeGradientSmall")
            )
                return RelativeGradientSmall;
            else if(m==Matlab::enumToNatural(
                "StoppingCondition","RelativeStepSmall")
            )
                return RelativeStepSmall;
            else if(m==Matlab::enumToNatural(
                "StoppingCondition","MaxItersExceeded")
            )
                return MaxItersExceeded;
            else if(m==Matlab::enumToNatural(
                "StoppingCondition","InteriorPointInstability")
            )
                return InteriorPointInstability;
            else if(m==Matlab::enumToNatural("StoppingCondition","UserDefined"))
                return UserDefined;
        }
    }
    
    namespace KrylovStop { 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & krylov_stop) {
            // Do the conversion
            switch(krylov_stop){
            case NegativeCurvature:
                return Matlab::enumToMxArray("KrylovStop","NegativeCurvature");
            case RelativeErrorSmall:
                return Matlab::enumToMxArray(
                    "KrylovStop","RelativeErrorSmall");
            case MaxItersExceeded:
                return Matlab::enumToMxArray("KrylovStop","MaxItersExceeded");
            case TrustRegionViolated:
                return Matlab::enumToMxArray(
                    "KrylovStop","TrustRegionViolated");
            case Instability:
                return Matlab::enumToMxArray("KrylovStop","Instability");
            case InvalidTrustRegionCenter:
                return Matlab::enumToMxArray(
                    "KrylovStop","InvalidTrustRegionCenter");
            default:
                throw;
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member) {
            // Convert the member to a Natural 
            Natural m(*mxGetPr(member));

            if(m==Matlab::enumToNatural("KrylovStop","NegativeCurvature"))
                return NegativeCurvature;
            else if(m==Matlab::enumToNatural("KrylovStop","RelativeErrorSmall"))
                return RelativeErrorSmall;
            else if(m==Matlab::enumToNatural("KrylovStop","MaxItersExceeded"))
                return MaxItersExceeded;
            else if(m==Matlab::enumToNatural(
                "KrylovStop","TrustRegionViolated")
            )
                return TrustRegionViolated;
            else if(m==Matlab::enumToNatural("KrylovStop","Instability"))
                return Instability;
            else if(m==Matlab::enumToNatural(
                "KrylovStop","InvalidTrustRegionCenter")
            )
                return InvalidTrustRegionCenter;
        }
    }
    
    namespace KrylovSolverTruncated { 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & truncated_krylov) {
            // Do the conversion
            switch(truncated_krylov){
            case ConjugateDirection:
                return Matlab::enumToMxArray(
                    "KrylovSolverTruncated","ConjugateDirection");
            case MINRES:
                return Matlab::enumToMxArray(
                    "KrylovSolverTruncated","MINRES");
            default:
                throw;
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member) {
            // Convert the member to a Natural 
            Natural m(*mxGetPr(member));

            if(m==Matlab::enumToNatural(
                "KrylovSolverTruncated","ConjugateDirection")
            )
                return ConjugateDirection;
            else if(m==Matlab::enumToNatural("KrylovSolverTruncated","MINRES"))
                return MINRES;
        }
    }

    namespace AlgorithmClass { 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & algorithm_class) {
            // Do the conversion
            switch(algorithm_class){
            case TrustRegion:
                return Matlab::enumToMxArray("AlgorithmClass","TrustRegion");
            case LineSearch:
                return Matlab::enumToMxArray("AlgorithmClass","LineSearch");
            case UserDefined:
                return Matlab::enumToMxArray("AlgorithmClass","UserDefined");
            default:
                throw;
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member) {
            // Convert the member to a Natural 
            Natural m(*mxGetPr(member));

            if(m==Matlab::enumToNatural("AlgorithmClass","TrustRegion"))
                return TrustRegion;
            else if(m==Matlab::enumToNatural("AlgorithmClass","LineSearch"))
                return LineSearch;
            else if(m==Matlab::enumToNatural("AlgorithmClass","UserDefined"))
                return UserDefined;
        }
    }

    namespace Operators { 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & op) {
            // Do the conversion
            switch(op){
            case Identity:
                return Matlab::enumToMxArray("Operators","Identity");
            case ScaledIdentity:
                return Matlab::enumToMxArray("Operators","ScaledIdentity");
            case BFGS:
                return Matlab::enumToMxArray("Operators","BFGS");
            case InvBFGS:
                return Matlab::enumToMxArray("Operators","InvBFGS");
            case SR1:
                return Matlab::enumToMxArray("Operators","SR1");
            case InvSR1:
                return Matlab::enumToMxArray("Operators","InvSR1");
            case UserDefined:
                return Matlab::enumToMxArray("Operators","UserDefined");
            default:
                throw;
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member) {
            // Convert the member to a Natural 
            Natural m(*mxGetPr(member));

            if(m==Matlab::enumToNatural("Operators","Identity"))
                return Identity;
            else if(m==Matlab::enumToNatural("Operators","ScaledIdentity"))
                return ScaledIdentity;
            else if(m==Matlab::enumToNatural("Operators","BFGS"))
                return BFGS;
            else if(m==Matlab::enumToNatural("Operators","InvBFGS"))
                return InvBFGS;
            else if(m==Matlab::enumToNatural("Operators","SR1"))
                return SR1;
            else if(m==Matlab::enumToNatural("Operators","InvSR1"))
                return InvSR1;
            else if(m==Matlab::enumToNatural("Operators","UserDefined"))
                return UserDefined;
        }
    }

    namespace LineSearchDirection {
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & dir) {
            // Do the conversion
            switch(dir){
            case SteepestDescent:
                return Matlab::enumToMxArray("LineSearchDirection",
                    "SteepestDescent");
            case FletcherReeves:
                return Matlab::enumToMxArray("LineSearchDirection",
                    "FletcherReeves");
            case PolakRibiere:
                return Matlab::enumToMxArray("LineSearchDirection",
                    "PolakRibiere");
            case HestenesStiefel:
                return Matlab::enumToMxArray("LineSearchDirection",
                    "HestenesStiefel");
            case BFGS:
                return Matlab::enumToMxArray("LineSearchDirection","BFGS");
            case NewtonCG:
                return Matlab::enumToMxArray("LineSearchDirection","NewtonCG");
            default:
                throw;
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member) {
            // Convert the member to a Natural 
            Natural m(*mxGetPr(member));

            if(m==Matlab::enumToNatural("LineSearchDirection",
                "SteepestDescent")
            )
                return SteepestDescent;
            else if(m==Matlab::enumToNatural("LineSearchDirection",
                "FletcherReeves")
            )
                return FletcherReeves;
            else if(m==Matlab::enumToNatural("LineSearchDirection",
                "PolakRibiere")
            )
                return PolakRibiere;
            else if(m==Matlab::enumToNatural("LineSearchDirection",
                "HestenesStiefel")
            )
                return HestenesStiefel;
            else if(m==Matlab::enumToNatural("LineSearchDirection","BFGS"))
                return BFGS;
            else if(m==Matlab::enumToNatural("LineSearchDirection","NewtonCG"))
                return NewtonCG;
        }
    }

    namespace LineSearchKind { 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & kind) {
            // Do the conversion
            switch(kind){
            case Brents:
                return Matlab::enumToMxArray("LineSearchKind","Brents");
            case GoldenSection:
                return Matlab::enumToMxArray("LineSearchKind","GoldenSection");
            case BackTracking:
                return Matlab::enumToMxArray("LineSearchKind","BackTracking");
            case TwoPointA:
                return Matlab::enumToMxArray("LineSearchKind","TwoPointA");
            case TwoPointB:
                return Matlab::enumToMxArray("LineSearchKind","TwoPointB");
            default:
                throw;
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member) {
            // Convert the member to a Natural 
            Natural m(*mxGetPr(member));

            if(m==Matlab::enumToNatural("LineSearchKind","Brents"))
                return Brents;
            else if(m==Matlab::enumToNatural("LineSearchKind","GoldenSection"))
                return GoldenSection;
            else if(m==Matlab::enumToNatural("LineSearchKind","BackTracking"))
                return BackTracking;
            else if(m==Matlab::enumToNatural("LineSearchKind","TwoPointA"))
                return TwoPointA;
            else if(m==Matlab::enumToNatural("LineSearchKind","TwoPointB"))
                return TwoPointB;
        }
    }

    namespace OptimizationLocation { 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & loc) {
            // Do the conversion
            switch(loc){
            case BeginningOfOptimization:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","BeginningOfOptimization");
            case BeforeInitialFuncAndGrad:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","BeforeInitialFuncAndGrad");
            case AfterInitialFuncAndGrad:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","AfterInitialFuncAndGrad");
            case BeforeOptimizationLoop:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","BeforeOptimizationLoop");
            case BeginningOfOptimizationLoop:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","BeginningOfOptimizationLoop");
            case BeforeSaveOld:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","BeforeSaveOld");
            case BeforeStep:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","BeforeStep");
            case BeforeGetStep:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","BeforeGetStep");
            case GetStep:
                return Matlab::enumToMxArray("OptimizationLocation","GetStep");
            case AfterStepBeforeGradient:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","AfterStepBeforeGradient");
            case AfterGradient:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","AfterGradient");
            case BeforeQuasi:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","BeforeQuasi");
            case AfterQuasi:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","AfterQuasi");
            case EndOfOptimizationIteration:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","EndOfOptimizationIteration");
            case BeforeLineSearch:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","BeforeLineSearch");
            case AfterRejectedTrustRegion:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","AfterRejectedTrustRegion");
            case AfterRejectedLineSearch:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","AfterRejectedLineSearch");
            case BeforeActualVersusPredicted:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","BeforeActualVersusPredicted");
            case EndOfKrylovIteration:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","EndOfKrylovIteration");
            case EndOfOptimization:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","EndOfOptimization");
            default:
                throw;
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member) {
            // Convert the member to a Natural 
            Natural m(*mxGetPr(member));

            if(m==Matlab::enumToNatural(
                "OptimizationLocation","BeginningOfOptimization"))
                return BeginningOfOptimization;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","BeforeInitialFuncAndGrad"))
                return BeforeInitialFuncAndGrad;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","AfterInitialFuncAndGrad"))
                return AfterInitialFuncAndGrad;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","BeforeOptimizationLoop"))
                return BeforeOptimizationLoop;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","BeginningOfOptimizationLoop"))
                return BeginningOfOptimizationLoop;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","BeforeSaveOld"))
                return BeforeSaveOld;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","BeforeStep"))
                return BeforeStep;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","BeforeGetStep"))
                return BeforeGetStep;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","GetStep"))
                return GetStep;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","AfterStepBeforeGradient"))
                return AfterStepBeforeGradient;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","AfterGradient"))
                return AfterGradient;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","BeforeQuasi"))
                return BeforeQuasi;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","AfterQuasi"))
                return AfterQuasi;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","EndOfOptimizationIteration"))
                return EndOfOptimizationIteration;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","BeforeLineSearch"))
                return BeforeLineSearch;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","AfterRejectedTrustRegion"))
                return AfterRejectedTrustRegion;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","AfterRejectedLineSearch"))
                return AfterRejectedLineSearch;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","BeforeActualVersusPredicted"))
                return BeforeActualVersusPredicted;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","EndOfKrylovIteration"))
                return EndOfKrylovIteration;
            else if(m==Matlab::enumToNatural(
                "OptimizationLocation","EndOfOptimization"))
                return EndOfOptimization;
        }
    }

    namespace InteriorPointMethod { 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & ipm) {
            // Do the conversion
            switch(ipm){
            case PrimalDual:
                return Matlab::enumToMxArray("InteriorPointMethod",
                    "PrimalDual");
            case PrimalDualLinked:
                return Matlab::enumToMxArray("InteriorPointMethod",
                    "PrimalDualLinked");
            case LogBarrier:
                return Matlab::enumToMxArray("InteriorPointMethod",
                    "LogBarrier");
            default:
                throw;
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member) {
            // Convert the member to a Natural 
            Natural m(*mxGetPr(member));

            if(m==Matlab::enumToNatural("InteriorPointMethod","PrimalDual"))
                return PrimalDual;
            else if(m==Matlab::enumToNatural("InteriorPointMethod",
                "PrimalDualLinked")
            )
                return PrimalDualLinked;
            else if(m==Matlab::enumToNatural("InteriorPointMethod",
                "LogBarrier")
            )
                return LogBarrier;
        }
    }

    namespace CentralityStrategy { 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & cstrat) {
            // Do the conversion
            switch(cstrat){
            case Constant:
                return Matlab::enumToMxArray("CentralityStrategy","Constant");
            case StairStep:
                return Matlab::enumToMxArray("CentralityStrategy","StairStep");
            case PredictorCorrector:
                return Matlab::enumToMxArray("CentralityStrategy",
                    "PredictorCorrector");
            default:
                throw;
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member) {
            // Convert the member to a Natural 
            Natural m(*mxGetPr(member));

            if(m==Matlab::enumToNatural("CentralityStrategy","Constant"))
                return Constant;
            else if(m==Matlab::enumToNatural("CentralityStrategy","StairStep"))
                return StairStep;
            else if(m==Matlab::enumToNatural("CentralityStrategy",
                "PredictorCorrector")
            )
                return PredictorCorrector;
        }
    }

    namespace FunctionDiagnostics { 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & diag) {
            // Do the conversion
            switch(diag){
            case NoDiagnostics:
                return Matlab::enumToMxArray("FunctionDiagnostics",
                    "NoDiagnostics");
            case FirstOrder:
                return Matlab::enumToMxArray("FunctionDiagnostics",
                    "FirstOrder");
            case SecondOrder:
                return Matlab::enumToMxArray("FunctionDiagnostics",
                    "SecondOrder");
            default:
                throw;
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member) {
            // Convert the member to a Natural 
            Natural m(*mxGetPr(member));

            if(m==Matlab::enumToNatural("FunctionDiagnostics","NoDiagnostics"))
                return NoDiagnostics;
            else if(m==Matlab::enumToNatural("FunctionDiagnostics",
                "FirstOrder")
            )
                return FirstOrder;
            else if(m==Matlab::enumToNatural("FunctionDiagnostics",
                "SecondOrder")
            )
                return SecondOrder;
        }
    }

    namespace DiagnosticScheme { 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & dscheme) {
            // Do the conversion
            switch(dscheme){
            case Never:
                return Matlab::enumToMxArray("DiagnosticScheme","Never");
            case DiagnosticsOnly:
                return Matlab::enumToMxArray("DiagnosticScheme",
                    "DiagnosticsOnly");
            case EveryIteration:
                return Matlab::enumToMxArray("DiagnosticScheme",
                    "EveryIteration");
            default:
                throw;
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member) {
            // Convert the member to a Natural 
            Natural m(*mxGetPr(member));

            if(m==Matlab::enumToNatural("DiagnosticScheme","Never"))
                return Never;
            else if(m==Matlab::enumToNatural("DiagnosticScheme",
                "DiagnosticsOnly")
            )
                return DiagnosticsOnly;
            else if(m==Matlab::enumToNatural("DiagnosticScheme",
                "EveryIteration")
            )
                return EveryIteration;
        }
    }

    namespace Matlab {
        // Calls a Matlab function with one argument 
        std::pair <mxArray *,int> mxArray_CallObject1(
            mxArray * const fn,
            mxArray * const arg1
        ) {
            mxArray* input[2]={fn,arg1};
            mxArray* output[1];
            int err=mexCallMATLAB(1,output,2,input,"feval");
            return std::pair <mxArray*,int> (output[0],err);
        }
        
        // Calls a Matlab function with two arguments
        std::pair <mxArray *,int> mxArray_CallObject2(
            mxArray * const fn,
            mxArray * const arg1,
            mxArray * const arg2
        ) {
            mxArray* input[3]={fn,arg1,arg2};
            mxArray* output[1];
            int err=mexCallMATLAB(1,output,3,input,"feval");
            return std::pair <mxArray*,int> (output[0],err);
        }
        
        // Calls a Matlab function with three arguments
        std::pair <mxArray *,int> mxArray_CallObject3(
            mxArray * const fn,
            mxArray * const arg1,
            mxArray * const arg2,
            mxArray * const arg3
        ) {
            mxArray* input[4]={fn,arg1,arg2,arg3};
            mxArray* output[1];
            int err=mexCallMATLAB(1,output,3,input,"feval");
            return std::pair <mxArray*,int> (output[0],err);
        }
        
        // Calls a Matlab function with four arguments
        std::pair <mxArray *,int> mxArray_CallObject4(
            mxArray * const fn,
            mxArray * const arg1,
            mxArray * const arg2,
            mxArray * const arg3,
            mxArray * const arg4
        ) {
            mxArray* input[5]={fn,arg1,arg2,arg3,arg4};
            mxArray* output[1];
            int err=mexCallMATLAB(1,output,5,input,"feval");
            return std::pair <mxArray*,int> (output[0],err);
        }

        // Creates a Matlab double from a C++ double
        mxArray * mxArray_FromDouble(double const x_) {
            mxArray * x(mxCreateDoubleMatrix(1,1,mxREAL));
            mxGetPr(x)[0]=x_;
        }

        // Imports the Optizelle structure
        mxArray * importOptizelle() {
            mxArray* output[1];
            int err=mexCallMATLAB(1,output,0,nullptr,"setupOptizelle");
            if(err)
                mexErrMsgTxt("Unable to import the Optizelle structure");
            return output[0];
        }

        // Converts an Optizelle enumerated type to a mxArray * 
        mxArray * enumToMxArray(
            std::string const & type,
            std::string const & member 
        ) {
            // Grab the Optizelle module
            mxArrayPtr optizelle(importOptizelle());

            // Grab the type enumerated type
            mxArrayPtr matenum(mxGetField(optizelle.get(),0,type.c_str()));

            // Finally, grab the actual specified member in this enumerated type
            return mxGetField(matenum.get(),0,member.c_str());
        }
       
        // Converts an Optizelle enumerated type to a Natural
        Natural enumToNatural(
            std::string const & type,
            std::string const & member 
        ) {
            // Grab the mxArray * for the type and member requested
            mxArrayPtr obj(enumToMxArray(type,member));

            // Convert and return the member
            return Natural(*mxGetPr(obj.get()));
        }

        // On construction, initialize the pointer and figure out if
        // we're capturing the pointer or attaching to it
        mxArrayPtr::mxArrayPtr(
            mxArray * const ptr_,
            mxArrayPtrMode::t const mode_ 
        ) : ptr(ptr_), mode(mode_) { }
            
        // Move constructor
        mxArrayPtr::mxArrayPtr(mxArrayPtr&& ptr_) noexcept
            : ptr(ptr_.release()), mode(ptr_.mode) {}
        
        // Move assignment operator
        mxArrayPtr const & mxArrayPtr::operator=(mxArrayPtr&& ptr_) noexcept {
            ptr=ptr_.release();
            mode=ptr_.mode;
            return *this;
        }

        // For a reset, we destroy the pointer and then assign a new
        // value.
        void mxArrayPtr::reset(mxArray * const ptr_) {
            if(ptr!=nullptr && mode!=mxArrayPtrMode::Attach)
                mxDestroyArray(ptr);
            ptr=ptr_;
            mode = mxArrayPtrMode::Capture;
        }

        // For an attach, we destroy the pointer then assign a new value
        void mxArrayPtr::attach(mxArray * const ptr_) {
            if(ptr!=nullptr && mode!=mxArrayPtrMode::Attach)
                mxDestroyArray(ptr);
            ptr=ptr_;
            mode = mxArrayPtrMode::Attach;
        }

        // On a get, we simply return the pointer.
        mxArray * mxArrayPtr::get() {
            return ptr;
        }
    
        // On a release, we return the underlying pointer and then clear
        // the vector.  This will prevent destruction later. 
        mxArray * mxArrayPtr::release() {
            mxArray * ptr_=ptr;
            ptr=nullptr;
            return ptr_;
        }

        // On destruction, destroy the pointer. 
        mxArrayPtr::~mxArrayPtr() {
            if(ptr!=nullptr && mode!=mxArrayPtrMode::Attach)
                mxDestroyArray(ptr);
            ptr=nullptr;
        }
            
        // On construction, we just grab the pointer to the messaging object
        Messaging::Messaging(
            mxArray * const ptr_,
            mxArrayPtrMode::t const mode
        ) : mxArrayPtr(ptr_,mode) {}
            
        // Move constructor
        Messaging::Messaging(Messaging && msg) noexcept
            : mxArrayPtr(msg.release(),msg.mode) {}

        // Move assignment operator
        Messaging const & Messaging::operator = (Messaging && msg) noexcept {
            ptr = msg.release();
            mode = msg.mode;
            return *this;
        }
            
        // Prints a message
        void Messaging::print(std::string const & msg_) const {
            // Call the print function on msg
            mxArrayPtr print(mxGetField(ptr,0,"print"));
            mxArrayPtr msg(mxCreateString(msg_.c_str()));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject1(
                print.get(),
                msg.get()));

            // Check errors
            if(ret_err.second !=0)
                error("Evaluation of the print function in the Messaging "
                    "object failed.");
        }

        // Prints an error
        void Messaging::error(std::string const & msg_) const {
            // Call the error function on msg
            mxArrayPtr error(mxGetField(ptr,0,"error"));
            mxArrayPtr msg(mxCreateString(msg_.c_str()));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject1(
                error.get(),
                msg.get()));

            // Check errors
            if(ret_err.second !=0) {
                std::string msg2="Evaluation of the error function in the "
                    "Messaging object failed.\n";
                mexErrMsgTxt(msg2.c_str());
            }
        }

        // Create a vector with the appropriate messaging and vector space 
        Vector::Vector(
            mxArray * const msg_,
            mxArray * const vs_,
            mxArray * const vec,
            mxArrayPtrMode::t mode
        ) : 
            mxArrayPtr(vec,mode),
            msg(msg_,mxArrayPtrMode::Attach),
            vs(vs_,mxArrayPtrMode::Attach)
        {}
            
        // Create a move constructor so we can interact with stl objects
        Vector::Vector(Vector && vec) noexcept :
            mxArrayPtr(std::move(vec)),
            msg(std::move(vec.msg)),
            vs(std::move(vec.vs))
        { }
            
        // Move assignment operator
        Vector const & Vector::operator = (Vector && vec) noexcept {
            ptr = vec.release(); 
            mode = vec.mode;
            msg = std::move(vec.msg);
            vs = std::move(vec.vs);
            return *this;
        }

        // Memory allocation and size setting 
        Vector Vector::init() { 
            // Call the init function on the internal and store in y 
            mxArrayPtr init(mxGetField(vs.get(),0,"init"));
            std::pair <mxArrayPtr,int> y_err(mxArray_CallObject1(
                init.get(),
                get()));

            // Check errors
            if(y_err.second!=0)
                msg.error(
                    "Evaluation of the vector space function init failed.");

            // Create and return a new vector based on y
            return std::move(Vector(msg.get(),vs.get(),y_err.first.release()));
        } 
        
        // y <- x (Shallow.  No memory allocation.)  Internal is y.
        void Vector::copy(Vector & x) { 
            // Call the copy function on x and the internal 
            mxArrayPtr copy(mxGetField(vs.get(),0,"copy"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject1(
                copy.get(),
                x.get()));

            // Check errors
            if(ret_err.second!=0)
                msg.error(
                    "Evaluation of the vector space function copy failed.");

            // Assign y
            reset(ret_err.first.release());
        } 

        // x <- alpha * x.  Internal is x.
        void Vector::scal(double const & alpha_) { 
            // Call the scal function on alpha and the internal storage 
            mxArrayPtr scal(mxGetField(vs.get(),0,"scal"));
            mxArrayPtr alpha(mxArray_FromDouble(alpha_));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject1(
                scal.get(),
                alpha.get()));

            // Check errors
            if(ret_err.second!=0)
                msg.error(
                    "Evaluation of the vector space function scal failed.");

            // Assign x
            reset(ret_err.first.release());
        } 

        // x <- 0.  Internal is x. 
        void Vector::zero() { 
            // Call the zero function on this vector.
            mxArrayPtr zero(mxGetField(vs.get(),0,"zero"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject1(
                zero.get(),
                get()));

            // Check errors
            if(ret_err.second!=0)
                msg.error(
                    "Evaluation of the vector space function zero failed.");

            // Assign x
            reset(ret_err.first.release());
        } 

        // y <- alpha * x + y.   Internal is y.
        void Vector::axpy(double const & alpha_,Vector & x) { 
            // Call the axpy function on alpha, x, and the internal storage.
            mxArrayPtr axpy(mxGetField(vs.get(),0,"axpy"));
            mxArrayPtr alpha(mxArray_FromDouble(alpha_));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject3(
                axpy.get(),
                alpha.get(),
                x.get(),
                get()));
           
            // Check errors
            if(ret_err.second!=0)
                msg.error(
                    "Evaluation of the vector space function axpy failed.");

            // Assign y
            reset(ret_err.first.release());
        } 

        // innr <- <x,y>.  Internal is y.
        double Vector::innr(Vector & x) { 
            // Call the innr function on x and the internal.  Store in z. 
            mxArrayPtr innr(mxGetField(vs.get(),0,"innr"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject2(
                innr.get(),
                x.get(),
                get()));

            // Check errors
            if(ret_err.second!=0)
                msg.error(
                    "Evaluation of the vector space function innr failed.");

            // Return the result 
            return mxGetPr(ret_err.first.get())[0]; 
        } 

        // x <- random.  Internal is x. 
        void Vector::rand() { 
            // Call the rand function on this vector.
            mxArrayPtr rand(mxGetField(vs.get(),0,"rand"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject1(
                rand.get(),
                get()));

            // Check errors
            if(ret_err.second!=0)
                msg.error(
                    "Evaluation of the vector space function rand failed.");
            
            // Assign x
            reset(ret_err.first.release());
        } 

        // Jordan product, z <- x o y.  Internal is z.
        void Vector::prod(Vector & x,Vector & y) { 
            // Call the prod function on x, y, and the internal 
            mxArrayPtr prod(mxGetField(vs.get(),0,"prod"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject2(
                prod.get(),
                x.get(),
                y.get()));

            // Check errors
            if(ret_err.second!=0)
                msg.error(
                    "Evaluation of the vector space function prod failed.");
            
            // Assign z
            reset(ret_err.first.release());
        } 

        // Identity element, x <- e such that x o e = x .  Internal is x.
        void Vector::id() { 
            // Call the id function on the internal.
            mxArrayPtr id(mxGetField(vs.get(),0,"id"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject1(
                id.get(),
                get()));

            // Check errors
            if(ret_err.second!=0)
                msg.error(
                    "Evaluation of the vector space function id failed.");
            
            // Assign x
            reset(ret_err.first.release());
        } 

        // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y.
        // Internal is z.
        void Vector::linv(Vector& x, Vector& y) { 
            // Call the linv function on x, y, and the internal
            mxArrayPtr linv(mxGetField(vs.get(),0,"linv"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject2(
                linv.get(),
                x.get(),
                y.get()));

            // Check errors
            if(ret_err.second!=0)
                msg.error(
                    "Evaluation of the vector space function linv failed.");
            
            // Assign z
            reset(ret_err.first.release());
        } 

        // Barrier function, barr <- barr(x) where x o grad barr(x) = e.
        // Internal is x.
        double Vector::barr() { 
            // Call the barr function on the internal.  Store in z.
            mxArrayPtr barr(mxGetField(vs.get(),0,"barr"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject1(
                barr.get(),
                get()));

            // Check errors
            if(ret_err.second!=0)
                msg.error(
                    "Evaluation of the vector space function barr failed.");

            // Return the result 
            return mxGetPr(ret_err.first.get())[0]; 
        } 

        // Line search, srch <- argmax {alpha in Real >= 0 : alpha x + y >= 0} 
        // where y > 0.  Internal is y.
        double Vector::srch(Vector& x) {  
            // Call the srch function on x and the internal.  Store in z.
            mxArrayPtr srch(mxGetField(vs.get(),0,"srch"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject2(
                srch.get(),
                x.get(),
                get()));

            // Check errors
            if(ret_err.second!=0)
                msg.error(
                    "Evaluation of the vector space function srch failed.");

            // Return the result 
            return mxGetPr(ret_err.first.get())[0]; 
        } 

        // Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
        // operator.  Internal is x.
        void Vector::symm() { 
            // Call the symm function on the internal.
            mxArrayPtr symm(mxGetField(vs.get(),0,"symm"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject1(
                symm.get(),
                get()));

            // Check errors
            if(ret_err.second!=0)
                msg.error(
                    "Evaluation of the vector space function symm failed.");
            
            // Assign x
            reset(ret_err.first.release());

        } 
        
        // Converts (copies) a value into Matlab.  
        mxArray * Vector::toMatlab() {
            // Call the copy function on the internal and x
            mxArrayPtr copy(mxGetField(vs.get(),0,"copy"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject1(
                copy.get(),
                get()));

            // Check errors
            if(ret_err.second!=0)
                msg.error(
                    "Evaluation of the vector space function copy failed.");
            
            // Return the pointer 
            return ret_err.first.release();
        } 
        
        // Converts (copies) a value from Matlab.  This assumes that the
        // vector space functions have already been properly assigned.
        void Vector::fromMatlab(mxArray * const ptr) {
            // Call the copy function on ptr and the internal 
            mxArrayPtr copy(mxGetField(vs.get(),0,"copy"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject1(
                copy.get(),
                ptr));

            // Check errors
            if(ret_err.second!=0)
                msg.error(
                    "Evaluation of the vector space function copy failed.");

            // Copy in the value
            reset(ret_err.first.release());
        } 

#if 0
        // Convert a C++ state to a Matlab state 
        template <>
        void State <MxUnconstrained>::toMatlab(
            typename MxUnconstrained::State::t const & state
        ) {
            Unconstrained::State::toMatlab(state,ptr);
        }
        template <>
        void State <MxEqualityConstrained>::toMatlab(
            typename MxEqualityConstrained::State::t const & state
        ) {
            EqualityConstrained::State::toMatlab(state,ptr);
        }
        template <>
        void State <MxInequalityConstrained>::toMatlab(
            typename MxInequalityConstrained::State::t const & state
        ) {
            InequalityConstrained::State::toMatlab(state,ptr);
        }
        template <>
        void State <MxConstrained>::toMatlab(
            typename MxConstrained::State::t const & state
        ) {
            Constrained::State::toMatlab(state,ptr);
        }

        // Convert a Matlab state to C++ 
        template <>
        void State <MxUnconstrained>::fromMatlab(
            typename MxUnconstrained::State::t & state
        ) {
            Unconstrained::State::fromMatlab(ptr,state);
        }
        template <>
        void State <MxEqualityConstrained>::fromMatlab(
            typename MxEqualityConstrained::State::t & state
        ) {
            EqualityConstrained::State::fromMatlab(ptr,state);
        }
        template <>
        void State <MxInequalityConstrained>::fromMatlab(
            typename MxInequalityConstrained::State::t & state
        ) {
            InequalityConstrained::State::fromMatlab(ptr,state);
        }
        template <>
        void State <MxConstrained>::fromMatlab(
            typename MxConstrained::State::t & state
        ) {
            Constrained::State::fromMatlab(ptr,state);
        }
        
        // Convert a Matlab bundle to C++ 
        template <>
        void Functions <MxUnconstrained>::fromMatlab(
            typename MxUnconstrained::Functions::t & fns 
        ) {
            Unconstrained::Functions::fromMatlab(
                msg.get(),ptr,mxstate.get(),state,fns);
        }
        template <>
        void Functions <MxEqualityConstrained>::fromMatlab(
            typename MxEqualityConstrained::Functions::t & fns 
        ) {
            EqualityConstrained::Functions::fromMatlab(
                msg.get(),ptr,mxstate.get(),state,fns);
        }
        template <>
        void Functions <MxInequalityConstrained>::fromMatlab(
            typename MxInequalityConstrained::Functions::t & fns 
        ) {
            InequalityConstrained::Functions::fromMatlab(
                msg.get(),ptr,mxstate.get(),state,fns);
        }
        template <>
        void Functions <MxConstrained>::fromMatlab(
            typename MxConstrained::Functions::t & fns 
        ) {
            Constrained::Functions::fromMatlab(
                msg.get(),ptr,mxstate.get(),state,fns);
        }
#endif

        // Create a function 
        ScalarValuedFunction::ScalarValuedFunction(
            mxArray * const msg_,
            mxArray * const f,
            mxArrayPtrMode::t mode
        ) :
            mxArrayPtr(f,mode),
            msg(msg_,mxArrayPtrMode::Attach)
        { }

        // <- f(x) 
        double ScalarValuedFunction::eval(Vector const & x) const { 
            // Call the objective function on x.
            mxArrayPtr eval(mxGetField(ptr,0,"eval"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject1(
                eval.get(),
                const_cast <Vector &> (x).get()));

            // Check errors
            if(ret_err.second!=0)
                msg.error("Evaluation of the objective f failed.");

            // Return the result 
            return mxGetPr(ret_err.first.get())[0]; 
        }

        // grad = grad f(x) 
        void ScalarValuedFunction::grad(
            Vector const & x,
            Vector & grad
        ) const { 
            // Call the gradient function on x
            mxArrayPtr mxgrad(mxGetField(ptr,0,"grad"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject1(
                mxgrad.get(),
                const_cast <Vector &>(x).get()));

            // Check errors
            if(ret_err.second!=0)
                msg.error("Evaluation of the gradient of f failed.");

            // Assign grad 
            grad.reset(ret_err.first.release());
        }

        // H_dx = hess f(x) dx 
        void ScalarValuedFunction::hessvec(
            Vector const & x,
            Vector const & dx,
            Vector & H_dx
        ) const {
            // Call the hessvec function on x and dx,
            mxArrayPtr hessvec(mxGetField(ptr,0,"hessvec"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject2(
                hessvec.get(),
                const_cast <Vector &> (x).get(),
                const_cast <Vector &> (dx).get()));

            // Check errors
            if(ret_err.second!=0)
                msg.error("Evaluation of the Hessian-vector product"
                    " of f failed.");
            
            // Assign H_dx
            H_dx.reset(ret_err.first.release());
        }

        // Create a function 
        VectorValuedFunction::VectorValuedFunction(
            std::string const & name_,
            mxArray * const msg_,
            mxArray * const f,
            mxArrayPtrMode::t mode
        ) :
            mxArrayPtr(f,mode),
            msg(msg_,mxArrayPtrMode::Attach),
            name(name_)
        {}

        // y=f(x)
        void VectorValuedFunction::eval(
            Vector const & x,
            VectorValuedFunction::Y_Vector& y
        ) const {
            // Call the objective function on x.
            mxArrayPtr eval(mxGetField(ptr,0,"eval"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject1(
                eval.get(),
                const_cast <Vector &> (x).get()));

            // Check errors
            if(ret_err.second!=0) {
                std::stringstream ss;
                ss << "Evaluation of the constraint " << name << " failed.";
                msg.error(ss.str());
            }
            
            // Assign y 
            y.reset(ret_err.first.release());
        }

        // y=f'(x)dx 
        void VectorValuedFunction::p(
            Vector const & x,
            Vector const & dx,
            VectorValuedFunction::Y_Vector& y
        ) const {
            // Call the prime function on x and dx
            mxArrayPtr p(mxGetField(ptr,0,"p"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject2(
                p.get(),
                const_cast <Vector &> (x).get(),
                const_cast <Vector &> (dx).get()));
           
            // Check errors
            if(ret_err.second!=0) {
                std::stringstream ss;
                ss << "Evaluation of the derivative of the constraint "
                    << name << " failed.";
                msg.error(ss.str());
            }
            
            // Assign y 
            y.reset(ret_err.first.release());
        }

        // z=f'(x)*dy
        void VectorValuedFunction::ps(
            Vector const & x,
            Vector const & dy,
            VectorValuedFunction::X_Vector& z
        ) const {
            // Call the prime-adjoint function on x and dy
            mxArrayPtr ps(mxGetField(ptr,0,"ps"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject2(
                ps.get(),
                const_cast <Vector &> (x).get(),
                const_cast <Vector &> (dy).get()));

            // Check errors
            if(ret_err.second!=0) {
                std::stringstream ss;
                ss << "Evaluation of the derivative-adjoint of the constraint "
                    << name << " failed.";
                msg.error(ss.str());
            }
            
            // Assign z
            z.reset(ret_err.first.release());
        }
             
        // z=(f''(x)dx)*dy
        void VectorValuedFunction::pps(
            Vector const & x,
            Vector const & dx,
            Vector const & dy,
            X_Vector& z
        ) const { 
            // Call the prime-adjoint function on x, dx, and dy
            mxArrayPtr pps(mxGetField(ptr,0,"pps"));
            std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject3(
                pps.get(),
                const_cast <Vector &> (x).get(),
                const_cast <Vector &> (dx).get(),
                const_cast <Vector &> (dy).get()));

            // Check errors
            if(ret_err.second!=0) {
                std::stringstream ss;
                ss << "Evaluation of the second derivative-adjoint of the "
                    "constraint " << name << " failed.";
                msg.error(ss.str());
            }
            
            // Assign z 
            z.reset(ret_err.first.release());
        }
    }
}
