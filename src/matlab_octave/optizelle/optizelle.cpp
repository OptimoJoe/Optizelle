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

#include "optizelle.h"

#define CATCH_MATLAB_ERRORS \
    catch(Matlab::Exception::t const & e) { \
        return; \
    } catch(std::exception const & e) { \
        mexErrMsgTxt( \
            Optizelle::Exception::exception_to_string(e).c_str()); \
        return; \
    }

namespace Optizelle {
    // In theory, I'd like to keep this variable local, but for some strange
    // reason Octave keeps moving around this memory.  As such, rather than
    // having it be static, I made it global and then grab the correct pointer
    // each time we enter these routines.
    static std::list<mxArray *> optizelle;

    // Enumerate types
    namespace OptimizationStop { 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & opt_stop) {
            // Do the conversion
            switch(opt_stop){
            case NotConverged:
                return Matlab::enumToMxArray(
                    "OptimizationStop","NotConverged");
            case GradientSmall:
                return Matlab::enumToMxArray(
                    "OptimizationStop","GradientSmall");
            case StepSmall:
                return Matlab::enumToMxArray(
                    "OptimizationStop","StepSmall");
            case MaxItersExceeded:
                return Matlab::enumToMxArray(
                    "OptimizationStop","MaxItersExceeded");
            case InteriorPointInstability:
                return Matlab::enumToMxArray(
                    "OptimizationStop","InteriorPointInstability");
            case GlobalizationFailure:
                return Matlab::enumToMxArray(
                    "OptimizationStop","GlobalizationFailure");
            case UserDefined:
                return Matlab::enumToMxArray(
                    "OptimizationStop","UserDefined");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = *mxGetPr(member.get());

            if(m==Matlab::enumToNatural("OptimizationStop","NotConverged"))
                return NotConverged;
            else if(m==Matlab::enumToNatural(
                "OptimizationStop","GradientSmall")
            )
                return GradientSmall;
            else if(m==Matlab::enumToNatural(
                "OptimizationStop","StepSmall")
            )
                return StepSmall;
            else if(m==Matlab::enumToNatural(
                "OptimizationStop","MaxItersExceeded")
            )
                return MaxItersExceeded;
            else if(m==Matlab::enumToNatural(
                "OptimizationStop","InteriorPointInstability")
            )
                return InteriorPointInstability;
            else if(m==Matlab::enumToNatural(
                "OptimizationStop","GlobalizationFailure")
            )
                return GlobalizationFailure;
            else if(m==Matlab::enumToNatural("OptimizationStop","UserDefined"))
                return UserDefined;
            else
                throw Matlab::Exception::t( __LOC__
                    + ", unknown OptimizationStop");
        }
    }
    
    namespace TruncatedStop { 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & trunc_stop) {
            // Do the conversion
            switch(trunc_stop){
            case NotConverged:
                return Matlab::enumToMxArray(
                    "TruncatedStop","NotConverged");
            case NegativeCurvature:
                return Matlab::enumToMxArray(
                    "TruncatedStop","NegativeCurvature");
            case RelativeErrorSmall:
                return Matlab::enumToMxArray(
                    "TruncatedStop","RelativeErrorSmall");
            case MaxItersExceeded:
                return Matlab::enumToMxArray(
                    "TruncatedStop","MaxItersExceeded");
            case TrustRegionViolated:
                return Matlab::enumToMxArray(
                    "TruncatedStop","TrustRegionViolated");
            case NanOperator:
                return Matlab::enumToMxArray("TruncatedStop","NanOperator");
            case NanPreconditioner:
                return Matlab::enumToMxArray(
                    "TruncatedStop","NanPreconditioner");
            case NonProjectorPreconditioner:
                return Matlab::enumToMxArray(
                    "TruncatedStop","NonProjectorPreconditioner");
            case NonSymmetricPreconditioner:
                return Matlab::enumToMxArray(
                    "TruncatedStop","NonSymmetricPreconditioner");
            case NonSymmetricOperator:
                return Matlab::enumToMxArray(
                    "TruncatedStop","NonSymmetricOperator");
            case LossOfOrthogonality:
                return Matlab::enumToMxArray(
                    "TruncatedStop","LossOfOrthogonality");
            case OffsetViolatesTrustRegion:
                return Matlab::enumToMxArray(
                    "TruncatedStop","OffsetViolatesTrustRegion");
            case OffsetViolatesSafeguard:
                return Matlab::enumToMxArray(
                    "TruncatedStop","OffsetViolatesSafeguard");
            case TooManyFailedSafeguard:
                return Matlab::enumToMxArray(
                    "TruncatedStop","TooManyFailedSafeguard");
            case ObjectiveIncrease:
                return Matlab::enumToMxArray(
                    "TruncatedStop","ObjectiveIncrease");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = *mxGetPr(member.get());

            if(m==Matlab::enumToNatural("TruncatedStop","NotConverged"))
                return NotConverged;
            else if(m==Matlab::enumToNatural(
                "TruncatedStop","NegativeCurvature")
            )
                return NegativeCurvature;
            else if(m==Matlab::enumToNatural(
                "TruncatedStop","RelativeErrorSmall")
            )
                return RelativeErrorSmall;
            else if(m==Matlab::enumToNatural(
                "TruncatedStop","MaxItersExceeded")
            )
                return MaxItersExceeded;
            else if(m==Matlab::enumToNatural(
                "TruncatedStop","TrustRegionViolated")
            )
                return TrustRegionViolated;
            else if(m==Matlab::enumToNatural("TruncatedStop","NanOperator"))
                return NanOperator;
            else if(m==Matlab::enumToNatural(
                "TruncatedStop","NanPreconditioner")
            )
                return NanPreconditioner;
            else if(m==Matlab::enumToNatural(
                "TruncatedStop","NonProjectorPreconditioner")
            )
                return NonProjectorPreconditioner;
            else if(m==Matlab::enumToNatural(
                "TruncatedStop","NonSymmetricPreconditioner")
            )
                return NonSymmetricPreconditioner;
            else if(m==Matlab::enumToNatural(
                "TruncatedStop","NonSymmetricOperator")
            )
                return NonSymmetricOperator;
            else if(m==Matlab::enumToNatural(
                "TruncatedStop","LossOfOrthogonality")
            )
                return LossOfOrthogonality;
            else if(m==Matlab::enumToNatural(
                "TruncatedStop","OffsetViolatesTrustRegion")
            )
                return OffsetViolatesTrustRegion;
            else if(m==Matlab::enumToNatural(
                "TruncatedStop","OffsetViolatesSafeguard")
            )
                return OffsetViolatesSafeguard;
            else if(m==Matlab::enumToNatural(
                "TruncatedStop","TooManyFailedSafeguard")
            )
                return TooManyFailedSafeguard;
            else if(m==Matlab::enumToNatural(
                "TruncatedStop","ObjectiveIncrease")
            )
                return ObjectiveIncrease;
            else
                throw Matlab::Exception::t( __LOC__
                    + ", unknown TruncatedStop");
        }
    }

    namespace AlgorithmClass { 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & algorithm_class) {
            // Do the conversion
            switch(algorithm_class){
            case TrustRegion:
                return Matlab::enumToMxArray("AlgorithmClass","TrustRegion");
            case LineSearch:
                return Matlab::enumToMxArray("AlgorithmClass","LineSearch");
            case UserDefined:
                return Matlab::enumToMxArray("AlgorithmClass","UserDefined");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = *mxGetPr(member.get());

            if(m==Matlab::enumToNatural("AlgorithmClass","TrustRegion"))
                return TrustRegion;
            else if(m==Matlab::enumToNatural("AlgorithmClass","LineSearch"))
                return LineSearch;
            else if(m==Matlab::enumToNatural("AlgorithmClass","UserDefined"))
                return UserDefined;
            else
                throw Matlab::Exception::t( __LOC__
                    + ", unknown AlgorithmClass");
        }
    }

    namespace Operators { 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & op) {
            // Do the conversion
            switch(op){
            case Identity:
                return Matlab::enumToMxArray("Operators","Identity");
            case Zero:
                return Matlab::enumToMxArray("Operators","Zero");
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
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = *mxGetPr(member.get());

            if(m==Matlab::enumToNatural("Operators","Identity"))
                return Identity;
            else if(m==Matlab::enumToNatural("Operators","Zero"))
                return Zero;
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
            else
                throw Matlab::Exception::t( __LOC__
                    + ", unknown Operators");
        }
    }

    namespace LineSearchDirection {
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & dir) {
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
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = *mxGetPr(member.get());

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
            else
                throw Matlab::Exception::t( __LOC__
                    + ", unknown LineSearchDirection");
        }
    }

    namespace LineSearchKind { 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & kind) {
            // Do the conversion
            switch(kind){
            case GoldenSection:
                return Matlab::enumToMxArray("LineSearchKind","GoldenSection");
            case BackTracking:
                return Matlab::enumToMxArray("LineSearchKind","BackTracking");
            case TwoPointA:
                return Matlab::enumToMxArray("LineSearchKind","TwoPointA");
            case TwoPointB:
                return Matlab::enumToMxArray("LineSearchKind","TwoPointB");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = *mxGetPr(member.get());

            if(m==Matlab::enumToNatural("LineSearchKind","GoldenSection"))
                return GoldenSection;
            else if(m==Matlab::enumToNatural("LineSearchKind","BackTracking"))
                return BackTracking;
            else if(m==Matlab::enumToNatural("LineSearchKind","TwoPointA"))
                return TwoPointA;
            else if(m==Matlab::enumToNatural("LineSearchKind","TwoPointB"))
                return TwoPointB;
            else
                throw Matlab::Exception::t( __LOC__
                    + ", unknown LineSearchKind");
        }
    }

    namespace OptimizationLocation { 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & loc) {
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
            case AfterCheckStop:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","AfterCheckStop");
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
            case EndOfOptimization:
                return Matlab::enumToMxArray(
                    "OptimizationLocation","EndOfOptimization");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = *mxGetPr(member.get());

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
                "OptimizationLocation","AfterCheckStop"))
                return AfterCheckStop;
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
                "OptimizationLocation","EndOfOptimization"))
                return EndOfOptimization;
            else
                throw Matlab::Exception::t( __LOC__
                    + ", unknown OptimizationLocation");
        }
    }

    namespace FunctionDiagnostics { 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & diag) {
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
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = *mxGetPr(member.get());

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
            else
                throw Matlab::Exception::t( __LOC__
                    + ", unknown FunctionDiagnostics");
        }
    }

    namespace VectorSpaceDiagnostics { 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & diag) {
            // Do the conversion
            switch(diag){
            case NoDiagnostics:
                return Matlab::enumToMxArray("VectorSpaceDiagnostics",
                    "NoDiagnostics");
            case Basic:
                return Matlab::enumToMxArray("VectorSpaceDiagnostics",
                    "Basic");
            case EuclideanJordan:
                return Matlab::enumToMxArray("VectorSpaceDiagnostics",
                    "EuclideanJordan");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = *mxGetPr(member.get());

            if(m==Matlab::enumToNatural("VectorSpaceDiagnostics",
                "NoDiagnostics")
            )
                return NoDiagnostics;
            else if(m==Matlab::enumToNatural("VectorSpaceDiagnostics",
                "Basic")
            )
                return Basic;
            else if(m==Matlab::enumToNatural("VectorSpaceDiagnostics",
                "EuclideanJordan")
            )
                return EuclideanJordan;
            else
                throw Matlab::Exception::t( __LOC__
                    + ", unknown VectorSpaceDiagnostics");
        }
    }

    namespace DiagnosticScheme { 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & dscheme) {
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
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = *mxGetPr(member.get());

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
            else
                throw Matlab::Exception::t( __LOC__
                    + ", unknown DiagnosticScheme");
        }
    }

    namespace ToleranceKind { 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & eps_kind) {
            // Do the conversion
            switch(eps_kind){
            case Relative:
                return Matlab::enumToMxArray("ToleranceKind",
                    "Relative");
            case Absolute:
                return Matlab::enumToMxArray("ToleranceKind",
                    "Absolute");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = *mxGetPr(member.get());

            if(m==Matlab::enumToNatural("ToleranceKind",
                "Relative")
            )
                return Relative;
            else if(m==Matlab::enumToNatural("ToleranceKind",
                "Absolute")
            )
                return Absolute;
            else
                throw Matlab::Exception::t( __LOC__
                    + ", unknown ToleranceKind");
        }
    }

    namespace QuasinormalStop{ 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & qn_stop) {
            // Do the conversion
            switch(qn_stop){
            case Newton:
                return Matlab::enumToMxArray("QuasinormalStop",
                    "Newton");
            case CauchyTrustRegion:
                return Matlab::enumToMxArray("QuasinormalStop",
                    "CauchyTrustRegion");
            case CauchySafeguard:
                return Matlab::enumToMxArray("QuasinormalStop",
                    "CauchySafeguard");
            case DoglegTrustRegion:
                return Matlab::enumToMxArray("QuasinormalStop",
                    "DoglegTrustRegion");
            case DoglegSafeguard:
                return Matlab::enumToMxArray("QuasinormalStop",
                    "DoglegSafeguard");
            case NewtonTrustRegion:
                return Matlab::enumToMxArray("QuasinormalStop",
                    "NewtonTrustRegion");
            case NewtonSafeguard:
                return Matlab::enumToMxArray("QuasinormalStop",
                    "NewtonSafeguard");
            case Feasible:
                return Matlab::enumToMxArray("QuasinormalStop",
                    "Feasible");
            case CauchySolved:
                return Matlab::enumToMxArray("QuasinormalStop",
                    "CauchySolved");
            case LocalMin:
                return Matlab::enumToMxArray("QuasinormalStop",
                    "LocalMin");
            case NewtonFailed:
                return Matlab::enumToMxArray("QuasinormalStop",
                    "NewtonFailed");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = *mxGetPr(member.get());

            if(m==Matlab::enumToNatural("QuasinormalStop",
                "Newton")
            )
                return Newton;
            else if(m==Matlab::enumToNatural("QuasinormalStop",
                "CauchyTrustRegion")
            )
                return CauchyTrustRegion;
            else if(m==Matlab::enumToNatural("QuasinormalStop",
                "CauchySafeguard")
            )
                return CauchySafeguard;
            else if(m==Matlab::enumToNatural("QuasinormalStop",
                "DoglegTrustRegion")
            )
                return DoglegTrustRegion;
            else if(m==Matlab::enumToNatural("QuasinormalStop",
                "DoglegSafeguard")
            )
                return DoglegSafeguard;
            else if(m==Matlab::enumToNatural("QuasinormalStop",
                "NewtonTrustRegion")
            )
                return NewtonTrustRegion;
            else if(m==Matlab::enumToNatural("QuasinormalStop",
                "NewtonSafeguard")
            )
                return NewtonSafeguard;
            else if(m==Matlab::enumToNatural("QuasinormalStop",
                "Feasible")
            )
                return Feasible;
            else if(m==Matlab::enumToNatural("QuasinormalStop",
                "CauchySolved")
            )
                return CauchySolved;
            else if(m==Matlab::enumToNatural("QuasinormalStop",
                "LocalMin")
            )
                return LocalMin;
            else if(m==Matlab::enumToNatural("QuasinormalStop",
                "NewtonFailed")
            )
                return NewtonFailed;
            else
                throw Matlab::Exception::t( __LOC__
                    + ", unknown QuasiNormalStep");
        }
    }

    namespace json {
        // Serialization utility for the Rm vector space
        template <>
        struct Serialization <double,Matlab::MatlabVS> {
            static std::string serialize (
                Matlab::Vector const & x,
                std::string const & name_,
                Natural const & iter_
            ) {
                // Convert the name and iteration to MATLAB
                auto name = Matlab::mxManaged(mxCreateString(name_.c_str()));
                auto iter = Matlab::mxArray_FromSize_t(iter_);

                // Call the serialize routine on the vector
                mxArray * input[3]={
                    x.data->get(),
                    name.get(),
                    iter.get()};
                mxArray * output[1] = {nullptr};
                auto err=mexCallMATLAB(1,output,3,input,"serialize");
                auto x_json = Matlab::mxManaged(output[0]);

                // Check errors
                if(err)
                    throw Matlab::Exception::t( __LOC__
                        + ", evaluation of the serialize function failed");

                // Convert the serialized vector to a string and return it 
                return std::string(mxArrayToString(x_json.get()));
            }

            static Matlab::Vector deserialize (
                Matlab::Vector const & x_,
                std::string const & x_json_
            ) {
                // Convert the inputed string into Matlab
                auto x_json =Matlab::mxManaged(mxCreateString(x_json_.c_str()));

                // Allocate memory for a new Matlab vector
                auto x = x_.init();

                // Call the deserialize routine on the reference vector and the
                // json vector
                mxArray * input[2] = {x_.data->get(),x_json.get()};
                mxArray * output[1] = {nullptr};
                auto err=mexCallMATLAB(1,output,2,input,"deserialize");
                x.data.reset(new Matlab::mxManaged(output[0]));

                // Check errors
                if(err)
                    throw Matlab::Exception::t( __LOC__
                        + ", evaluation of the deserialize function failed");

                // Move out the new vector
                return x;
            }
        };
    }

    namespace Matlab {
        // Grab the pointer
        mxArrayPtr::mxArrayPtr(
            mxArray const * const & ptr_
        ) : ptr(const_cast <mxArray *> (ptr_)) {}

        // Allow move semantics
        mxArrayPtr::mxArrayPtr(mxArrayPtr && p) : ptr(p.ptr) {
            p.ptr = nullptr;
        }
        mxArrayPtr & mxArrayPtr::operator = (mxArrayPtr && p) {
            ptr = p.ptr;
            p.ptr = nullptr;
            return *this;
        }

        // Allow the deallocator to be overwritten 
        mxArrayPtr::~mxArrayPtr() {} 

        // Grab the internal pointer
        mxArray * mxArrayPtr::get() const {
            return ptr;
        }

        // Release ownership and management of the pointer
        mxArray * mxArrayPtr::release() {
            // Grab a copy of the pointer
            auto ptr_ = ptr;

            // Set the internal pointer to null
            ptr = nullptr;

            // Return the pointer
            return ptr_;
        }

        // Grab the pointer
        mxUnmanaged::mxUnmanaged(
            mxArray const * const & ptr
        ) : mxArrayPtr(ptr) {}

        // Allow copy semantics
        mxUnmanaged::mxUnmanaged(mxUnmanaged const & p) : mxArrayPtr(p.get()) {}
        mxUnmanaged & mxUnmanaged::operator = (mxUnmanaged const & p) {
            ptr = p.ptr;
        }

        // Grab the pointer
        mxManaged::mxManaged(
            mxArray const * const & ptr
        ) : mxArrayPtr(ptr) {}

        // Free the memory
        mxManaged::~mxManaged() {
            mxDestroyArray(ptr); 
        }

        // Reset the pointer to a given value
        void mxManaged::reset(mxArray * ptr_) {
            // Free the internal memory if necessary
            if(ptr)
                mxDestroyArray(ptr);

            // Grab the pointer
            ptr = ptr_;
        }

        // Calls a Matlab function with one argument and no returns 
        int mxArray_CallObject1_0(
            mxArray * const fn,
            mxArray * const arg1
        ) {
            mxArray * input[2]={fn,arg1};
            mxArray * output[0];
            auto err=mexCallMATLAB(0,output,2,input,"feval");
            return err;
        }

        // Calls a Matlab function with one argument 
        std::tuple <mxManaged,int> mxArray_CallObject1(
            mxArray * const fn,
            mxArray * const arg1
        ) {
            mxArray * input[2]={fn,arg1};
            mxArray * output[1] = {nullptr};
            auto err=mexCallMATLAB(1,output,2,input,"feval");
            return std::make_tuple(output[0],err);
        }
        
        // Calls a Matlab function with two arguments
        std::tuple <mxManaged,int> mxArray_CallObject2(
            mxArray * const fn,
            mxArray * const arg1,
            mxArray * const arg2
        ) {
            mxArray* input[3]={fn,arg1,arg2};
            mxArray* output[1] = {nullptr};
            auto err=mexCallMATLAB(1,output,3,input,"feval");
            return std::make_tuple(output[0],err);
        }
        
        // Calls a Matlab function with three arguments
        std::tuple <mxManaged,int> mxArray_CallObject3(
            mxArray * const fn,
            mxArray * const arg1,
            mxArray * const arg2,
            mxArray * const arg3
        ) {
            mxArray* input[4]={fn,arg1,arg2,arg3};
            mxArray* output[1] = {nullptr};
            auto err=mexCallMATLAB(1,output,4,input,"feval");
            return std::make_tuple(output[0],err);
        }
        
        // Calls a Matlab function with four arguments
        std::tuple <mxManaged,int> mxArray_CallObject4(
            mxArray * const fn,
            mxArray * const arg1,
            mxArray * const arg2,
            mxArray * const arg3,
            mxArray * const arg4
        ) {
            mxArray* input[5]={fn,arg1,arg2,arg3,arg4};
            mxArray* output[1] = {nullptr};
            auto err=mexCallMATLAB(1,output,5,input,"feval");
            return std::make_tuple(output[0],err);
        }

        // Creates a Matlab double from a C++ double
        mxManaged mxArray_FromDouble(double const x_) {
            auto x = mxManaged(mxCreateDoubleMatrix(1,1,mxREAL));
            mxGetPr(x.get())[0]=x_;
            return x;
        }

        // Creates a Matlab int from a C++ size_t 
        mxManaged mxArray_FromSize_t(Natural const x_) {
            auto x = mxManaged(mxCreateDoubleMatrix(1,1,mxREAL));
            mxGetPr(x.get())[0]=x_;
            return x;
        }

        // Converts an Optizelle enumerated type to a mxArray *
        mxManaged enumToMxArray(
            std::string const & type_,
            std::string const & member_
        ) {
            // Grab the type 
            mxArray * type = mxGetField(optizelle.back(),0,type_.c_str());

            // Grab the member
            mxArray * member = mxGetField(type,0,member_.c_str());
                
            // Return the member 
            return mxArray_FromSize_t(fromDouble(mxGetPr(member)[0]));; 
        }
       
        // Converts an Optizelle enumerated type to a Natural
        Natural enumToNatural(
            std::string const & type,
            std::string const & member 
        ) {
            // Grab the mxArray * for the type and member requested
            auto obj = enumToMxArray(type,member);

            // Convert and return the member
            return Natural(*mxGetPr(obj.get()));
        }
        
        // Converts a MATLAB double to an Optizelle Natural
        Natural fromDouble(double value) {
            // Check that we're in bounds.  We were hitting a precision
            // problem on 64-bit MATLAB/Octave that made this necessary.
            if(value >= double(std::numeric_limits<Optizelle::Natural>::max())
            )
                return std::numeric_limits <Optizelle::Natural>::max();
            else if(
                value <= double(std::numeric_limits <Optizelle::Natural>::min())
            )
                return std::numeric_limits <Optizelle::Natural>::min();
            else
                return Optizelle::Natural(value);
        }

        // Memory allocation and size setting 
        Vector Vector::init() const {
            // Call the init function on the internal and store in y 
            auto init = mxGetField(vs->get(),0,"init");
            auto ret_err(mxArray_CallObject1(
                init,
                data->get()));

            // Check errors
            if(std::get<1>(ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the vector space function init failed");

            // Return the vector
            return Vector(vs,std::move(std::get<0> (ret_err)));
        } 
        
        // y <- x (Shallow.  No memory allocation.)  Internal is y.
        void Vector::copy(Vector const & x) { 
            // Call the copy function on x and the internal 
            auto copy = mxGetField(vs->get(),0,"copy");
            auto ret_err = mxArray_CallObject1(copy,x.data->get());

            // Check errors
            if(std::get<1>(ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the vector space function copy failed");

            // Assign y
            data.reset(new mxManaged(std::get<0>(ret_err).release()));
        } 

        // x <- alpha * x.  Internal is x.
        void Vector::scal(double const & alpha_) { 
            // Call the scal function on alpha and the internal storage 
            auto scal = mxGetField(vs->get(),0,"scal");
            auto alpha = mxArray_FromDouble(alpha_);
            auto ret_err = mxArray_CallObject2(
                scal,
                alpha.get(),
                data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the vector space function scal failed");

            // Assign x
            data.reset(new mxManaged(std::get<0>(ret_err).release()));
        } 

        // x <- 0.  Internal is x. 
        void Vector::zero() { 
            // Call the zero function on this vector.
            auto zero = mxGetField(vs->get(),0,"zero");
            auto ret_err = mxArray_CallObject1(
                zero,
                data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the vector space function zero failed");

            // Assign x
            data.reset(new mxManaged(std::get<0>(ret_err).release()));
        } 

        // y <- alpha * x + y.   Internal is y.
        void Vector::axpy(double const & alpha_,Vector const & x) { 
            // Call the axpy function on alpha, x, and the internal storage.
            auto axpy = mxGetField(vs->get(),0,"axpy");
            auto alpha = mxArray_FromDouble(alpha_);
            auto ret_err = mxArray_CallObject3(
                axpy,
                alpha.get(),
                x.data->get(),
                data->get());
           
            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the vector space function axpy failed");

            // Assign y
            data.reset(new mxManaged(std::get<0>(ret_err).release()));
        } 

        // innr <- <x,y>.  Internal is y.
        double Vector::innr(Vector const & x) const {
            // Call the innr function on x and the internal.  Store in z. 
            auto innr = mxGetField(vs->get(),0,"innr");
            auto ret_err = mxArray_CallObject2(
                innr,
                x.data->get(),
                data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the vector space function innr failed.");

            // Return the result 
            return mxGetPr(std::get <0> (ret_err).get())[0]; 
        } 

        // x <- random.  Internal is x. 
        void Vector::rand() { 
            // Call the rand function on this vector.
            auto rand = mxGetField(vs->get(),0,"rand");
            auto ret_err = mxArray_CallObject1(rand,data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the vector space function rand failed");
            
            // Assign x
            data.reset(new mxManaged(std::get<0>(ret_err).release()));
        }

        // Jordan product, z <- x o y.  Internal is z.
        void Vector::prod(Vector const & x,Vector const & y) { 
            // Call the prod function on x, y, and the internal 
            auto prod = mxGetField(vs->get(),0,"prod");
            auto ret_err= mxArray_CallObject2(prod,x.data->get(),y.data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the vector space function prod failed");
            
            // Assign z
            data.reset(new mxManaged(std::get<0>(ret_err).release()));
        } 

        // Identity element, x <- e such that x o e = x .  Internal is x.
        void Vector::id() { 
            // Call the id function on the internal.
            auto id = mxGetField(vs->get(),0,"id");
            auto ret_err = mxArray_CallObject1(id,data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the vector space function id failed");
            
            // Assign x
            data.reset(new mxManaged(std::get<0>(ret_err).release()));
        } 

        // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y.
        // Internal is z.
        void Vector::linv(Vector const & x, Vector const & y) { 
            // Call the linv function on x, y, and the internal
            auto linv = mxGetField(vs->get(),0,"linv");
            auto ret_err= mxArray_CallObject2(linv,x.data->get(),y.data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the vector space function linv failed.");
            
            // Assign z
            data.reset(new mxManaged(std::get<0>(ret_err).release()));
        } 

        // Barrier function, barr <- barr(x) where x o grad barr(x) = e.
        // Internal is x.
        double Vector::barr() const {
            // Call the barr function on the internal.  Store in z.
            auto barr = mxGetField(vs->get(),0,"barr");
            auto ret_err = mxArray_CallObject1(barr,data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__ +
                    ", evaluation of the vector space function barr failed.");

            // Return the result 
            return mxGetPr(std::get <0> (ret_err).get())[0]; 
        } 

        // Line search, srch <- argmax {alpha in Real >= 0 : alpha x + y >= 0} 
        // where y > 0.  Internal is y.
        double Vector::srch(Vector const & x) const {
            // Call the srch function on x and the internal.  Store in z.
            auto srch = mxGetField(vs->get(),0,"srch");
            auto ret_err = mxArray_CallObject2(srch,x.data->get(),data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the vector space function srch failed");

            // Return the result 
            return mxGetPr(std::get <0> (ret_err).get())[0]; 
        } 

        // Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
        // operator.  Internal is x.
        void Vector::symm() { 
            // Call the symm function on the internal.
            auto symm = mxGetField(vs->get(),0,"symm");
            auto ret_err = mxArray_CallObject1(symm,data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the vector space function symm failed");
            
            // Assign x
            data.reset(new mxManaged(std::get<0>(ret_err).release()));
        } 
        
        // Converts (copies) a value into Matlab.  
        mxManaged Vector::toMatlab() const {
            // Call the copy function on the internal and x
            auto copy = mxGetField(vs->get(),0,"copy");
            auto ret_err = mxArray_CallObject1(
                copy,
                data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the vector space function copy failed");
            
            // Return the pointer 
            return std::get <0> (std::move(ret_err));
        } 
        
        // Converts (copies) a value from Matlab.  This assumes that the
        // vector space functions have already been properly assigned.
        void Vector::fromMatlab(mxArrayPtr const & ptr) {
            // Call the copy function on ptr and the internal 
            auto copy = mxGetField(vs->get(),0,"copy");
            auto ret_err = mxArray_CallObject1(copy,ptr.get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the vector space function copy failed");

            // Copy in the value
            data.reset(new mxManaged(std::get<0>(ret_err).release()));
        } 

        // <- f(x)-> 
        double ScalarValuedFunction::eval(Vector const & x) const { 
            // Call the objective function on x
            auto eval = mxGetField(data->get(),0,"eval");
            auto ret_err = mxArray_CallObject1(eval,x.data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the objective f failed");

            // Return the result 
            return mxGetPr(std::get <0> (ret_err).get())[0]; 
        }

        // grad = grad f(x) 
        void ScalarValuedFunction::grad(
            Vector const & x,
            Vector & grad
        ) const { 
            // Call the gradient function on x
            auto mxgrad = mxGetField(data->get(),0,"grad");
            auto ret_err = mxArray_CallObject1(mxgrad,x.data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the gradient of f failed.");

            // Assign grad 
            grad.data.reset(new mxManaged(std::move(std::get<0> (ret_err))));
        }

        // H_dx = hess f(x) dx 
        void ScalarValuedFunction::hessvec(
            Vector const & x,
            Vector const & dx,
            Vector & H_dx
        ) const {
            // Call the hessvec function on x and dx,
            auto hessvec = mxGetField(data->get(),0,"hessvec");
            auto ret_err = mxArray_CallObject2(
                hessvec,x.data->get(),dx.data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the Hessian-vector product of f failed");
            
            // Assign H_dx
            H_dx.data.reset(new mxManaged(std::move(std::get<0> (ret_err))));
        }

        // y=f(x)
        void VectorValuedFunction::eval(
            X_Vector const & x,
            Y_Vector & y
        ) const {
            // Call the objective function on x.
            auto eval = mxGetField(data->get(),0,"eval");
            auto ret_err = mxArray_CallObject1(eval,x.data->get());

            // Check errors
            if(std::get<1>(ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the constraint " + name + " failed");
            
            // Assign y 
            y.data.reset(new mxManaged(std::move(std::get<0> (ret_err))));
        }

        // y=f'(x)dx 
        void VectorValuedFunction::p(
            X_Vector const & x,
            X_Vector const & dx,
            Y_Vector& y
        ) const {
            // Call the prime function on x and dx
            auto p = mxGetField(data->get(),0,"p");
            auto ret_err = mxArray_CallObject2(p,x.data->get(),dx.data->get());
           
            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the derivative of the constraint "
                    + name + " failed");
            
            // Assign y 
            y.data.reset(new mxManaged(std::move(std::get<0> (ret_err))));
        }

        // xhat=f'(x)*dy
        void VectorValuedFunction::ps(
            X_Vector const & x,
            Y_Vector const & dy,
            X_Vector & xhat 
        ) const {
            // Call the prime-adjoint function on x and dy
            auto ps = mxGetField(data->get(),0,"ps");
            auto ret_err = mxArray_CallObject2(ps,x.data->get(),dy.data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    +", evaluation of the derivative-adjoint of the constraint "
                    + name + " failed");
            
            // Assign z
            xhat.data.reset(new mxManaged(std::move(std::get<0> (ret_err))));
        }
             
        // xhat=(f''(x)dx)*dy
        void VectorValuedFunction::pps(
            X_Vector const & x,
            X_Vector const & dx,
            Y_Vector const & dy,
            X_Vector & xhat 
        ) const { 
            // Call the prime-adjoint function on x, dx, and dy
            auto pps = mxGetField(data->get(),0,"pps");
            auto ret_err = mxArray_CallObject3(
                pps,x.data->get(),dx.data->get(),dy.data->get());

            // Check errors
            if(std::get <1> (ret_err))
                throw Matlab::Exception::t( __LOC__
                    + ", evaluation of the second derivative-adjoint of the "
                    + "constraint " + name + " failed");
            
            // Assign xhat 
            xhat.data.reset(new mxManaged(std::move(std::get<0> (ret_err))));
        }

        // Converts elements from C++ to Matlab 
        namespace toMatlab {
                    
            // Sets a real in a Matlab state 
            void Real(
                std::string const & name,
                double const & value,
                mxArrayPtr & mxstate 
            ) {
                // Get the field
                mxArray * field = mxGetField(mxstate.get(),0,name.c_str());

                // If it's empty, allocate new memory
                if(!field) {
                    auto item = mxArray_FromDouble(value);
                    mxSetField(mxstate.get(),0,name.c_str(),item.release());

                // If it's not empty, just modify the element
                } else
                    mxGetPr(field)[0]=value;
            }
        
            // Sets a natural in a Matlab state 
            void Natural(
                std::string const & name,
                Optizelle::Natural const & value,
                mxArrayPtr & mxstate 
            ) {
                // Get the field
                auto field = mxGetField(mxstate.get(),0,name.c_str());

                // If it's empty, allocate new memory
                if(!field) {
                    auto item = mxArray_FromSize_t(value);
                    mxSetField(mxstate.get(),0,name.c_str(),item.release());

                // If it's not empty, just modify the element
                } else
                    mxGetPr(field)[0]=value;
            }
        
            // Sets a vector in a Matlab state 
            void Vector(
                std::string const & name,
                Matlab::Vector const & value,
                mxArrayPtr & mxstate 
            ) {
                // Get the field
                auto field = mxGetField(mxstate.get(),0,name.c_str());

                // If it's not empty, we need to free the memory 
                if(field)
                    mxDestroyArray(field);

                // Now, set the element
                mxSetField(mxstate.get(),0,name.c_str(),
                    value.toMatlab().release());
            }
        
            // Sets a list of vectors in a Matlab state 
            void VectorList(
                std::string const & name,
                std::list <Matlab::Vector> const & vectors,
                mxArrayPtr & mxstate 
            ) {
                // Get the field
                auto field = mxGetField(mxstate.get(),0,name.c_str());

                // If it's not empty, we need to free the memory 
                if(field)
                    mxDestroyArray(field);

                // Create a new Matlab cell array that we insert elements into
                auto mxvectors = mxCreateCellMatrix(1,vectors.size());

                // Loop over all of the items inside vectors and then insert 
                // them into items 
                auto i = Optizelle::Natural(0);
                for(auto const & vector : vectors) {
                    // Allocate memory for a new vector
                    auto v = vector.init();

                    // Copy the information from the current iterator into this
                    // new vector
                    v.copy(vector);

                    // Release the pointer into the Matlab cell array
                    mxSetCell(mxvectors,i,v.data->release());

                    // Increment our counter
                    i++;
                }
                
                // Insert the items into mxstate 
                mxSetField(mxstate.get(),0,name.c_str(),mxvectors);
            }
        
            // Sets restart vectors in Matlab 
            void Vectors(
                Matlab::Vectors const & values,
                mxArrayPtr & mxvalues 
            ) {
            
                // Loop over all of the items inside values and then insert 
                // them into mxvalues 
                auto i = Optizelle::Natural(0);
                for(auto const & value : values) {
                    // Allocate memory for a new vector
                    auto mxvalue = value.second.init();

                    // Copy the information from the current iterator into this
                    // new vector
                    mxvalue.copy(value.second);

                    // Create a 2-element cell array with the name and value
                    auto tuple = mxCreateCellMatrix(1,2);
                    mxSetCell(tuple,0,
                        mxCreateString(value.first.c_str()));
                    mxSetCell(tuple,1,mxvalue.data->release());

                    // Release the tuple into the Matlab cell array
                    mxSetCell(mxvalues.get(),i,tuple);

                    // Increment our counter
                    i++;
                }
            }
        
            // Sets restart reals in Matlab 
            void Reals(
                Matlab::Reals const & values,
                mxArrayPtr & mxvalues 
            ) {
                // Loop over all of the items inside values and then insert 
                // them into mxvalues 
                auto i = Optizelle::Natural(0);
                for(auto const & value : values) {
                    // Create a 2-element cell array with the name and value
                    auto tuple = mxCreateCellMatrix(1,2);
                    mxSetCell(tuple,0,
                        mxCreateString(value.first.c_str()));
                    mxSetCell(tuple,1,
                        mxArray_FromDouble(value.second).release());

                    // Release the tuple into the Matlab cell array
                    mxSetCell(mxvalues.get(),i,tuple);

                    // Increment our counter
                    i++;
                }
            }
        
            // Converts a list of naturals to a Matlab list 
            void Naturals(
                Matlab::Naturals const & values,
                mxArrayPtr & mxvalues 
            ) {
                // Loop over all of the items inside values and then insert 
                // them into mxvalues 
                auto i = Optizelle::Natural(0);
                for(auto const & value : values) {
                    // Create a 2-element cell array with the name and value
                    auto tuple = mxCreateCellMatrix(1,2);
                    mxSetCell(tuple,0,
                        mxCreateString(value.first.c_str()));
                    mxSetCell(tuple,1,
                        mxArray_FromSize_t(value.second).release());

                    // Release the tuple into the Matlab cell array
                    mxSetCell(mxvalues.get(),i,tuple);

                    // Increment our counter
                    i++;
                }
            }
        
            // Sets restart parameters in Matlab 
            void Params(
                Matlab::Params const & values,
                mxArrayPtr & mxvalues 
            ) {
                // Loop over all of the items inside values and then insert 
                // them into mxvalues 
                auto i = Optizelle::Natural(0);
                for(auto const & value : values) {
                    // Create a 2-element cell array with the name and value
                    auto tuple = mxCreateCellMatrix(1,2);
                    mxSetCell(tuple,0,
                        mxCreateString(value.first.c_str()));
                    mxSetCell(tuple,1,
                        mxCreateString(value.second.c_str()));

                    // Release the tuple into the Matlab cell array
                    mxSetCell(mxvalues.get(),i,tuple);

                    // Increment our counter
                    i++;
                }
            }
        }
        
        // Converts elements from Matlab to C++ 
        namespace fromMatlab {
        
            // Sets a real in a C++ state 
            void Real(
                std::string const & name,
                mxArrayPtr const & mxstate,
                double & value
            ) {
                auto item = mxGetField(mxstate.get(),0,name.c_str());
                value=mxGetPr(item)[0];
            }
            
            // Sets a natural in a C++ state 
            void Natural(
                std::string const & name,
                mxArrayPtr const & obj,
                Optizelle::Natural & value
            ) {
                auto item = mxGetField(obj.get(),0,name.c_str());
                value = fromDouble(mxGetPr(item)[0]);
            }
            
            // Sets a vector in a C++ state 
            void Vector(
                std::string const & name,
                mxArrayPtr const & obj,
                Matlab::Vector & value
            ) {
                auto item = mxGetField(obj.get(),0,name.c_str());
                value.fromMatlab(mxUnmanaged(item));
            }
            
            // Sets a list of vectors in a C++ state 
            void VectorList(
                std::string const & name,
                mxArrayPtr const & mxstate,
                Matlab::Vector const & vec,
                std::list <Matlab::Vector> & values
            ) {
                // Grab the list of items
                auto items = mxGetField(mxstate.get(),0,name.c_str());

                // Loop over all the elements in items and insert them one
                // at a time into values
                values.clear();
                for(auto i=0;i<mxGetN(items);i++) {
                    // Grab the current item from Matlab
                    auto item = mxGetCell(items,i);

                    // Create a new vector in values 
                    values.emplace_back(vec.init());

                    // Copy the Matlab item into the new value
                    values.back().fromMatlab(mxUnmanaged(item));
                }
            }
            
            // Sets a scalar-valued function in a C++ function bundle 
            void ScalarValuedFunction(
                std::string const & name,
                mxArrayPtr const & fns,
                std::unique_ptr <MxScalarValuedFunction> & value
            ) {
                value.reset(new Matlab::ScalarValuedFunction(
                    mxUnmanaged(mxGetField(fns.get(),0,name.c_str()))));
            }
        
            // Sets restart vectors in C++ 
            void Vectors(
                Matlab::Vector const & vec,
                mxArrayPtr const & mxvalues,
                Matlab::Vectors & values
            ) {
                // Loop over all the elements in mxvalues and insert them one
                // at a time into values
                values.clear();
                for(auto i=0;i<mxGetN(mxvalues.get());i++) {
                    // Grab the current item from Matlab
                    auto mxvalue = mxGetCell(mxvalues.get(),i);

                    // Create the elements in values 
                    values.emplace_back(
                        mxArrayToString(mxGetCell(mxvalue,0)),
                        vec.init());

                    // Copy the Matlab value into the C++ value
                    values.back().second.fromMatlab(
                        mxUnmanaged(mxGetCell(mxvalue,1)));
                }
            }
            
            // Sets restart reals in C++ 
            void Reals(
                mxArrayPtr const & mxvalues,
                Matlab::Reals & values
            ) {
                // Loop over all the elements in mxvalues and insert them one
                // at a time into values
                values.clear();
                for(auto i=0;i<mxGetN(mxvalues.get());i++) {
                    // Grab the current item from Matlab
                    auto mxvalue = mxGetCell(mxvalues.get(),i);
                    
                    // Create the elements in values 
                    values.emplace_back(
                        mxArrayToString(mxGetCell(mxvalue,0)),
                        mxGetPr(mxGetCell(mxvalue,1))[0]);
                }
            }
            
            // Sets restart naturals in C++ 
            void Naturals(
                mxArrayPtr const & mxvalues,
                Matlab::Naturals & values
            ) {
                // Loop over all the elements in mxvalues and insert them one
                // at a time into values
                values.clear();
                for(auto i=0;i<mxGetN(mxvalues.get());i++) {
                    // Grab the current item from Matlab
                    auto mxvalue = mxGetCell(mxvalues.get(),i);
                    
                    // Create the elements in values 
                    values.emplace_back(
                        mxArrayToString(mxGetCell(mxvalue,0)),
                        fromDouble(mxGetPr(mxGetCell(mxvalue,1))[0]));
                }
            }
            
            // Sets restart parameters in C++ 
            void Params(
                mxArrayPtr const & mxvalues,
                Matlab::Params & values
            ) {
                // Loop over all the elements in mxvalues and insert them one
                // at a time into values
                values.clear();
                for(auto i=0;i<mxGetN(mxvalues.get());i++) {
                    // Grab the current item from Matlab
                    auto mxvalue = mxGetCell(mxvalues.get(),i);
                    
                    // Create the elements in values 
                    values.emplace_back(
                        mxArrayToString(mxGetCell(mxvalue,0)),
                        mxArrayToString(mxGetCell(mxvalue,1)));
                }
            }
        }

        // Convert a C++ state to a Matlab state 
        template <>
        void State <MxUnconstrained>::toMatlab(
            typename MxUnconstrained::State::t const & state
        ) {
            Unconstrained::State::toMatlab(state,*data);
        }
        template <>
        void State <MxEqualityConstrained>::toMatlab(
            typename MxEqualityConstrained::State::t const & state
        ) {
            EqualityConstrained::State::toMatlab(state,*data);
        }
        template <>
        void State <MxInequalityConstrained>::toMatlab(
            typename MxInequalityConstrained::State::t const & state
        ) {
            InequalityConstrained::State::toMatlab(state,*data);
        }
        template <>
        void State <MxConstrained>::toMatlab(
            typename MxConstrained::State::t const & state
        ) {
            Constrained::State::toMatlab(state,*data);
        }

        // Convert a Matlab state to C++ 
        template <>
        void State <MxUnconstrained>::fromMatlab(
            typename MxUnconstrained::State::t & state
        ) {
            Unconstrained::State::fromMatlab(*data,state);
        }
        template <>
        void State <MxEqualityConstrained>::fromMatlab(
            typename MxEqualityConstrained::State::t & state
        ) {
            EqualityConstrained::State::fromMatlab(*data,state);
        }
        template <>
        void State <MxInequalityConstrained>::fromMatlab(
            typename MxInequalityConstrained::State::t & state
        ) {
            InequalityConstrained::State::fromMatlab(*data,state);
        }
        template <>
        void State <MxConstrained>::fromMatlab(
            typename MxConstrained::State::t & state
        ) {
            Constrained::State::fromMatlab(*data,state);
        }
        
        // Convert a Matlab bundle to C++ 
        template <>
        void Functions <MxUnconstrained>::fromMatlab(
            typename MxUnconstrained::Functions::t & fns 
        ) {
            Unconstrained::Functions::fromMatlab(
                *this,mxstate,state,fns);
        }
        template <>
        void Functions <MxEqualityConstrained>::fromMatlab(
            typename MxEqualityConstrained::Functions::t & fns 
        ) {
            EqualityConstrained::Functions::fromMatlab(
                *this,mxstate,state,fns);
        }
        template <>
        void Functions <MxInequalityConstrained>::fromMatlab(
            typename MxInequalityConstrained::Functions::t & fns 
        ) {
            InequalityConstrained::Functions::fromMatlab(
                *this,mxstate,state,fns);
        }
        template <>
        void Functions <MxConstrained>::fromMatlab(
            typename MxConstrained::Functions::t & fns 
        ) {
            Constrained::Functions::fromMatlab(
                *this,mxstate,state,fns);
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
                // Returns the fields names of the state
                std::vector <char const *> fieldNames_() {
                    std::vector <const char*> names = {
                        "eps_grad",
                        "eps_dx",
                        "stored_history",
                        "iter",
                        "iter_max",
                        "glob_iter",
                        "glob_iter_max",
                        "glob_iter_total",
                        "opt_stop",
                        "trunc_iter",
                        "trunc_iter_max",
                        "trunc_iter_total",
                        "trunc_orthog_storage_max",
                        "trunc_orthog_iter_max",
                        "trunc_stop",
                        "trunc_err",
                        "eps_trunc",
                        "algorithm_class",
                        "PH_type",
                        "H_type",
                        "norm_gradtyp",
                        "norm_dxtyp",
                        "x",
                        "grad",
                        "dx",
                        "x_old",
                        "grad_old",
                        "dx_old",
                        "oldY",
                        "oldS",
                        "f_x",
                        "f_xpdx",
                        "msg_level",
                        "safeguard_failed_max",
                        "safeguard_failed",
                        "safeguard_failed_total",
                        "alpha_x",
                        "alpha_x_qn",
                        "delta",
                        "eta1",
                        "eta2",
                        "ared",
                        "pred",
                        "alpha0",
                        "alpha",
                        "c1",
                        "ls_iter",
                        "ls_iter_max",
                        "ls_iter_total",
                        "eps_ls",
                        "dir",
                        "kind",
                        "f_diag",
                        "L_diag",
                        "x_diag",
                        "dscheme",
                        "eps_kind"};

                    return std::move(names);
                }
                std::vector <char const *> fieldNames() {
                    return fieldNames_();
                }

                // Create the structure for a Matlab state
                mxManaged mxCreate() {
                    auto names = Unconstrained::State::fieldNames();
                    return mxCreateStructMatrix(1,1,names.size(),&(names[0]));
                }

                // Convert a C++ state to a Matlab state 
                void toMatlab_(
                    typename MxUnconstrained::State::t const & state,
                    mxArrayPtr & mxstate
                ){
                    // Set each of the required items in the Matlab state
                    toMatlab::Real("eps_grad",state.eps_grad,mxstate);
                    toMatlab::Real("eps_dx",state.eps_dx,mxstate);
                    toMatlab::Natural("stored_history",
                        state.stored_history,mxstate);
                    toMatlab::Natural("iter",state.iter,mxstate);
                    toMatlab::Natural("iter_max",state.iter_max,mxstate);
                    toMatlab::Natural("glob_iter",state.glob_iter,mxstate);
                    toMatlab::Natural("glob_iter_max",
                        state.glob_iter_max,mxstate);
                    toMatlab::Natural("glob_iter_total",
                        state.glob_iter_total,mxstate);
                    toMatlab::Param <OptimizationStop::t> (
                        "opt_stop",
                        OptimizationStop::toMatlab,
                        state.opt_stop,
                        mxstate);
                    toMatlab::Natural("trunc_iter",state.trunc_iter,mxstate);
                    toMatlab::Natural("trunc_iter_max",
                        state.trunc_iter_max,mxstate);
                    toMatlab::Natural("trunc_iter_total",
                        state.trunc_iter_total,mxstate);
                    toMatlab::Natural("trunc_orthog_storage_max",
                        state.trunc_orthog_storage_max,mxstate);
                    toMatlab::Natural("trunc_orthog_iter_max",
                        state.trunc_orthog_iter_max,mxstate);
                    toMatlab::Param <TruncatedStop::t> (
                        "trunc_stop",
                        TruncatedStop::toMatlab,
                        state.trunc_stop,
                        mxstate);
                    toMatlab::Real("trunc_err",
                        state.trunc_err,mxstate);
                    toMatlab::Real("eps_trunc",state.eps_trunc,mxstate);
                    toMatlab::Param <AlgorithmClass::t> (
                        "algorithm_class",
                        AlgorithmClass::toMatlab,
                        state.algorithm_class,
                        mxstate);
                    toMatlab::Param <Operators::t> (
                        "PH_type",
                        Operators::toMatlab,
                        state.PH_type,
                        mxstate);
                    toMatlab::Param <Operators::t> (
                        "H_type",
                        Operators::toMatlab,
                        state.H_type,
                        mxstate);
                    toMatlab::Real("norm_gradtyp",state.norm_gradtyp,mxstate);
                    toMatlab::Real("norm_dxtyp",state.norm_dxtyp,mxstate);
                    toMatlab::Vector("x",state.x,mxstate);
                    toMatlab::Vector("grad",state.grad,mxstate);
                    toMatlab::Vector("dx",state.dx,mxstate);
                    toMatlab::Vector("x_old",state.x_old,mxstate);
                    toMatlab::Vector("grad_old",state.grad_old,mxstate);
                    toMatlab::Vector("dx_old",state.dx_old,mxstate);
                    toMatlab::VectorList("oldY",state.oldY,mxstate);
                    toMatlab::VectorList("oldS",state.oldS,mxstate);
                    toMatlab::Real("f_x",state.f_x,mxstate);
                    toMatlab::Real("f_xpdx",state.f_xpdx,mxstate);
                    toMatlab::Natural("msg_level",state.msg_level,mxstate);
                    toMatlab::Natural("safeguard_failed_max",
                        state.safeguard_failed_max,mxstate);
                    toMatlab::Natural("safeguard_failed",
                        state.safeguard_failed,mxstate);
                    toMatlab::Natural("safeguard_failed_total",
                        state.safeguard_failed_total,mxstate);
                    toMatlab::Real("alpha_x",state.alpha_x,mxstate);
                    toMatlab::Real("alpha_x_qn",state.alpha_x_qn,mxstate);
                    toMatlab::Real("delta",state.delta,mxstate);
                    toMatlab::Real("eta1",state.eta1,mxstate);
                    toMatlab::Real("eta2",state.eta2,mxstate);
                    toMatlab::Real("ared",state.ared,mxstate);
                    toMatlab::Real("pred",state.pred,mxstate);
                    toMatlab::Real("alpha0",state.alpha0,mxstate);
                    toMatlab::Real("alpha",state.alpha,mxstate);
                    toMatlab::Real("c1",state.c1,mxstate);
                    toMatlab::Natural("ls_iter",
                        state.ls_iter,mxstate);
                    toMatlab::Natural("ls_iter_max",
                        state.ls_iter_max,mxstate);
                    toMatlab::Natural("ls_iter_total",
                        state.ls_iter_total,mxstate);
                    toMatlab::Real("eps_ls",state.eps_ls,mxstate);
                    toMatlab::Param <LineSearchDirection::t> (
                        "dir",
                        LineSearchDirection::toMatlab,
                        state.dir,
                        mxstate);
                    toMatlab::Param <LineSearchKind::t> (
                        "kind",
                        LineSearchKind::toMatlab,
                        state.kind,
                        mxstate);
                    toMatlab::Param <FunctionDiagnostics::t> (
                        "f_diag",
                        FunctionDiagnostics::toMatlab,
                        state.f_diag,
                        mxstate);
                    toMatlab::Param <FunctionDiagnostics::t> (
                        "L_diag",
                        FunctionDiagnostics::toMatlab,
                        state.L_diag,
                        mxstate);
                    toMatlab::Param <VectorSpaceDiagnostics::t> (
                        "x_diag",
                        VectorSpaceDiagnostics::toMatlab,
                        state.x_diag,
                        mxstate);
                    toMatlab::Param <DiagnosticScheme::t> (
                        "dscheme",
                        DiagnosticScheme::toMatlab,
                        state.dscheme,
                        mxstate);
                    toMatlab::Param <ToleranceKind::t> (
                        "eps_kind",
                        ToleranceKind::toMatlab,
                        state.eps_kind,
                        mxstate);
                }
                void toMatlab(
                    typename MxUnconstrained::State::t const & state,
                    mxArrayPtr & mxstate
                ){
                    Unconstrained::State::toMatlab_(state,mxstate);
                }
                
                // Convert a Matlab state to C++ 
                void fromMatlab_(
                    mxArrayPtr const & mxstate,
                    typename MxUnconstrained::State::t & state
                ){
                    // Set each of the required items in the Matlab state
                    fromMatlab::Real("eps_grad",mxstate,state.eps_grad);
                    fromMatlab::Real("eps_dx",mxstate,state.eps_dx);
                    fromMatlab::Natural("stored_history",
                        mxstate,state.stored_history);
                    fromMatlab::Natural("iter",mxstate,state.iter);
                    fromMatlab::Natural("iter_max",mxstate,state.iter_max);
                    fromMatlab::Natural("glob_iter",
                        mxstate,state.glob_iter);
                    fromMatlab::Natural("glob_iter_max",
                        mxstate,state.glob_iter_max);
                    fromMatlab::Natural("glob_iter_total",
                        mxstate,state.glob_iter_total);
                    fromMatlab::Param <OptimizationStop::t> (
                        "opt_stop",
                        OptimizationStop::fromMatlab,
                        mxstate,
                        state.opt_stop);
                    fromMatlab::Natural("trunc_iter",
                        mxstate,state.trunc_iter);
                    fromMatlab::Natural("trunc_iter_max",
                        mxstate,state.trunc_iter_max);
                    fromMatlab::Natural("trunc_iter_total",
                        mxstate,state.trunc_iter_total);
                    fromMatlab::Natural("trunc_orthog_storage_max",
                        mxstate,state.trunc_orthog_storage_max);
                    fromMatlab::Natural("trunc_orthog_iter_max",
                        mxstate,state.trunc_orthog_iter_max);
                    fromMatlab::Param <TruncatedStop::t> (
                        "trunc_stop",
                        TruncatedStop::fromMatlab,
                        mxstate,
                        state.trunc_stop);
                    fromMatlab::Real("trunc_err",
                        mxstate,state.trunc_err);
                    fromMatlab::Real("eps_trunc",mxstate,state.eps_trunc);
                    fromMatlab::Param <AlgorithmClass::t> (
                        "algorithm_class",
                        AlgorithmClass::fromMatlab,
                        mxstate,
                        state.algorithm_class);
                    fromMatlab::Param <Operators::t> (
                        "PH_type",
                        Operators::fromMatlab,
                        mxstate,
                        state.PH_type);
                    fromMatlab::Param <Operators::t> (
                        "H_type",
                        Operators::fromMatlab,
                        mxstate,
                        state.H_type);
                    fromMatlab::Real("norm_gradtyp",
                        mxstate,state.norm_gradtyp);
                    fromMatlab::Real("norm_dxtyp",mxstate,state.norm_dxtyp);
                    fromMatlab::Vector("x",mxstate,state.x);
                    fromMatlab::Vector("grad",mxstate,state.grad);
                    fromMatlab::Vector("dx",mxstate,state.dx);
                    fromMatlab::Vector("x_old",mxstate,state.x_old);
                    fromMatlab::Vector("grad_old",mxstate,state.grad_old);
                    fromMatlab::Vector("dx_old",mxstate,state.dx_old);
                    fromMatlab::VectorList("oldY",mxstate,state.x,state.oldY);
                    fromMatlab::VectorList("oldS",mxstate,state.x,state.oldS);
                    fromMatlab::Real("f_x",mxstate,state.f_x);
                    fromMatlab::Real("f_xpdx",mxstate,state.f_xpdx);
                    fromMatlab::Natural("msg_level",mxstate,state.msg_level);
                    fromMatlab::Natural("safeguard_failed_max",
                        mxstate,state.safeguard_failed_max);
                    fromMatlab::Natural("safeguard_failed",
                        mxstate,state.safeguard_failed);
                    fromMatlab::Natural("safeguard_failed_total",
                        mxstate,state.safeguard_failed_total);
                    fromMatlab::Real("alpha_x",mxstate,state.alpha_x);
                    fromMatlab::Real("alpha_x_qn",mxstate,state.alpha_x_qn);
                    fromMatlab::Real("delta",mxstate,state.delta);
                    fromMatlab::Real("eta1",mxstate,state.eta1);
                    fromMatlab::Real("eta2",mxstate,state.eta2);
                    fromMatlab::Real("ared",mxstate,state.ared);
                    fromMatlab::Real("pred",mxstate,state.pred);
                    fromMatlab::Real("alpha0",mxstate,state.alpha0);
                    fromMatlab::Real("alpha",mxstate,state.alpha);
                    fromMatlab::Real("c1",mxstate,state.c1);
                    fromMatlab::Natural("ls_iter",
                        mxstate,state.ls_iter);
                    fromMatlab::Natural("ls_iter_max",
                        mxstate,state.ls_iter_max);
                    fromMatlab::Natural("ls_iter_total",mxstate,
                        state.ls_iter_total);
                    fromMatlab::Real("eps_ls",mxstate,state.eps_ls);
                    fromMatlab::Param <LineSearchDirection::t> (
                        "dir",
                        LineSearchDirection::fromMatlab,
                        mxstate,
                        state.dir);
                    fromMatlab::Param <LineSearchKind::t> (
                        "kind",
                        LineSearchKind::fromMatlab,
                        mxstate,
                        state.kind);
                    fromMatlab::Param <FunctionDiagnostics::t> (
                        "f_diag",
                        FunctionDiagnostics::fromMatlab,
                        mxstate,
                        state.f_diag);
                    fromMatlab::Param <FunctionDiagnostics::t> (
                        "L_diag",
                        FunctionDiagnostics::fromMatlab,
                        mxstate,
                        state.L_diag);
                    fromMatlab::Param <VectorSpaceDiagnostics::t> (
                        "x_diag",
                        VectorSpaceDiagnostics::fromMatlab,
                        mxstate,
                        state.x_diag);
                    fromMatlab::Param <DiagnosticScheme::t> (
                        "dscheme",
                        DiagnosticScheme::fromMatlab,
                        mxstate,
                        state.dscheme);
                    fromMatlab::Param <ToleranceKind::t> (
                        "eps_kind",
                        ToleranceKind::fromMatlab,
                        mxstate,
                        state.eps_kind);
                }
                void fromMatlab(
                    mxArrayPtr const & mxstate,
                    typename MxUnconstrained::State::t & state
                ){
                    Unconstrained::State::fromMatlab_(mxstate,state);
                }

                // Creates a state and inserts the elements into mxstate 
                void create(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be (X,x) -> (mxstate_out)
                    mxUnmanaged X=pInput[0],
                                x_=pInput[1];

                    // Create a vector from the user input
                    auto x = Vector(std::move(X),std::move(x_));

                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxUnconstrained> (
                        Unconstrained::State::mxCreate());

                    // Create a new C++ state
                    typename MxUnconstrained::State::t state(x);

                    // Convert the state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS; 
        
                // Read json parameters from file
                void readJson(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be (X,fname,mxstate)
                    // -> (mxstate_out)
                    mxUnmanaged X=pInput[0],
                                fname_=pInput[1],
                                mxstate_=pInput[2];

                    // Grab the file name
                    auto fname = mxArrayToString(fname_.get());

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxUnconstrained> (
                        std::move(mxstate_));
                    auto mxstate_out = Matlab::State <MxUnconstrained> (
                        Unconstrained::State::mxCreate());
                
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));

                    // Create a new C++ state
                    typename MxUnconstrained::State::t state(x);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Read the JSON file into the C++ state
                    MxJsonUnconstrained::read(fname,state);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();
                            
                    // Return nothing 
                    return; 
                } CATCH_MATLAB_ERRORS; 
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Matlab bundle to C++ 
                void fromMatlab(
                    Matlab::Functions <MxUnconstrained> const & mxfns,
                    Matlab::State <MxUnconstrained> & mxstate,
                    typename MxUnconstrained::State::t const & state,
                    typename MxUnconstrained::Functions::t & fns 
                ) {
                    Unconstrained::Functions::fromMatlab_
                        <MxUnconstrained> (mxfns,mxstate,state,fns);
                }
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem 
                void getMin(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    // (X,msg,mxfns,mxstate,smanip) -> (mxstate_out)
                    mxUnmanaged X=pInput[0],
                                msg_=pInput[1],
                                mxfns_=pInput[2],
                                mxstate_=pInput[3],
                                smanip_=pInput[4];

                    // Create a messaging object
                    auto msg = Optizelle::Matlab::Messaging::matlab(
                        std::move(msg_));

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxUnconstrained> (
                        std::move(mxstate_));
                    auto mxstate_out = Matlab::State <MxUnconstrained> (
                        Unconstrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));

                    // Create a C++ state
                    typename MxUnconstrained::State::t state(x);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Create a Matlab bundle of functions
                    Matlab::Functions <MxUnconstrained> mxfns(
                        mxstate_out,
                        state,
                        std::move(mxfns_));

                    // Create a C++ bundle of functions
                    typename MxUnconstrained::Functions::t fns;
                    
                    // Convert the Matlab bundle of functions to C++ 
                    mxfns.fromMatlab(fns);

                    // Create a state manipulator 
                    Matlab::StateManipulator <MxUnconstrained> smanip(
                        mxstate_out,
                        mxfns,
                        std::move(smanip_));
                   
                    // Minimize
                    MxUnconstrained::Algorithms::getMin(msg,fns,state,smanip);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();
                    
                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return; 
                } CATCH_MATLAB_ERRORS;
            }
            
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                void release(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    // (X,mxstate) -> (xs,reals,nats,params)
                    mxUnmanaged X=pInput[0],
                                mxstate_=pInput[1];

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxUnconstrained> (
                        std::move(mxstate_));
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));

                    // Create a C++ state
                    typename MxUnconstrained::State::t state(x);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Do a release 
                    MxUnconstrained::Restart::X_Vectors xs;
                    MxUnconstrained::Restart::Reals reals;
                    MxUnconstrained::Restart::Naturals nats;
                    MxUnconstrained::Restart::Params params;
                    MxUnconstrained::Restart
                        ::release(state,xs,reals,nats,params);

                    // Allocate memory for the Matlab versions
                    auto mxxs = mxManaged(mxCreateCellMatrix(1,xs.size()));
                    auto mxreals= mxManaged(mxCreateCellMatrix(1,reals.size()));
                    auto mxnats = mxManaged(mxCreateCellMatrix(1,nats.size()));
                    auto mxparams
                        = mxManaged(mxCreateCellMatrix(1,params.size()));

                    // Convert the restart information to Matlab 
                    toMatlab::Vectors(xs,mxxs);
                    toMatlab::Reals(reals,mxreals);
                    toMatlab::Naturals(nats,mxnats);
                    toMatlab::Params(params,mxparams);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the ouptuts
                    pOutput[0]=mxxs.release();
                    pOutput[1]=mxreals.release();
                    pOutput[2]=mxnats.release();
                    pOutput[3]=mxparams.release();

                    // Return nothing 
                    return; 
                
                } CATCH_MATLAB_ERRORS;

                // Capture data from structures controlled by the user.  
                void capture(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    // (X,mxstate,mxxs,mxreals,mxnats,mxparams)
                    // -> mxstate_out
                    mxUnmanaged X=pInput[0],
                                mxstate_=pInput[1],
                                mxxs=pInput[2],
                                mxreals=pInput[3],
                                mxnats=pInput[4],
                                mxparams=pInput[5];

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxUnconstrained> (
                        std::move(mxstate_));
                    auto mxstate_out = Matlab::State <MxUnconstrained> (
                        Unconstrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));

                    // Create a C++ state
                    typename MxUnconstrained::State::t state(x);
                   
                    // Allocate memory for the released vectors
                    MxUnconstrained::Restart::X_Vectors xs;
                    MxUnconstrained::Restart::Reals reals;
                    MxUnconstrained::Restart::Naturals nats;
                    MxUnconstrained::Restart::Params params;
                    
                    // Convert the restart information from Matlab 
                    fromMatlab::Vectors(x,mxxs,xs);
                    fromMatlab::Reals(mxreals,reals);
                    fromMatlab::Naturals(mxnats,nats);
                    fromMatlab::Params(mxparams,params);

                    // Do a capture 
                    MxUnconstrained::Restart
                        ::capture(state,xs,reals,nats,params);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS;
                
                // Writes a json restart file
                void write_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be (X,fname,state) -> ()
                    mxUnmanaged X=pInput[0],
                                fname_=pInput[1],
                                mxstate_=pInput[2];

                    // Grab the file name
                    auto fname = mxArrayToString(fname_.get());

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxUnconstrained> (
                        std::move(mxstate_));
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    
                    // Create a C++ state
                    typename MxUnconstrained::State::t state(x);
                    
                    // Convert Matlab state to C++ 
                    mxstate.fromMatlab(state);

                    // Write the restart file
                    MxJsonUnconstrained::write_restart(fname,state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();
                    
                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS;
                
                // Reads a json restart file
                void read_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be (X,fname,x) -> (mxstate)
                    mxUnmanaged X=pInput[0],
                                fname_=pInput[1],
                                x_=pInput[2];

                    // Grab the file name
                    auto fname = mxArrayToString(fname_.get());

                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxUnconstrained> (
                        Unconstrained::State::mxCreate());
                    
                    // Grab the reference vector 
                    auto x = Vector(std::move(X),std::move(x_));
                    
                    // Create a C++ state
                    typename MxUnconstrained::State::t state(x);

                    // Read the restart file into the C++ state 
                    MxJsonUnconstrained::read_restart(fname,x,state);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return; 
                } CATCH_MATLAB_ERRORS;
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
                // Returns the fields names of the state
                std::vector <char const *> fieldNames_() {
                    std::vector <const char*> names = {
                        "y",
                        "dy",
                        "zeta",
                        "eta0",
                        "rho",
                        "rho_old",
                        "rho_bar",
                        "eps_constr",
                        "xi_qn",
                        "xi_pg",
                        "xi_proj",
                        "xi_tang",
                        "xi_lmh",
                        "xi_lmg",
                        "xi_4",
                        "rpred",
                        "PSchur_left_type",
                        "PSchur_right_type",
                        "augsys_iter_max",
                        "augsys_rst_freq",
                        "augsys_qn_iter",
                        "augsys_pg_iter",
                        "augsys_proj_iter",
                        "augsys_tang_iter",
                        "augsys_lmh_iter",
                        "augsys_qn_iter_total",
                        "augsys_pg_iter_total",
                        "augsys_proj_iter_total",
                        "augsys_tang_iter_total",
                        "augsys_lmh_iter_total",
                        "augsys_qn_err",
                        "augsys_pg_err",
                        "augsys_proj_err",
                        "augsys_tang_err",
                        "augsys_lmh_err",
                        "augsys_qn_err_target",
                        "augsys_pg_err_target",
                        "augsys_proj_err_target",
                        "augsys_tang_err_target",
                        "augsys_lmh_err_target",
                        "augsys_iter_total",
                        "augsys_qn_failed",
                        "augsys_pg_failed",
                        "augsys_proj_failed",
                        "augsys_tang_failed",
                        "augsys_lmh_failed",
                        "augsys_failed_total",
                        "g_x",
                        "norm_gxtyp",
                        "norm_gpsgxtyp",
                        "gpxdxn_p_gx",
                        "gpxdxt",
                        "norm_gpxdxnpgx",
                        "dx_n",
                        "dx_ncp",
                        "dx_t",
                        "dx_t_uncorrected",
                        "dx_tcp_uncorrected",
                        "H_dxn",
                        "W_gradpHdxn",
                        "H_dxtuncorrected",
                        "g_diag",
                        "y_diag",
                        "qn_stop"};

                    return std::move(names);
                }
                std::vector <char const *> fieldNames() {
                    std::vector <char const*> un
                        = Unconstrained::State::fieldNames_();
                    std::vector <char const*> eq
                        = EqualityConstrained::State::fieldNames_();
                    un.reserve(un.size()+eq.size());
                    un.insert(un.end(),eq.begin(),eq.end());
                    return std::move(un); 
                }

                // Create the structure for a Matlab state
                mxManaged mxCreate() {
                    auto names = EqualityConstrained::State::fieldNames();
                    return mxCreateStructMatrix(1,1,names.size(),&(names[0]));
                }

                // Convert a C++ state to a Matlab state 
                void toMatlab_(
                    typename MxEqualityConstrained::State::t const & state,
                    mxArrayPtr & mxstate
                ){
                    toMatlab::Vector("y",state.y,mxstate);
                    toMatlab::Vector("dy",state.dy,mxstate);
                    toMatlab::Real("zeta",state.zeta,mxstate);
                    toMatlab::Real("eta0",state.eta0,mxstate);
                    toMatlab::Real("rho",state.rho,mxstate);
                    toMatlab::Real("rho_old",state.rho_old,mxstate);
                    toMatlab::Real("rho_bar",state.rho_bar,mxstate);
                    toMatlab::Real("eps_constr",state.eps_constr,mxstate);
                    toMatlab::Real("xi_qn",state.xi_qn,mxstate);
                    toMatlab::Real("xi_pg",state.xi_pg,mxstate);
                    toMatlab::Real("xi_proj",state.xi_proj,mxstate);
                    toMatlab::Real("xi_tang",state.xi_tang,mxstate);
                    toMatlab::Real("xi_lmh",state.xi_lmh,mxstate);
                    toMatlab::Real("xi_lmg",state.xi_lmg,mxstate);
                    toMatlab::Real("xi_4",state.xi_4,mxstate);
                    toMatlab::Real("rpred",state.rpred,mxstate);
                    toMatlab::Param <Operators::t> (
                        "PSchur_left_type",
                        Operators::toMatlab,
                        state.PSchur_left_type,
                        mxstate);
                    toMatlab::Param <Operators::t> (
                        "PSchur_right_type",
                        Operators::toMatlab,
                        state.PSchur_right_type,
                        mxstate);
                    toMatlab::Natural("augsys_iter_max",
                        state.augsys_iter_max,mxstate);
                    toMatlab::Natural("augsys_rst_freq",
                        state.augsys_rst_freq,mxstate);
                    toMatlab::Natural("augsys_qn_iter",
                        state.augsys_qn_iter,mxstate);
                    toMatlab::Natural("augsys_pg_iter",
                        state.augsys_pg_iter,mxstate);
                    toMatlab::Natural("augsys_proj_iter",
                        state.augsys_proj_iter,mxstate);
                    toMatlab::Natural("augsys_tang_iter",
                        state.augsys_tang_iter,mxstate);
                    toMatlab::Natural("augsys_lmh_iter",
                        state.augsys_lmh_iter,mxstate);
                    toMatlab::Natural("augsys_qn_iter_total",
                        state.augsys_qn_iter_total,mxstate);
                    toMatlab::Natural("augsys_pg_iter_total",
                        state.augsys_pg_iter_total,mxstate);
                    toMatlab::Natural("augsys_proj_iter_total",
                        state.augsys_proj_iter_total,mxstate);
                    toMatlab::Natural("augsys_tang_iter_total",
                        state.augsys_tang_iter_total,mxstate);
                    toMatlab::Natural("augsys_lmh_iter_total",
                        state.augsys_lmh_iter_total,mxstate);
                    toMatlab::Real("augsys_qn_err",
                        state.augsys_qn_err,mxstate);
                    toMatlab::Real("augsys_pg_err",
                        state.augsys_pg_err,mxstate);
                    toMatlab::Real("augsys_proj_err",
                        state.augsys_proj_err,mxstate);
                    toMatlab::Real("augsys_tang_err",
                        state.augsys_tang_err,mxstate);
                    toMatlab::Real("augsys_lmh_err",
                        state.augsys_lmh_err,mxstate);
                    toMatlab::Real("augsys_qn_err_target",
                        state.augsys_qn_err_target,mxstate);
                    toMatlab::Real("augsys_pg_err_target",
                        state.augsys_pg_err_target,mxstate);
                    toMatlab::Real("augsys_proj_err_target",
                        state.augsys_proj_err_target,mxstate);
                    toMatlab::Real("augsys_tang_err_target",
                        state.augsys_tang_err_target,mxstate);
                    toMatlab::Real("augsys_lmh_err_target",
                        state.augsys_lmh_err_target,mxstate);
                    toMatlab::Natural("augsys_iter_total",
                        state.augsys_iter_total,mxstate);
                    toMatlab::Natural("augsys_qn_failed",
                        state.augsys_qn_failed,mxstate);
                    toMatlab::Natural("augsys_pg_failed",
                        state.augsys_pg_failed,mxstate);
                    toMatlab::Natural("augsys_proj_failed",
                        state.augsys_proj_failed,mxstate);
                    toMatlab::Natural("augsys_tang_failed",
                        state.augsys_tang_failed,mxstate);
                    toMatlab::Natural("augsys_lmh_failed",
                        state.augsys_lmh_failed,mxstate);
                    toMatlab::Natural("augsys_failed_total",
                        state.augsys_failed_total,mxstate);
                    toMatlab::Vector("g_x",state.g_x,mxstate);
                    toMatlab::Real("norm_gxtyp",state.norm_gxtyp,mxstate);
                    toMatlab::Real("norm_gpsgxtyp",state.norm_gpsgxtyp,mxstate);
                    toMatlab::Vector("gpxdxn_p_gx",state.gpxdxn_p_gx,mxstate);
                    toMatlab::Vector("gpxdxt",state.gpxdxt,mxstate);
                    toMatlab::Real("norm_gpxdxnpgx",
                        state.norm_gpxdxnpgx,mxstate);
                    toMatlab::Vector("dx_n",state.dx_n,mxstate);
                    toMatlab::Vector("dx_ncp",state.dx_ncp,mxstate);
                    toMatlab::Vector("dx_t",state.dx_t,mxstate);
                    toMatlab::Vector("dx_t_uncorrected",
                        state.dx_t_uncorrected,mxstate);
                    toMatlab::Vector("dx_tcp_uncorrected",
                        state.dx_tcp_uncorrected,mxstate);
                    toMatlab::Vector("H_dxn",state.H_dxn,mxstate);
                    toMatlab::Vector("W_gradpHdxn",state.W_gradpHdxn,mxstate);
                    toMatlab::Vector("H_dxtuncorrected",
                        state.H_dxtuncorrected,mxstate);
                    toMatlab::Param <FunctionDiagnostics::t> (
                        "g_diag",
                        FunctionDiagnostics::toMatlab,
                        state.g_diag,
                        mxstate);
                    toMatlab::Param <VectorSpaceDiagnostics::t> (
                        "y_diag",
                        VectorSpaceDiagnostics::toMatlab,
                        state.y_diag,
                        mxstate);
                    toMatlab::Param <QuasinormalStop::t> (
                        "qn_stop",
                        QuasinormalStop::toMatlab,
                        state.qn_stop,
                        mxstate);
                }
                void toMatlab(
                    typename MxEqualityConstrained::State::t const & state,
                    mxArrayPtr & mxstate
                ){
                    Unconstrained::State::toMatlab_(state,mxstate);
                    EqualityConstrained::State::toMatlab_(state,mxstate);
                }
                
                // Convert a Matlab state to C++ 
                void fromMatlab_(
                    mxArrayPtr const & mxstate,
                    typename MxEqualityConstrained::State::t & state
                ){
                    fromMatlab::Vector("y",mxstate,state.y);
                    fromMatlab::Vector("dy",mxstate,state.dy);
                    fromMatlab::Real("zeta",mxstate,state.zeta);
                    fromMatlab::Real("eta0",mxstate,state.eta0);
                    fromMatlab::Real("rho",mxstate,state.rho);
                    fromMatlab::Real("rho_old",mxstate,state.rho_old);
                    fromMatlab::Real("rho_bar",mxstate,state.rho_bar);
                    fromMatlab::Real("eps_constr",mxstate,state.eps_constr);
                    fromMatlab::Real("xi_qn",mxstate,state.xi_qn);
                    fromMatlab::Real("xi_pg",mxstate,state.xi_pg);
                    fromMatlab::Real("xi_proj",mxstate,state.xi_proj);
                    fromMatlab::Real("xi_tang",mxstate,state.xi_tang);
                    fromMatlab::Real("xi_lmh",mxstate,state.xi_lmh);
                    fromMatlab::Real("xi_lmg",mxstate,state.xi_lmg);
                    fromMatlab::Real("xi_4",mxstate,state.xi_4);
                    fromMatlab::Real("rpred",mxstate,state.rpred);
                    fromMatlab::Param <Operators::t> (
                        "PSchur_left_type",
                        Operators::fromMatlab,
                        mxstate,
                        state.PSchur_left_type);
                    fromMatlab::Param <Operators::t> (
                        "PSchur_right_type",
                        Operators::fromMatlab,
                        mxstate,
                        state.PSchur_right_type);
                    fromMatlab::Natural("augsys_iter_max",
                        mxstate,state.augsys_iter_max);
                    fromMatlab::Natural("augsys_rst_freq",
                        mxstate,state.augsys_rst_freq);
                    fromMatlab::Natural("augsys_qn_iter",
                        mxstate,state.augsys_qn_iter);
                    fromMatlab::Natural("augsys_pg_iter",
                        mxstate,state.augsys_pg_iter);
                    fromMatlab::Natural("augsys_proj_iter",
                        mxstate,state.augsys_proj_iter);
                    fromMatlab::Natural("augsys_tang_iter",
                        mxstate,state.augsys_tang_iter);
                    fromMatlab::Natural("augsys_lmh_iter",
                        mxstate,state.augsys_lmh_iter);
                    fromMatlab::Natural("augsys_qn_iter_total",
                        mxstate,state.augsys_qn_iter_total);
                    fromMatlab::Natural("augsys_pg_iter_total",
                        mxstate,state.augsys_pg_iter_total);
                    fromMatlab::Natural("augsys_proj_iter_total",
                        mxstate,state.augsys_proj_iter_total);
                    fromMatlab::Natural("augsys_tang_iter_total",
                        mxstate,state.augsys_tang_iter_total);
                    fromMatlab::Natural("augsys_lmh_iter_total",
                        mxstate,state.augsys_lmh_iter_total);
                    fromMatlab::Real("augsys_qn_err",
                        mxstate,state.augsys_qn_err);
                    fromMatlab::Real("augsys_pg_err",
                        mxstate,state.augsys_pg_err);
                    fromMatlab::Real("augsys_proj_err",
                        mxstate,state.augsys_proj_err);
                    fromMatlab::Real("augsys_tang_err",
                        mxstate,state.augsys_tang_err);
                    fromMatlab::Real("augsys_lmh_err",
                        mxstate,state.augsys_lmh_err);
                    fromMatlab::Real("augsys_qn_err_target",
                        mxstate,state.augsys_qn_err_target);
                    fromMatlab::Real("augsys_pg_err_target",
                        mxstate,state.augsys_pg_err_target);
                    fromMatlab::Real("augsys_proj_err_target",
                        mxstate,state.augsys_proj_err_target);
                    fromMatlab::Real("augsys_tang_err_target",
                        mxstate,state.augsys_tang_err_target);
                    fromMatlab::Real("augsys_lmh_err_target",
                        mxstate,state.augsys_lmh_err_target);
                    fromMatlab::Natural("augsys_iter_total",
                        mxstate,state.augsys_iter_total);
                    fromMatlab::Natural("augsys_qn_failed",
                        mxstate,state.augsys_qn_failed);
                    fromMatlab::Natural("augsys_pg_failed",
                        mxstate,state.augsys_pg_failed);
                    fromMatlab::Natural("augsys_proj_failed",
                        mxstate,state.augsys_proj_failed);
                    fromMatlab::Natural("augsys_tang_failed",
                        mxstate,state.augsys_tang_failed);
                    fromMatlab::Natural("augsys_lmh_failed",
                        mxstate,state.augsys_lmh_failed);
                    fromMatlab::Natural("augsys_failed_total",
                        mxstate,state.augsys_failed_total);
                    fromMatlab::Vector("g_x",mxstate,state.g_x);
                    fromMatlab::Real("norm_gxtyp",mxstate,state.norm_gxtyp);
                    fromMatlab::Real("norm_gpsgxtyp",
                        mxstate,state.norm_gpsgxtyp);
                    fromMatlab::Vector("gpxdxn_p_gx",mxstate,state.gpxdxn_p_gx);
                    fromMatlab::Vector("gpxdxt",mxstate,state.gpxdxt);
                    fromMatlab::Real("norm_gpxdxnpgx",
                        mxstate,state.norm_gpxdxnpgx);
                    fromMatlab::Vector("dx_n",mxstate,state.dx_n);
                    fromMatlab::Vector("dx_ncp",mxstate,state.dx_ncp);
                    fromMatlab::Vector("dx_t",mxstate,state.dx_t);
                    fromMatlab::Vector("dx_t_uncorrected",
                        mxstate,state.dx_t_uncorrected);
                    fromMatlab::Vector("dx_tcp_uncorrected",
                        mxstate,state.dx_tcp_uncorrected);
                    fromMatlab::Vector("H_dxn",mxstate,state.H_dxn);
                    fromMatlab::Vector("W_gradpHdxn",mxstate,state.W_gradpHdxn);
                    fromMatlab::Vector("H_dxtuncorrected",
                        mxstate,state.H_dxtuncorrected);
                    fromMatlab::Param <FunctionDiagnostics::t> (
                        "g_diag",
                        FunctionDiagnostics::fromMatlab,
                        mxstate,
                        state.g_diag);
                    fromMatlab::Param <VectorSpaceDiagnostics::t> (
                        "y_diag",
                        VectorSpaceDiagnostics::fromMatlab,
                        mxstate,
                        state.y_diag);
                    fromMatlab::Param <QuasinormalStop::t> (
                        "qn_stop",
                        QuasinormalStop::fromMatlab,
                        mxstate,
                        state.qn_stop);
                }
                void fromMatlab(
                    mxArrayPtr const & mxstate,
                    typename MxEqualityConstrained::State::t & state
                ){
                    Unconstrained::State::fromMatlab_(mxstate,state);
                    EqualityConstrained::State::fromMatlab_(mxstate,state);
                }

                // Creates a state and inserts the elements into mxstate 
                void create(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be (X,Y,x,y)->(mxstate_out)
                    mxUnmanaged X=pInput[0],
                                Y=pInput[1],
                                x_=pInput[2],
                                y_=pInput[3];

                    // Create a vector from the user input
                    auto x = Vector(std::move(X),std::move(x_));
                    auto y = Vector(std::move(Y),std::move(y_));

                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxEqualityConstrained> (
                        EqualityConstrained::State::mxCreate());

                    // Create a new C++ state
                    typename MxEqualityConstrained::State::t state(x,y);

                    // Convert the state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS;
        
                // Read json parameters from file
                void readJson(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be (X,Y,fname,mxstate)
                    // -> (mxstate_out)
                    mxUnmanaged X=pInput[0],
                                Y=pInput[1],
                                fname_=pInput[2],
                                mxstate_=pInput[3];

                    // Grab the file name
                    auto fname = mxArrayToString(fname_.get());

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxEqualityConstrained> (
                        std::move(mxstate_));
                    auto mxstate_out = Matlab::State <MxEqualityConstrained> (
                        EqualityConstrained::State::mxCreate());
                
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    auto y_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"y"));
                    auto y = Vector(std::move(Y),std::move(y_));

                    // Create a new C++ state
                    typename MxEqualityConstrained::State::t state(x,y);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Read the JSON file into the C++ state
                    MxJsonEqualityConstrained::read(fname,state);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();
                            
                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS;
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Matlab bundle to C++ 
                void fromMatlab(
                    Matlab::Functions <MxEqualityConstrained> const & mxfns,
                    Matlab::State <MxEqualityConstrained> & mxstate,
                    typename MxEqualityConstrained::State::t const & state,
                    typename MxEqualityConstrained::Functions::t & fns 
                ) {
                    Unconstrained::Functions::fromMatlab_
                        <MxEqualityConstrained> (mxfns,mxstate,state,fns);
                    EqualityConstrained::Functions::fromMatlab_
                        <MxEqualityConstrained> (mxfns,mxstate,state,fns);
                }
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem 
                void getMin(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    // (X,Y,msg,mxfns,mxstate,smanip) -> (mxstate_out)
                    mxUnmanaged X=pInput[0],
                                Y=pInput[1],
                                msg_=pInput[2],
                                mxfns_=pInput[3],
                                mxstate_=pInput[4],
                                smanip_=pInput[5];

                    // Create a messaging object
                    auto msg = Optizelle::Matlab::Messaging::matlab(
                        std::move(msg_));

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxEqualityConstrained> (
                        std::move(mxstate_));
                    auto mxstate_out = Matlab::State <MxEqualityConstrained> (
                        EqualityConstrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    auto y_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"y"));
                    auto y = Vector(std::move(Y),std::move(y_));

                    // Create a C++ state
                    typename MxEqualityConstrained::State::t state(x,y);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Create a Matlab bundle of functions
                    Matlab::Functions <MxEqualityConstrained> mxfns(
                        mxstate_out,
                        state,
                        std::move(mxfns_));

                    // Create a C++ bundle of functions
                    typename MxEqualityConstrained::Functions::t fns;
                    
                    // Convert the Matlab bundle of functions to C++ 
                    mxfns.fromMatlab(fns);
                    
                    // Create a state manipulator 
                    Matlab::StateManipulator <MxEqualityConstrained> smanip(
                        mxstate_out,
                        mxfns,
                        std::move(smanip_));
                   
                    // Minimize
                    MxEqualityConstrained::Algorithms::getMin(
                        msg,fns,state,smanip);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS;
            }
            
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                void release(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    // (X,Y,mxstate) -> (xs,ys,reals,nats,params)
                    mxUnmanaged X=pInput[0],
                                Y=pInput[1],
                                mxstate_=pInput[2];

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxEqualityConstrained> (
                        std::move(mxstate_));
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    auto y_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"y"));
                    auto y = Vector(std::move(Y),std::move(y_));

                    // Create a C++ state
                    typename MxEqualityConstrained::State::t state(x,y);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Do a release 
                    MxEqualityConstrained::Restart::X_Vectors xs;
                    MxEqualityConstrained::Restart::Y_Vectors ys;
                    MxEqualityConstrained::Restart::Reals reals;
                    MxEqualityConstrained::Restart::Naturals nats;
                    MxEqualityConstrained::Restart::Params params;
                    MxEqualityConstrained::Restart
                        ::release(state,xs,ys,reals,nats,params);

                    // Allocate memory for the Matlab versions
                    auto mxxs = mxManaged(mxCreateCellMatrix(1,xs.size()));
                    auto mxys = mxManaged(mxCreateCellMatrix(1,ys.size()));
                    auto mxreals= mxManaged(mxCreateCellMatrix(1,reals.size()));
                    auto mxnats = mxManaged(mxCreateCellMatrix(1,nats.size()));
                    auto mxparams
                        = mxManaged(mxCreateCellMatrix(1,params.size()));

                    // Convert the restart information to Matlab 
                    toMatlab::Vectors(xs,mxxs);
                    toMatlab::Vectors(ys,mxys);
                    toMatlab::Reals(reals,mxreals);
                    toMatlab::Naturals(nats,mxnats);
                    toMatlab::Params(params,mxparams);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the ouptuts
                    pOutput[0]=mxxs.release();
                    pOutput[1]=mxys.release();
                    pOutput[2]=mxreals.release();
                    pOutput[3]=mxnats.release();
                    pOutput[4]=mxparams.release();

                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS;

                // Capture data from structures controlled by the user.  
                void capture(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    // (X,Y,mxstate,mxxs,mxys,mxreals,mxnats,mxparams)
                    // -> mxstate_out
                    mxUnmanaged X=pInput[0],
                                Y=pInput[1],
                                mxstate_=pInput[2],
                                mxxs=pInput[3],
                                mxys=pInput[4],
                                mxreals=pInput[5],
                                mxnats=pInput[6],
                                mxparams=pInput[7];

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxEqualityConstrained> (
                        std::move(mxstate_));
                    auto mxstate_out = Matlab::State <MxEqualityConstrained> (
                        EqualityConstrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    auto y_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"y"));
                    auto y = Vector(std::move(Y),std::move(y_));

                    // Create a C++ state
                    typename MxEqualityConstrained::State::t state(x,y);
                   
                    // Allocate memory for the released vectors
                    MxEqualityConstrained::Restart::X_Vectors xs;
                    MxEqualityConstrained::Restart::Y_Vectors ys;
                    MxEqualityConstrained::Restart::Reals reals;
                    MxEqualityConstrained::Restart::Naturals nats;
                    MxEqualityConstrained::Restart::Params params;
                    
                    // Convert the restart information from Matlab 
                    fromMatlab::Vectors(x,mxxs,xs);
                    fromMatlab::Vectors(y,mxys,ys);
                    fromMatlab::Reals(mxreals,reals);
                    fromMatlab::Naturals(mxnats,nats);
                    fromMatlab::Params(mxparams,params);

                    // Do a capture 
                    MxEqualityConstrained::Restart
                        ::capture(state,xs,ys,reals,nats,params);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS;
                
                // Writes a json restart file
                void write_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be (X,Y,msg,fname,state) -> ()
                    mxUnmanaged X=pInput[0],
                                Y=pInput[1],
                                fname_=pInput[2],
                                mxstate_=pInput[3];

                    // Grab the file name
                    auto fname = mxArrayToString(fname_.get());

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxEqualityConstrained> (
                        std::move(mxstate_));
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    auto y_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"y"));
                    auto y = Vector(std::move(Y),std::move(y_));

                    // Create a C++ state
                    typename MxEqualityConstrained::State::t state(x,y);
                    
                    // Convert Matlab state to C++ 
                    mxstate.fromMatlab(state);

                    // Write the restart file
                    MxJsonEqualityConstrained::write_restart(fname,state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();
                    
                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS;
                
                // Reads a json restart file
                void read_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    // (X,Y,fname,x,y) -> (mxstate)
                    mxUnmanaged X=pInput[0],
                                Y=pInput[1],
                                fname_=pInput[2],
                                x_=pInput[3],
                                y_=pInput[4];

                    // Grab the file name
                    auto fname = mxArrayToString(fname_.get());
                    
                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxEqualityConstrained> (
                        EqualityConstrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x = Vector(std::move(X),std::move(x_));
                    auto y = Vector(std::move(Y),std::move(y_));

                    // Create a C++ state
                    typename MxEqualityConstrained::State::t state(x,y);

                    // Read the restart file into the C++ state 
                    MxJsonEqualityConstrained::read_restart(fname,x,y,state);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return;

                } CATCH_MATLAB_ERRORS;
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
                // Returns the fields names of the state
                std::vector <char const *> fieldNames_() {
                    std::vector <const char*> names = {
                        "z",
                        "dz",
                        "h_x",
                        "mu",
                        "mu_est",
                        "mu_typ",
                        "eps_mu",
                        "sigma",
                        "gamma",
                        "alpha_z",
                        "h_diag",
                        "z_diag"};

                    return std::move(names);
                }
                std::vector <char const *> fieldNames() {
                    std::vector <char const*> un
                        = Unconstrained::State::fieldNames_();
                    std::vector <char const*> iq
                        = InequalityConstrained::State::fieldNames_();
                    un.reserve(un.size()+iq.size());
                    un.insert(un.end(),iq.begin(),iq.end());
                    return std::move(un); 
                }

                // Create the structure for a Matlab state
                mxManaged mxCreate() {
                    auto names = InequalityConstrained::State::fieldNames();
                    return mxCreateStructMatrix(1,1,names.size(),&(names[0]));
                }

                // Convert a C++ state to a Matlab state 
                void toMatlab_(
                    typename MxInequalityConstrained::State::t const & state,
                    mxArrayPtr & mxstate
                ){
                    toMatlab::Vector("z",state.z,mxstate);
                    toMatlab::Vector("dz",state.dz,mxstate);
                    toMatlab::Vector("h_x",state.h_x,mxstate);
                    toMatlab::Real("mu",state.mu,mxstate);
                    toMatlab::Real("mu_est",state.mu_est,mxstate);
                    toMatlab::Real("mu_typ",state.mu_typ,mxstate);
                    toMatlab::Real("eps_mu",state.eps_mu,mxstate);
                    toMatlab::Real("sigma",state.sigma,mxstate);
                    toMatlab::Real("gamma",state.gamma,mxstate);
                    toMatlab::Real("alpha_z",state.alpha_z,mxstate);
                    toMatlab::Param <FunctionDiagnostics::t> (
                        "h_diag",
                        FunctionDiagnostics::toMatlab,
                        state.h_diag,
                        mxstate);
                    toMatlab::Param <VectorSpaceDiagnostics::t> (
                        "z_diag",
                        VectorSpaceDiagnostics::toMatlab,
                        state.z_diag,
                        mxstate);
                }
                void toMatlab(
                    typename MxInequalityConstrained::State::t const & state,
                    mxArrayPtr & mxstate
                ){
                    Unconstrained::State::toMatlab_(state,mxstate);
                    InequalityConstrained::State::toMatlab_(state,mxstate);
                }
                
                // Convert a Matlab state to C++ 
                void fromMatlab_(
                    mxArrayPtr const & mxstate,
                    typename MxInequalityConstrained::State::t & state
                ){
                    fromMatlab::Vector("z",mxstate,state.z);
                    fromMatlab::Vector("dz",mxstate,state.dz);
                    fromMatlab::Vector("h_x",mxstate,state.h_x);
                    fromMatlab::Real("mu",mxstate,state.mu);
                    fromMatlab::Real("mu_est",mxstate,state.mu_est);
                    fromMatlab::Real("mu_typ",mxstate,state.mu_typ);
                    fromMatlab::Real("eps_mu",mxstate,state.eps_mu);
                    fromMatlab::Real("sigma",mxstate,state.sigma);
                    fromMatlab::Real("gamma",mxstate,state.gamma);
                    fromMatlab::Real("alpha_z",mxstate,state.alpha_z);
                    fromMatlab::Param <FunctionDiagnostics::t> (
                        "h_diag",
                        FunctionDiagnostics::fromMatlab,
                        mxstate,
                        state.h_diag);
                    fromMatlab::Param <VectorSpaceDiagnostics::t> (
                        "z_diag",
                        VectorSpaceDiagnostics::fromMatlab,
                        mxstate,
                        state.z_diag);
                }
                void fromMatlab(
                    mxArrayPtr const & mxstate,
                    typename MxInequalityConstrained::State::t & state
                ){
                    Unconstrained::State::fromMatlab_(mxstate,state);
                    InequalityConstrained::State::fromMatlab_(mxstate,state);
                }

                // Creates a state and inserts the elements into mxstate 
                void create(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be (X,Z,x,z)->(mxstate_out)
                    mxUnmanaged X=pInput[0],
                                Z=pInput[1],
                                x_=pInput[2],
                                z_=pInput[3];

                    // Create a vector from the user input
                    auto x = Vector(std::move(X),std::move(x_));
                    auto z = Vector(std::move(Z),std::move(z_));

                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxInequalityConstrained> (
                        InequalityConstrained::State::mxCreate());

                    // Create a new C++ state
                    typename MxInequalityConstrained::State::t state(x,z);

                    // Convert the state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return; 
                } CATCH_MATLAB_ERRORS;
        
                // Read json parameters from file
                void readJson(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be (X,Z,fname,mxstate)
                    // -> (mxstate_out)
                    mxUnmanaged X=pInput[0],
                                Z=pInput[1],
                                fname_=pInput[2],
                                mxstate_=pInput[3];

                    // Grab the file name
                    auto fname = mxArrayToString(fname_.get());

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxInequalityConstrained> (
                        std::move(mxstate_));
                    auto mxstate_out = Matlab::State <MxInequalityConstrained> (
                        InequalityConstrained::State::mxCreate());
                
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    auto z_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"z"));
                    auto z = Vector(std::move(Z),std::move(z_));

                    // Create a new C++ state
                    typename MxInequalityConstrained::State::t state(x,z);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Read the JSON file into the C++ state
                    MxJsonInequalityConstrained::read(fname,state);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();
                            
                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS;
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Matlab bundle to C++ 
                void fromMatlab(
                    Matlab::Functions <MxInequalityConstrained> const & mxfns,
                    Matlab::State <MxInequalityConstrained> & mxstate,
                    typename MxInequalityConstrained::State::t const & state,
                    typename MxInequalityConstrained::Functions::t & fns 
                ) {
                    Unconstrained::Functions::fromMatlab_
                        <MxInequalityConstrained> (mxfns,mxstate,state,fns);
                    InequalityConstrained::Functions::fromMatlab_
                        <MxInequalityConstrained> (mxfns,mxstate,state,fns);
                }
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem 
                void getMin(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    // (X,Z,msg,mxfns,mxstate,smanip) -> (mxstate_out)
                    mxUnmanaged X=pInput[0],
                                Z=pInput[1],
                                msg_=pInput[2],
                                mxfns_=pInput[3],
                                mxstate_=pInput[4],
                                smanip_=pInput[5];
                    
                    // Create a messaging object
                    auto msg = Optizelle::Matlab::Messaging::matlab(
                        std::move(msg_));
                        
                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxInequalityConstrained> (
                        std::move(mxstate_));
                    auto mxstate_out = Matlab::State <MxInequalityConstrained> (
                        InequalityConstrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    auto z_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"z"));
                    auto z = Vector(std::move(Z),std::move(z_));

                    // Create a C++ state
                    typename MxInequalityConstrained::State::t state(x,z);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Create a Matlab bundle of functions
                    Matlab::Functions <MxInequalityConstrained> mxfns(
                        mxstate_out,
                        state,
                        std::move(mxfns_));

                    // Create a C++ bundle of functions
                    typename MxInequalityConstrained::Functions::t fns;
                    
                    // Convert the Matlab bundle of functions to C++ 
                    mxfns.fromMatlab(fns);
                    
                    // Create a state manipulator 
                    Matlab::StateManipulator <MxInequalityConstrained> smanip(
                        mxstate_out,
                        mxfns,
                        std::move(smanip_));
                   
                    // Minimize
                    MxInequalityConstrained::Algorithms::getMin(
                        msg,fns,state,smanip);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return; 
                } CATCH_MATLAB_ERRORS;
            }
            
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                void release(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    // (X,Z,mxstate) -> (xs,zs,reals,nats,params)
                    mxUnmanaged X=pInput[0],
                                Z=pInput[1],
                                mxstate_=pInput[2];

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxInequalityConstrained> (
                        std::move(mxstate_));
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    auto z_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"z"));
                    auto z = Vector(std::move(Z),std::move(z_));

                    // Create a C++ state
                    typename MxInequalityConstrained::State::t state(x,z);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Do a release 
                    MxInequalityConstrained::Restart::X_Vectors xs;
                    MxInequalityConstrained::Restart::Z_Vectors zs;
                    MxInequalityConstrained::Restart::Reals reals;
                    MxInequalityConstrained::Restart::Naturals nats;
                    MxInequalityConstrained::Restart::Params params;
                    MxInequalityConstrained::Restart
                        ::release(state,xs,zs,reals,nats,params);

                    // Allocate memory for the Matlab versions
                    auto mxxs = mxManaged(mxCreateCellMatrix(1,xs.size()));
                    auto mxzs = mxManaged(mxCreateCellMatrix(1,zs.size()));
                    auto mxreals= mxManaged(mxCreateCellMatrix(1,reals.size()));
                    auto mxnats = mxManaged(mxCreateCellMatrix(1,nats.size()));
                    auto mxparams
                        = mxManaged(mxCreateCellMatrix(1,params.size()));

                    // Convert the restart information to Matlab 
                    toMatlab::Vectors(xs,mxxs);
                    toMatlab::Vectors(zs,mxzs);
                    toMatlab::Reals(reals,mxreals);
                    toMatlab::Naturals(nats,mxnats);
                    toMatlab::Params(params,mxparams);

                    // Set the ouptuts
                    pOutput[0]=mxxs.release();
                    pOutput[1]=mxzs.release();
                    pOutput[2]=mxreals.release();
                    pOutput[3]=mxnats.release();
                    pOutput[4]=mxparams.release();
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Return nothing 
                    return; 
                } CATCH_MATLAB_ERRORS;

                // Capture data from structures controlled by the user.  
                void capture(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    // (X,Z,mxstate,mxxs,mxzs,mxreals,mxnats,mxparams)
                    // -> mxstate_out
                    mxUnmanaged X=pInput[0],
                                Z=pInput[1],
                                mxstate_=pInput[2],
                                mxxs=pInput[3],
                                mxzs=pInput[4],
                                mxreals=pInput[5],
                                mxnats=pInput[6],
                                mxparams=pInput[7];

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxInequalityConstrained> (
                        std::move(mxstate_));
                    auto mxstate_out = Matlab::State <MxInequalityConstrained> (
                        InequalityConstrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    auto z_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"z"));
                    auto z = Vector(std::move(Z),std::move(z_));

                    // Create a C++ state
                    typename MxInequalityConstrained::State::t state(x,z);
                   
                    // Allocate memory for the released vectors
                    MxInequalityConstrained::Restart::X_Vectors xs;
                    MxInequalityConstrained::Restart::Z_Vectors zs;
                    MxInequalityConstrained::Restart::Reals reals;
                    MxInequalityConstrained::Restart::Naturals nats;
                    MxInequalityConstrained::Restart::Params params;
                    
                    // Convert the restart information from Matlab 
                    fromMatlab::Vectors(x,mxxs,xs);
                    fromMatlab::Vectors(z,mxzs,zs);
                    fromMatlab::Reals(mxreals,reals);
                    fromMatlab::Naturals(mxnats,nats);
                    fromMatlab::Params(mxparams,params);

                    // Do a capture 
                    MxInequalityConstrained::Restart
                        ::capture(state,xs,zs,reals,nats,params);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS;
                
                // Writes a json restart file
                void write_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be (X,Z,fname,state) -> ()
                    mxUnmanaged X=pInput[0],
                                Z=pInput[1],
                                fname_=pInput[2],
                                mxstate_=pInput[3];

                    // Grab the file name
                    auto fname = mxArrayToString(fname_.get());

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxInequalityConstrained> (
                        std::move(mxstate_));
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    auto z_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"z"));
                    auto z = Vector(std::move(Z),std::move(z_));
                    
                    // Create a C++ state
                    typename MxInequalityConstrained::State::t state(x,z);
                    
                    // Convert Matlab state to C++ 
                    mxstate.fromMatlab(state);

                    // Write the restart file
                    MxJsonInequalityConstrained::write_restart(fname,state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();
                    
                    // Return nothing 
                    return; 
                } CATCH_MATLAB_ERRORS;
                
                // Reads a json restart file
                void read_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    // (X,Z,fname,x,z) -> (mxstate)
                    mxUnmanaged X=pInput[0],
                                Z=pInput[1],
                                fname_=pInput[2],
                                x_=pInput[3],
                                z_=pInput[4];

                    // Grab the file name
                    auto fname = mxArrayToString(fname_.get());

                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxInequalityConstrained> (
                        InequalityConstrained::State::mxCreate());
                    
                    // Grab the reference vector 
                    auto x = Vector(std::move(X),std::move(x_));
                    auto z = Vector(std::move(Z),std::move(z_));
                    
                    // Create a C++ state
                    typename MxInequalityConstrained::State::t state(x,z);

                    // Read the restart file into the C++ state 
                    MxJsonInequalityConstrained::read_restart(fname,x,z,state);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return; 
                } CATCH_MATLAB_ERRORS;
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
                // Returns the fields names of the state
                std::vector <char const *> fieldNames() {
                    std::vector <char const*> un
                        = Unconstrained::State::fieldNames_();
                    std::vector <char const*> eq
                        = EqualityConstrained::State::fieldNames_();
                    std::vector <char const*> iq
                        = InequalityConstrained::State::fieldNames_();
                    un.reserve(un.size()+eq.size()+iq.size());
                    un.insert(un.end(),eq.begin(),eq.end());
                    un.insert(un.end(),iq.begin(),iq.end());
                    return std::move(un); 
                }

                // Create the structure for a Matlab state
                mxManaged mxCreate() {
                    auto names = Constrained::State::fieldNames();
                    return mxCreateStructMatrix(1,1,names.size(),&(names[0]));
                }

                // Convert a C++ state to a Matlab state 
                void toMatlab(
                    typename MxConstrained::State::t const & state,
                    mxArrayPtr & mxstate
                ){
                    Unconstrained::State::toMatlab_(state,mxstate);
                    EqualityConstrained::State::toMatlab_(state,mxstate);
                    InequalityConstrained::State::toMatlab_(state,mxstate);
                }
                
                // Convert a Matlab state to C++ 
                void fromMatlab(
                    mxArrayPtr const & mxstate,
                    typename MxConstrained::State::t & state
                ){
                    Unconstrained::State::fromMatlab_(mxstate,state);
                    EqualityConstrained::State::fromMatlab_(mxstate,state);
                    InequalityConstrained::State::fromMatlab_(mxstate,state);
                }

                // Creates a state and inserts the elements into mxstate 
                void create(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    // (X,Y,Z,x,y,z)->(mxstate_out)
                    mxUnmanaged X=pInput[0],
                                Y=pInput[1],
                                Z=pInput[2],
                                x_=pInput[3],
                                y_=pInput[4],
                                z_=pInput[5];

                    // Create a vector from the user input
                    auto x = Vector(std::move(X),std::move(x_));
                    auto y = Vector(std::move(Y),std::move(y_));
                    auto z = Vector(std::move(Z),std::move(z_));

                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxConstrained> (
                        Constrained::State::mxCreate());

                    // Create a new C++ state
                    typename MxConstrained::State::t state(x,y,z);

                    // Convert the state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS;
        
                // Read json parameters from file
                void readJson(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be (X,Y,Z,fname,mxstate)
                    // -> (mxstate_out)
                    mxUnmanaged X=pInput[0],
                                Y=pInput[1],
                                Z=pInput[2],
                                fname_=pInput[3],
                                mxstate_=pInput[4];

                    // Grab the file name
                    auto fname = mxArrayToString(fname_.get());

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxConstrained> (
                        std::move(mxstate_));
                    auto mxstate_out = Matlab::State <MxConstrained> (
                        Constrained::State::mxCreate());
                
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    auto y_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"y"));
                    auto y = Vector(std::move(Y),std::move(y_));
                    auto z_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"z"));
                    auto z = Vector(std::move(Z),std::move(z_));

                    // Create a new C++ state
                    typename MxConstrained::State::t state(x,y,z);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Read the JSON file into the C++ state
                    MxJsonConstrained::read(fname,state);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();
                            
                    // Return nothing 
                    return; 
                } CATCH_MATLAB_ERRORS;
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Matlab bundle to C++ 
                void fromMatlab(
                    Matlab::Functions <MxConstrained> const & mxfns,
                    Matlab::State <MxConstrained> & mxstate,
                    typename MxConstrained::State::t const & state,
                    typename MxConstrained::Functions::t & fns 
                ) {
                    Unconstrained::Functions::fromMatlab_
                        <MxConstrained> (mxfns,mxstate,state,fns);
                    EqualityConstrained::Functions::fromMatlab_
                        <MxConstrained> (mxfns,mxstate,state,fns);
                    InequalityConstrained::Functions::fromMatlab_
                        <MxConstrained> (mxfns,mxstate,state,fns);
                }
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem 
                void getMin(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    // (X,Z,msg,mxfns,mxstate,smanip) -> (mxstate_out)
                    mxUnmanaged X=pInput[0],
                                Y=pInput[1],
                                Z=pInput[2],
                                msg_=pInput[3],
                                mxfns_=pInput[4],
                                mxstate_=pInput[5],
                                smanip_=pInput[6];

                    // Create a messaging object
                    auto msg = Optizelle::Matlab::Messaging::matlab(
                        std::move(msg_));

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxConstrained> (
                        std::move(mxstate_));
                    auto mxstate_out = Matlab::State <MxConstrained> (
                        Constrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    auto y_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"y"));
                    auto y = Vector(std::move(Y),std::move(y_));
                    auto z_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"z"));
                    auto z = Vector(std::move(Z),std::move(z_));

                    // Create a C++ state
                    typename MxConstrained::State::t state(x,y,z);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Create a Matlab bundle of functions
                    Matlab::Functions <MxConstrained> mxfns(
                        mxstate_out,
                        state,
                        std::move(mxfns_));

                    // Create a C++ bundle of functions
                    typename MxConstrained::Functions::t fns;
                    
                    // Convert the Matlab bundle of functions to C++ 
                    mxfns.fromMatlab(fns);
                    
                    // Create a state manipulator 
                    Matlab::StateManipulator <MxConstrained> smanip(
                        mxstate_out,
                        mxfns,
                        std::move(smanip_));
                   
                    // Minimize
                    MxConstrained::Algorithms::getMin(
                        msg,fns,state,smanip);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS;
            }
            
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                void release(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    // (X,Y,Z,mxstate) -> (xs,ys,zs,reals,nats,params)
                    mxUnmanaged X=pInput[0],
                                Y=pInput[1],
                                Z=pInput[2],
                                mxstate_=pInput[3];

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxEqualityConstrained> (
                        std::move(mxstate_));
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    auto y_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"y"));
                    auto y = Vector(std::move(Y),std::move(y_));
                    auto z_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"z"));
                    auto z = Vector(std::move(Z),std::move(z_));

                    // Create a C++ state
                    typename MxConstrained::State::t state(x,y,z);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Do a release 
                    MxConstrained::Restart::X_Vectors xs;
                    MxConstrained::Restart::Y_Vectors ys;
                    MxConstrained::Restart::Z_Vectors zs;
                    MxConstrained::Restart::Reals reals;
                    MxConstrained::Restart::Naturals nats;
                    MxConstrained::Restart::Params params;
                    MxConstrained::Restart
                        ::release(state,xs,ys,zs,reals,nats,params);

                    // Allocate memory for the Matlab versions
                    auto mxxs = mxManaged(mxCreateCellMatrix(1,xs.size()));
                    auto mxys = mxManaged(mxCreateCellMatrix(1,ys.size()));
                    auto mxzs = mxManaged(mxCreateCellMatrix(1,zs.size()));
                    auto mxreals= mxManaged(mxCreateCellMatrix(1,reals.size()));
                    auto mxnats = mxManaged(mxCreateCellMatrix(1,nats.size()));
                    auto mxparams
                        = mxManaged(mxCreateCellMatrix(1,params.size()));

                    // Convert the restart information to Matlab 
                    toMatlab::Vectors(xs,mxxs);
                    toMatlab::Vectors(ys,mxys);
                    toMatlab::Vectors(zs,mxzs);
                    toMatlab::Reals(reals,mxreals);
                    toMatlab::Naturals(nats,mxnats);
                    toMatlab::Params(params,mxparams);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the ouptuts
                    pOutput[0]=mxxs.release();
                    pOutput[1]=mxys.release();
                    pOutput[2]=mxzs.release();
                    pOutput[3]=mxreals.release();
                    pOutput[4]=mxnats.release();
                    pOutput[5]=mxparams.release();

                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS;

                // Capture data from structures controlled by the user.  
                void capture(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    //(X,Y,Z,mxstate,mxxs,mxys,mxzs,mxreals,mxnats,mxparams)
                    // -> mxstate_out
                    mxUnmanaged X=pInput[0],
                                Y=pInput[1],
                                Z=pInput[2],
                                mxstate_=pInput[3],
                                mxxs=pInput[4],
                                mxys=pInput[5],
                                mxzs=pInput[6],
                                mxreals=pInput[7],
                                mxnats=pInput[8],
                                mxparams=pInput[9];

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxConstrained> (
                        std::move(mxstate_));
                    auto mxstate_out = Matlab::State <MxConstrained> (
                        Constrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    auto y_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"y"));
                    auto y = Vector(std::move(Y),std::move(y_));
                    auto z_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"z"));
                    auto z = Vector(std::move(Z),std::move(z_));

                    // Create a C++ state
                    typename MxConstrained::State::t state(x,y,z);
                   
                    // Allocate memory for the released vectors
                    MxConstrained::Restart::X_Vectors xs;
                    MxConstrained::Restart::Y_Vectors ys;
                    MxConstrained::Restart::Z_Vectors zs;
                    MxConstrained::Restart::Reals reals;
                    MxConstrained::Restart::Naturals nats;
                    MxConstrained::Restart::Params params;
                    
                    // Convert the restart information from Matlab 
                    fromMatlab::Vectors(x,mxxs,xs);
                    fromMatlab::Vectors(y,mxys,ys);
                    fromMatlab::Vectors(z,mxzs,zs);
                    fromMatlab::Reals(mxreals,reals);
                    fromMatlab::Naturals(mxnats,nats);
                    fromMatlab::Params(mxparams,params);

                    // Do a capture 
                    MxConstrained::Restart
                        ::capture(state,xs,ys,zs,reals,nats,params);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS;
                
                // Writes a json restart file
                void write_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be (X,Y,Z,fname,state)-> ()
                    mxUnmanaged X=pInput[0],
                                Y=pInput[1],
                                Z=pInput[2],
                                fname_=pInput[3],
                                mxstate_=pInput[4];

                    // Grab the file name
                    auto fname = mxArrayToString(fname_.get());

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxConstrained> (
                        std::move(mxstate_));
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"x"));
                    auto x = Vector(std::move(X),std::move(x_));
                    auto y_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"y"));
                    auto y = Vector(std::move(Y),std::move(y_));
                    auto z_= mxUnmanaged(mxGetField(mxstate.data->get(),0,"z"));
                    auto z = Vector(std::move(Z),std::move(z_));
                    
                    // Create a C++ state
                    typename MxConstrained::State::t state(x,y,z);
                    
                    // Convert Matlab state to C++ 
                    mxstate.fromMatlab(state);

                    // Write the restart file
                    MxJsonConstrained::write_restart(fname,state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();
                    
                    // Return nothing 
                    return; 

                } CATCH_MATLAB_ERRORS;
                
                // Reads a json restart file
                void read_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab Optizelle
                    optizelle.push_back(mexGetVariable("global","Optizelle"));

                    // Calling convention should be
                    // (X,Y,Z,fname,x,y,z) -> (mxstate)
                    mxUnmanaged X=pInput[0],
                                Y=pInput[1],
                                Z=pInput[2],
                                fname_=pInput[3],
                                x_=pInput[4],
                                y_=pInput[5],
                                z_=pInput[6];

                    // Grab the file name
                    auto fname = mxArrayToString(fname_.get());

                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxConstrained> (
                        Constrained::State::mxCreate());
                    
                    // Grab the reference vector 
                    auto x = Vector(std::move(X),std::move(x_));
                    auto y = Vector(std::move(Y),std::move(y_));
                    auto z = Vector(std::move(Z),std::move(z_));
                    
                    // Create a C++ state
                    typename MxConstrained::State::t state(x,y,z);

                    // Read the restart file into the C++ state 
                    MxJsonConstrained::read_restart(fname,x,y,z,state);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Cleanup Optizelle
                    optizelle.pop_back();

                    // Set the output
                    pOutput[0] = mxstate_out.data->release();

                    // Return nothing 
                    return; 
                } CATCH_MATLAB_ERRORS;
            }
        }
    }
}
