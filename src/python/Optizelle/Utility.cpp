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

// Handle Python and Optizelle errors
#define CATCH_PYTHON_ERRORS \
    catch(Python::Exception::t const & e) { \
        return nullptr; \
    } catch(std::exception const & e) { \
        capi::PyErr_SetString_Optizelle( \
            Optizelle::Exception::exception_to_string(e)); \
        return nullptr; \
    }

// Define a raw Python variable that we'll parse into
#define PY_RAW_VAR(name) PyObject *name##__;

// Parse a raw Python variable
#define PY_PARSE_VAR(name) &name##__

// Convert a raw Python variable into our memory management scheme
#define PY_INPUT_VAR(name) PyObjectPtr name##_(name##__,PyObjectPtr::Borrowed);

// Grab and convert the Python inputs
#define PY_VAR_3(v1,v2,v3) \
    PY_RAW_VAR(v1); \
    PY_RAW_VAR(v2); \
    PY_RAW_VAR(v3); \
    if(!PyArg_ParseTuple(args,"OOO", \
        PY_PARSE_VAR(v1), \
        PY_PARSE_VAR(v2), \
        PY_PARSE_VAR(v3) \
    )) \
        return nullptr; \
    PY_INPUT_VAR(v1); \
    PY_INPUT_VAR(v2); \
    PY_INPUT_VAR(v3);
#define PY_VAR_4(v1,v2,v3,v4) \
    PY_RAW_VAR(v1); \
    PY_RAW_VAR(v2); \
    PY_RAW_VAR(v3); \
    PY_RAW_VAR(v4); \
    if(!PyArg_ParseTuple(args,"OOOO", \
        PY_PARSE_VAR(v1), \
        PY_PARSE_VAR(v2), \
        PY_PARSE_VAR(v3), \
        PY_PARSE_VAR(v4) \
    )) \
        return nullptr; \
    PY_INPUT_VAR(v1); \
    PY_INPUT_VAR(v2); \
    PY_INPUT_VAR(v3); \
    PY_INPUT_VAR(v4);
#define PY_VAR_5(v1,v2,v3,v4,v5) \
    PY_RAW_VAR(v1); \
    PY_RAW_VAR(v2); \
    PY_RAW_VAR(v3); \
    PY_RAW_VAR(v4); \
    PY_RAW_VAR(v5); \
    if(!PyArg_ParseTuple(args,"OOOOO", \
        PY_PARSE_VAR(v1), \
        PY_PARSE_VAR(v2), \
        PY_PARSE_VAR(v3), \
        PY_PARSE_VAR(v4), \
        PY_PARSE_VAR(v5) \
    )) \
        return nullptr; \
    PY_INPUT_VAR(v1); \
    PY_INPUT_VAR(v2); \
    PY_INPUT_VAR(v3); \
    PY_INPUT_VAR(v4); \
    PY_INPUT_VAR(v5);
#define PY_VAR_6(v1,v2,v3,v4,v5,v6) \
    PY_RAW_VAR(v1); \
    PY_RAW_VAR(v2); \
    PY_RAW_VAR(v3); \
    PY_RAW_VAR(v4); \
    PY_RAW_VAR(v5); \
    PY_RAW_VAR(v6); \
    if(!PyArg_ParseTuple(args,"OOOOOO", \
        PY_PARSE_VAR(v1), \
        PY_PARSE_VAR(v2), \
        PY_PARSE_VAR(v3), \
        PY_PARSE_VAR(v4), \
        PY_PARSE_VAR(v5), \
        PY_PARSE_VAR(v6) \
    )) \
        return nullptr; \
    PY_INPUT_VAR(v1); \
    PY_INPUT_VAR(v2); \
    PY_INPUT_VAR(v3); \
    PY_INPUT_VAR(v4); \
    PY_INPUT_VAR(v5); \
    PY_INPUT_VAR(v6);
#define PY_VAR_7(v1,v2,v3,v4,v5,v6,v7) \
    PY_RAW_VAR(v1); \
    PY_RAW_VAR(v2); \
    PY_RAW_VAR(v3); \
    PY_RAW_VAR(v4); \
    PY_RAW_VAR(v5); \
    PY_RAW_VAR(v6); \
    PY_RAW_VAR(v7); \
    if(!PyArg_ParseTuple(args,"OOOOOOO", \
        PY_PARSE_VAR(v1), \
        PY_PARSE_VAR(v2), \
        PY_PARSE_VAR(v3), \
        PY_PARSE_VAR(v4), \
        PY_PARSE_VAR(v5), \
        PY_PARSE_VAR(v6), \
        PY_PARSE_VAR(v7) \
    )) \
        return nullptr; \
    PY_INPUT_VAR(v1); \
    PY_INPUT_VAR(v2); \
    PY_INPUT_VAR(v3); \
    PY_INPUT_VAR(v4); \
    PY_INPUT_VAR(v5); \
    PY_INPUT_VAR(v6); \
    PY_INPUT_VAR(v7);
#define PY_VAR_8(v1,v2,v3,v4,v5,v6,v7,v8) \
    PY_RAW_VAR(v1); \
    PY_RAW_VAR(v2); \
    PY_RAW_VAR(v3); \
    PY_RAW_VAR(v4); \
    PY_RAW_VAR(v5); \
    PY_RAW_VAR(v6); \
    PY_RAW_VAR(v7); \
    PY_RAW_VAR(v8); \
    if(!PyArg_ParseTuple(args,"OOOOOOOO", \
        PY_PARSE_VAR(v1), \
        PY_PARSE_VAR(v2), \
        PY_PARSE_VAR(v3), \
        PY_PARSE_VAR(v4), \
        PY_PARSE_VAR(v5), \
        PY_PARSE_VAR(v6), \
        PY_PARSE_VAR(v7), \
        PY_PARSE_VAR(v8) \
    )) \
        return nullptr; \
    PY_INPUT_VAR(v1); \
    PY_INPUT_VAR(v2); \
    PY_INPUT_VAR(v3); \
    PY_INPUT_VAR(v4); \
    PY_INPUT_VAR(v5); \
    PY_INPUT_VAR(v6); \
    PY_INPUT_VAR(v7); \
    PY_INPUT_VAR(v8);
#define PY_VAR_10(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10) \
    PY_RAW_VAR(v1); \
    PY_RAW_VAR(v2); \
    PY_RAW_VAR(v3); \
    PY_RAW_VAR(v4); \
    PY_RAW_VAR(v5); \
    PY_RAW_VAR(v6); \
    PY_RAW_VAR(v7); \
    PY_RAW_VAR(v8); \
    PY_RAW_VAR(v9); \
    PY_RAW_VAR(v10); \
    if(!PyArg_ParseTuple(args,"OOOOOOOOOO", \
        PY_PARSE_VAR(v1), \
        PY_PARSE_VAR(v2), \
        PY_PARSE_VAR(v3), \
        PY_PARSE_VAR(v4), \
        PY_PARSE_VAR(v5), \
        PY_PARSE_VAR(v6), \
        PY_PARSE_VAR(v7), \
        PY_PARSE_VAR(v8), \
        PY_PARSE_VAR(v9), \
        PY_PARSE_VAR(v10) \
    )) \
        return nullptr; \
    PY_INPUT_VAR(v1); \
    PY_INPUT_VAR(v2); \
    PY_INPUT_VAR(v3); \
    PY_INPUT_VAR(v4); \
    PY_INPUT_VAR(v5); \
    PY_INPUT_VAR(v6); \
    PY_INPUT_VAR(v7); \
    PY_INPUT_VAR(v8); \
    PY_INPUT_VAR(v9); \
    PY_INPUT_VAR(v10);

namespace Optizelle {
    namespace OptimizationStop { 
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & opt_stop) {
            // Do the conversion
            switch(opt_stop){
            case NotConverged:
                return Python::capi::enumToPyObject(
                    "OptimizationStop","NotConverged");
            case GradientSmall:
                return Python::capi::enumToPyObject(
                    "OptimizationStop","GradientSmall");
            case StepSmall:
                return Python::capi::enumToPyObject(
                    "OptimizationStop","StepSmall");
            case MaxItersExceeded:
                return Python::capi::enumToPyObject(
                    "OptimizationStop","MaxItersExceeded");
            case InteriorPointInstability:
                return Python::capi::enumToPyObject(
                    "OptimizationStop","InteriorPointInstability");
            case GlobalizationFailure:
                return Python::capi::enumToPyObject(
                    "OptimizationStop","GlobalizationFailure");
            case UserDefined:
                return Python::capi::enumToPyObject(
                    "OptimizationStop","UserDefined");
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(Python::PyObjectPtr const & member) {
            // Convert the member to a Natural 
            auto m=Python::capi::PyInt_AsNatural(member);

            if(m==Python::capi::enumToNatural(
                "OptimizationStop","NotConverged")
            )
                return NotConverged;
            else if(m==Python::capi::enumToNatural(
                "OptimizationStop","GradientSmall")
            )
                return GradientSmall;
            else if(m==Python::capi::enumToNatural(
                "OptimizationStop","StepSmall")
            )
                return StepSmall;
            else if(m==Python::capi::enumToNatural(
                "OptimizationStop","MaxItersExceeded")
            )
                return MaxItersExceeded;
            else if(m==Python::capi::enumToNatural(
                "OptimizationStop","InteriorPointInstability")
            )
                return InteriorPointInstability;
            else if(m==Python::capi::enumToNatural(
                "OptimizationStop","GlobalizationFailure")
            )
                return GlobalizationFailure;
            else if(m==Python::capi::enumToNatural(
                "OptimizationStop","UserDefined")
            )
                return UserDefined;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown OptimizationStop");
        }
    }
    
    namespace TruncatedStop { 
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & trunc_stop) {
            // Do the conversion
            switch(trunc_stop){
            case NotConverged:
                return Python::capi::enumToPyObject(
                    "TruncatedStop","NotConverged");
            case NegativeCurvature:
                return Python::capi::enumToPyObject(
                    "TruncatedStop","NegativeCurvature");
            case RelativeErrorSmall:
                return Python::capi::enumToPyObject(
                    "TruncatedStop","RelativeErrorSmall");
            case MaxItersExceeded:
                return Python::capi::enumToPyObject(
                    "TruncatedStop","MaxItersExceeded");
            case TrustRegionViolated:
                return Python::capi::enumToPyObject(
                    "TruncatedStop","TrustRegionViolated");
            case NanOperator:
                return Python::capi::enumToPyObject(
                    "TruncatedStop","NanOperator");
            case NanPreconditioner:
                return Python::capi::enumToPyObject(
                    "TruncatedStop","NanPreconditioner");
            case NonProjectorPreconditioner:
                return Python::capi::enumToPyObject(
                    "TruncatedStop","NonProjectorPreconditioner");
            case NonSymmetricPreconditioner:
                return Python::capi::enumToPyObject(
                    "TruncatedStop","NonSymmetricPreconditioner");
            case NonSymmetricOperator:
                return Python::capi::enumToPyObject(
                    "TruncatedStop","NonSymmetricOperator");
            case LossOfOrthogonality:
                return Python::capi::enumToPyObject(
                    "TruncatedStop","LossOfOrthogonality");
            case OffsetViolatesTrustRegion:
                return Python::capi::enumToPyObject(
                    "TruncatedStop","OffsetViolatesTrustRegion");
            case OffsetViolatesSafeguard:
                return Python::capi::enumToPyObject(
                    "TruncatedStop","OffsetViolatesSafeguard");
            case TooManyFailedSafeguard:
                return Python::capi::enumToPyObject(
                    "TruncatedStop","TooManyFailedSafeguard");
            case ObjectiveIncrease:
                return Python::capi::enumToPyObject(
                    "TruncatedStop","ObjectiveIncrease");
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(Python::PyObjectPtr const & member) {
            // Convert the member to a Natural 
            auto m=Python::capi::PyInt_AsNatural(member);

            if(m==Python::capi::enumToNatural("TruncatedStop","NotConverged"))
                return NotConverged;
            else if(m==Python::capi::enumToNatural(
                "TruncatedStop","NegativeCurvature")
            )
                return NegativeCurvature;
            else if(m==Python::capi::enumToNatural(
                "TruncatedStop","RelativeErrorSmall")
            )
                return RelativeErrorSmall;
            else if(m==Python::capi::enumToNatural(
                "TruncatedStop","MaxItersExceeded")
            )
                return MaxItersExceeded;
            else if(m==Python::capi::enumToNatural(
                "TruncatedStop","TrustRegionViolated")
            )
                return TrustRegionViolated;
            else if(m==Python::capi::enumToNatural(
                "TruncatedStop","NanOperator")
            )
                return NanOperator;
            else if(m==Python::capi::enumToNatural(
                "TruncatedStop","NanPreconditioner")
            )
                return NanPreconditioner;
            else if(m==Python::capi::enumToNatural(
                "TruncatedStop","NonProjectorPreconditioner")
            )
                return NonProjectorPreconditioner;
            else if(m==Python::capi::enumToNatural(
                "TruncatedStop","NonSymmetricPreconditioner")
            )
                return NonSymmetricPreconditioner;
            else if(m==Python::capi::enumToNatural(
                "TruncatedStop","NonSymmetricOperator")
            )
                return NonSymmetricOperator;
            else if(m==Python::capi::enumToNatural(
                "TruncatedStop","LossOfOrthogonality")
            )
                return LossOfOrthogonality;
            else if(m==Python::capi::enumToNatural(
                "TruncatedStop","OffsetViolatesTrustRegion")
            )
                return OffsetViolatesTrustRegion;
            else if(m==Python::capi::enumToNatural(
                "TruncatedStop","OffsetViolatesSafeguard")
            )
                return OffsetViolatesSafeguard;
            else if(m==Python::capi::enumToNatural(
                "TruncatedStop","TooManyFailedSafeguard")
            )
                return TooManyFailedSafeguard;
            else if(m==Python::capi::enumToNatural(
                "TruncatedStop","ObjectiveIncrease")
            )
                return ObjectiveIncrease;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown TruncatedStop");
        }
    }

    namespace AlgorithmClass { 
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & algorithm_class) {
            // Do the conversion
            switch(algorithm_class){
            case TrustRegion:
                return Python::capi::enumToPyObject(
                    "AlgorithmClass","TrustRegion");
            case LineSearch:
                return Python::capi::enumToPyObject(
                    "AlgorithmClass","LineSearch");
            case UserDefined:
                return Python::capi::enumToPyObject(
                    "AlgorithmClass","UserDefined");
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(Python::PyObjectPtr const & member) {
            // Convert the member to a Natural 
            auto m=Python::capi::PyInt_AsNatural(member);

            if(m==Python::capi::enumToNatural(
                "AlgorithmClass","TrustRegion")
            )
                return TrustRegion;
            else if(m==Python::capi::enumToNatural(
                "AlgorithmClass","LineSearch")
            )
                return LineSearch;
            else if(m==Python::capi::enumToNatural(
                "AlgorithmClass","UserDefined")
            )
                return UserDefined;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown AlgorithmClass");
        }
    }

    namespace Operators { 
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & op) {
            // Do the conversion
            switch(op){
            case Identity:
                return Python::capi::enumToPyObject("Operators","Identity");
            case Zero:
                return Python::capi::enumToPyObject("Operators","Zero");
            case ScaledIdentity:
                return Python::capi::enumToPyObject(
                    "Operators","ScaledIdentity");
            case BFGS:
                return Python::capi::enumToPyObject("Operators","BFGS");
            case InvBFGS:
                return Python::capi::enumToPyObject("Operators","InvBFGS");
            case SR1:
                return Python::capi::enumToPyObject("Operators","SR1");
            case InvSR1:
                return Python::capi::enumToPyObject("Operators","InvSR1");
            case UserDefined:
                return Python::capi::enumToPyObject("Operators","UserDefined");
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(Python::PyObjectPtr const & member) {
            // Convert the member to a Natural 
            auto m=Python::capi::PyInt_AsNatural(member);

            if(m==Python::capi::enumToNatural("Operators","Identity"))
                return Identity;
            else if(m==Python::capi::enumToNatural("Operators","Zero"))
                return Zero;
            else if(m==Python::capi::enumToNatural(
                    "Operators","ScaledIdentity")
            )
                return ScaledIdentity;
            else if(m==Python::capi::enumToNatural("Operators","BFGS"))
                return BFGS;
            else if(m==Python::capi::enumToNatural("Operators","InvBFGS"))
                return InvBFGS;
            else if(m==Python::capi::enumToNatural("Operators","SR1"))
                return SR1;
            else if(m==Python::capi::enumToNatural("Operators","InvSR1"))
                return InvSR1;
            else if(m==Python::capi::enumToNatural("Operators","UserDefined"))
                return UserDefined;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown Operators");
        }
    }

    namespace LineSearchDirection {
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & dir) {
            // Do the conversion
            switch(dir){
            case SteepestDescent:
                return Python::capi::enumToPyObject("LineSearchDirection",
                    "SteepestDescent");
            case FletcherReeves:
                return Python::capi::enumToPyObject("LineSearchDirection",
                    "FletcherReeves");
            case PolakRibiere:
                return Python::capi::enumToPyObject("LineSearchDirection",
                    "PolakRibiere");
            case HestenesStiefel:
                return Python::capi::enumToPyObject("LineSearchDirection",
                    "HestenesStiefel");
            case BFGS:
                return Python::capi::enumToPyObject(
                    "LineSearchDirection","BFGS");
            case NewtonCG:
                return Python::capi::enumToPyObject(
                    "LineSearchDirection","NewtonCG");
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(Python::PyObjectPtr const & member) {
            // Convert the member to a Natural 
            auto m=Python::capi::PyInt_AsNatural(member);

            if(m==Python::capi::enumToNatural("LineSearchDirection",
                "SteepestDescent")
            )
                return SteepestDescent;
            else if(m==Python::capi::enumToNatural("LineSearchDirection",
                "FletcherReeves")
            )
                return FletcherReeves;
            else if(m==Python::capi::enumToNatural("LineSearchDirection",
                "PolakRibiere")
            )
                return PolakRibiere;
            else if(m==Python::capi::enumToNatural("LineSearchDirection",
                "HestenesStiefel")
            )
                return HestenesStiefel;
            else if(m==Python::capi::enumToNatural(
                "LineSearchDirection","BFGS")
            )
                return BFGS;
            else if(m==Python::capi::enumToNatural(
                "LineSearchDirection","NewtonCG")
            )
                return NewtonCG;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown LineSearchDirection");
        }
    }

    namespace LineSearchKind { 
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & kind) {
            // Do the conversion
            switch(kind){
            case GoldenSection:
                return Python::capi::enumToPyObject(
                    "LineSearchKind","GoldenSection");
            case BackTracking:
                return Python::capi::enumToPyObject(
                    "LineSearchKind","BackTracking");
            case TwoPointA:
                return Python::capi::enumToPyObject(
                    "LineSearchKind","TwoPointA");
            case TwoPointB:
                return Python::capi::enumToPyObject(
                    "LineSearchKind","TwoPointB");
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(Python::PyObjectPtr const & member) {
            // Convert the member to a Natural 
            auto m=Python::capi::PyInt_AsNatural(member);

            if(m==Python::capi::enumToNatural(
                "LineSearchKind","GoldenSection")
            )
                return GoldenSection;
            else if(m==Python::capi::enumToNatural(
                "LineSearchKind","BackTracking")
            )
                return BackTracking;
            else if(m==Python::capi::enumToNatural(
                "LineSearchKind","TwoPointA")
            )
                return TwoPointA;
            else if(m==Python::capi::enumToNatural(
                "LineSearchKind","TwoPointB")
            )
                return TwoPointB;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown LineSearchKind");
        }
    }

    namespace OptimizationLocation { 
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & loc) {
            // Do the conversion
            switch(loc){
            case BeginningOfOptimization:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","BeginningOfOptimization");
            case BeforeInitialFuncAndGrad:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","BeforeInitialFuncAndGrad");
            case AfterInitialFuncAndGrad:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","AfterInitialFuncAndGrad");
            case BeforeOptimizationLoop:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","BeforeOptimizationLoop");
            case BeginningOfOptimizationLoop:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","BeginningOfOptimizationLoop");
            case BeforeSaveOld:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","BeforeSaveOld");
            case BeforeStep:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","BeforeStep");
            case BeforeGetStep:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","BeforeGetStep");
            case GetStep:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","GetStep");
            case AfterStepBeforeGradient:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","AfterStepBeforeGradient");
            case AfterGradient:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","AfterGradient");
            case BeforeQuasi:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","BeforeQuasi");
            case AfterQuasi:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","AfterQuasi");
            case AfterCheckStop:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","AfterCheckStop");
            case EndOfOptimizationIteration:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","EndOfOptimizationIteration");
            case BeforeLineSearch:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","BeforeLineSearch");
            case AfterRejectedTrustRegion:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","AfterRejectedTrustRegion");
            case AfterRejectedLineSearch:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","AfterRejectedLineSearch");
            case BeforeActualVersusPredicted:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","BeforeActualVersusPredicted");
            case EndOfOptimization:
                return Python::capi::enumToPyObject(
                    "OptimizationLocation","EndOfOptimization");
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(Python::PyObjectPtr const & member) {
            // Convert the member to a Natural 
            auto m=Python::capi::PyInt_AsNatural(member);

            if(m==Python::capi::enumToNatural(
                "OptimizationLocation","BeginningOfOptimization"))
                return BeginningOfOptimization;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","BeforeInitialFuncAndGrad"))
                return BeforeInitialFuncAndGrad;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","AfterInitialFuncAndGrad"))
                return AfterInitialFuncAndGrad;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","BeforeOptimizationLoop"))
                return BeforeOptimizationLoop;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","BeginningOfOptimizationLoop"))
                return BeginningOfOptimizationLoop;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","BeforeSaveOld"))
                return BeforeSaveOld;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","BeforeStep"))
                return BeforeStep;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","BeforeGetStep"))
                return BeforeGetStep;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","GetStep"))
                return GetStep;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","AfterStepBeforeGradient"))
                return AfterStepBeforeGradient;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","AfterGradient"))
                return AfterGradient;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","BeforeQuasi"))
                return BeforeQuasi;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","AfterQuasi"))
                return AfterQuasi;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","AfterCheckStop"))
                return AfterCheckStop;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","EndOfOptimizationIteration"))
                return EndOfOptimizationIteration;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","BeforeLineSearch"))
                return BeforeLineSearch;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","AfterRejectedTrustRegion"))
                return AfterRejectedTrustRegion;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","AfterRejectedLineSearch"))
                return AfterRejectedLineSearch;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","BeforeActualVersusPredicted"))
                return BeforeActualVersusPredicted;
            else if(m==Python::capi::enumToNatural(
                "OptimizationLocation","EndOfOptimization"))
                return EndOfOptimization;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown OptimizationLocation");
        }
    }

    namespace FunctionDiagnostics { 
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & diag) {
            // Do the conversion
            switch(diag){
            case NoDiagnostics:
                return Python::capi::enumToPyObject("FunctionDiagnostics",
                    "NoDiagnostics");
            case FirstOrder:
                return Python::capi::enumToPyObject("FunctionDiagnostics",
                    "FirstOrder");
            case SecondOrder:
                return Python::capi::enumToPyObject("FunctionDiagnostics",
                    "SecondOrder");
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(Python::PyObjectPtr const & member) {
            // Convert the member to a Natural 
            auto m=Python::capi::PyInt_AsNatural(member);

            if(m==Python::capi::enumToNatural(
                "FunctionDiagnostics","NoDiagnostics")
            )
                return NoDiagnostics;
            else if(m==Python::capi::enumToNatural("FunctionDiagnostics",
                "FirstOrder")
            )
                return FirstOrder;
            else if(m==Python::capi::enumToNatural("FunctionDiagnostics",
                "SecondOrder")
            )
                return SecondOrder;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown FunctionDiagnostics");
        }
    }

    namespace VectorSpaceDiagnostics { 
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & diag) {
            // Do the conversion
            switch(diag){
            case NoDiagnostics:
                return Python::capi::enumToPyObject("VectorSpaceDiagnostics",
                    "NoDiagnostics");
            case Basic:
                return Python::capi::enumToPyObject("VectorSpaceDiagnostics",
                    "Basic");
            case EuclideanJordan:
                return Python::capi::enumToPyObject("VectorSpaceDiagnostics",
                    "EuclideanJordan");
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(Python::PyObjectPtr const & member) {
            // Convert the member to a Natural 
            auto m=Python::capi::PyInt_AsNatural(member);

            if(m==Python::capi::enumToNatural("VectorSpaceDiagnostics",
                "NoDiagnostics")
            )
                return NoDiagnostics;
            else if(m==Python::capi::enumToNatural("VectorSpaceDiagnostics",
                "Basic")
            )
                return Basic;
            else if(m==Python::capi::enumToNatural("VectorSpaceDiagnostics",
                "EuclideanJordan")
            )
                return EuclideanJordan;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown VectorSpaceDiagnostics");
        }
    }

    namespace DiagnosticScheme { 
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & dscheme) {
            // Do the conversion
            switch(dscheme){
            case Never:
                return Python::capi::enumToPyObject("DiagnosticScheme","Never");
            case DiagnosticsOnly:
                return Python::capi::enumToPyObject("DiagnosticScheme",
                    "DiagnosticsOnly");
            case EveryIteration:
                return Python::capi::enumToPyObject("DiagnosticScheme",
                    "EveryIteration");
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(Python::PyObjectPtr const & member) {
            // Convert the member to a Natural 
            auto m=Python::capi::PyInt_AsNatural(member);

            if(m==Python::capi::enumToNatural("DiagnosticScheme","Never"))
                return Never;
            else if(m==Python::capi::enumToNatural("DiagnosticScheme",
                "DiagnosticsOnly")
            )
                return DiagnosticsOnly;
            else if(m==Python::capi::enumToNatural("DiagnosticScheme",
                "EveryIteration")
            )
                return EveryIteration;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown DiagnosticScheme");
        }
    }

    namespace ToleranceKind { 
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & eps_kind) {
            // Do the conversion
            switch(eps_kind){
            case Relative:
                return Python::capi::enumToPyObject("ToleranceKind",
                    "Relative");
            case Absolute:
                return Python::capi::enumToPyObject("ToleranceKind",
                    "Absolute");
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(Python::PyObjectPtr const & member) {
            // Convert the member to a Natural 
            auto m=Python::capi::PyInt_AsNatural(member);

            if(m==Python::capi::enumToNatural("ToleranceKind",
                "Relative")
            )
                return Relative;
            else if(m==Python::capi::enumToNatural("ToleranceKind",
                "Absolute")
            )
                return Absolute;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown ToleranceKind");
        }
    }

    namespace QuasinormalStop{ 
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & qn_stop) {
            // Do the conversion
            switch(qn_stop){
            case Newton:
                return Python::capi::enumToPyObject("QuasinormalStop",
                    "Newton");
            case CauchyTrustRegion:
                return Python::capi::enumToPyObject("QuasinormalStop",
                    "CauchyTrustRegion");
            case CauchySafeguard:
                return Python::capi::enumToPyObject("QuasinormalStop",
                    "CauchySafeguard");
            case DoglegTrustRegion:
                return Python::capi::enumToPyObject("QuasinormalStop",
                    "DoglegTrustRegion");
            case DoglegSafeguard:
                return Python::capi::enumToPyObject("QuasinormalStop",
                    "DoglegSafeguard");
            case NewtonTrustRegion:
                return Python::capi::enumToPyObject("QuasinormalStop",
                    "NewtonTrustRegion");
            case NewtonSafeguard:
                return Python::capi::enumToPyObject("QuasinormalStop",
                    "NewtonSafeguard");
            case Feasible:
                return Python::capi::enumToPyObject("QuasinormalStop",
                    "Feasible");
            case CauchySolved:
                return Python::capi::enumToPyObject("QuasinormalStop",
                    "CauchySolved");
            case LocalMin:
                return Python::capi::enumToPyObject("QuasinormalStop",
                    "LocalMin");
            case NewtonFailed:
                return Python::capi::enumToPyObject("QuasinormalStop",
                    "NewtonFailed");
            }
        }

        // Converts a Python enumerated type to t 
        t fromPython(Python::PyObjectPtr const & member) {
            // Convert the member to a Natural 
            auto m=Python::capi::PyInt_AsNatural(member);

            if(m==Python::capi::enumToNatural("QuasinormalStop",
                "Newton")
            )
                return Newton;
            else if(m==Python::capi::enumToNatural("QuasinormalStop",
                "CauchyTrustRegion")
            )
                return CauchyTrustRegion;
            else if(m==Python::capi::enumToNatural("QuasinormalStop",
                "CauchySafeguard")
            )
                return CauchySafeguard;
            else if(m==Python::capi::enumToNatural("QuasinormalStop",
                "DoglegTrustRegion")
            )
                return DoglegTrustRegion;
            else if(m==Python::capi::enumToNatural("QuasinormalStop",
                "DoglegSafeguard")
            )
                return DoglegSafeguard;
            else if(m==Python::capi::enumToNatural("QuasinormalStop",
                "NewtonTrustRegion")
            )
                return NewtonTrustRegion;
            else if(m==Python::capi::enumToNatural("QuasinormalStop",
                "NewtonSafeguard")
            )
                return NewtonSafeguard;
            else if(m==Python::capi::enumToNatural("QuasinormalStop",
                "Feasible")
            )
                return Feasible;
            else if(m==Python::capi::enumToNatural("QuasinormalStop",
                "CauchySolved")
            )
                return CauchySolved;
            else if(m==Python::capi::enumToNatural("QuasinormalStop",
                "LocalMin")
            )
                return LocalMin;
            else if(m==Python::capi::enumToNatural("QuasinormalStop",
                "NewtonFailed")
            )
                return NewtonFailed;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown QuasinormalStep");
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
                auto module = Python::capi::PyImport_ImportModule(
                    "Optizelle.json.Serialization"); 

                // Now, get the serialize routine
                auto serialize =
                    Python::capi::PyObject_GetAttrString(module,"serialize");
                
                // Make a Python object of the name and iteration 
                auto name = Python::capi::PyString_FromString(name_.c_str());
                auto iter = Python::capi::PyInt_FromNatural(iter_);

                // Call the serialize routine on the vector
                auto x_json = Python::capi::PyObject_CallObject3(
                    serialize,
                    x.data,
                    name,
                    iter,
                    __LOC__
                        + ", evaluation of the serialize function failed");

                // Convert the serialized vector to a string and return it 
                return Python::capi::PyString_AsString(x_json);
            }

            static Python::Vector deserialize (
                Python::Vector const & x_,
                std::string const & x_json_
            ) {
                // Grab the serialization module 
                auto module = Python::capi::PyImport_ImportModule(
                    "Optizelle.json.Serialization"); 

                // Now, get the deserialize routine
                auto deserialize = Python::capi::PyObject_GetAttrString(
                    module,"deserialize");

                // Convert the inputed string into Python
                auto x_json= Python::capi::PyString_FromString(x_json_.c_str());

                // Call the deserialize routine on the reference vector and the
                // json vector
                auto x_raw = Python::capi::PyObject_CallObject2(
                    deserialize,
                    x_.init().data,
                    x_json,
                    __LOC__
                        + ", evaluation of the deserialize function failed");

                // Return a vector based on this information 
                return Python::Vector(x_.vs,x_raw);
            }
        };
    }

    namespace Python {
        // Grab the pointer 
        PyObjectPtr::PyObjectPtr(
            PyObject * const ptr_,
            PyObjectPtr::Mode const & mode
        ) : ptr(ptr_) {
            // If we have a borrowed reference, increase the reference count
            if(mode==PyObjectPtr::Borrowed)
                Py_XINCREF(ptr);
        }

        // Copy semantics 
        PyObjectPtr::PyObjectPtr(PyObjectPtr const & p) : ptr(p.ptr) {
            Py_XINCREF(ptr);
        }
        PyObjectPtr & PyObjectPtr::operator = (PyObjectPtr& p) {
            // Decrease the reference count on this object first
            if(ptr)
                Py_XDECREF(ptr);

            // Then, grab the new pointer and increase the reference to that
            ptr =p.ptr;
            Py_XINCREF(ptr);
        }

        // Move semantics
        PyObjectPtr::PyObjectPtr(PyObjectPtr && p) : ptr(p.ptr) {
            p.ptr = nullptr;
        }
        PyObjectPtr & PyObjectPtr::operator = (PyObjectPtr && p) { 
            // Decrease the reference count on this object first
            if(ptr)
                Py_XDECREF(ptr);

            // Then, grab the new pointer
            ptr=p.ptr;
            p.ptr = nullptr;
        }

        // On a get, we simply return the pointer.
        PyObject * PyObjectPtr::get() const {
            return ptr;
        }

        // On destruction, decrement the reference count 
        PyObjectPtr::~PyObjectPtr() {
            if(ptr)
                Py_XDECREF(ptr);
        }
        namespace capi {
            PyObjectPtr PyImport_ImportModule(const char *name) {
                auto ret = ::PyImport_ImportModule(name); 
                if(!ret)
                    throw Python::Exception::t(__LOC__
                        + ", unable to open the module " + name);
                else
                    return ret;
            }

            PyObjectPtr PyString_FromString(const char *v) {
                if(!v)
                    throw Python::Exception::t(__LOC__
                        + ", can't convert a nullstring into a Python string");
                auto ret = ::PyString_FromString(v); 
                if(!ret)
                    throw Python::Exception::t(__LOC__
                        + ", unable to convert the string " + v 
                        + " into a Python string") ;
                else
                    return ret;
            }
            std::string PyString_AsString(PyObjectPtr const & string) {
                auto ret = ::PyString_AsString(string.get()); 
                if(!ret)
                    throw Python::Exception::t(__LOC__
                        + ", unable to convert a Python object into a string");
                else
                    return std::string(ret);
            }

            Natural PyInt_AsNatural(PyObjectPtr const & io) {
                return PyInt_AsSsize_t(io.get());
            }
            PyObjectPtr PyInt_FromNatural(Natural const & ival) {
                return PyInt_FromSize_t(ival);
            }

            PyObjectPtr PyFloat_FromDouble(double v) {
                auto ret = ::PyFloat_FromDouble(v);
                if(!ret)
                    throw Python::Exception::t(__LOC__
                        + ", unable to convert the float " + std::to_string(v)
                        + " in a Python object");
                else
                    return ret;
            }
            double PyFloat_AsDouble(PyObjectPtr const & pyfloat) {
                auto x = ::PyFloat_AsDouble(pyfloat.get());
                if(::PyErr_Occurred())
                    throw Python::Exception::t(__LOC__
                        + ", unable to convert object into a float");
                return x;
            }

            PyObjectPtr PyObject_GetAttrString(
                PyObjectPtr const & o,
                const char *attr_name
            ) {
                auto ret = ::PyObject_GetAttrString(o.get(),attr_name);
                if(!ret)
                    throw Python::Exception::t(__LOC__
                        + ", unable to get the attribute " + attr_name
                        + " in a Python object");
                else
                    return ret;
            }
            void PyObject_SetAttrString(
                PyObjectPtr & o,
                const char * attr_name,
                PyObjectPtr const & v
            ) {
                auto ret = PyObject_SetAttrString(o.get(),attr_name,v.get());
                if(ret==-1)
                    throw Python::Exception::t(__LOC__
                        + ", unable set the attribute " + attr_name 
                        + " in a Python object");
            }
            PyObjectPtr PyObject_CallObject1(
                PyObjectPtr const & fn,
                PyObjectPtr const & arg1,
                std::string const & errmsg
            ) {
                auto args = capi::PyTuple_New(1);
                capi::PyTuple_SetItem(args,0,arg1);
                auto ret = ::PyObject_CallObject(fn.get(),args.get()); 
                if(!ret)
                    throw Python::Exception::t(errmsg);
                else
                    return ret;
            }
            PyObjectPtr PyObject_CallObject2(
                PyObjectPtr const & fn,
                PyObjectPtr const & arg1,
                PyObjectPtr const & arg2,
                std::string const & errmsg
            ) {
                auto args = capi::PyTuple_New(2);
                capi::PyTuple_SetItem(args,0,arg1);
                capi::PyTuple_SetItem(args,1,arg2);
                auto ret = ::PyObject_CallObject(fn.get(),args.get()); 
                if(!ret)
                    throw Python::Exception::t(errmsg);
                else
                    return ret;
            }
            PyObjectPtr PyObject_CallObject3(
                PyObjectPtr const & fn,
                PyObjectPtr const & arg1,
                PyObjectPtr const & arg2,
                PyObjectPtr const & arg3,
                std::string const & errmsg
            ) {
                auto args = capi::PyTuple_New(3);
                capi::PyTuple_SetItem(args,0,arg1);
                capi::PyTuple_SetItem(args,1,arg2);
                capi::PyTuple_SetItem(args,2,arg3);
                auto ret = ::PyObject_CallObject(fn.get(),args.get()); 
                if(!ret)
                    throw Python::Exception::t(errmsg);
                else
                    return ret;
            }
            PyObjectPtr PyObject_CallObject4(
                PyObjectPtr const & fn,
                PyObjectPtr const & arg1,
                PyObjectPtr const & arg2,
                PyObjectPtr const & arg3,
                PyObjectPtr const & arg4,
                std::string const & errmsg
            ) {
                auto args = capi::PyTuple_New(4);
                capi::PyTuple_SetItem(args,0,arg1);
                capi::PyTuple_SetItem(args,1,arg2);
                capi::PyTuple_SetItem(args,2,arg3);
                capi::PyTuple_SetItem(args,3,arg4);
                auto ret = ::PyObject_CallObject(fn.get(),args.get()); 
                if(!ret)
                    throw Python::Exception::t(errmsg);
                else
                    return ret;
            }

            PyObjectPtr PyTuple_New(Py_ssize_t const & len) {
                auto ret = ::PyTuple_New(len);
                if(!ret)
                    throw Python::Exception::t(__LOC__
                        + ", unable to create a tuple of size "
                        + std::to_string(len) + " in a Python tuple");
                else
                    return ret;
            }
            void PyTuple_SetItem(
                PyObjectPtr const & p,
                Py_ssize_t const & pos,
                PyObjectPtr const & o
            ) {
                Py_INCREF(o.get());
                auto ret = ::PyTuple_SetItem(p.get(),pos,o.get());
                if(ret)
                    throw Python::Exception::t(__LOC__
                        + ", error setting an item at position "
                        + std::to_string(pos) + " in a Python tuple");
            }
            PyObjectPtr PyTuple_GetItem(
                PyObjectPtr const & p,
                Py_ssize_t const & pos
            ) {
                auto ret = ::PyTuple_GetItem(p.get(),pos);
                if(!ret)
                    throw Python::Exception::t(__LOC__
                        + ", error getting an item at position "
                        + std::to_string(pos) + " in a Python tuple");
                return PyObjectPtr(ret,PyObjectPtr::Borrowed);
            }
            PyObjectPtr PyTuple_Pack_2(
                PyObjectPtr const & item1,
                PyObjectPtr const & item2
            ) {
                auto ret = ::PyTuple_Pack(2,item1.get(),item2.get());
                if(!ret)
                    throw Python::Exception::t(__LOC__
                        + ", unable to create a tuple with two items"); 
                return ret;
            }

            PyObjectPtr PyList_New(Py_ssize_t const & len) {
                auto ret = ::PyList_New(len);
                if(!ret)
                    throw Python::Exception::t(__LOC__
                        + ", error creating a list of size"
                        + std::to_string(len));
                return ret;
            }
            void PyList_Append(PyObjectPtr & list, PyObjectPtr const & item) {
                auto ret = ::PyList_Append(list.get(),item.get());
                if(ret == -1)
                    throw Python::Exception::t(__LOC__
                        + ", error appending an item to a list");
            }
            Natural PyList_Size(PyObjectPtr const & list) {
                auto ret = PyList_Size(list.get());
                return ret  < 0 ? 0 : ret;
            }
            PyObjectPtr PyList_GetItem(
                PyObjectPtr const & list,
                Py_ssize_t const & index
            ) {
                auto ret = PyList_GetItem(list.get(),index);
                if(!ret)
                    throw Python::Exception::t(__LOC__
                        + ", unable to get item " + std::to_string(index)
                        + " in a list");
                return PyObjectPtr(ret,PyObjectPtr::Borrowed);
            }

            // Calls the Optizelle exception with a string
            void PyErr_SetString_Optizelle(std::string const & msg) {
                auto module = capi::PyImport_ImportModule("Optizelle"); 
                auto exception = capi::PyObject_GetAttrString(
                    module,
                    "Exception");
                ::PyErr_SetString(exception.get(),msg.c_str());
            }

            // Deep copy of a Python object and return the result
            PyObjectPtr deepcopy(PyObjectPtr const & in) {
                // Grab the deepcopy function from the copy module 
                auto module = capi::PyImport_ImportModule("copy"); 
                auto deepcopy = capi::PyObject_GetAttrString(module,"deepcopy");

                // Call deepcopy on in and return the result
                return capi::PyObject_CallObject1(
                    deepcopy,
                    in,
                    __LOC__ + ", failed to deep copy an object");
            }

            // Converts an Optizelle enumerated type to a PyObject * 
            PyObjectPtr enumToPyObject(
                std::string const & type,
                std::string const & member 
            ) {
                // Grab the enumerated type object from the Optizelle module.
                // We just use simple classes in Python to represent the
                // enumerated type
                auto module = capi::PyImport_ImportModule("Optizelle"); 
                auto pyclass= capi::PyObject_GetAttrString(module,type.c_str());

                // Grab and return the member
                return PyObject_GetAttrString(pyclass.get(),member.c_str());
            }
           
            // Converts an Optizelle enumerated type to a Natural based on
            // the scheme in the Python enumerated type
            Natural enumToNatural(
                std::string const & type,
                std::string const & member 
            ) {
                // Grab the PyObject * for the type and member requested
                auto obj = capi::enumToPyObject(type,member);

                // Convert and return the member
                return capi::PyInt_AsNatural(obj);
            }
        }

        // A messaging utility that hooks directly into Python
        namespace Messaging {
            Optizelle::Messaging::t python(PyObjectPtr const & print) {
                return [print](std::string const & msg_) {
                    // Call the print function
                    auto msg = capi::PyString_FromString(msg_.c_str());
                    auto ret = capi::PyObject_CallObject1(
                        print,
                        msg,
                         __LOC__
                            + ", evaluation of the Messaging function failed");
                };
            }
        }

        // Create a vector with the appropriate vector space 
        Vector::Vector(PyObjectPtr const & vs_,PyObjectPtr const & vec_) :
            vs(vs_), data(vec_) {}
            
        // Memory allocation and size setting 
        Vector Vector::init() const {
            // Call the init function on the internal and store in y 
            auto init = capi::PyObject_GetAttrString(vs,"init");
            auto y = capi::PyObject_CallObject1(
                init,
                data,
                __LOC__
                    + ", evaluation of the vector space function init failed");

            // Create and return a new vector based on y
            return Vector(vs,y);
        }
        
        // y <- x (Shallow.  No memory allocation.)  Internal is y.
        void Vector::copy(Vector const & x) { 
            // Call the copy function on x and the internal 
            auto copy = capi::PyObject_GetAttrString(vs,"copy");
            capi::PyObject_CallObject2(
                copy,
                x.data,
                data,
                __LOC__
                    + ", evaluation of the vector space function copy failed");
        }

        // x <- alpha * x.  Internal is x.
        void Vector::scal(double const & alpha_) { 
            // Call the scal function on alpha and the internal storage 
            auto scal = capi::PyObject_GetAttrString(vs,"scal");
            auto alpha = capi::PyFloat_FromDouble(alpha_);
            capi::PyObject_CallObject2(
                scal,
                alpha,
                data,
                __LOC__
                    + ", evaluation of the vector space function scal failed");
        } 

        // x <- 0.  Internal is x. 
        void Vector::zero() { 
            // Call the zero function on this vector.
            auto zero = capi::PyObject_GetAttrString(vs,"zero");
            capi::PyObject_CallObject1(
                zero,
                data,
                __LOC__
                    + ", evaluation of the vector space function zero failed");
        }

        // y <- alpha * x + y.   Internal is y.
        void Vector::axpy(double const & alpha_,Vector const & x) { 
            // Call the axpy function on alpha, x, and the internal storage.
            auto axpy = capi::PyObject_GetAttrString(vs,"axpy");
            auto alpha = capi::PyFloat_FromDouble(alpha_);
            capi::PyObject_CallObject3(
                axpy,
                alpha,
                x.data,
                data,
                __LOC__
                    + ", evaluation of the vector space function axpy failed");
        }

        // innr <- <x,y>.  Internal is y.
        double Vector::innr(Vector const & x) const {
            // Call the innr function on x and the internal.  Store in z. 
            auto innr = capi::PyObject_GetAttrString(vs,"innr");
            auto z = capi::PyObject_CallObject2(
                innr,
                x.data,
                data,
                __LOC__
                    + ", evaluation of the vector space function innr failed");

            // Return the result 
            return capi::PyFloat_AsDouble(z); 
        }

        // x <- random.  Internal is x. 
        void Vector::rand() { 
            // Call the rand function on this vector.
            auto rand = capi::PyObject_GetAttrString(vs,"rand");
            capi::PyObject_CallObject1(
                rand,
                data,
                __LOC__
                    + ", evaluation of the vector space function rand failed");
        }

        // Jordan product, z <- x o y.  Internal is z.
        void Vector::prod(Vector const & x,Vector const & y) { 
            // Call the prod function on x, y, and the internal 
            auto prod = capi::PyObject_GetAttrString(vs,"prod");
            capi::PyObject_CallObject3(
                prod,
                x.data,
                y.data,
                data,
                __LOC__
                    + ", evaluation of the vector space function prod failed");
        }

        // Identity element, x <- e such that x o e = x .  Internal is x.
        void Vector::id() {
            // Call the id function on the internal.
            auto id = capi::PyObject_GetAttrString(vs,"id");
            capi::PyObject_CallObject1(
                id,
                data,
                __LOC__
                    + ", evaluation of the vector space function id failed");
        } 

        // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y.
        // Internal is z.
        void Vector::linv(Vector const & x, Vector const & y) { 
            // Call the linv function on x, y, and the internal
            auto linv = capi::PyObject_GetAttrString(vs,"linv");
            capi::PyObject_CallObject3(
                linv,
                x.data,
                y.data,
                data,
                __LOC__
                    + ", evaluation of the vector space function linv failed");
        }

        // Barrier function, barr <- barr(x) where x o grad barr(x) = e.
        // Internal is x.
        double Vector::barr() const {
            // Call the barr function on the internal.  Store in z.
            auto barr = capi::PyObject_GetAttrString(vs,"barr");
            auto z = capi::PyObject_CallObject1(
                barr,
                data,
                __LOC__
                    + ", evaluation of the vector space function barr failed");

            // Return the result 
            return capi::PyFloat_AsDouble(z);
        }

        // Line search, srch <- argmax {alpha in Real >= 0 : alpha x + y >= 0} 
        // where y > 0.  Internal is y.
        double Vector::srch(Vector const & x) const {
            // Call the srch function on x and the internal.  Store in z.
            auto srch = capi::PyObject_GetAttrString(vs,"srch");
            auto z = capi::PyObject_CallObject2(
                srch,
                x.data,
                data,
                __LOC__
                    + ", evaluation of the vector space function srch failed");

            // Return the result 
            return capi::PyFloat_AsDouble(z);
        } 

        // Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
        // operator.  Internal is x.
        void Vector::symm() { 
            // Call the symm function on the internal.
            auto symm = capi::PyObject_GetAttrString(vs,"symm");
            capi::PyObject_CallObject1(
                symm,
                data,
                __LOC__
                    + ", evaluation of the vector space function symm failed");
        }
        
        // Converts (copies) a value into Python.  This assumes memory
        // has been allocated both in the vector as well as Python.
        void Vector::toPython(PyObjectPtr const & ptr) const {
            // Call the copy function on the internal and x
            auto copy = capi::PyObject_GetAttrString(vs,"copy");
            capi::PyObject_CallObject2(
                copy,
                data,
                ptr,
                __LOC__
                    + ", evaluation of the vector space function copy failed");
        }
        
        // Converts (copies) a value from Python.  This assumes memory
        // has been allocated both in the vector as well as Python.
        void Vector::fromPython(PyObjectPtr const & ptr) {
            // Call the copy function on ptr and the internal 
            auto copy = capi::PyObject_GetAttrString(vs,"copy");
            capi::PyObject_CallObject2(
                copy,
                ptr,
                data,
                __LOC__
                    + ", evaluation of the vector space function copy failed");
        }
            
        // Convert a C++ state to a Python state 
        template <>
        void State <PyUnconstrained>::toPython(
            typename PyUnconstrained::State::t const & state
        ) {
            Unconstrained::State::toPython(state,*this);
        }
        template <>
        void State <PyEqualityConstrained>::toPython(
            typename PyEqualityConstrained::State::t const & state
        ) {
            EqualityConstrained::State::toPython(state,*this);
        }
        template <>
        void State <PyInequalityConstrained>::toPython(
            typename PyInequalityConstrained::State::t const & state
        ) {
            InequalityConstrained::State::toPython(state,*this);
        }
        template <>
        void State <PyConstrained>::toPython(
            typename PyConstrained::State::t const & state
        ) {
            Constrained::State::toPython(state,*this);
        }

        // Convert a Python state to C++ 
        template <>
        void State <PyUnconstrained>::fromPython(
            typename PyUnconstrained::State::t & state
        ) {
            Unconstrained::State::fromPython(data,state);
        }
        template <>
        void State <PyEqualityConstrained>::fromPython(
            typename PyEqualityConstrained::State::t & state
        ) {
            EqualityConstrained::State::fromPython(data,state);
        }
        template <>
        void State <PyInequalityConstrained>::fromPython(
            typename PyInequalityConstrained::State::t & state
        ) {
            InequalityConstrained::State::fromPython(data,state);
        }
        template <>
        void State <PyConstrained>::fromPython(
            typename PyConstrained::State::t & state
        ) {
            Constrained::State::fromPython(data,state);
        }
        
        // Convert a Python bundle to C++ 
        template <>
        void Functions <PyUnconstrained>::fromPython(
            typename PyUnconstrained::Functions::t & fns 
        ) {
            Unconstrained::Functions::fromPython(
                *this,pystate,state,fns);
        }
        template <>
        void Functions <PyEqualityConstrained>::fromPython(
            typename PyEqualityConstrained::Functions::t & fns 
        ) {
            EqualityConstrained::Functions::fromPython(
                *this,pystate,state,fns);
        }
        template <>
        void Functions <PyInequalityConstrained>::fromPython(
            typename PyInequalityConstrained::Functions::t & fns 
        ) {
            InequalityConstrained::Functions::fromPython(
                *this,pystate,state,fns);
        }
        template <>
        void Functions <PyConstrained>::fromPython(
            typename PyConstrained::Functions::t & fns 
        ) {
            Constrained::Functions::fromPython(
                *this,pystate,state,fns);
        }

        // Create a function 
        ScalarValuedFunction::ScalarValuedFunction(PyObjectPtr const & data_)
            : data(data_) {}

        // <- f(x) 
        double ScalarValuedFunction::eval(Vector const & x) const { 
            // Call the objective function on x.  Store in z.
            auto eval = capi::PyObject_GetAttrString(data,"eval");
            auto ret = capi::PyObject_CallObject1(
                eval,
                x.data,
                __LOC__
                    + ", evaluation of the objective f failed");

            // Return the result
            return capi::PyFloat_AsDouble(ret);
        }

        // grad = grad f(x) 
        void ScalarValuedFunction::grad(
            Vector const & x,
            Vector & grad
        ) const { 
            // Call the gradient function on x and grad. 
            auto pygrad = capi::PyObject_GetAttrString(data,"grad");
            capi::PyObject_CallObject2(
                pygrad,
                x.data,
                grad.data,
                __LOC__
                    + ", evaluation of the gradient of f failed");
        }

        // H_dx = hess f(x) dx 
        void ScalarValuedFunction::hessvec(
            Vector const & x,
            Vector const & dx,
            Vector & H_dx
        ) const {
            // Call the hessvec function on x, dx, and H_dx.
            auto hessvec = capi::PyObject_GetAttrString(data,"hessvec");
            capi::PyObject_CallObject3(
                hessvec,
                x.data,
                dx.data,
                H_dx.data,
                __LOC__
                    + ", evaluation of the Hessian-vector product of f failed");
        }

        // Create a function 
        VectorValuedFunction::VectorValuedFunction(
            std::string const & name_,
            PyObjectPtr const & data_
        ) :
            name(name_),
            data(data_)
        {}

        // y=f(x)
        void VectorValuedFunction::eval(
            X_Vector const & x,
            VectorValuedFunction::Y_Vector& y
        ) const {
            // Call the evaluate function on x and y.
            auto eval = capi::PyObject_GetAttrString(data,"eval");
            capi::PyObject_CallObject2(
                eval,
                x.data,
                y.data,
                __LOC__
                    + ", evaluation of the constraint " + name + " failed");
        }

        // y=f'(x)dx 
        void VectorValuedFunction::p(
            X_Vector const & x,
            X_Vector const & dx,
            Y_Vector & y
        ) const {
            // Call the prime function on x, dx, and y
            auto p = capi::PyObject_GetAttrString(data,"p");
            capi::PyObject_CallObject3(
                p,
                x.data,
                dx.data,
                y.data,
                __LOC__
                    + ", evaluation of the derivative of the constraint "
                    + name + " failed");
        }

        // xhat=f'(x)*dy
        void VectorValuedFunction::ps(
            X_Vector const & x,
            Y_Vector const & dy,
            X_Vector & xhat 
        ) const {
            // Call the prime-adjoint function on x, dy, and z
            auto ps = capi::PyObject_GetAttrString(data,"ps");
            capi::PyObject_CallObject3(
                ps,
                x.data,
                dy.data,
                xhat.data,
                __LOC__
                    +", evaluation of the derivative-adjoint of the constraint "
                    + name + " failed");
        }
             
        // xhat=(f''(x)dx)*dy
        void VectorValuedFunction::pps(
            X_Vector const & x,
            X_Vector const & dx,
            Y_Vector const & dy,
            X_Vector & xhat 
        ) const { 
            // Call the prime-adjoint function on x, dx, dy, and z
            auto pps = capi::PyObject_GetAttrString(data,"pps");
            capi::PyObject_CallObject4(
                pps,
                x.data,
                dx.data,
                dy.data,
                xhat.data,
                __LOC__
                    + ", evaluation of the second derivative-adjoint of the "
                    "constraint " + name + " failed");
        }
        
        // Converts elements from C++ to Python 
        namespace toPython {
                    
            // Sets a real in a Python state 
            void Real(
                std::string const & name,
                double const & value,
                PyObjectPtr & pystate 
            ) {
                auto item = capi::PyFloat_FromDouble(value);
                capi::PyObject_SetAttrString(pystate,name.c_str(),item);
            }
        
            // Sets a natural in a Python state 
            void Natural(
                std::string const & name,
                Optizelle::Natural const & value,
                PyObjectPtr & pystate 
            ) {
                auto item = capi::PyInt_FromNatural(value);
                capi::PyObject_SetAttrString(pystate,name.c_str(),item);
            }
        
            // Sets a vector in a Python state 
            void Vector(
                std::string const & name,
                Python::Vector const & value,
                PyObjectPtr & pystate 
            ) {
                auto item = capi::PyObject_GetAttrString(pystate,name.c_str());
                value.toPython(item);
            }
        
            // Sets a list of vectors in a Python state 
            void VectorList(
                std::string const & name,
                std::list <Python::Vector> const & values,
                PyObjectPtr & pystate 
            ) {
                // Create a new Python list that we insert elements into
                auto items = capi::PyList_New(0);

                // Loop over all of the items inside values and then insert 
                // them into items 
                for(auto const & value : values) {
                    // Allocate memory for a new vector
                    auto item = value.init();

                    // Copy the information from the current iterator into this
                    // new vector
                    item.copy(value);

                    // Append this item to the list 
                    capi::PyList_Append(items,item.data);
                }
                
                // Insert the items into pystate 
                capi::PyObject_SetAttrString(pystate,name.c_str(),items);
            }

            // Sets restart vectors in Python 
            void Vectors(
                Python::Vectors const & values,
                PyObjectPtr & pyvalues 
            ) {
            
                // Loop over all of the items inside values and then insert 
                // them into pyvalues 
                for(auto const & value : values) {
                    // Allocate memory for a new vector
                    auto pyvalue = value.second.init();

                    // Copy the information from the current iterator into this
                    // new vector
                    pyvalue.copy(value.second);

                    // Release the pointer into the Python list
                    capi::PyList_Append(
                        pyvalues,
                        capi::PyTuple_Pack_2(
                            capi::PyString_FromString(value.first.c_str()),
                            pyvalue.data));
                }
            }
        
            // Sets restart reals in Python 
            void Reals(
                Python::Reals const & values,
                PyObjectPtr & pyvalues 
            ) {
                // Loop over all of the items inside values and then insert 
                // them into pyvalues 
                for(auto const & value : values) {
                    // Insert the double into the Python list 
                    capi::PyList_Append(
                        pyvalues,
                        capi::PyTuple_Pack_2(
                            capi::PyString_FromString(value.first.c_str()),
                            capi::PyFloat_FromDouble(value.second)));
                }
            }
        
            // Converts a list of naturals to a Python list 
            void Naturals(
                Python::Naturals const & values,
                PyObjectPtr & pyvalues 
            ) {
                // Loop over all of the items inside values and then insert 
                // them into pyvalues 
                for(auto const & value : values) {
                    // Insert the double into the Python list 
                    capi::PyList_Append(
                        pyvalues,
                        capi::PyTuple_Pack_2(
                            capi::PyString_FromString(value.first.c_str()),
                            capi::PyInt_FromNatural(value.second)));
                }
            }
        
            // Sets restart parameters in Python 
            void Params(
                Python::Params const & values,
                PyObjectPtr & pyvalues 
            ) {
                // Loop over all of the items inside values and then insert 
                // them into pyvalues 
                for(auto const & value : values) {
                    // Insert the double into the Python list 
                    capi::PyList_Append(
                        pyvalues,
                        capi::PyTuple_Pack_2(
                            capi::PyString_FromString(value.first.c_str()),
                            capi::PyString_FromString(value.second.c_str())));
                }
            }
        }
        
        // Converts elements from Python to C++ 
        namespace fromPython {
        
            // Sets a real in a C++ state 
            void Real(
                std::string const & name,
                PyObjectPtr const & pystate,
                double & value
            ) {
                auto item = capi::PyObject_GetAttrString(pystate,name.c_str());
                value=capi::PyFloat_AsDouble(item);
            }
            
            // Sets a natural in a C++ state 
            void Natural(
                std::string const & name,
                PyObjectPtr const & pystate,
                Optizelle::Natural & value
            ) {
                auto item = capi::PyObject_GetAttrString(pystate,name.c_str());
                value=capi::PyInt_AsNatural(item);
            }
            
            // Sets a vector in a C++ state 
            void Vector(
                std::string const & name,
                PyObjectPtr const & pystate,
                Python::Vector & value
            ) {
                auto item = capi::PyObject_GetAttrString(pystate,name.c_str());
                value.fromPython(item);
            }
            
            // Sets a list of vectors in a C++ state 
            void VectorList(
                std::string const & name,
                PyObjectPtr const & pystate,
                Python::Vector const & vec,
                std::list <Python::Vector> & values
            ) {
                // Grab the list of items
                auto items = capi::PyObject_GetAttrString(pystate,name.c_str());

                // Loop over all the elements in items and insert them one
                // at a time into values
                values.clear();
                for(auto i=Optizelle::Natural(0);
                    i<capi::PyList_Size(items);
                    i++
                ) {
                    // Grab the current item from Python
                    auto item = capi::PyList_GetItem(items,i);

                    // Create a new vector in values 
                    values.emplace_back(vec.init());

                    // Copy the Python item into the new value
                    values.back().fromPython(item);
                }
            }
        
            // Sets a scalar-valued function in a C++ function bundle 
            void ScalarValuedFunction(
                std::string const & name,
                PyObjectPtr const & pyfns,
                std::unique_ptr <PyScalarValuedFunction> & value
            ) {
                value.reset(
                    new Python::ScalarValuedFunction(
                        capi::PyObject_GetAttrString(
                            pyfns,
                            name.c_str())));
            }
            
            // Sets a vector-valued function in a C++ function bundle 
            void VectorValuedFunction(
                std::string const & name,
                PyObjectPtr const & pyfns,
                std::unique_ptr <PyVectorValuedFunction> & value
            ) {
                value.reset(
                    new Python::VectorValuedFunction(name,
                        capi::PyObject_GetAttrString(
                            pyfns,
                            name.c_str())));
            }
        
            // Sets restart vectors in C++ 
            void Vectors(
                Python::Vector const & vec,
                PyObjectPtr const & pyvalues,
                Python::Vectors & values
            ) {
                // Loop over all the elements in pyvalues and insert them one
                // at a time into values
                values.clear();
                for(auto i = Optizelle::Natural(0);
                    i<capi::PyList_Size(pyvalues);
                    i++
                ) {
                    // Grab the current item from Python
                    auto pyvalue = capi::PyList_GetItem(pyvalues,i);

                    // Grab the values
                    auto name = capi::PyTuple_GetItem(pyvalue,0);
                    auto val = capi::PyTuple_GetItem(pyvalue,1);

                    // Create the elements in values 
                    values.emplace_back(
                        capi::PyString_AsString(name),
                        vec.init());

                    // Copy the Python value into the C++ value
                    values.back().second.fromPython(val);
                }
            }
            
            // Sets restart reals in C++ 
            void Reals(
                PyObjectPtr const & pyvalues,
                Python::Reals & values
            ) {
                // Loop over all the elements in pyvalues and insert them one
                // at a time into values
                values.clear();
                for(auto i = Optizelle::Natural(0);
                    i<capi::PyList_Size(pyvalues);
                    i++
                ) {
                    // Grab the current item from Python
                    auto pyvalue = capi::PyList_GetItem(pyvalues,i);

                    // Grab the values
                    auto name = capi::PyTuple_GetItem(pyvalue,0);
                    auto val = capi::PyTuple_GetItem(pyvalue,1);
                    
                    // Create the elements in values 
                    values.emplace_back(
                        capi::PyString_AsString(name),
                        capi::PyFloat_AsDouble(val));
                }
            }
            
            // Sets restart naturals in C++ 
            void Naturals(
                PyObjectPtr const & pyvalues,
                Python::Naturals & values
            ) {
                // Loop over all the elements in pyvalues and insert them one
                // at a time into values
                values.clear();
                for(auto i = Optizelle::Natural(0);
                    i<capi::PyList_Size(pyvalues);
                    i++
                ) {
                    // Grab the current item from Python
                    auto pyvalue = capi::PyList_GetItem(pyvalues,i);

                    // Grab the values
                    auto name = capi::PyTuple_GetItem(pyvalue,0);
                    auto val = capi::PyTuple_GetItem(pyvalue,1);
                    
                    // Create the elements in values 
                    values.emplace_back(
                        capi::PyString_AsString(name),
                        capi::PyInt_AsNatural(val));
                }
            }
            
            // Sets restart parameters in C++ 
            void Params(
                PyObjectPtr const & pyvalues,
                Python::Params & values
            ) {
                // Loop over all the elements in pyvalues and insert them one
                // at a time into values
                values.clear();
                for(auto i = Optizelle::Natural(0);
                    i<capi::PyList_Size(pyvalues);
                    i++
                ) {
                    // Grab the current item from Python
                    auto pyvalue = capi::PyList_GetItem(pyvalues,i);

                    // Grab the values
                    auto name = capi::PyTuple_GetItem(pyvalue,0);
                    auto val = capi::PyTuple_GetItem(pyvalue,1);
                    
                    // Create the elements in values 
                    values.emplace_back(
                        capi::PyString_AsString(name),
                        capi::PyString_AsString(val));
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
                    PyObjectPtr & pystate
                ){
                    // Set each of the required items in the Python state
                    toPython::Real("eps_grad",state.eps_grad,pystate);
                    toPython::Real("eps_dx",state.eps_dx,pystate);
                    toPython::Natural("stored_history",
                        state.stored_history,pystate);
                    toPython::Natural("iter",state.iter,pystate);
                    toPython::Natural("iter_max",state.iter_max,pystate);
                    toPython::Natural("glob_iter",state.glob_iter,pystate);
                    toPython::Natural("glob_iter_max",
                        state.glob_iter_max,pystate);
                    toPython::Natural("glob_iter_total",
                        state.glob_iter_total,pystate);
                    toPython::Param <OptimizationStop::t> (
                        "opt_stop",
                        OptimizationStop::toPython,
                        state.opt_stop,
                        pystate);
                    toPython::Natural("trunc_iter",state.trunc_iter,pystate);
                    toPython::Natural("trunc_iter_max",
                        state.trunc_iter_max,pystate);
                    toPython::Natural("trunc_iter_total",
                        state.trunc_iter_total,pystate);
                    toPython::Natural("trunc_orthog_storage_max",
                        state.trunc_orthog_storage_max,pystate);
                    toPython::Natural("trunc_orthog_iter_max",
                        state.trunc_orthog_iter_max,pystate);
                    toPython::Param <TruncatedStop::t> (
                        "trunc_stop",
                        TruncatedStop::toPython,
                        state.trunc_stop,
                        pystate);
                    toPython::Real("trunc_err",
                        state.trunc_err,pystate);
                    toPython::Real("eps_trunc",state.eps_trunc,pystate);
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
                    toPython::Natural("safeguard_failed_max",
                        state.safeguard_failed_max,pystate);
                    toPython::Natural("safeguard_failed",
                        state.safeguard_failed,pystate);
                    toPython::Natural("safeguard_failed_total",
                        state.safeguard_failed_total,pystate);
                    toPython::Real("alpha_x",state.alpha_x,pystate);
                    toPython::Real("alpha_x_qn",state.alpha_x_qn,pystate);
                    toPython::Real("delta",state.delta,pystate);
                    toPython::Real("eta1",state.eta1,pystate);
                    toPython::Real("eta2",state.eta2,pystate);
                    toPython::Real("ared",state.ared,pystate);
                    toPython::Real("pred",state.pred,pystate);
                    toPython::Real("alpha0",state.alpha0,pystate);
                    toPython::Real("alpha",state.alpha,pystate);
                    toPython::Real("c1",state.c1,pystate);
                    toPython::Natural("ls_iter",
                        state.ls_iter,pystate);
                    toPython::Natural("ls_iter_max",
                        state.ls_iter_max,pystate);
                    toPython::Natural("ls_iter_total",
                        state.ls_iter_total,pystate);
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
                    toPython::Param <ToleranceKind::t> (
                        "eps_kind",
                        ToleranceKind::toPython,
                        state.eps_kind,
                        pystate);
                }
                void toPython(
                    typename PyUnconstrained::State::t const & state,
                    Python::State <PyUnconstrained> & pystate
                ){
                    Unconstrained::State::toPython_(state,pystate.data);
                }
                
                // Convert a Python state to C++ 
                void fromPython_(
                    PyObjectPtr const & pystate,
                    typename PyUnconstrained::State::t & state
                ){
                    // Set each of the required items in the Python state
                    fromPython::Real("eps_grad",pystate,state.eps_grad);
                    fromPython::Real("eps_dx",pystate,state.eps_dx);
                    fromPython::Natural("stored_history",
                        pystate,state.stored_history);
                    fromPython::Natural("iter",pystate,state.iter);
                    fromPython::Natural("iter_max",pystate,state.iter_max);
                    fromPython::Natural("glob_iter",
                        pystate,state.glob_iter);
                    fromPython::Natural("glob_iter_max",
                        pystate,state.glob_iter_max);
                    fromPython::Natural("glob_iter_total",
                        pystate,state.glob_iter_total);
                    fromPython::Param <OptimizationStop::t> (
                        "opt_stop",
                        OptimizationStop::fromPython,
                        pystate,
                        state.opt_stop);
                    fromPython::Natural("trunc_iter",
                        pystate,state.trunc_iter);
                    fromPython::Natural("trunc_iter_max",
                        pystate,state.trunc_iter_max);
                    fromPython::Natural("trunc_iter_total",
                        pystate,state.trunc_iter_total);
                    fromPython::Natural("trunc_orthog_storage_max",
                        pystate,state.trunc_orthog_storage_max);
                    fromPython::Natural("trunc_orthog_iter_max",
                        pystate,state.trunc_orthog_iter_max);
                    fromPython::Param <TruncatedStop::t> (
                        "trunc_stop",
                        TruncatedStop::fromPython,
                        pystate,
                        state.trunc_stop);
                    fromPython::Real("trunc_err",
                        pystate,state.trunc_err);
                    fromPython::Real("eps_trunc",pystate,state.eps_trunc);
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
                    fromPython::Natural("safeguard_failed_max",
                        pystate,state.safeguard_failed_max);
                    fromPython::Natural("safeguard_failed",
                        pystate,state.safeguard_failed);
                    fromPython::Natural("safeguard_failed_total",
                        pystate,state.safeguard_failed_total);
                    fromPython::Real("alpha_x",pystate,state.alpha_x);
                    fromPython::Real("alpha_x_qn",pystate,state.alpha_x_qn);
                    fromPython::Real("delta",pystate,state.delta);
                    fromPython::Real("eta1",pystate,state.eta1);
                    fromPython::Real("eta2",pystate,state.eta2);
                    fromPython::Real("ared",pystate,state.ared);
                    fromPython::Real("pred",pystate,state.pred);
                    fromPython::Real("alpha0",pystate,state.alpha0);
                    fromPython::Real("alpha",pystate,state.alpha);
                    fromPython::Real("c1",pystate,state.c1);
                    fromPython::Natural("ls_iter",
                        pystate,state.ls_iter);
                    fromPython::Natural("ls_iter_max",
                        pystate,state.ls_iter_max);
                    fromPython::Natural("ls_iter_total",pystate,
                        state.ls_iter_total);
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
                    fromPython::Param <ToleranceKind::t> (
                        "eps_kind",
                        ToleranceKind::fromPython,
                        pystate,
                        state.eps_kind);
                }
                void fromPython(
                    Python::State <PyUnconstrained> const & pystate,
                    typename PyUnconstrained::State::t & state
                ){
                    Unconstrained::State::fromPython_(pystate.data,state);
                }

                // Creates a state and inserts the elements into pystate 
                PyObject * create(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables 
                    PY_VAR_3(pystate,X,x);

                    // Create a vector from the user input
                    auto x = Vector(X_,x_);

                    // Create a Python state 
                    auto pystate = Python::State <PyUnconstrained> (pystate_);

                    // Create a new C++ state
                    typename PyUnconstrained::State::t state(x);

                    // Convert the state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE;

                // Catch errors 
                } CATCH_PYTHON_ERRORS;

                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables 
                    PY_VAR_3(X,fname,pystate);

                    // Grab the file name
                    auto fname = capi::PyString_AsString(fname_);

                    // Create a Python state 
                    auto pystate = Python::State <PyUnconstrained> (pystate_);
                
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);

                    // Create a new C++ state
                    typename PyUnconstrained::State::t state(x);
                    
                    // Convert the Python state to a C++ state
                    pystate.fromPython(state);

                    // Read the JSON file into the C++ state
                    PyJsonUnconstrained::read(fname,state);

                    // Convert the C++ state to a Python state
                    pystate.toPython(state);
                            
                    // Return nothing 
                    Py_RETURN_NONE;
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Python bundle to C++ 
                void fromPython(
                    Python::Functions <PyUnconstrained> const & pyfns,
                    Python::State <PyUnconstrained> & pystate,
                    typename PyUnconstrained::State::t const & state,
                    typename PyUnconstrained::Functions::t & fns 
                ) {
                    Unconstrained::Functions::fromPython_
                        <PyUnconstrained> (pyfns,pystate,state,fns);
                }
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem 
                PyObject * getMin(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables 
                    PY_VAR_5(X,msg,pyfns,pystate,smanip);

                    // Create a messaging object
                    auto msg = Optizelle::Python::Messaging::python(msg_);
                        
                    // Create a Python state 
                    auto pystate = Python::State <PyUnconstrained> (pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);

                    // Create a C++ state
                    typename PyUnconstrained::State::t state(x);
                    
                    // Convert the Python state to a C++ state
                    pystate.fromPython(state);

                    // Create a Python bundle of functions
                    auto pyfns = Python::Functions <PyUnconstrained>(
                        pystate,
                        state,
                        pyfns_);

                    // Create a C++ bundle of functions
                    typename PyUnconstrained::Functions::t fns;
                    
                    // Convert the Python bundle of functions to C++ 
                    pyfns.fromPython(fns);
                    
                    // Create a state manipulator 
                    auto smanip = Python::StateManipulator <PyUnconstrained>(
                        pystate,
                        pyfns,
                        smanip_);
                   
                    // Minimize
                    PyUnconstrained::Algorithms::getMin(
                        msg,fns,state,smanip);
                    
                    // Convert the C++ state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE;
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;
            }
            
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                PyObject * release(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables
                    PY_VAR_6(X,pystate,pyxs,pyreals,pynats,pyparams);

                    // Create a Python state 
                    auto pystate = Python::State <PyUnconstrained>(pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);

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
                    toPython::Vectors(xs,pyxs_);
                    toPython::Reals(reals,pyreals_);
                    toPython::Naturals(nats,pynats_);
                    toPython::Params(params,pyparams_);

                    // Return nothing 
                    Py_RETURN_NONE;
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;

                // Capture data from structures controlled by the user.  
                PyObject * capture(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables
                    PY_VAR_6(X,pystate,pyxs,pyreals,pynats,pyparams);

                    // Create a Python state 
                    auto pystate = Python::State <PyUnconstrained> (pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);

                    // Create a C++ state
                    typename PyUnconstrained::State::t state(x);
                   
                    // Allocate memory for the released vectors
                    PyUnconstrained::Restart::X_Vectors xs;
                    PyUnconstrained::Restart::Reals reals;
                    PyUnconstrained::Restart::Naturals nats;
                    PyUnconstrained::Restart::Params params;
                    
                    // Convert the restart information from Python 
                    fromPython::Vectors(x,pyxs_,xs);
                    fromPython::Reals(pyreals_,reals);
                    fromPython::Naturals(pynats_,nats);
                    fromPython::Params(pyparams_,params);

                    // Do a capture 
                    PyUnconstrained::Restart
                        ::capture(state,xs,reals,nats,params);

                    // Convert the C++ state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE;
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;

                // Writes a json restart file
                PyObject * write_restart(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables 
                    PY_VAR_3(X,fname,pystate);

                    // Grab the file name
                    auto fname = capi::PyString_AsString(fname_);

                    // Create a Python state 
                    auto pystate = Python::State <PyUnconstrained> (pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    
                    // Create a C++ state
                    typename PyUnconstrained::State::t state(x);
                    
                    // Convert Python state to C++ 
                    pystate.fromPython(state);

                    // Write the restart file
                    PyJsonUnconstrained::write_restart(fname,state);
                    
                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;

                // Reads a json restart file
                PyObject * read_restart(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables
                    PY_VAR_4(X,fname,x,pystate);

                    // Grab the file name
                    auto fname = capi::PyString_AsString(fname_);

                    // Create a Python state 
                    auto pystate = Python::State <PyUnconstrained> (pystate_);
                    
                    // Grab the reference vector 
                    auto x = Vector(X_,x_);
                    
                    // Create a C++ state
                    typename PyUnconstrained::State::t state(x);

                    // Read the restart file into the C++ state 
                    PyJsonUnconstrained::read_restart(fname,x,state);
                    
                    // Convert the C++ state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;
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
                    PyObjectPtr & pystate
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
                    toPython::Natural("augsys_qn_failed",
                        state.augsys_qn_failed,pystate);
                    toPython::Natural("augsys_pg_failed",
                        state.augsys_pg_failed,pystate);
                    toPython::Natural("augsys_proj_failed",
                        state.augsys_proj_failed,pystate);
                    toPython::Natural("augsys_tang_failed",
                        state.augsys_tang_failed,pystate);
                    toPython::Natural("augsys_lmh_failed",
                        state.augsys_lmh_failed,pystate);
                    toPython::Natural("augsys_failed_total",
                        state.augsys_failed_total,pystate);
                    toPython::Vector("g_x",state.g_x,pystate);
                    toPython::Real("norm_gxtyp",state.norm_gxtyp,pystate);
                    toPython::Real("norm_gpsgxtyp",state.norm_gpsgxtyp,pystate);
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
                    toPython::Param <QuasinormalStop::t> (
                        "qn_stop",
                        QuasinormalStop::toPython,
                        state.qn_stop,
                        pystate);
                }
                void toPython(
                    typename PyEqualityConstrained::State::t const & state,
                    Python::State <PyEqualityConstrained> & pystate
                ){
                    Unconstrained::State::toPython_(state,pystate.data);
                    EqualityConstrained::State::toPython_(state,pystate.data);
                }
                
                // Convert a Python state to C++ 
                void fromPython_(
                    PyObjectPtr const & pystate,
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
                    fromPython::Natural("augsys_qn_failed",
                        pystate,state.augsys_qn_failed);
                    fromPython::Natural("augsys_pg_failed",
                        pystate,state.augsys_pg_failed);
                    fromPython::Natural("augsys_proj_failed",
                        pystate,state.augsys_proj_failed);
                    fromPython::Natural("augsys_tang_failed",
                        pystate,state.augsys_tang_failed);
                    fromPython::Natural("augsys_lmh_failed",
                        pystate,state.augsys_lmh_failed);
                    fromPython::Natural("augsys_failed_total",
                        pystate,state.augsys_failed_total);
                    fromPython::Vector("g_x",pystate,state.g_x);
                    fromPython::Real("norm_gxtyp",pystate,state.norm_gxtyp);
                    fromPython::Real("norm_gpsgxtyp",
                        pystate,state.norm_gpsgxtyp);
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
                    fromPython::Param <QuasinormalStop::t> (
                        "qn_stop",
                        QuasinormalStop::fromPython,
                        pystate,
                        state.qn_stop);
                }
                void fromPython(
                    PyObjectPtr const & pystate,
                    typename PyEqualityConstrained::State::t & state
                ){
                    Unconstrained::State::fromPython_(pystate,state);
                    EqualityConstrained::State::fromPython_(pystate,state);
                }
        
                // Creates a state and inserts the elements into pystate 
                PyObject * create(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables 
                    PY_VAR_5(pystate,X,Y,x,y);

                    // Create vectors from the user input
                    auto x = Vector(X_,x_);
                    auto y = Vector(Y_,y_);

                    // Create a Python state 
                    auto pystate = Python::State <PyEqualityConstrained> (
                        pystate_);

                    // Create a new C++ state
                    typename PyEqualityConstrained::State::t state(x,y);

                    // Convert the state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;
        
                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables 
                    PY_VAR_4(X,Y,fname,pystate);

                    // Grab the file name
                    auto fname = capi::PyString_AsString(fname_);

                    // Create a Python state 
                    auto pystate = Python::State <PyEqualityConstrained> (
                        pystate_);
                
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    auto y_ = capi::PyObject_GetAttrString(pystate.data,"y");
                    auto y = Vector(Y_,y_);

                    // Create a new C++ state
                    typename PyEqualityConstrained::State::t state(x,y);
                    
                    // Convert the Python state to a C++ state
                    pystate.fromPython(state);

                    // Read the JSON file into the C++ state
                    PyJsonEqualityConstrained::read(fname,state);

                    // Convert the C++ state to a Python state
                    pystate.toPython(state);
                            
                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Python bundle to C++ 
                void fromPython(
                    Python::Functions <PyEqualityConstrained> const & pyfns,
                    Python::State <PyEqualityConstrained> & pystate,
                    typename PyEqualityConstrained::State::t const & state,
                    typename PyEqualityConstrained::Functions::t & fns 
                ) {
                    Unconstrained::Functions::fromPython_
                        <PyEqualityConstrained> (pyfns,pystate,state,fns);
                    EqualityConstrained::Functions::fromPython_
                        <PyEqualityConstrained> (pyfns,pystate,state,fns);
                }
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem 
                PyObject * getMin(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables 
                    PY_VAR_6(X,Y,msg,pyfns,pystate,smanip);

                    // Create a messaging object
                    auto msg = Optizelle::Python::Messaging::python(msg_);
                        
                    // Create a Python state 
                    auto pystate = Python::State <PyEqualityConstrained> (
                        pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    auto y_ = capi::PyObject_GetAttrString(pystate.data,"y");
                    auto y = Vector(Y_,y_);

                    // Create a C++ state
                    typename PyEqualityConstrained::State::t state(x,y);
                    
                    // Convert the Python state to a C++ state
                    pystate.fromPython(state);

                    // Create a Python bundle of functions
                    auto pyfns = Python::Functions <PyEqualityConstrained>(
                        pystate,
                        state,
                        pyfns_);

                    // Create a C++ bundle of functions
                    typename PyEqualityConstrained::Functions::t fns;
                    
                    // Convert the Python bundle of functions to C++ 
                    pyfns.fromPython(fns);
                    
                    // Create a state manipulator 
                    auto smanip =
                        Python::StateManipulator <PyEqualityConstrained>(
                        pystate,
                        pyfns,
                        smanip_);
                   
                    // Minimize
                    PyEqualityConstrained::Algorithms::getMin(
                        msg,fns,state,smanip);
                    
                    // Convert the C++ state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;
            }
            
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                PyObject * release(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables
                    PY_VAR_8(X,Y,pystate,pyxs,pyys,pyreals,pynats,pyparams);

                    // Create a Python state 
                    auto pystate = Python::State <PyEqualityConstrained> (
                        pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    auto y_ = capi::PyObject_GetAttrString(pystate.data,"y");
                    auto y = Vector(Y_,y_);

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
                    toPython::Vectors(xs,pyxs_);
                    toPython::Vectors(ys,pyys_);
                    toPython::Reals(reals,pyreals_);
                    toPython::Naturals(nats,pynats_);
                    toPython::Params(params,pyparams_);

                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;

                // Capture data from structures controlled by the user.  
                PyObject * capture(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables
                    PY_VAR_8(X,Y,pystate,pyxs,pyys,pyreals,pynats,pyparams);

                    // Create a Python state 
                    auto pystate = Python::State <PyEqualityConstrained> (
                        pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    auto y_ = capi::PyObject_GetAttrString(pystate.data,"y");
                    auto y = Vector(Y_,y_);

                    // Create a C++ state
                    typename PyEqualityConstrained::State::t state(x,y);
                   
                    // Allocate memory for the released vectors
                    PyEqualityConstrained::Restart::X_Vectors xs;
                    PyEqualityConstrained::Restart::Y_Vectors ys;
                    PyEqualityConstrained::Restart::Reals reals;
                    PyEqualityConstrained::Restart::Naturals nats;
                    PyEqualityConstrained::Restart::Params params;
                    
                    // Convert the restart information from Python 
                    fromPython::Vectors(x,pyxs_,xs);
                    fromPython::Vectors(y,pyys_,ys);
                    fromPython::Reals(pyreals_,reals);
                    fromPython::Naturals(pynats_,nats);
                    fromPython::Params(pyparams_,params);

                    // Do a capture 
                    PyEqualityConstrained::Restart
                        ::capture(state,xs,ys,reals,nats,params);

                    // Convert the C++ state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;

                // Writes a json restart file
                PyObject * write_restart(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables 
                    PY_VAR_4(X,Y,fname,pystate);

                    // Grab the file name
                    auto fname = capi::PyString_AsString(fname_);

                    // Create a Python state 
                    auto pystate = Python::State <PyEqualityConstrained> (
                        pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    auto y_ = capi::PyObject_GetAttrString(pystate.data,"y");
                    auto y = Vector(Y_,y_);
                    
                    // Create a C++ state
                    typename PyEqualityConstrained::State::t state(x,y);
                    
                    // Convert Python state to C++ 
                    pystate.fromPython(state);

                    // Write the restart file
                    PyJsonEqualityConstrained::write_restart(fname,state);
                    
                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;

                // Reads a json restart file
                PyObject * read_restart(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables
                    PY_VAR_6(X,Y,fname,x,y,pystate);

                    // Grab the file name
                    auto fname = capi::PyString_AsString(fname_);

                    // Create a Python state 
                    auto pystate = Python::State <PyEqualityConstrained> (
                        pystate_);
                    
                    // Grab the reference vector 
                    auto x = Vector(X_,x_);
                    auto y = Vector(Y_,y_);
                    
                    // Create a C++ state
                    typename PyEqualityConstrained::State::t state(x,y);

                    // Read the restart file into the C++ state 
                    PyJsonEqualityConstrained::read_restart(
                        fname,x,y,state);
                    
                    // Convert the C++ state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;
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
                    PyObjectPtr & pystate
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
                    Python::State <PyInequalityConstrained> & pystate
                ){
                    Unconstrained::State::toPython_(state,pystate.data);
                    InequalityConstrained::State::toPython_(state,pystate.data);
                }
                
                // Convert a Python state to C++ 
                void fromPython_(
                    PyObjectPtr const & pystate,
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
                    PyObjectPtr const & pystate,
                    typename PyInequalityConstrained::State::t & state
                ){
                    Unconstrained::State::fromPython_(pystate,state);
                    InequalityConstrained::State::fromPython_(pystate,state);
                }

                // Creates a state and inserts the elements into pystate 
                PyObject * create( 
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables 
                    PY_VAR_5(pystate,X,Z,x,z);

                    // Create vectors from the user input
                    auto x = Vector(X_,x_);
                    auto z = Vector(Z_,z_);

                    // Create a Python state 
                    auto pystate = Python::State <PyInequalityConstrained> (
                        pystate_);

                    // Create a new C++ state
                    typename PyInequalityConstrained::State::t state(x,z);

                    // Convert the state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;

                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables 
                    PY_VAR_4(X,Z,fname,pystate);

                    // Grab the file name
                    auto fname = capi::PyString_AsString(fname_);

                    // Create a Python state 
                    auto pystate = Python::State <PyInequalityConstrained> (
                        pystate_);
                
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    auto z_ = capi::PyObject_GetAttrString(pystate.data,"z");
                    auto z = Vector(Z_,z_);

                    // Create a new C++ state
                    typename PyInequalityConstrained::State::t state(x,z);
                    
                    // Convert the Python state to a C++ state
                    pystate.fromPython(state);

                    // Read the JSON file into the C++ state
                    PyJsonInequalityConstrained::read(fname,state);

                    // Convert the C++ state to a Python state
                    pystate.toPython(state);
                            
                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;
            }
            
            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Python bundle to C++ 
                void fromPython(
                    Python::Functions <PyInequalityConstrained> const & pyfns,
                    Python::State <PyInequalityConstrained> & pystate,
                    typename PyInequalityConstrained::State::t const & state,
                    typename PyInequalityConstrained::Functions::t & fns 
                ) {
                    Unconstrained::Functions::fromPython_
                        <PyInequalityConstrained> (pyfns,pystate,state,fns);
                    InequalityConstrained::Functions::fromPython_
                        <PyInequalityConstrained> (pyfns,pystate,state,fns);
                }
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem 
                PyObject * getMin(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables 
                    PY_VAR_6(X,Z,msg,pyfns,pystate,smanip);

                    // Create a messaging object
                    auto msg = Optizelle::Python::Messaging::python(msg_);
                        
                    // Create a Python state 
                    auto pystate = Python::State <PyInequalityConstrained> (
                        pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    auto z_ = capi::PyObject_GetAttrString(pystate.data,"z");
                    auto z = Vector(Z_,z_);

                    // Create a C++ state
                    typename PyInequalityConstrained::State::t state(x,z);
                    
                    // Convert the Python state to a C++ state
                    pystate.fromPython(state);

                    // Create a Python bundle of functions
                    auto pyfns = Python::Functions <PyInequalityConstrained>(
                        pystate,
                        state,
                        pyfns_);

                    // Create a C++ bundle of functions
                    typename PyInequalityConstrained::Functions::t fns;
                    
                    // Convert the Python bundle of functions to C++ 
                    pyfns.fromPython(fns);
                    
                    // Create a state manipulator 
                    auto smanip =
                        Python::StateManipulator <PyInequalityConstrained>(
                        pystate,
                        pyfns,
                        smanip_);
                   
                    // Minimize
                    PyInequalityConstrained::Algorithms::getMin(
                        msg,fns,state,smanip);
                    
                    // Convert the C++ state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;
            }
            
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                PyObject * release(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables
                    PY_VAR_8(X,Z,pystate,pyxs,pyzs,pyreals,pynats,pyparams);

                    // Create a Python state 
                    auto pystate = Python::State <PyInequalityConstrained> (
                        pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    auto z_ = capi::PyObject_GetAttrString(pystate.data,"z");
                    auto z = Vector(Z_,z_);

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
                    toPython::Vectors(xs,pyxs_);
                    toPython::Vectors(zs,pyzs_);
                    toPython::Reals(reals,pyreals_);
                    toPython::Naturals(nats,pynats_);
                    toPython::Params(params,pyparams_);

                    // Return nothing 
                    Py_RETURN_NONE; 

                // Catch errors 
                } CATCH_PYTHON_ERRORS;

                // Capture data from structures controlled by the user.  
                PyObject * capture(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables
                    PY_VAR_8(X,Z,pystate,pyxs,pyzs,pyreals,pynats,pyparams);

                    // Create a Python state 
                    auto pystate = Python::State <PyInequalityConstrained> (
                        pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    auto z_ = capi::PyObject_GetAttrString(pystate.data,"z");
                    auto z = Vector(Z_,z_);

                    // Create a C++ state
                    typename PyInequalityConstrained::State::t state(x,z);
                   
                    // Allocate memory for the released vectors
                    PyInequalityConstrained::Restart::X_Vectors xs;
                    PyInequalityConstrained::Restart::Z_Vectors zs;
                    PyInequalityConstrained::Restart::Reals reals;
                    PyInequalityConstrained::Restart::Naturals nats;
                    PyInequalityConstrained::Restart::Params params;
                    
                    // Convert the restart information from Python 
                    fromPython::Vectors(x,pyxs_,xs);
                    fromPython::Vectors(z,pyzs_,zs);
                    fromPython::Reals(pyreals_,reals);
                    fromPython::Naturals(pynats_,nats);
                    fromPython::Params(pyparams_,params);

                    // Do a capture 
                    PyInequalityConstrained::Restart
                        ::capture(state,xs,zs,reals,nats,params);

                    // Convert the C++ state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE; 

                // Python error
                } CATCH_PYTHON_ERRORS;

                // Writes a json restart file
                PyObject * write_restart(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables 
                    PY_VAR_4(X,Z,fname,pystate);

                    // Grab the file name
                    auto fname = capi::PyString_AsString(fname_);

                    // Create a Python state 
                    auto pystate = Python::State <PyInequalityConstrained> (
                        pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    auto z_ = capi::PyObject_GetAttrString(pystate.data,"z");
                    auto z = Vector(Z_,z_);
                    
                    // Create a C++ state
                    typename PyInequalityConstrained::State::t state(x,z);
                    
                    // Convert Python state to C++ 
                    pystate.fromPython(state);

                    // Write the restart file
                    PyJsonInequalityConstrained::write_restart(fname,state);
                    
                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;
                
                // Reads a json restart file
                PyObject * read_restart(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables
                    PY_VAR_6(X,Z,fname,x,z,pystate);

                    // Grab the file name
                    auto fname = capi::PyString_AsString(fname_);

                    // Create a Python state 
                    auto pystate = Python::State <PyInequalityConstrained> (
                        pystate_);
                    
                    // Grab the reference vector 
                    auto x = Vector(X_,x_);
                    auto z = Vector(Z_,z_);
                    
                    // Create a C++ state
                    typename PyInequalityConstrained::State::t state(x,z);

                    // Read the restart file into the C++ state 
                    PyJsonInequalityConstrained::read_restart(
                        fname,x,z,state);
                    
                    // Convert the C++ state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE; 

                // Catch errors 
                } CATCH_PYTHON_ERRORS;
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
                    Python::State <PyConstrained> & pystate
                ){
                    Unconstrained::State::toPython_(state,pystate.data);
                    EqualityConstrained::State::toPython_(state,pystate.data);
                    InequalityConstrained::State::toPython_(state,pystate.data);
                }
                
                // Convert a Python state to C++ 
                void fromPython(
                    PyObjectPtr const & pystate,
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
                ) try {
                    // Grab variables 
                    PY_VAR_7(pystate,X,Y,Z,x,y,z);

                    // Create vectors from the user input
                    auto x = Vector(X_,x_);
                    auto y = Vector(Y_,y_);
                    auto z = Vector(Z_,z_);

                    // Create a Python state 
                    auto pystate = Python::State <PyConstrained> (pystate_);

                    // Create a new C++ state
                    typename PyConstrained::State::t state(x,y,z);

                    // Convert the state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;

                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables 
                    PY_VAR_5(X,Y,Z,fname,pystate);

                    // Grab the file name
                    auto fname = capi::PyString_AsString(fname_);

                    // Create a Python state 
                    auto pystate = Python::State <PyConstrained> (pystate_);
                
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    auto y_ = capi::PyObject_GetAttrString(pystate.data,"y");
                    auto y = Vector(Y_,y_);
                    auto z_ = capi::PyObject_GetAttrString(pystate.data,"z");
                    auto z = Vector(Z_,z_);

                    // Create a new C++ state
                    typename PyConstrained::State::t state(x,y,z);
                    
                    // Convert the Python state to a C++ state
                    pystate.fromPython(state);

                    // Read the JSON file into the C++ state
                    PyJsonConstrained::read(fname,state);

                    // Convert the C++ state to a Python state
                    pystate.toPython(state);
                            
                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;
            }
            
            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Python bundle to C++ 
                void fromPython(
                    Python::Functions <PyConstrained> const & pyfns,
                    Python::State <PyConstrained> & pystate,
                    typename PyConstrained::State::t const & state,
                    typename PyConstrained::Functions::t & fns 
                ) {
                    Unconstrained::Functions::fromPython_
                        <PyConstrained> (pyfns,pystate,state,fns);
                    EqualityConstrained::Functions::fromPython_
                        <PyConstrained> (pyfns,pystate,state,fns);
                    InequalityConstrained::Functions::fromPython_
                        <PyConstrained> (pyfns,pystate,state,fns);
                }
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem 
                PyObject * getMin(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables 
                    PY_VAR_7(X,Y,Z,msg,pyfns,pystate,smanip);

                    // Create a messaging object
                    auto msg = Optizelle::Python::Messaging::python(msg_);
                        
                    // Create a Python state 
                    auto pystate = Python::State <PyConstrained> (pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    auto y_ = capi::PyObject_GetAttrString(pystate.data,"y");
                    auto y = Vector(Y_,y_);
                    auto z_ = capi::PyObject_GetAttrString(pystate.data,"z");
                    auto z = Vector(Z_,z_);

                    // Create a C++ state
                    typename PyConstrained::State::t state(x,y,z);
                    
                    // Convert the Python state to a C++ state
                    pystate.fromPython(state);

                    // Create a Python bundle of functions
                    auto pyfns = Python::Functions <PyConstrained>(
                        pystate,
                        state,
                        pyfns_);

                    // Create a C++ bundle of functions
                    typename PyConstrained::Functions::t fns;
                    
                    // Convert the Python bundle of functions to C++ 
                    pyfns.fromPython(fns);
                    
                    // Create a state manipulator 
                    auto smanip = Python::StateManipulator <PyConstrained>(
                        pystate,
                        pyfns,
                        smanip_);
                   
                    // Minimize
                    PyConstrained::Algorithms::getMin(
                        msg,fns,state,smanip);
                    
                    // Convert the C++ state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;
            }

            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                PyObject * release(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables
                    PY_VAR_10(X,Y,Z,pystate,pyxs,pyys,pyzs,
                        pyreals,pynats,pyparams);

                    // Create a Python state 
                    auto pystate = Python::State <PyConstrained> (pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    auto y_ = capi::PyObject_GetAttrString(pystate.data,"y");
                    auto y = Vector(Y_,y_);
                    auto z_ = capi::PyObject_GetAttrString(pystate.data,"z");
                    auto z = Vector(Z_,z_);

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
                    toPython::Vectors(xs,pyxs_);
                    toPython::Vectors(ys,pyys_);
                    toPython::Vectors(zs,pyzs_);
                    toPython::Reals(reals,pyreals_);
                    toPython::Naturals(nats,pynats_);
                    toPython::Params(params,pyparams_);

                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;

                // Capture data from structures controlled by the user.  
                PyObject * capture(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables
                    PY_VAR_10(X,Y,Z,pystate,pyxs,pyys,pyzs,
                        pyreals,pynats,pyparams);

                    // Create a Python state 
                    auto pystate = Python::State <PyConstrained> (pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    auto y_ = capi::PyObject_GetAttrString(pystate.data,"y");
                    auto y = Vector(Y_,y_);
                    auto z_ = capi::PyObject_GetAttrString(pystate.data,"z");
                    auto z = Vector(Z_,z_);

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
                    fromPython::Vectors(x,pyxs_,xs);
                    fromPython::Vectors(y,pyys_,ys);
                    fromPython::Vectors(z,pyzs_,zs);
                    fromPython::Reals(pyreals_,reals);
                    fromPython::Naturals(pynats_,nats);
                    fromPython::Params(pyparams_,params);

                    // Do a capture 
                    PyConstrained::Restart
                        ::capture(state,xs,ys,zs,reals,nats,params);

                    // Convert the C++ state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;
                
                // Writes a json restart file
                PyObject * write_restart(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables 
                    PY_VAR_5(X,Y,Z,fname,pystate);

                    // Grab the file name
                    auto fname = capi::PyString_AsString(fname_);

                    // Create a Python state 
                    auto pystate = Python::State <PyConstrained> (pystate_);
                    
                    // Grab the base vectors from the Python state
                    auto x_ = capi::PyObject_GetAttrString(pystate.data,"x");
                    auto x = Vector(X_,x_);
                    auto y_ = capi::PyObject_GetAttrString(pystate.data,"y");
                    auto y = Vector(Y_,y_);
                    auto z_ = capi::PyObject_GetAttrString(pystate.data,"z");
                    auto z = Vector(Z_,z_);
                    
                    // Create a C++ state
                    typename PyConstrained::State::t state(x,y,z);
                    
                    // Convert Python state to C++ 
                    pystate.fromPython(state);

                    // Write the restart file
                    PyJsonConstrained::write_restart(fname,state);
                    
                    // Return nothing 
                    Py_RETURN_NONE; 
                
                // Catch errors 
                } CATCH_PYTHON_ERRORS;

                // Reads a json restart file
                PyObject * read_restart(
                    PyObject * self,
                    PyObject * args
                ) try {
                    // Grab variables
                    PY_VAR_8(X,Y,Z,fname,x,y,z,pystate);

                    // Grab the file name
                    auto fname = capi::PyString_AsString(fname_);

                    // Create a Python state 
                    auto pystate = Python::State <PyConstrained> (pystate_);
                    
                    // Grab the reference vector 
                    auto x = Vector(X_,x_);
                    auto y = Vector(Y_,y_);
                    auto z = Vector(Z_,z_);
                    
                    // Create a C++ state
                    typename PyConstrained::State::t state(x,y,z);

                    // Read the restart file into the C++ state 
                    PyJsonConstrained::read_restart(fname,x,y,z,state);
                    
                    // Convert the C++ state to a Python state
                    pystate.toPython(state);

                    // Return nothing 
                    Py_RETURN_NONE; 

                // Catch errors 
                } CATCH_PYTHON_ERRORS;
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
