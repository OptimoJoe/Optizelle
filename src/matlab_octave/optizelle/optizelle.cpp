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

// Catch and handle MATLAB/Octave and Optizelle errors.  In theory, we should
// never have to handle Matlab::Exception::t since MATLAB/Octave should exit
// the mex file immediately.  MATLAB does this correctly, but Octave does not,
// so we require this safeguard.  Though, it swallows the actual error, which
// is bad.  I've this reported as bug #48185 with Octave's bugtracker at
// savannah.gnu.org.
#define CATCH_MATLAB_ERRORS \
    catch(Matlab::Exception::t const & e) { \
        mexErrMsgTxt( \
            Optizelle::Exception::to_string(e).c_str()); \
    } catch(std::exception const & e) { \
        mexErrMsgTxt( \
            Optizelle::Exception::to_string(e).c_str()); \
        return; \
    }

// Grab and globalize the Optizelle module
#define MX_OPT_BEGIN \
    optizelle.emplace_back(capi::mexGetVariablePtr("global","Optizelle"));

// Clean up the Optizelle module
#define MX_OPT_END optizelle.pop_back(); 

// Return arguments to MATLAB/Octave
#define MX_RETURN_0 \
    MX_OPT_END; \
    return;
#define MX_RETURN_1(v1) \
    MX_OPT_END; \
    pOutput[0] = v1.leak(); \
    return;
#define MX_RETURN_4(v1,v2,v3,v4) \
    MX_OPT_END; \
    pOutput[0] = v1.leak(); \
    pOutput[1] = v2.leak(); \
    pOutput[2] = v3.leak(); \
    pOutput[3] = v4.leak(); \
    return;
#define MX_RETURN_5(v1,v2,v3,v4,v5) \
    MX_OPT_END; \
    pOutput[0] = v1.leak(); \
    pOutput[1] = v2.leak(); \
    pOutput[2] = v3.leak(); \
    pOutput[3] = v4.leak(); \
    pOutput[4] = v5.leak(); \
    return;
#define MX_RETURN_6(v1,v2,v3,v4,v5,v6) \
    MX_OPT_END; \
    pOutput[0] = v1.leak(); \
    pOutput[1] = v2.leak(); \
    pOutput[2] = v3.leak(); \
    pOutput[3] = v4.leak(); \
    pOutput[4] = v5.leak(); \
    pOutput[5] = v6.leak(); \
    return;

// Grab arguments from MATLAB/Octave
#define MX_VAR_2(v1,v2) \
    MX_OPT_BEGIN; \
    auto v1##_ = mxArrayPtr(pInput[0],mxArrayPtr::Unmanaged); \
    auto v2##_ = mxArrayPtr(pInput[1],mxArrayPtr::Unmanaged);
#define MX_VAR_3(v1,v2,v3) \
    MX_OPT_BEGIN; \
    auto v1##_ = mxArrayPtr(pInput[0],mxArrayPtr::Unmanaged); \
    auto v2##_ = mxArrayPtr(pInput[1],mxArrayPtr::Unmanaged); \
    auto v3##_ = mxArrayPtr(pInput[2],mxArrayPtr::Unmanaged);
#define MX_VAR_4(v1,v2,v3,v4) \
    MX_OPT_BEGIN; \
    auto v1##_ = mxArrayPtr(pInput[0],mxArrayPtr::Unmanaged); \
    auto v2##_ = mxArrayPtr(pInput[1],mxArrayPtr::Unmanaged); \
    auto v3##_ = mxArrayPtr(pInput[2],mxArrayPtr::Unmanaged); \
    auto v4##_ = mxArrayPtr(pInput[3],mxArrayPtr::Unmanaged);
#define MX_VAR_5(v1,v2,v3,v4,v5) \
    MX_OPT_BEGIN; \
    auto v1##_ = mxArrayPtr(pInput[0],mxArrayPtr::Unmanaged); \
    auto v2##_ = mxArrayPtr(pInput[1],mxArrayPtr::Unmanaged); \
    auto v3##_ = mxArrayPtr(pInput[2],mxArrayPtr::Unmanaged); \
    auto v4##_ = mxArrayPtr(pInput[3],mxArrayPtr::Unmanaged); \
    auto v5##_ = mxArrayPtr(pInput[4],mxArrayPtr::Unmanaged);
#define MX_VAR_6(v1,v2,v3,v4,v5,v6) \
    MX_OPT_BEGIN; \
    auto v1##_ = mxArrayPtr(pInput[0],mxArrayPtr::Unmanaged); \
    auto v2##_ = mxArrayPtr(pInput[1],mxArrayPtr::Unmanaged); \
    auto v3##_ = mxArrayPtr(pInput[2],mxArrayPtr::Unmanaged); \
    auto v4##_ = mxArrayPtr(pInput[3],mxArrayPtr::Unmanaged); \
    auto v5##_ = mxArrayPtr(pInput[4],mxArrayPtr::Unmanaged); \
    auto v6##_ = mxArrayPtr(pInput[5],mxArrayPtr::Unmanaged);
#define MX_VAR_7(v1,v2,v3,v4,v5,v6,v7) \
    MX_OPT_BEGIN; \
    auto v1##_ = mxArrayPtr(pInput[0],mxArrayPtr::Unmanaged); \
    auto v2##_ = mxArrayPtr(pInput[1],mxArrayPtr::Unmanaged); \
    auto v3##_ = mxArrayPtr(pInput[2],mxArrayPtr::Unmanaged); \
    auto v4##_ = mxArrayPtr(pInput[3],mxArrayPtr::Unmanaged); \
    auto v5##_ = mxArrayPtr(pInput[4],mxArrayPtr::Unmanaged); \
    auto v6##_ = mxArrayPtr(pInput[5],mxArrayPtr::Unmanaged); \
    auto v7##_ = mxArrayPtr(pInput[6],mxArrayPtr::Unmanaged);
#define MX_VAR_8(v1,v2,v3,v4,v5,v6,v7,v8) \
    MX_OPT_BEGIN; \
    auto v1##_ = mxArrayPtr(pInput[0],mxArrayPtr::Unmanaged); \
    auto v2##_ = mxArrayPtr(pInput[1],mxArrayPtr::Unmanaged); \
    auto v3##_ = mxArrayPtr(pInput[2],mxArrayPtr::Unmanaged); \
    auto v4##_ = mxArrayPtr(pInput[3],mxArrayPtr::Unmanaged); \
    auto v5##_ = mxArrayPtr(pInput[4],mxArrayPtr::Unmanaged); \
    auto v6##_ = mxArrayPtr(pInput[5],mxArrayPtr::Unmanaged); \
    auto v7##_ = mxArrayPtr(pInput[6],mxArrayPtr::Unmanaged); \
    auto v8##_ = mxArrayPtr(pInput[7],mxArrayPtr::Unmanaged);
#define MX_VAR_10(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10) \
    MX_OPT_BEGIN; \
    auto v1##_ = mxArrayPtr(pInput[0],mxArrayPtr::Unmanaged); \
    auto v2##_ = mxArrayPtr(pInput[1],mxArrayPtr::Unmanaged); \
    auto v3##_ = mxArrayPtr(pInput[2],mxArrayPtr::Unmanaged); \
    auto v4##_ = mxArrayPtr(pInput[3],mxArrayPtr::Unmanaged); \
    auto v5##_ = mxArrayPtr(pInput[4],mxArrayPtr::Unmanaged); \
    auto v6##_ = mxArrayPtr(pInput[5],mxArrayPtr::Unmanaged); \
    auto v7##_ = mxArrayPtr(pInput[6],mxArrayPtr::Unmanaged); \
    auto v8##_ = mxArrayPtr(pInput[7],mxArrayPtr::Unmanaged); \
    auto v9##_ = mxArrayPtr(pInput[8],mxArrayPtr::Unmanaged); \
    auto v10##_ = mxArrayPtr(pInput[9],mxArrayPtr::Unmanaged);

namespace Optizelle {
    // In theory, I'd like to keep this variable local, but for some strange
    // reason Octave keeps moving around this memory.  As such, rather than
    // having it be static, I made it global and then grab the correct pointer
    // each time we enter these routines.
    static std::list<Matlab::mxArrayPtr> optizelle;

    // Enumerate types
    namespace OptimizationStop { 
        // Converts t to a Matlab enumerated type
        Matlab::mxArrayPtr toMatlab(t const & opt_stop) {
            // Do the conversion
            switch(opt_stop){
            case NotConverged:
                return Matlab::capi::enumToMxArray(
                    "OptimizationStop","NotConverged");
            case GradientSmall:
                return Matlab::capi::enumToMxArray(
                    "OptimizationStop","GradientSmall");
            case StepSmall:
                return Matlab::capi::enumToMxArray(
                    "OptimizationStop","StepSmall");
            case MaxItersExceeded:
                return Matlab::capi::enumToMxArray(
                    "OptimizationStop","MaxItersExceeded");
            case InteriorPointInstability:
                return Matlab::capi::enumToMxArray(
                    "OptimizationStop","InteriorPointInstability");
            case GlobalizationFailure:
                return Matlab::capi::enumToMxArray(
                    "OptimizationStop","GlobalizationFailure");
            case UserDefined:
                return Matlab::capi::enumToMxArray(
                    "OptimizationStop","UserDefined");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = Matlab::capi::mxArrayToNatural(member);

            if(m==Matlab::capi::enumToNatural(
                "OptimizationStop","NotConverged")
            )
                return NotConverged;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationStop","GradientSmall")
            )
                return GradientSmall;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationStop","StepSmall")
            )
                return StepSmall;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationStop","MaxItersExceeded")
            )
                return MaxItersExceeded;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationStop","InteriorPointInstability")
            )
                return InteriorPointInstability;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationStop","GlobalizationFailure")
            )
                return GlobalizationFailure;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationStop","UserDefined")
            )
                return UserDefined;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown OptimizationStop");
        }
    }
    
    namespace TruncatedStop { 
        // Converts t to a Matlab enumerated type
        Matlab::mxArrayPtr toMatlab(t const & trunc_stop) {
            // Do the conversion
            switch(trunc_stop){
            case NotConverged:
                return Matlab::capi::enumToMxArray(
                    "TruncatedStop","NotConverged");
            case NegativeCurvature:
                return Matlab::capi::enumToMxArray(
                    "TruncatedStop","NegativeCurvature");
            case RelativeErrorSmall:
                return Matlab::capi::enumToMxArray(
                    "TruncatedStop","RelativeErrorSmall");
            case MaxItersExceeded:
                return Matlab::capi::enumToMxArray(
                    "TruncatedStop","MaxItersExceeded");
            case TrustRegionViolated:
                return Matlab::capi::enumToMxArray(
                    "TruncatedStop","TrustRegionViolated");
            case NanOperator:
                return Matlab::capi::enumToMxArray("TruncatedStop","NanOperator");
            case NanPreconditioner:
                return Matlab::capi::enumToMxArray(
                    "TruncatedStop","NanPreconditioner");
            case NonProjectorPreconditioner:
                return Matlab::capi::enumToMxArray(
                    "TruncatedStop","NonProjectorPreconditioner");
            case NonSymmetricPreconditioner:
                return Matlab::capi::enumToMxArray(
                    "TruncatedStop","NonSymmetricPreconditioner");
            case NonSymmetricOperator:
                return Matlab::capi::enumToMxArray(
                    "TruncatedStop","NonSymmetricOperator");
            case LossOfOrthogonality:
                return Matlab::capi::enumToMxArray(
                    "TruncatedStop","LossOfOrthogonality");
            case OffsetViolatesTrustRegion:
                return Matlab::capi::enumToMxArray(
                    "TruncatedStop","OffsetViolatesTrustRegion");
            case OffsetViolatesSafeguard:
                return Matlab::capi::enumToMxArray(
                    "TruncatedStop","OffsetViolatesSafeguard");
            case TooManyFailedSafeguard:
                return Matlab::capi::enumToMxArray(
                    "TruncatedStop","TooManyFailedSafeguard");
            case ObjectiveIncrease:
                return Matlab::capi::enumToMxArray(
                    "TruncatedStop","ObjectiveIncrease");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = Matlab::capi::mxArrayToNatural(member);

            if(m==Matlab::capi::enumToNatural("TruncatedStop","NotConverged"))
                return NotConverged;
            else if(m==Matlab::capi::enumToNatural(
                "TruncatedStop","NegativeCurvature")
            )
                return NegativeCurvature;
            else if(m==Matlab::capi::enumToNatural(
                "TruncatedStop","RelativeErrorSmall")
            )
                return RelativeErrorSmall;
            else if(m==Matlab::capi::enumToNatural(
                "TruncatedStop","MaxItersExceeded")
            )
                return MaxItersExceeded;
            else if(m==Matlab::capi::enumToNatural(
                "TruncatedStop","TrustRegionViolated")
            )
                return TrustRegionViolated;
            else if(m==Matlab::capi::enumToNatural(
                "TruncatedStop","NanOperator")
            )
                return NanOperator;
            else if(m==Matlab::capi::enumToNatural(
                "TruncatedStop","NanPreconditioner")
            )
                return NanPreconditioner;
            else if(m==Matlab::capi::enumToNatural(
                "TruncatedStop","NonProjectorPreconditioner")
            )
                return NonProjectorPreconditioner;
            else if(m==Matlab::capi::enumToNatural(
                "TruncatedStop","NonSymmetricPreconditioner")
            )
                return NonSymmetricPreconditioner;
            else if(m==Matlab::capi::enumToNatural(
                "TruncatedStop","NonSymmetricOperator")
            )
                return NonSymmetricOperator;
            else if(m==Matlab::capi::enumToNatural(
                "TruncatedStop","LossOfOrthogonality")
            )
                return LossOfOrthogonality;
            else if(m==Matlab::capi::enumToNatural(
                "TruncatedStop","OffsetViolatesTrustRegion")
            )
                return OffsetViolatesTrustRegion;
            else if(m==Matlab::capi::enumToNatural(
                "TruncatedStop","OffsetViolatesSafeguard")
            )
                return OffsetViolatesSafeguard;
            else if(m==Matlab::capi::enumToNatural(
                "TruncatedStop","TooManyFailedSafeguard")
            )
                return TooManyFailedSafeguard;
            else if(m==Matlab::capi::enumToNatural(
                "TruncatedStop","ObjectiveIncrease")
            )
                return ObjectiveIncrease;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown TruncatedStop");
        }
    }

    namespace AlgorithmClass { 
        // Converts t to a Matlab enumerated type
        Matlab::mxArrayPtr toMatlab(t const & algorithm_class) {
            // Do the conversion
            switch(algorithm_class){
            case TrustRegion:
                return Matlab::capi::enumToMxArray("AlgorithmClass","TrustRegion");
            case LineSearch:
                return Matlab::capi::enumToMxArray("AlgorithmClass","LineSearch");
            case UserDefined:
                return Matlab::capi::enumToMxArray("AlgorithmClass","UserDefined");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = Matlab::capi::mxArrayToNatural(member);

            if(m==Matlab::capi::enumToNatural("AlgorithmClass","TrustRegion"))
                return TrustRegion;
            else if(m==Matlab::capi::enumToNatural(
                "AlgorithmClass","LineSearch")
            )
                return LineSearch;
            else if(m==Matlab::capi::enumToNatural(
                "AlgorithmClass","UserDefined")
            )
                return UserDefined;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown AlgorithmClass");
        }
    }

    namespace Operators { 
        // Converts t to a Matlab enumerated type
        Matlab::mxArrayPtr toMatlab(t const & op) {
            // Do the conversion
            switch(op){
            case Identity:
                return Matlab::capi::enumToMxArray("Operators","Identity");
            case Zero:
                return Matlab::capi::enumToMxArray("Operators","Zero");
            case ScaledIdentity:
                return Matlab::capi::enumToMxArray("Operators","ScaledIdentity");
            case BFGS:
                return Matlab::capi::enumToMxArray("Operators","BFGS");
            case InvBFGS:
                return Matlab::capi::enumToMxArray("Operators","InvBFGS");
            case SR1:
                return Matlab::capi::enumToMxArray("Operators","SR1");
            case InvSR1:
                return Matlab::capi::enumToMxArray("Operators","InvSR1");
            case UserDefined:
                return Matlab::capi::enumToMxArray("Operators","UserDefined");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = Matlab::capi::mxArrayToNatural(member);

            if(m==Matlab::capi::enumToNatural("Operators","Identity"))
                return Identity;
            else if(m==Matlab::capi::enumToNatural("Operators","Zero"))
                return Zero;
            else if(m==Matlab::capi::enumToNatural(
                "Operators","ScaledIdentity")
            )
                return ScaledIdentity;
            else if(m==Matlab::capi::enumToNatural("Operators","BFGS"))
                return BFGS;
            else if(m==Matlab::capi::enumToNatural("Operators","InvBFGS"))
                return InvBFGS;
            else if(m==Matlab::capi::enumToNatural("Operators","SR1"))
                return SR1;
            else if(m==Matlab::capi::enumToNatural("Operators","InvSR1"))
                return InvSR1;
            else if(m==Matlab::capi::enumToNatural("Operators","UserDefined"))
                return UserDefined;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown Operators");
        }
    }

    namespace LineSearchDirection {
        // Converts t to a Matlab enumerated type
        Matlab::mxArrayPtr toMatlab(t const & dir) {
            // Do the conversion
            switch(dir){
            case SteepestDescent:
                return Matlab::capi::enumToMxArray("LineSearchDirection",
                    "SteepestDescent");
            case FletcherReeves:
                return Matlab::capi::enumToMxArray("LineSearchDirection",
                    "FletcherReeves");
            case PolakRibiere:
                return Matlab::capi::enumToMxArray("LineSearchDirection",
                    "PolakRibiere");
            case HestenesStiefel:
                return Matlab::capi::enumToMxArray("LineSearchDirection",
                    "HestenesStiefel");
            case BFGS:
                return Matlab::capi::enumToMxArray(
                    "LineSearchDirection","BFGS");
            case NewtonCG:
                return Matlab::capi::enumToMxArray(
                    "LineSearchDirection","NewtonCG");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = Matlab::capi::mxArrayToNatural(member);

            if(m==Matlab::capi::enumToNatural(
                "LineSearchDirection","SteepestDescent")
            )
                return SteepestDescent;
            else if(m==Matlab::capi::enumToNatural(
                "LineSearchDirection","FletcherReeves")
            )
                return FletcherReeves;
            else if(m==Matlab::capi::enumToNatural(
                "LineSearchDirection","PolakRibiere")
            )
                return PolakRibiere;
            else if(m==Matlab::capi::enumToNatural(
                "LineSearchDirection","HestenesStiefel")
            )
                return HestenesStiefel;
            else if(m==Matlab::capi::enumToNatural(
                "LineSearchDirection","BFGS")
            )
                return BFGS;
            else if(m==Matlab::capi::enumToNatural(
                "LineSearchDirection","NewtonCG")
            )
                return NewtonCG;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown LineSearchDirection");
        }
    }

    namespace LineSearchKind { 
        // Converts t to a Matlab enumerated type
        Matlab::mxArrayPtr toMatlab(t const & kind) {
            // Do the conversion
            switch(kind){
            case GoldenSection:
                return Matlab::capi::enumToMxArray(
                    "LineSearchKind","GoldenSection");
            case BackTracking:
                return Matlab::capi::enumToMxArray(
                    "LineSearchKind","BackTracking");
            case TwoPointA:
                return Matlab::capi::enumToMxArray(
                    "LineSearchKind","TwoPointA");
            case TwoPointB:
                return Matlab::capi::enumToMxArray(
                    "LineSearchKind","TwoPointB");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = Matlab::capi::mxArrayToNatural(member);

            if(m==Matlab::capi::enumToNatural("LineSearchKind","GoldenSection"))
                return GoldenSection;
            else if(m==Matlab::capi::enumToNatural(
                "LineSearchKind","BackTracking")
            )
                return BackTracking;
            else if(m==Matlab::capi::enumToNatural(
                "LineSearchKind","TwoPointA")
            )
                return TwoPointA;
            else if(m==Matlab::capi::enumToNatural(
                "LineSearchKind","TwoPointB")
            )
                return TwoPointB;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown LineSearchKind");
        }
    }

    namespace OptimizationLocation { 
        // Converts t to a Matlab enumerated type
        Matlab::mxArrayPtr toMatlab(t const & loc) {
            // Do the conversion
            switch(loc){
            case BeginningOfOptimization:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","BeginningOfOptimization");
            case BeforeInitialFuncAndGrad:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","BeforeInitialFuncAndGrad");
            case AfterInitialFuncAndGrad:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","AfterInitialFuncAndGrad");
            case BeforeOptimizationLoop:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","BeforeOptimizationLoop");
            case BeginningOfOptimizationLoop:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","BeginningOfOptimizationLoop");
            case BeforeSaveOld:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","BeforeSaveOld");
            case BeforeStep:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","BeforeStep");
            case BeforeGetStep:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","BeforeGetStep");
            case GetStep:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","GetStep");
            case AfterStepBeforeGradient:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","AfterStepBeforeGradient");
            case AfterGradient:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","AfterGradient");
            case BeforeQuasi:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","BeforeQuasi");
            case AfterQuasi:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","AfterQuasi");
            case AfterCheckStop:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","AfterCheckStop");
            case EndOfOptimizationIteration:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","EndOfOptimizationIteration");
            case BeforeLineSearch:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","BeforeLineSearch");
            case AfterRejectedTrustRegion:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","AfterRejectedTrustRegion");
            case AfterRejectedLineSearch:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","AfterRejectedLineSearch");
            case BeforeActualVersusPredicted:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","BeforeActualVersusPredicted");
            case EndOfOptimization:
                return Matlab::capi::enumToMxArray(
                    "OptimizationLocation","EndOfOptimization");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = Matlab::capi::mxArrayToNatural(member);

            if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","BeginningOfOptimization"))
                return BeginningOfOptimization;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","BeforeInitialFuncAndGrad"))
                return BeforeInitialFuncAndGrad;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","AfterInitialFuncAndGrad"))
                return AfterInitialFuncAndGrad;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","BeforeOptimizationLoop"))
                return BeforeOptimizationLoop;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","BeginningOfOptimizationLoop"))
                return BeginningOfOptimizationLoop;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","BeforeSaveOld"))
                return BeforeSaveOld;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","BeforeStep"))
                return BeforeStep;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","BeforeGetStep"))
                return BeforeGetStep;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","GetStep"))
                return GetStep;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","AfterStepBeforeGradient"))
                return AfterStepBeforeGradient;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","AfterGradient"))
                return AfterGradient;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","BeforeQuasi"))
                return BeforeQuasi;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","AfterQuasi"))
                return AfterQuasi;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","AfterCheckStop"))
                return AfterCheckStop;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","EndOfOptimizationIteration"))
                return EndOfOptimizationIteration;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","BeforeLineSearch"))
                return BeforeLineSearch;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","AfterRejectedTrustRegion"))
                return AfterRejectedTrustRegion;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","AfterRejectedLineSearch"))
                return AfterRejectedLineSearch;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","BeforeActualVersusPredicted"))
                return BeforeActualVersusPredicted;
            else if(m==Matlab::capi::enumToNatural(
                "OptimizationLocation","EndOfOptimization"))
                return EndOfOptimization;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown OptimizationLocation");
        }
    }

    namespace FunctionDiagnostics { 
        // Converts t to a Matlab enumerated type
        Matlab::mxArrayPtr toMatlab(t const & diag) {
            // Do the conversion
            switch(diag){
            case NoDiagnostics:
                return Matlab::capi::enumToMxArray(
                    "FunctionDiagnostics","NoDiagnostics");
            case FirstOrder:
                return Matlab::capi::enumToMxArray(
                    "FunctionDiagnostics","FirstOrder");
            case SecondOrder:
                return Matlab::capi::enumToMxArray(
                    "FunctionDiagnostics","SecondOrder");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = Matlab::capi::mxArrayToNatural(member);

            if(m==Matlab::capi::enumToNatural(
                "FunctionDiagnostics","NoDiagnostics")
            )
                return NoDiagnostics;
            else if(m==Matlab::capi::enumToNatural(
                "FunctionDiagnostics","FirstOrder")
            )
                return FirstOrder;
            else if(m==Matlab::capi::enumToNatural(
                "FunctionDiagnostics","SecondOrder")
            )
                return SecondOrder;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown FunctionDiagnostics");
        }
    }

    namespace VectorSpaceDiagnostics { 
        // Converts t to a Matlab enumerated type
        Matlab::mxArrayPtr toMatlab(t const & diag) {
            // Do the conversion
            switch(diag){
            case NoDiagnostics:
                return Matlab::capi::enumToMxArray(
                    "VectorSpaceDiagnostics","NoDiagnostics");
            case Basic:
                return Matlab::capi::enumToMxArray(
                    "VectorSpaceDiagnostics","Basic");
            case EuclideanJordan:
                return Matlab::capi::enumToMxArray(
                    "VectorSpaceDiagnostics","EuclideanJordan");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = Matlab::capi::mxArrayToNatural(member);

            if(m==Matlab::capi::enumToNatural(
                "VectorSpaceDiagnostics","NoDiagnostics")
            )
                return NoDiagnostics;
            else if(m==Matlab::capi::enumToNatural(
                "VectorSpaceDiagnostics","Basic")
            )
                return Basic;
            else if(m==Matlab::capi::enumToNatural(
                "VectorSpaceDiagnostics","EuclideanJordan")
            )
                return EuclideanJordan;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown VectorSpaceDiagnostics");
        }
    }

    namespace DiagnosticScheme { 
        // Converts t to a Matlab enumerated type
        Matlab::mxArrayPtr toMatlab(t const & dscheme) {
            // Do the conversion
            switch(dscheme){
            case Never:
                return Matlab::capi::enumToMxArray("DiagnosticScheme","Never");
            case DiagnosticsOnly:
                return Matlab::capi::enumToMxArray(
                    "DiagnosticScheme","DiagnosticsOnly");
            case EveryIteration:
                return Matlab::capi::enumToMxArray(
                    "DiagnosticScheme","EveryIteration");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = Matlab::capi::mxArrayToNatural(member);

            if(m==Matlab::capi::enumToNatural("DiagnosticScheme","Never"))
                return Never;
            else if(m==Matlab::capi::enumToNatural(
                "DiagnosticScheme","DiagnosticsOnly")
            )
                return DiagnosticsOnly;
            else if(m==Matlab::capi::enumToNatural(
                "DiagnosticScheme","EveryIteration")
            )
                return EveryIteration;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown DiagnosticScheme");
        }
    }

    namespace ToleranceKind { 
        // Converts t to a Matlab enumerated type
        Matlab::mxArrayPtr toMatlab(t const & eps_kind) {
            // Do the conversion
            switch(eps_kind){
            case Relative:
                return Matlab::capi::enumToMxArray(
                    "ToleranceKind","Relative");
            case Absolute:
                return Matlab::capi::enumToMxArray(
                    "ToleranceKind","Absolute");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = Matlab::capi::mxArrayToNatural(member);

            if(m==Matlab::capi::enumToNatural(
                "ToleranceKind","Relative")
            )
                return Relative;
            else if(m==Matlab::capi::enumToNatural(
                "ToleranceKind","Absolute")
            )
                return Absolute;
            else
                throw Optizelle::Exception::t( __LOC__
                    + ", unknown ToleranceKind");
        }
    }

    namespace QuasinormalStop{ 
        // Converts t to a Matlab enumerated type
        Matlab::mxArrayPtr toMatlab(t const & qn_stop) {
            // Do the conversion
            switch(qn_stop){
            case Newton:
                return Matlab::capi::enumToMxArray(
                    "QuasinormalStop","Newton");
            case CauchyTrustRegion:
                return Matlab::capi::enumToMxArray(
                    "QuasinormalStop","CauchyTrustRegion");
            case CauchySafeguard:
                return Matlab::capi::enumToMxArray(
                    "QuasinormalStop","CauchySafeguard");
            case DoglegTrustRegion:
                return Matlab::capi::enumToMxArray(
                    "QuasinormalStop","DoglegTrustRegion");
            case DoglegSafeguard:
                return Matlab::capi::enumToMxArray(
                    "QuasinormalStop","DoglegSafeguard");
            case NewtonTrustRegion:
                return Matlab::capi::enumToMxArray(
                    "QuasinormalStop","NewtonTrustRegion");
            case NewtonSafeguard:
                return Matlab::capi::enumToMxArray(
                    "QuasinormalStop","NewtonSafeguard");
            case Feasible:
                return Matlab::capi::enumToMxArray(
                    "QuasinormalStop","Feasible");
            case CauchySolved:
                return Matlab::capi::enumToMxArray(
                    "QuasinormalStop","CauchySolved");
            case LocalMin:
                return Matlab::capi::enumToMxArray(
                    "QuasinormalStop","LocalMin");
            case NewtonFailed:
                return Matlab::capi::enumToMxArray(
                    "QuasinormalStop","NewtonFailed");
            }
        }

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member) {
            // Convert the member to a Natural 
            auto m = Matlab::capi::mxArrayToNatural(member);

            if(m==Matlab::capi::enumToNatural(
                "QuasinormalStop","Newton")
            )
                return Newton;
            else if(m==Matlab::capi::enumToNatural(
                "QuasinormalStop","CauchyTrustRegion")
            )
                return CauchyTrustRegion;
            else if(m==Matlab::capi::enumToNatural(
                "QuasinormalStop","CauchySafeguard")
            )
                return CauchySafeguard;
            else if(m==Matlab::capi::enumToNatural(
                "QuasinormalStop","DoglegTrustRegion")
            )
                return DoglegTrustRegion;
            else if(m==Matlab::capi::enumToNatural(
                "QuasinormalStop","DoglegSafeguard")
            )
                return DoglegSafeguard;
            else if(m==Matlab::capi::enumToNatural(
                "QuasinormalStop","NewtonTrustRegion")
            )
                return NewtonTrustRegion;
            else if(m==Matlab::capi::enumToNatural(
                "QuasinormalStop","NewtonSafeguard")
            )
                return NewtonSafeguard;
            else if(m==Matlab::capi::enumToNatural(
                "QuasinormalStop","Feasible")
            )
                return Feasible;
            else if(m==Matlab::capi::enumToNatural(
                "QuasinormalStop","CauchySolved")
            )
                return CauchySolved;
            else if(m==Matlab::capi::enumToNatural(
                "QuasinormalStop","LocalMin")
            )
                return LocalMin;
            else if(m==Matlab::capi::enumToNatural(
                "QuasinormalStop","NewtonFailed")
            )
                return NewtonFailed;
            else
                throw Optizelle::Exception::t( __LOC__
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
                auto name = Matlab::capi::mxCreateString(name_);
                auto iter = Matlab::capi::mxArrayFromNatural(iter_);

                // Call the serialize routine on the vector
                auto x_json = Matlab::capi::mexCallMATLAB3(
                    std::string("serialize"),
                    x.data,
                    name,
                    iter,
                    __LOC__
                        + ", evaluation of the serialize function failed");

                // Convert the serialized vector to a string and return it 
                return Matlab::capi::mxArrayToString(x_json);
            }

            static Matlab::Vector deserialize (
                Matlab::Vector const & x_,
                std::string const & x_json_
            ) {
                // Convert the inputed string into Matlab
                auto x_json =Matlab::capi::mxCreateString(x_json_);

                // Allocate memory for a new Matlab vector
                auto x = x_.init();

                // Call the deserialize routine on the reference vector and the
                // json vector
                x.data = Matlab::capi::mexCallMATLAB2(
                    std::string("deserialize"),
                    x_.data,
                    x_json,
                    __LOC__
                        + ", evaluation of the deserialize function failed");

                // Move out the new vector
                return x;
            }
        };
    }

    namespace Matlab {
        // Grab the pointer and the method we use to destruct it
        mxArrayPtr::mxArrayPtrData::mxArrayPtrData(
            mxArray const * const & ptr_,
            std::function <void(mxArray *)> const & destructor_
        ) : ptr(const_cast <mxArray * const> (ptr_)),
            destructor(destructor_)
        {}

        // Call our custom destructor on the data 
        mxArrayPtr::mxArrayPtrData::~mxArrayPtrData() { 
            if(ptr)
                destructor(ptr);
        }

        // Grab the pointer
        mxArrayPtr::mxArrayPtr(
            mxArray const * const & ptr,
            Mode const & mode
        ) : data(std::make_shared <mxArrayPtrData> (
                ptr,
                mode==Managed?mxDestroyArray:[](mxArray *){}))
        {}

        // Grab the internal pointer
        mxArray * mxArrayPtr::get() const {
            return data->ptr;
        }

        // Grab the pointer and convert the destructor to no longer free
        // the memory
        mxArray * mxArrayPtr::leak() const {
            data->destructor = [](auto){}; 
            return get();
        }

        namespace capi {
            mxArrayPtr mexCallMATLAB1(
                mxArrayPtr const & fn,
                mxArrayPtr const & arg1,
                std::string const & errmsg
            ) {
                mxArray * input[2] = {fn.get(),arg1.get()};
                mxArray * output[1] = {nullptr}; 
                auto err=::mexCallMATLAB(1,output,2,input,"feval");
                if(err)
                    throw Matlab::Exception::t(errmsg); 
                return output[0];
            }
            void mexCallMATLAB1_0(
                mxArrayPtr const & fn,
                mxArrayPtr const & arg1,
                std::string const & errmsg
            ) {
                mxArray * input[2] = {fn.get(),arg1.get()};
                auto err=::mexCallMATLAB(0,nullptr,2,input,"feval");
                if(err)
                    throw Matlab::Exception::t(errmsg); 
            }
            mxArrayPtr mexCallMATLAB2(
                mxArrayPtr const & fn,
                mxArrayPtr const & arg1,
                mxArrayPtr const & arg2,
                std::string const & errmsg
            ) {
                mxArray * input[3] = {fn.get(),arg1.get(),arg2.get()};
                mxArray * output[1] = {nullptr}; 
                auto err=::mexCallMATLAB(1,output,3,input,"feval");
                if(err)
                    throw Matlab::Exception::t(errmsg); 
                return output[0];
            }
            mxArrayPtr mexCallMATLAB2(
                std::string const & fn,
                mxArrayPtr const & arg1,
                mxArrayPtr const & arg2,
                std::string const & errmsg
            ) {
                mxArray * input[2] = {arg1.get(),arg2.get()};
                mxArray * output[1] = {nullptr}; 
                auto err=::mexCallMATLAB(1,output,2,input,fn.c_str());
                if(err)
                    throw Matlab::Exception::t(errmsg); 
                return output[0];
            }
            mxArrayPtr mexCallMATLAB3(
                mxArrayPtr const & fn,
                mxArrayPtr const & arg1,
                mxArrayPtr const & arg2,
                mxArrayPtr const & arg3,
                std::string const & errmsg
            ) {
                mxArray * input[4] = {fn.get(),arg1.get(),arg2.get(),
                    arg3.get()};
                mxArray * output[1] = {nullptr}; 
                auto err=::mexCallMATLAB(1,output,4,input,"feval");
                if(err)
                    throw Matlab::Exception::t(errmsg); 
                return output[0];
            }
            mxArrayPtr mexCallMATLAB3(
                std::string const & fn,
                mxArrayPtr const & arg1,
                mxArrayPtr const & arg2,
                mxArrayPtr const & arg3,
                std::string const & errmsg
            ) {
                mxArray * input[3] = {arg1.get(),arg2.get(),arg3.get()};
                mxArray * output[1] = {nullptr}; 
                auto err=::mexCallMATLAB(1,output,3,input,fn.c_str());
                if(err)
                    throw Matlab::Exception::t(errmsg); 
                return output[0];
            }
            mxArrayPtr mexCallMATLAB4(
                mxArrayPtr const & fn,
                mxArrayPtr const & arg1,
                mxArrayPtr const & arg2,
                mxArrayPtr const & arg3,
                mxArrayPtr const & arg4,
                std::string const & errmsg
            ) {
                mxArray * input[5] = {fn.get(),arg1.get(),arg2.get(),
                    arg3.get(),arg4.get()};
                mxArray * output[1] = {nullptr}; 
                auto err=::mexCallMATLAB(1,output,5,input,"feval");
                if(err)
                    throw Matlab::Exception::t(errmsg); 
                return output[0];
            }
            mxArrayPtr mxDuplicateArray(mxArrayPtr const & x) {
                auto ret = ::mxDuplicateArray(x.get());
                if(!ret)
                    throw Matlab::Exception::t(__LOC__
                        + ", error when duplicating an array"); 
                return ret;
            }
            mxArrayPtr mxGetField(
                mxArrayPtr const & pm,
                mwIndex const & index,
                std::string const & fieldname
            ) {
                auto ret = ::mxGetField(pm.get(),index,fieldname.c_str());
                if(!ret)
                    throw Matlab::Exception::t(__LOC__
                        + ", unable to get the field " + fieldname
                        + " at position " + std::to_string(index));
                return {ret,mxArrayPtr::Unmanaged};
            }
            void mxSetField(
                mxArrayPtr & pm,
                mwIndex const & index,
                std::string const & fieldname,
                mxArrayPtr const & pvalue
            ) {
                // First, check if there's memory in that location
                auto field = ::mxGetField(pm.get(),index,fieldname.c_str());

                // If we have memory, free it
                if(field)
                    ::mxDestroyArray(field);

                // Now, set the field with our value.  Make sure to duplicate
                // the memory first, since this method steals the memory
                ::mxSetField(
                    pm.get(),
                    index,
                    fieldname.c_str(),
                    capi::mxDuplicateArray(pvalue).leak());
            }
            mxArrayPtr mxCreateString(std::string const & str) {
                auto ret = ::mxCreateString(str.c_str());
                if(!ret)
                    throw Matlab::Exception::t(__LOC__
                        + ", unable to convert the string " + str
                        + " into an mxArray");
                return ret;
            }
            std::string mxArrayToString(mxArrayPtr const & array_ptr) {
                auto ret = ::mxArrayToString(array_ptr.get());
                if(!ret)
                    throw Matlab::Exception::t(__LOC__
                        + ", unable to convert an mxArray into a string");
                return ret;
            }
            mxArrayPtr mxCreateCellMatrix(mwSize const & m, mwSize const & n) {
                auto ret = ::mxCreateCellMatrix(m,n);
                if(!ret)
                    throw Matlab::Exception::t(__LOC__
                        + ", unable to create a cell array of size "
                        + std::to_string(m) + " x " + std::to_string(n));
                return ret;
            }
            mxArrayPtr mxGetCell(mxArrayPtr const & pm, mwIndex const & index) {
                auto ret = ::mxGetCell(pm.get(),index);
                if(!ret)
                    throw Matlab::Exception::t(__LOC__
                        + ", unable to get the cell item at position "
                        + std::to_string(index));
                return {ret,mxArrayPtr::Unmanaged};
            }
            void mxSetCell(
                mxArrayPtr & pm,
                mwIndex const & index,
                mxArrayPtr const & value
            ) {
                // First, check if there's memory in that location
                auto item = ::mxGetCell(pm.get(),index);

                // If we have memory, free it
                if(item)
                    ::mxDestroyArray(item);

                // Now, set the field with our value.  Make sure to duplicate
                // the memory first, since this method steals the memory
                ::mxSetCell(
                    pm.get(),
                    index,
                    capi::mxDuplicateArray(value).leak());
            }
            size_t mxGetN(mxArrayPtr const & pm) {
                return ::mxGetN(pm.get());
            }
            mxArrayPtr mexGetVariablePtr(
                std::string const & workspace,
                std::string const & varname
            ) {
                auto ret=::mexGetVariablePtr(workspace.c_str(),varname.c_str());
                if(!ret)
                    throw Matlab::Exception::t(__LOC__
                        + ", unable to get the variable " + varname
                        + " in the workspace " + workspace);
                return {ret,mxArrayPtr::Unmanaged};
            }
            mxArrayPtr mxCreateStructMatrix(
                mwSize const & m,
                mwSize const & n,
                int const & nfields, 
                const char **fieldnames
            ) {
                auto ret = ::mxCreateStructMatrix(m,n,nfields,fieldnames);
                if(!ret)
                    throw Matlab::Exception::t(__LOC__
                        + ", unable to create a struct array of size "
                        + std::to_string(m) + " x " + std::to_string(n)
                        + " with " + std::to_string(nfields) + " fields");
                return ret;
            }

            // Creates a MATLAB/Octave double from a C++ double
            mxArrayPtr mxArrayFromDouble(double const x_) {
                auto x = mxArrayPtr(::mxCreateDoubleMatrix(1,1,mxREAL));
                ::mxGetPr(x.get())[0]=x_;
                return x;
            }

            // Creates a C++ double from a MATLAB/Octave double
            double mxArrayToDouble(mxArrayPtr const & x_) {
                return ::mxGetPr(x_.get())[0]; 
            }

            // Creates a MATLAB/Octave int from a C++ Natural 
            mxArrayPtr mxArrayFromNatural(Natural const x_) {
                auto x = mxArrayPtr(mxCreateDoubleMatrix(1,1,mxREAL));
                ::mxGetPr(x.get())[0]=x_;
                return x;
            }

            // Creates a C++ Natural from a MATLAB/Octave integer 
            Natural mxArrayToNatural(mxArrayPtr const & x_) {
                return ::mxGetPr(x_.get())[0]; 
            }

            // Converts an Optizelle enumerated type to a mxArray *
            mxArrayPtr enumToMxArray(
                std::string const & type_,
                std::string const & member_
            ) {
                // Grab the type 
                auto type = capi::mxGetField(optizelle.back(),0,type_);

                // Grab the member
                auto member = capi::mxGetField(type,0,member_);
                    
                // Return the member 
                return mxArrayFromNatural(fromDouble(mxGetPr(member.get())[0]));
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
                if(value>=double(std::numeric_limits<Optizelle::Natural>::max())
                )
                    return std::numeric_limits <Optizelle::Natural>::max();
                else if(
                    value<=double(
                        std::numeric_limits <Optizelle::Natural>::min())
                )
                    return std::numeric_limits <Optizelle::Natural>::min();
                else
                    return Optizelle::Natural(value);
            }
        }

        // A messaging utility that hooks directly into MATLAB/Octave 
        namespace Messaging {
            Optizelle::Messaging::t matlab(mxArrayPtr const & print) {
                return [print](std::string const & msg_) {
                    // Call the print function
                    auto msg = capi::mxCreateString(msg_);
                    capi::mexCallMATLAB1_0(
                        print,
                        msg,
                        __LOC__
                            + ", evaluation of the Messaging function failed");
                };
            }
        }

        // Grab the vector space and data 
        Vector::Vector(mxArrayPtr const & vs_, mxArrayPtr const & data_) :
            vs(vs_),
            data(data_)
        {}

        // Memory allocation and size setting 
        Vector Vector::init() const {
            // Call the init function on the internal and store in y 
            auto init = capi::mxGetField(vs,0,"init");
            auto ret = capi::mexCallMATLAB1(
                init,
                data,
                 __LOC__
                    + ", evaluation of the vector space function init failed");

            // Return the vector
            return Vector(vs,ret);
        }
        
        // y <- x (Shallow.  No memory allocation.)  Internal is y.
        void Vector::copy(Vector const & x) { 
            // Call the copy function on x and the internal 
            auto copy = capi::mxGetField(vs,0,"copy");
            data = capi::mexCallMATLAB1(
                copy,
                x.data,
                 __LOC__
                    + ", evaluation of the vector space function copy failed");
        }

        // x <- alpha * x.  Internal is x.
        void Vector::scal(double const & alpha_) { 
            // Call the scal function on alpha and the internal storage 
            auto scal = capi::mxGetField(vs,0,"scal");
            auto alpha = capi::mxArrayFromDouble(alpha_);
            data = capi::mexCallMATLAB2(
                scal,
                alpha,
                data,
                __LOC__
                    + ", evaluation of the vector space function scal failed");
        } 

        // x <- 0.  Internal is x. 
        void Vector::zero() { 
            // Call the zero function on this vector.
            auto zero = capi::mxGetField(vs,0,"zero");
            data = capi::mexCallMATLAB1(
                zero,
                data,
                __LOC__
                    + ", evaluation of the vector space function zero failed");
        } 

        // y <- alpha * x + y.   Internal is y.
        void Vector::axpy(double const & alpha_,Vector const & x) { 
            // Call the axpy function on alpha, x, and the internal storage.
            auto axpy = capi::mxGetField(vs,0,"axpy");
            auto alpha = capi::mxArrayFromDouble(alpha_);
            data = capi::mexCallMATLAB3(
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
            auto innr = capi::mxGetField(vs,0,"innr");
            return capi::mxArrayToDouble(capi::mexCallMATLAB2(
                innr,
                x.data,
                data,
                __LOC__
                    + ", evaluation of the vector space function innr failed"));
        } 

        // x <- random.  Internal is x. 
        void Vector::rand() { 
            // Call the rand function on this vector.
            auto rand = capi::mxGetField(vs,0,"rand");
            data = capi::mexCallMATLAB1(
                rand,
                data,
                __LOC__
                    + ", evaluation of the vector space function rand failed");
        }

        // Jordan product, z <- x o y.  Internal is z.
        void Vector::prod(Vector const & x,Vector const & y) { 
            // Call the prod function on x, y, and the internal 
            auto prod = capi::mxGetField(vs,0,"prod");
            data = capi::mexCallMATLAB2(
                prod,
                x.data,
                y.data,
                __LOC__
                    + ", evaluation of the vector space function prod failed");
        } 

        // Identity element, x <- e such that x o e = x .  Internal is x.
        void Vector::id() { 
            // Call the id function on the internal.
            auto id = capi::mxGetField(vs,0,"id");
            data = capi::mexCallMATLAB1(
                id,
                data,
                __LOC__
                    + ", evaluation of the vector space function id failed");
        } 

        // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y.
        // Internal is z.
        void Vector::linv(Vector const & x, Vector const & y) { 
            // Call the linv function on x, y, and the internal
            auto linv = capi::mxGetField(vs,0,"linv");
            data = capi::mexCallMATLAB2(
                linv,
                x.data,
                y.data,
                __LOC__
                    + ", evaluation of the vector space function linv failed");
        } 

        // Barrier function, barr <- barr(x) where x o grad barr(x) = e.
        // Internal is x.
        double Vector::barr() const {
            // Call the barr function on the internal.  Store in z.
            auto barr = capi::mxGetField(vs,0,"barr");
            return capi::mxArrayToDouble(capi::mexCallMATLAB1(
                barr,
                data,
                __LOC__ +
                    ", evaluation of the vector space function barr failed."));
        } 

        // Line search, srch <- argmax {alpha in Real >= 0 : alpha x + y >= 0} 
        // where y > 0.  Internal is y.
        double Vector::srch(Vector const & x) const {
            // Call the srch function on x and the internal.  Store in z.
            auto srch = capi::mxGetField(vs,0,"srch");
            return capi::mxArrayToDouble(capi::mexCallMATLAB2(
                srch,
                x.data,
                data,
                __LOC__
                    + ", evaluation of the vector space function srch failed"));
        } 

        // Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
        // operator.  Internal is x.
        void Vector::symm() { 
            // Call the symm function on the internal.
            auto symm = capi::mxGetField(vs,0,"symm");
            data = capi::mexCallMATLAB1(
                symm,
                data,
                __LOC__
                    + ", evaluation of the vector space function symm failed");
        } 
        
        // Converts (copies) a value into Matlab.  
        mxArrayPtr Vector::toMatlab() const {
            // Call the copy function on the internal and x
            auto copy = capi::mxGetField(vs,0,"copy");
            return capi::mexCallMATLAB1(
                copy,
                data,
                __LOC__
                    + ", evaluation of the vector space function copy failed");
        } 
        
        // Converts (copies) a value from Matlab.  This assumes that the
        // vector space functions have already been properly assigned.
        void Vector::fromMatlab(mxArrayPtr const & ptr) {
            // Call the copy function on ptr and the internal 
            auto copy = capi::mxGetField(vs,0,"copy");
            data = capi::mexCallMATLAB1(
                copy,
                ptr,
                __LOC__
                    + ", evaluation of the vector space function copy failed");
        } 

        // Grab the pointer to the function information 
        ScalarValuedFunction::ScalarValuedFunction(mxArrayPtr const & data_)
            : data(data_) {}

        // <- f(x)-> 
        double ScalarValuedFunction::eval(Vector const & x) const { 
            // Call the objective function on x
            auto eval = capi::mxGetField(data,0,"eval");
            return capi::mxArrayToDouble(capi::mexCallMATLAB1(
                eval,
                x.data,
                __LOC__
                    + ", evaluation of the objective f failed"));
        }

        // grad = grad f(x) 
        void ScalarValuedFunction::grad(
            Vector const & x,
            Vector & grad
        ) const { 
            // Call the gradient function on x
            auto mxgrad = capi::mxGetField(data,0,"grad");
            grad.data = capi::mexCallMATLAB1(
                mxgrad,
                x.data,
                __LOC__
                    + ", evaluation of the gradient of f failed.");
        }

        // H_dx = hess f(x) dx 
        void ScalarValuedFunction::hessvec(
            Vector const & x,
            Vector const & dx,
            Vector & H_dx
        ) const {
            // Call the hessvec function on x and dx,
            auto hessvec = capi::mxGetField(data,0,"hessvec");
            H_dx.data = capi::mexCallMATLAB2(
                hessvec,
                x.data,
                dx.data,
                __LOC__
                    + ", evaluation of the Hessian-vector product of f failed");
        }

        // Grab the function's name and a pointer to the underlying data 
        VectorValuedFunction::VectorValuedFunction(
            std::string const & name_,
            mxArrayPtr const & data_ 
        ) :
            name(name_), data(data_)
        {}

        // y=f(x)
        void VectorValuedFunction::eval(
            X_Vector const & x,
            Y_Vector & y
        ) const {
            // Call the objective function on x.
            auto eval = capi::mxGetField(data,0,"eval");
            y.data = capi::mexCallMATLAB1(
                eval,
                x.data,
                __LOC__
                    + ", evaluation of the constraint " + name + " failed");
        }

        // y=f'(x)dx 
        void VectorValuedFunction::p(
            X_Vector const & x,
            X_Vector const & dx,
            Y_Vector & y
        ) const {
            // Call the prime function on x and dx
            auto p = capi::mxGetField(data,0,"p");
            y.data = capi::mexCallMATLAB2(
                p,
                x.data,
                dx.data,
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
            // Call the prime-adjoint function on x and dy
            auto ps = capi::mxGetField(data,0,"ps");
            xhat.data = capi::mexCallMATLAB2(
                ps,
                x.data,
                dy.data,
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
            // Call the prime-adjoint function on x, dx, and dy
            auto pps = capi::mxGetField(data,0,"pps");
            xhat.data = capi::mexCallMATLAB3(
                pps,
                x.data,
                dx.data,
                dy.data,
                __LOC__
                    + ", evaluation of the second derivative-adjoint of the "
                    + "constraint " + name + " failed");
        }

        // Converts elements from C++ to Matlab 
        namespace toMatlab {
                    
            // Sets a real in a Matlab state 
            void Real(
                std::string const & name,
                double const & value,
                mxArrayPtr & mxstate 
            ) {
                auto item = capi::mxArrayFromDouble(value);
                capi::mxSetField(mxstate,0,name,item);
            }
        
            // Sets a natural in a Matlab state 
            void Natural(
                std::string const & name,
                Optizelle::Natural const & value,
                mxArrayPtr & mxstate 
            ) {
                auto item = capi::mxArrayFromNatural(value);
                capi::mxSetField(mxstate,0,name,item);
            }
        
            // Sets a vector in a Matlab state 
            void Vector(
                std::string const & name,
                Matlab::Vector const & value,
                mxArrayPtr & mxstate 
            ) {
                capi::mxSetField(mxstate,0,name,value.toMatlab());
            }
        
            // Sets a list of vectors in a Matlab state 
            void VectorList(
                std::string const & name,
                std::list <Matlab::Vector> const & vectors,
                mxArrayPtr & mxstate 
            ) {
                // Create a new Matlab cell array that we insert elements into
                auto mxvectors = capi::mxCreateCellMatrix(1,vectors.size());

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
                    capi::mxSetCell(mxvectors,i,v.data);

                    // Increment our counter
                    i++;
                }
                
                // Insert the items into mxstate 
                capi::mxSetField(mxstate,0,name,mxvectors);
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
                    auto tuple = capi::mxCreateCellMatrix(1,2);
                    capi::mxSetCell(tuple,0,
                        capi::mxCreateString(value.first));
                    capi::mxSetCell(tuple,1,mxvalue.data);

                    // Put tuple into the cell array
                    capi::mxSetCell(mxvalues,i,tuple);

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
                    auto tuple = capi::mxCreateCellMatrix(1,2);
                    capi::mxSetCell(tuple,0,capi::mxCreateString(value.first));
                    capi::mxSetCell(tuple,1,
                        capi::mxArrayFromDouble(value.second));

                    // Put tuple into the cell array
                    capi::mxSetCell(mxvalues,i,tuple);

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
                    auto tuple = capi::mxCreateCellMatrix(1,2);
                    capi::mxSetCell(tuple,0,capi::mxCreateString(value.first));
                    capi::mxSetCell(tuple,1,
                        capi::mxArrayFromNatural(value.second));

                    // Release the tuple into the Matlab cell array
                    capi::mxSetCell(mxvalues,i,tuple);

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
                    auto tuple = capi::mxCreateCellMatrix(1,2);
                    capi::mxSetCell(tuple,0,capi::mxCreateString(value.first));
                    capi::mxSetCell(tuple,1,capi::mxCreateString(value.second));

                    // Release the tuple into the Matlab cell array
                    capi::mxSetCell(mxvalues,i,tuple);

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
                auto item = capi::mxGetField(mxstate,0,name);
                value=capi::mxArrayToDouble(item);
            }
            
            // Sets a natural in a C++ state 
            void Natural(
                std::string const & name,
                mxArrayPtr const & mxstate,
                Optizelle::Natural & value
            ) {
                auto item = capi::mxGetField(mxstate,0,name);
                value=capi::mxArrayToNatural(item);
            }
            
            // Sets a vector in a C++ state 
            void Vector(
                std::string const & name,
                mxArrayPtr const & mxstate,
                Matlab::Vector & value
            ) {
                auto item = capi::mxGetField(mxstate,0,name);
                value.fromMatlab(item);
            }
            
            // Sets a list of vectors in a C++ state 
            void VectorList(
                std::string const & name,
                mxArrayPtr const & mxstate,
                Matlab::Vector const & vec,
                std::list <Matlab::Vector> & values
            ) {
                // Grab the list of items
                auto items = capi::mxGetField(mxstate,0,name);

                // Loop over all the elements in items and insert them one
                // at a time into values
                values.clear();
                for(auto i=0;i<capi::mxGetN(items);i++) {
                    // Grab the current item from Matlab
                    auto item = capi::mxGetCell(items,i);

                    // Create a new vector in values 
                    values.emplace_back(vec.init());

                    // Copy the Matlab item into the new value
                    values.back().fromMatlab(item);
                }
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
                for(auto i=0;i<capi::mxGetN(mxvalues);i++) {
                    // Grab the current item from Matlab
                    auto mxvalue = capi::mxGetCell(mxvalues,i);

                    // Create the elements in values 
                    values.emplace_back(
                        capi::mxArrayToString(capi::mxGetCell(mxvalue,0)),
                        vec.init());

                    // Copy the Matlab value into the C++ value
                    values.back().second.fromMatlab(
                        capi::mxGetCell(mxvalue,1));
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
                for(auto i=0;i<capi::mxGetN(mxvalues);i++) {
                    // Grab the current item from Matlab
                    auto mxvalue = capi::mxGetCell(mxvalues,i);
                    
                    // Create the elements in values 
                    values.emplace_back(
                        capi::mxArrayToString(capi::mxGetCell(mxvalue,0)),
                        capi::mxArrayToDouble(capi::mxGetCell(mxvalue,1)));
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
                for(auto i=0;i<capi::mxGetN(mxvalues);i++) {
                    // Grab the current item from Matlab
                    auto mxvalue = capi::mxGetCell(mxvalues,i);
                    
                    // Create the elements in values 
                    values.emplace_back(
                        capi::mxArrayToString(capi::mxGetCell(mxvalue,0)),
                        capi::mxArrayToNatural(capi::mxGetCell(mxvalue,1)));
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
                for(auto i=0;i<capi::mxGetN(mxvalues);i++) {
                    // Grab the current item from Matlab
                    auto mxvalue = capi::mxGetCell(mxvalues,i);
                    
                    // Create the elements in values 
                    values.emplace_back(
                        capi::mxArrayToString(capi::mxGetCell(mxvalue,0)),
                        capi::mxArrayToString(capi::mxGetCell(mxvalue,1)));
                }
            }
        }

        // Convert a C++ state to a Matlab state 
        template <>
        void State <MxUnconstrained>::toMatlab(
            typename MxUnconstrained::State::t const & state
        ) {
            Unconstrained::State::toMatlab(state,data);
        }
        template <>
        void State <MxEqualityConstrained>::toMatlab(
            typename MxEqualityConstrained::State::t const & state
        ) {
            EqualityConstrained::State::toMatlab(state,data);
        }
        template <>
        void State <MxInequalityConstrained>::toMatlab(
            typename MxInequalityConstrained::State::t const & state
        ) {
            InequalityConstrained::State::toMatlab(state,data);
        }
        template <>
        void State <MxConstrained>::toMatlab(
            typename MxConstrained::State::t const & state
        ) {
            Constrained::State::toMatlab(state,data);
        }

        // Convert a Matlab state to C++ 
        template <>
        void State <MxUnconstrained>::fromMatlab(
            typename MxUnconstrained::State::t & state
        ) {
            Unconstrained::State::fromMatlab(data,state);
        }
        template <>
        void State <MxEqualityConstrained>::fromMatlab(
            typename MxEqualityConstrained::State::t & state
        ) {
            EqualityConstrained::State::fromMatlab(data,state);
        }
        template <>
        void State <MxInequalityConstrained>::fromMatlab(
            typename MxInequalityConstrained::State::t & state
        ) {
            InequalityConstrained::State::fromMatlab(data,state);
        }
        template <>
        void State <MxConstrained>::fromMatlab(
            typename MxConstrained::State::t & state
        ) {
            Constrained::State::fromMatlab(data,state);
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

                    return names;
                }
                std::vector <char const *> fieldNames() {
                    return fieldNames_();
                }

                // Create the structure for a Matlab state
                mxArrayPtr mxCreate() {
                    auto names = Unconstrained::State::fieldNames();
                    return capi::mxCreateStructMatrix(
                        1,1,names.size(),&(names[0]));
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
                    // Grab the variables
                    MX_VAR_2(X,x);

                    // Create a vector from the user input
                    auto x = Vector(X_,x_);
                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxUnconstrained> (
                        Unconstrained::State::mxCreate());

                    // Create a new C++ state
                    typename MxUnconstrained::State::t state(x);

                    // Convert the state to a Matlab state
                    mxstate_out.toMatlab(state);
                   
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

                } CATCH_MATLAB_ERRORS; 
        
                // Read json parameters from file
                void readJson(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_3(X,fname,mxstate);

                    // Grab the file name
                    auto fname = capi::mxArrayToString(fname_);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxUnconstrained> (mxstate_);
                    auto mxstate_out = Matlab::State <MxUnconstrained> (
                        Unconstrained::State::mxCreate());
                
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);

                    // Create a new C++ state
                    typename MxUnconstrained::State::t state(x);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Read the JSON file into the C++ state
                    MxJsonUnconstrained::read(fname,state);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);
                   
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
                    // Grab the variables
                    MX_VAR_5(X,msg,mxfns,mxstate,smanip);

                    // Create a messaging object
                    auto msg = Optizelle::Matlab::Messaging::matlab(msg_);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxUnconstrained> (mxstate_);
                    auto mxstate_out = Matlab::State <MxUnconstrained> (
                        Unconstrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);

                    // Create a C++ state
                    typename MxUnconstrained::State::t state(x);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Create a Matlab bundle of functions
                    Matlab::Functions <MxUnconstrained> mxfns(
                        mxstate_out,
                        state,
                        mxfns_);

                    // Create a C++ bundle of functions
                    typename MxUnconstrained::Functions::t fns;
                    
                    // Convert the Matlab bundle of functions to C++ 
                    mxfns.fromMatlab(fns);

                    // Create a state manipulator 
                    Matlab::StateManipulator <MxUnconstrained> smanip(
                        mxstate_out,
                        mxfns,
                        smanip_);
                   
                    // Minimize
                    MxUnconstrained::Algorithms::getMin(msg,fns,state,smanip);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

                } CATCH_MATLAB_ERRORS;
            }
            
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                void release(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_2(X,mxstate);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxUnconstrained> (mxstate_);
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);

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
                    auto mxxs = capi::mxCreateCellMatrix(1,xs.size());
                    auto mxreals = capi::mxCreateCellMatrix(1,reals.size());
                    auto mxnats = capi::mxCreateCellMatrix(1,nats.size());
                    auto mxparams = capi::mxCreateCellMatrix(1,params.size());

                    // Convert the restart information to Matlab 
                    toMatlab::Vectors(xs,mxxs);
                    toMatlab::Reals(reals,mxreals);
                    toMatlab::Naturals(nats,mxnats);
                    toMatlab::Params(params,mxparams);
                    
                    // Return the result 
                    MX_RETURN_4(mxxs,mxreals,mxnats,mxparams);
                
                } CATCH_MATLAB_ERRORS;

                // Capture data from structures controlled by the user.  
                void capture(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_6(X,mxstate,mxxs,mxreals,mxnats,mxparams);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxUnconstrained> (mxstate_);
                    auto mxstate_out = Matlab::State <MxUnconstrained> (
                        Unconstrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);

                    // Create a C++ state
                    typename MxUnconstrained::State::t state(x);
                   
                    // Allocate memory for the released vectors
                    MxUnconstrained::Restart::X_Vectors xs;
                    MxUnconstrained::Restart::Reals reals;
                    MxUnconstrained::Restart::Naturals nats;
                    MxUnconstrained::Restart::Params params;
                    
                    // Convert the restart information from Matlab 
                    fromMatlab::Vectors(x,mxxs_,xs);
                    fromMatlab::Reals(mxreals_,reals);
                    fromMatlab::Naturals(mxnats_,nats);
                    fromMatlab::Params(mxparams_,params);

                    // Do a capture 
                    MxUnconstrained::Restart
                        ::capture(state,xs,reals,nats,params);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

                } CATCH_MATLAB_ERRORS;
                
                // Writes a json restart file
                void write_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_3(X,fname,mxstate);

                    // Grab the file name
                    auto fname = capi::mxArrayToString(fname_);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxUnconstrained> (mxstate_);
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    
                    // Create a C++ state
                    typename MxUnconstrained::State::t state(x);
                    
                    // Convert Matlab state to C++ 
                    mxstate.fromMatlab(state);

                    // Write the restart file
                    MxJsonUnconstrained::write_restart(fname,state);

                    // Return nothing
                    MX_RETURN_0;

                } CATCH_MATLAB_ERRORS;
                
                // Reads a json restart file
                void read_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_3(X,fname,x);

                    // Grab the file name
                    auto fname = capi::mxArrayToString(fname_);

                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxUnconstrained> (
                        Unconstrained::State::mxCreate());
                    
                    // Grab the reference vector 
                    auto x = Vector(X_,x_);
                    
                    // Create a C++ state
                    typename MxUnconstrained::State::t state(x);

                    // Read the restart file into the C++ state 
                    MxJsonUnconstrained::read_restart(fname,x,state);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

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

                    return names;
                }
                std::vector <char const *> fieldNames() {
                    std::vector <char const*> un
                        = Unconstrained::State::fieldNames_();
                    std::vector <char const*> eq
                        = EqualityConstrained::State::fieldNames_();
                    un.reserve(un.size()+eq.size());
                    un.insert(un.end(),eq.begin(),eq.end());
                    return un; 
                }

                // Create the structure for a Matlab state
                mxArrayPtr mxCreate() {
                    auto names = EqualityConstrained::State::fieldNames();
                    return capi::mxCreateStructMatrix(
                        1,1,names.size(),&(names[0]));
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
                    // Grab the variables
                    MX_VAR_4(X,Y,x,y);
                    
                    // Create a vector from the user input
                    auto x = Vector(X_,x_);
                    auto y = Vector(Y_,y_);

                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxEqualityConstrained> (
                        EqualityConstrained::State::mxCreate());

                    // Create a new C++ state
                    typename MxEqualityConstrained::State::t state(x,y);

                    // Convert the state to a Matlab state
                    mxstate_out.toMatlab(state);

                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

                } CATCH_MATLAB_ERRORS;
        
                // Read json parameters from file
                void readJson(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_4(X,Y,fname,mxstate);

                    // Grab the file name
                    auto fname = capi::mxArrayToString(fname_);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxEqualityConstrained> (
                        mxstate_);
                    auto mxstate_out = Matlab::State <MxEqualityConstrained> (
                        EqualityConstrained::State::mxCreate());
                
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    auto y_= capi::mxGetField(mxstate.data,0,"y");
                    auto y = Vector(Y_,y_);

                    // Create a new C++ state
                    typename MxEqualityConstrained::State::t state(x,y);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Read the JSON file into the C++ state
                    MxJsonEqualityConstrained::read(fname,state);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

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
                    // Grab the variables
                    MX_VAR_6(X,Y,msg,mxfns,mxstate,smanip);

                    // Create a messaging object
                    auto msg = Optizelle::Matlab::Messaging::matlab(msg_);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxEqualityConstrained> (
                        mxstate_);
                    auto mxstate_out = Matlab::State <MxEqualityConstrained> (
                        EqualityConstrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    auto y_= capi::mxGetField(mxstate.data,0,"y");
                    auto y = Vector(Y_,y_);

                    // Create a C++ state
                    typename MxEqualityConstrained::State::t state(x,y);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Create a Matlab bundle of functions
                    Matlab::Functions <MxEqualityConstrained> mxfns(
                        mxstate_out,
                        state,
                        mxfns_);

                    // Create a C++ bundle of functions
                    typename MxEqualityConstrained::Functions::t fns;
                    
                    // Convert the Matlab bundle of functions to C++ 
                    mxfns.fromMatlab(fns);
                    
                    // Create a state manipulator 
                    Matlab::StateManipulator <MxEqualityConstrained> smanip(
                        mxstate_out,
                        mxfns,
                        smanip_);
                   
                    // Minimize
                    MxEqualityConstrained::Algorithms::getMin(
                        msg,fns,state,smanip);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

                } CATCH_MATLAB_ERRORS;
            }
            
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                void release(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_3(X,Y,mxstate);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxEqualityConstrained> (
                        mxstate_);
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    auto y_= capi::mxGetField(mxstate.data,0,"y");
                    auto y = Vector(Y_,y_);

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
                    auto mxxs = capi::mxCreateCellMatrix(1,xs.size());
                    auto mxys = capi::mxCreateCellMatrix(1,ys.size());
                    auto mxreals = capi::mxCreateCellMatrix(1,reals.size());
                    auto mxnats = capi::mxCreateCellMatrix(1,nats.size());
                    auto mxparams = capi::mxCreateCellMatrix(1,params.size());

                    // Convert the restart information to Matlab 
                    toMatlab::Vectors(xs,mxxs);
                    toMatlab::Vectors(ys,mxys);
                    toMatlab::Reals(reals,mxreals);
                    toMatlab::Naturals(nats,mxnats);
                    toMatlab::Params(params,mxparams);
                    
                    // Return the result 
                    MX_RETURN_5(mxxs,mxys,mxreals,mxnats,mxparams);

                } CATCH_MATLAB_ERRORS;

                // Capture data from structures controlled by the user.  
                void capture(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_8(X,Y,mxstate,mxxs,mxys,mxreals,mxnats,mxparams);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxEqualityConstrained> (
                        mxstate_);
                    auto mxstate_out = Matlab::State <MxEqualityConstrained> (
                        EqualityConstrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    auto y_= capi::mxGetField(mxstate.data,0,"y");
                    auto y = Vector(Y_,y_);

                    // Create a C++ state
                    typename MxEqualityConstrained::State::t state(x,y);
                   
                    // Allocate memory for the released vectors
                    MxEqualityConstrained::Restart::X_Vectors xs;
                    MxEqualityConstrained::Restart::Y_Vectors ys;
                    MxEqualityConstrained::Restart::Reals reals;
                    MxEqualityConstrained::Restart::Naturals nats;
                    MxEqualityConstrained::Restart::Params params;
                    
                    // Convert the restart information from Matlab 
                    fromMatlab::Vectors(x,mxxs_,xs);
                    fromMatlab::Vectors(y,mxys_,ys);
                    fromMatlab::Reals(mxreals_,reals);
                    fromMatlab::Naturals(mxnats_,nats);
                    fromMatlab::Params(mxparams_,params);

                    // Do a capture 
                    MxEqualityConstrained::Restart
                        ::capture(state,xs,ys,reals,nats,params);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);
                   
                } CATCH_MATLAB_ERRORS;
                
                // Writes a json restart file
                void write_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_4(X,Y,fname,mxstate);

                    // Grab the file name
                    auto fname = capi::mxArrayToString(fname_);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxEqualityConstrained> (
                        mxstate_);
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    auto y_= capi::mxGetField(mxstate.data,0,"y");
                    auto y = Vector(Y_,y_);

                    // Create a C++ state
                    typename MxEqualityConstrained::State::t state(x,y);
                    
                    // Convert Matlab state to C++ 
                    mxstate.fromMatlab(state);

                    // Write the restart file
                    MxJsonEqualityConstrained::write_restart(fname,state);

                    // Return nothing
                    MX_RETURN_0;
                   
                } CATCH_MATLAB_ERRORS;
                
                // Reads a json restart file
                void read_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_5(X,Y,fname,x,y);

                    // Grab the file name
                    auto fname = capi::mxArrayToString(fname_);
                    
                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxEqualityConstrained> (
                        EqualityConstrained::State::mxCreate());
                    
                    // Grab the reference vectors
                    auto x = Vector(X_,x_);
                    auto y = Vector(Y_,y_);

                    // Create a C++ state
                    typename MxEqualityConstrained::State::t state(x,y);

                    // Read the restart file into the C++ state 
                    MxJsonEqualityConstrained::read_restart(fname,x,y,state);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

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

                    return names;
                }
                std::vector <char const *> fieldNames() {
                    std::vector <char const*> un
                        = Unconstrained::State::fieldNames_();
                    std::vector <char const*> iq
                        = InequalityConstrained::State::fieldNames_();
                    un.reserve(un.size()+iq.size());
                    un.insert(un.end(),iq.begin(),iq.end());
                    return un; 
                }

                // Create the structure for a Matlab state
                mxArrayPtr mxCreate() {
                    auto names = InequalityConstrained::State::fieldNames();
                    return capi::mxCreateStructMatrix(
                        1,1,names.size(),&(names[0]));
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
                    // Grab the variables
                    MX_VAR_4(X,Z,x,z);

                    // Create a vector from the user input
                    auto x = Vector(X_,x_);
                    auto z = Vector(Z_,z_);

                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxInequalityConstrained> (
                        InequalityConstrained::State::mxCreate());

                    // Create a new C++ state
                    typename MxInequalityConstrained::State::t state(x,z);

                    // Convert the state to a Matlab state
                    mxstate_out.toMatlab(state);

                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

                } CATCH_MATLAB_ERRORS;
        
                // Read json parameters from file
                void readJson(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_4(X,Z,fname,mxstate);

                    // Grab the file name
                    auto fname = capi::mxArrayToString(fname_);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxInequalityConstrained> (
                        mxstate_);
                    auto mxstate_out = Matlab::State <MxInequalityConstrained> (
                        InequalityConstrained::State::mxCreate());
                
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    auto z_= capi::mxGetField(mxstate.data,0,"z");
                    auto z = Vector(Z_,z_);

                    // Create a new C++ state
                    typename MxInequalityConstrained::State::t state(x,z);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Read the JSON file into the C++ state
                    MxJsonInequalityConstrained::read(fname,state);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);
                   
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
                    // Grab the variables
                    MX_VAR_6(X,Z,msg,mxfns,mxstate,smanip);
                    
                    // Create a messaging object
                    auto msg = Optizelle::Matlab::Messaging::matlab(msg_);
                        
                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxInequalityConstrained> (
                        mxstate_);
                    auto mxstate_out = Matlab::State <MxInequalityConstrained> (
                        InequalityConstrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    auto z_= capi::mxGetField(mxstate.data,0,"z");
                    auto z = Vector(Z_,z_);

                    // Create a C++ state
                    typename MxInequalityConstrained::State::t state(x,z);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Create a Matlab bundle of functions
                    Matlab::Functions <MxInequalityConstrained> mxfns(
                        mxstate_out,
                        state,
                        mxfns_);

                    // Create a C++ bundle of functions
                    typename MxInequalityConstrained::Functions::t fns;
                    
                    // Convert the Matlab bundle of functions to C++ 
                    mxfns.fromMatlab(fns);
                    
                    // Create a state manipulator 
                    Matlab::StateManipulator <MxInequalityConstrained> smanip(
                        mxstate_out,
                        mxfns,
                        smanip_);
                   
                    // Minimize
                    MxInequalityConstrained::Algorithms::getMin(
                        msg,fns,state,smanip);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

                } CATCH_MATLAB_ERRORS;
            }
            
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                void release(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_3(X,Z,mxstate);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxInequalityConstrained> (
                        mxstate_);
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    auto z_= capi::mxGetField(mxstate.data,0,"z");
                    auto z = Vector(Z_,z_);

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
                    auto mxxs = capi::mxCreateCellMatrix(1,xs.size());
                    auto mxzs = capi::mxCreateCellMatrix(1,zs.size());
                    auto mxreals = capi::mxCreateCellMatrix(1,reals.size());
                    auto mxnats = capi::mxCreateCellMatrix(1,nats.size());
                    auto mxparams = capi::mxCreateCellMatrix(1,params.size());

                    // Convert the restart information to Matlab 
                    toMatlab::Vectors(xs,mxxs);
                    toMatlab::Vectors(zs,mxzs);
                    toMatlab::Reals(reals,mxreals);
                    toMatlab::Naturals(nats,mxnats);
                    toMatlab::Params(params,mxparams);

                    // Return the result 
                    MX_RETURN_5(mxxs,mxzs,mxreals,mxnats,mxparams);
                   
                } CATCH_MATLAB_ERRORS;

                // Capture data from structures controlled by the user.  
                void capture(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_8(X,Z,mxstate,mxxs,mxzs,mxreals,mxnats,mxparams);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxInequalityConstrained> (
                        mxstate_);
                    auto mxstate_out = Matlab::State <MxInequalityConstrained> (
                        InequalityConstrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    auto z_= capi::mxGetField(mxstate.data,0,"z");
                    auto z = Vector(Z_,z_);

                    // Create a C++ state
                    typename MxInequalityConstrained::State::t state(x,z);
                   
                    // Allocate memory for the released vectors
                    MxInequalityConstrained::Restart::X_Vectors xs;
                    MxInequalityConstrained::Restart::Z_Vectors zs;
                    MxInequalityConstrained::Restart::Reals reals;
                    MxInequalityConstrained::Restart::Naturals nats;
                    MxInequalityConstrained::Restart::Params params;
                    
                    // Convert the restart information from Matlab 
                    fromMatlab::Vectors(x,mxxs_,xs);
                    fromMatlab::Vectors(z,mxzs_,zs);
                    fromMatlab::Reals(mxreals_,reals);
                    fromMatlab::Naturals(mxnats_,nats);
                    fromMatlab::Params(mxparams_,params);

                    // Do a capture 
                    MxInequalityConstrained::Restart
                        ::capture(state,xs,zs,reals,nats,params);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

                } CATCH_MATLAB_ERRORS;
                
                // Writes a json restart file
                void write_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_4(X,Z,fname,mxstate);

                    // Grab the file name
                    auto fname = capi::mxArrayToString(fname_);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxInequalityConstrained> (
                        mxstate_);
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    auto z_= capi::mxGetField(mxstate.data,0,"z");
                    auto z = Vector(Z_,z_);
                    
                    // Create a C++ state
                    typename MxInequalityConstrained::State::t state(x,z);
                    
                    // Convert Matlab state to C++ 
                    mxstate.fromMatlab(state);

                    // Write the restart file
                    MxJsonInequalityConstrained::write_restart(fname,state);

                    // Return nothing
                    MX_RETURN_0;
                   
                } CATCH_MATLAB_ERRORS;
                
                // Reads a json restart file
                void read_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_5(X,Z,fname,x,z);

                    // Grab the file name
                    auto fname = capi::mxArrayToString(fname_);

                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxInequalityConstrained> (
                        InequalityConstrained::State::mxCreate());
                    
                    // Grab the reference vectors
                    auto x = Vector(X_,x_);
                    auto z = Vector(Z_,z_);
                    
                    // Create a C++ state
                    typename MxInequalityConstrained::State::t state(x,z);

                    // Read the restart file into the C++ state 
                    MxJsonInequalityConstrained::read_restart(fname,x,z,state);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);
                   
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
                    return un; 
                }

                // Create the structure for a Matlab state
                mxArrayPtr mxCreate() {
                    auto names = Constrained::State::fieldNames();
                    return capi::mxCreateStructMatrix(
                        1,1,names.size(),&(names[0]));
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
                    // Grab the variables
                    MX_VAR_6(X,Y,Z,x,y,z);

                    // Create a vector from the user input
                    auto x = Vector(X_,x_);
                    auto y = Vector(Y_,y_);
                    auto z = Vector(Z_,z_);

                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxConstrained> (
                        Constrained::State::mxCreate());

                    // Create a new C++ state
                    typename MxConstrained::State::t state(x,y,z);

                    // Convert the state to a Matlab state
                    mxstate_out.toMatlab(state);

                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

                } CATCH_MATLAB_ERRORS;
        
                // Read json parameters from file
                void readJson(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_5(X,Y,Z,fname,mxstate);

                    // Grab the file name
                    auto fname = capi::mxArrayToString(fname_);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxConstrained> (mxstate_);
                    auto mxstate_out = Matlab::State <MxConstrained> (
                        Constrained::State::mxCreate());
                
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    auto y_= capi::mxGetField(mxstate.data,0,"y");
                    auto y = Vector(Y_,y_);
                    auto z_= capi::mxGetField(mxstate.data,0,"z");
                    auto z = Vector(Z_,z_);

                    // Create a new C++ state
                    typename MxConstrained::State::t state(x,y,z);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Read the JSON file into the C++ state
                    MxJsonConstrained::read(fname,state);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

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
                    // Grab the variables
                    MX_VAR_7(X,Y,Z,msg,mxfns,mxstate,smanip);

                    // Create a messaging object
                    auto msg = Optizelle::Matlab::Messaging::matlab(msg_);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxConstrained> (
                        mxstate_);
                    auto mxstate_out = Matlab::State <MxConstrained> (
                        Constrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    auto y_= capi::mxGetField(mxstate.data,0,"y");
                    auto y = Vector(Y_,y_);
                    auto z_= capi::mxGetField(mxstate.data,0,"z");
                    auto z = Vector(Z_,z_);

                    // Create a C++ state
                    typename MxConstrained::State::t state(x,y,z);
                    
                    // Convert the Matlab state to a C++ state
                    mxstate.fromMatlab(state);

                    // Create a Matlab bundle of functions
                    Matlab::Functions <MxConstrained> mxfns(
                        mxstate_out,
                        state,
                        mxfns_);

                    // Create a C++ bundle of functions
                    typename MxConstrained::Functions::t fns;
                    
                    // Convert the Matlab bundle of functions to C++ 
                    mxfns.fromMatlab(fns);
                    
                    // Create a state manipulator 
                    Matlab::StateManipulator <MxConstrained> smanip(
                        mxstate_out,
                        mxfns,
                        smanip_);
                   
                    // Minimize
                    MxConstrained::Algorithms::getMin(msg,fns,state,smanip);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

                } CATCH_MATLAB_ERRORS;
            }
            
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                void release(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_4(X,Y,Z,mxstate);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxConstrained> (mxstate_);
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    auto y_= capi::mxGetField(mxstate.data,0,"y");
                    auto y = Vector(Y_,y_);
                    auto z_= capi::mxGetField(mxstate.data,0,"z");
                    auto z = Vector(Z_,z_);

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
                    auto mxxs = capi::mxCreateCellMatrix(1,xs.size());
                    auto mxys = capi::mxCreateCellMatrix(1,ys.size());
                    auto mxzs = capi::mxCreateCellMatrix(1,zs.size());
                    auto mxreals = capi::mxCreateCellMatrix(1,reals.size());
                    auto mxnats = capi::mxCreateCellMatrix(1,nats.size());
                    auto mxparams = capi::mxCreateCellMatrix(1,params.size());

                    // Convert the restart information to Matlab 
                    toMatlab::Vectors(xs,mxxs);
                    toMatlab::Vectors(ys,mxys);
                    toMatlab::Vectors(zs,mxzs);
                    toMatlab::Reals(reals,mxreals);
                    toMatlab::Naturals(nats,mxnats);
                    toMatlab::Params(params,mxparams);
                    
                    // Return the result 
                    MX_RETURN_6(mxxs,mxys,mxzs,mxreals,mxnats,mxparams);

                } CATCH_MATLAB_ERRORS;

                // Capture data from structures controlled by the user.  
                void capture(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_10(X,Y,Z,mxstate,mxxs,mxys,mxzs,
                        mxreals,mxnats,mxparams);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxConstrained> (mxstate_);
                    auto mxstate_out = Matlab::State <MxConstrained> (
                        Constrained::State::mxCreate());
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    auto y_= capi::mxGetField(mxstate.data,0,"y");
                    auto y = Vector(Y_,y_);
                    auto z_= capi::mxGetField(mxstate.data,0,"z");
                    auto z = Vector(Z_,z_);

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
                    fromMatlab::Vectors(x,mxxs_,xs);
                    fromMatlab::Vectors(y,mxys_,ys);
                    fromMatlab::Vectors(z,mxzs_,zs);
                    fromMatlab::Reals(mxreals_,reals);
                    fromMatlab::Naturals(mxnats_,nats);
                    fromMatlab::Params(mxparams_,params);

                    // Do a capture 
                    MxConstrained::Restart
                        ::capture(state,xs,ys,zs,reals,nats,params);

                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

                } CATCH_MATLAB_ERRORS;
                
                // Writes a json restart file
                void write_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_5(X,Y,Z,fname,mxstate);

                    // Grab the file name
                    auto fname = capi::mxArrayToString(fname_);

                    // Create a Matlab state 
                    auto mxstate = Matlab::State <MxConstrained> (mxstate_);
                    
                    // Grab the base vectors from the Matlab state
                    auto x_= capi::mxGetField(mxstate.data,0,"x");
                    auto x = Vector(X_,x_);
                    auto y_= capi::mxGetField(mxstate.data,0,"y");
                    auto y = Vector(Y_,y_);
                    auto z_= capi::mxGetField(mxstate.data,0,"z");
                    auto z = Vector(Z_,z_);
                    
                    // Create a C++ state
                    typename MxConstrained::State::t state(x,y,z);
                    
                    // Convert Matlab state to C++ 
                    mxstate.fromMatlab(state);

                    // Write the restart file
                    MxJsonConstrained::write_restart(fname,state);

                    // Return nothing
                    MX_RETURN_0;

                } CATCH_MATLAB_ERRORS;
                
                // Reads a json restart file
                void read_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                ) try {
                    // Grab the variables
                    MX_VAR_7(X,Y,Z,fname,x,y,z);

                    // Grab the file name
                    auto fname = capi::mxArrayToString(fname_);

                    // Create a Matlab state 
                    auto mxstate_out = Matlab::State <MxConstrained> (
                        Constrained::State::mxCreate());
                    
                    // Grab the reference vectors
                    auto x = Vector(X_,x_);
                    auto y = Vector(Y_,y_);
                    auto z = Vector(Z_,z_);
                    
                    // Create a C++ state
                    typename MxConstrained::State::t state(x,y,z);

                    // Read the restart file into the C++ state 
                    MxJsonConstrained::read_restart(fname,x,y,z,state);
                    
                    // Convert the C++ state to a Matlab state
                    mxstate_out.toMatlab(state);
                    
                    // Return the result 
                    MX_RETURN_1(mxstate_out.data);

                } CATCH_MATLAB_ERRORS;
            }
        }
    }
}
