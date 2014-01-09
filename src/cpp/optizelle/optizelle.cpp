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

#include<optizelle/optizelle.h>

namespace Optizelle{

    // Prints a message
    void Messaging::print(std::string const & msg) const {
        std::cout << msg << std::endl;
    }

    // Prints an error
    void Messaging::error(std::string const & msg) const {
        std::cerr << msg << std::endl;
        exit(EXIT_FAILURE);
    }

    // Allow a derived class to deallocate memory
    Messaging::~Messaging() {}

    // Converts the algorithm class to a string
    std::string AlgorithmClass::to_string(
        AlgorithmClass::t const & algorithm_class
    ){
        switch(algorithm_class){
        case TrustRegion:
            return "TrustRegion";
        case LineSearch:
            return "LineSearch";
        case UserDefined:
            return "UserDefined";
        default:
            throw;
        }
    }
    
    // Converts a string to an algorithm class 
    AlgorithmClass::t AlgorithmClass::from_string(
        std::string const & algorithm_class
    ){
        if(algorithm_class=="TrustRegion")
            return TrustRegion;
        else if(algorithm_class=="LineSearch")
            return LineSearch;
        else if(algorithm_class=="UserDefined")
            return UserDefined;
        else
            throw;
    }

    // Checks whether or not a string is valid
    bool AlgorithmClass::is_valid::operator () (std::string const & name)const {
        if( name=="TrustRegion" ||
            name=="LineSearch" ||
            name=="UserDefined"
        )
            return true;
        else
            return false;
    }

    // Converts the stopping condition to a string 
    std::string StoppingCondition::to_string(
        StoppingCondition::t const & opt_stop
    ){
        switch(opt_stop){
        case NotConverged:
            return "NotConverged";
        case RelativeGradientSmall:
            return "RelativeGradientSmall";
        case RelativeStepSmall:
            return "RelativeStepSmall";
        case MaxItersExceeded:
            return "MaxItersExceeded";
        case InteriorPointInstability:
            return "InteriorPointInstability";
        case UserDefined:
            return "UserDefined";
        default:
            throw;
        }
    }

    // Converts a string to a stopping condition
    StoppingCondition::t StoppingCondition::from_string(
        std::string const & opt_stop
    ){
        if(opt_stop=="NotConverged")
            return StoppingCondition::NotConverged;
        else if(opt_stop=="RelativeGradientSmall")
            return StoppingCondition::RelativeGradientSmall;
        else if(opt_stop=="RelativeStepSmall")
            return StoppingCondition::RelativeStepSmall;
        else if(opt_stop=="MaxItersExceeded")
            return StoppingCondition::MaxItersExceeded;
        else if(opt_stop=="InteriorPointInstability")
            return StoppingCondition::InteriorPointInstability; 
        else if(opt_stop=="UserDefined")
            return StoppingCondition::UserDefined;
        else
            throw;
    }

    // Checks whether or not a string is valid
    bool StoppingCondition::is_valid::operator () (
        std::string const & name
    ) const {
        if( name=="NotConverged" ||
            name=="RelativeGradientSmall" ||
            name=="RelativeStepSmall" ||
            name=="MaxItersExceeded" ||
            name=="InteriorPointInstability" ||
            name=="UserDefined"
        )
            return true;
        else
            return false;
    }

    // Converts the operator type to a string 
    std::string Operators::to_string(Operators::t const & op){
        switch(op){
        case Operators::Identity:
            return "Identity";
        case Operators::ScaledIdentity:
            return "ScaledIdentity";
        case Operators::BFGS:
            return "BFGS";
        case Operators::InvBFGS:
            return "InvBFGS";
        case Operators::SR1:
            return "SR1";
        case Operators::InvSR1:
            return "InvSR1";
        case Operators::UserDefined:
            return "UserDefined";
        default:
            throw;
        }
    }
    
    // Converts a string to a operator 
    Operators::t Operators::from_string(std::string const & op){
        if(op=="Identity")
            return Operators::Identity; 
        else if(op=="ScaledIdentity")
            return Operators::ScaledIdentity; 
        else if(op=="BFGS")
            return Operators::BFGS; 
        else if(op=="InvBFGS")
            return Operators::InvBFGS; 
        else if(op=="SR1")
            return Operators::SR1; 
        else if(op=="InvSR1")
            return Operators::InvSR1; 
        else if(op=="UserDefined")
            return Operators::UserDefined; 
        else
            throw;
    }

    // Checks whether or not a string is valid
    bool Operators::is_valid::operator () (std::string const & name) const {
        if( name=="Identity" ||
            name=="ScaledIdentity" ||
            name=="BFGS" ||
            name=="InvBFGS" ||
            name=="SR1" ||
            name=="InvSR1" ||
            name=="UserDefined" 
        )
            return true;
        else
            return false;
    }

    // Converts the line-search direction to a string 
    std::string LineSearchDirection::to_string(t const & dir){
        switch(dir){
        case LineSearchDirection::SteepestDescent:
            return "SteepestDescent";
        case LineSearchDirection::FletcherReeves:
            return "FletcherReeves";
        case LineSearchDirection::PolakRibiere:
            return "PolakRibiere";
        case LineSearchDirection::HestenesStiefel:
            return "HestenesStiefel";
        case LineSearchDirection::BFGS:
            return "BFGS";
        case LineSearchDirection::NewtonCG:
            return "NewtonCG";
        default:
            throw;
        }
    }
    
    // Converts a string to a line-search direction 
    LineSearchDirection::t LineSearchDirection::from_string(
        std::string const & dir
    ){
        if(dir=="SteepestDescent")
            return LineSearchDirection::SteepestDescent; 
        else if(dir=="FletcherReeves")
            return LineSearchDirection::FletcherReeves; 
        else if(dir=="PolakRibiere")
            return LineSearchDirection::PolakRibiere; 
        else if(dir=="HestenesStiefel")
            return LineSearchDirection::HestenesStiefel; 
        else if(dir=="BFGS")
            return LineSearchDirection::BFGS; 
        else if(dir=="NewtonCG")
            return LineSearchDirection::NewtonCG; 
        else
            throw;
    }

    // Checks whether or not a string is valid
    bool LineSearchDirection::is_valid::operator () (
        std::string const & name
    ) const {
        if( name=="SteepestDescent" ||
            name=="FletcherReeves" ||
            name=="PolakRibiere" ||
            name=="HestenesStiefel" ||
            name=="BFGS" ||
            name=="NewtonCG"
        )
            return true;
        else
            return false;
    }

    // Converts the line-search kind to a string 
    std::string LineSearchKind::to_string(t const & kind){
        switch(kind){
        case LineSearchKind::Brents:
            return "Brents";
        case LineSearchKind::GoldenSection:
            return "GoldenSection";
        case LineSearchKind::BackTracking:
            return "BackTracking";
        case LineSearchKind::TwoPointA:
            return "TwoPointA";
        case LineSearchKind::TwoPointB:
            return "TwoPointB";
        default:
            throw;
        }
    }
    
    // Converts a string to a line-search kind 
    LineSearchKind::t LineSearchKind::from_string(std::string const & kind){
        if(kind=="Brents")
            return LineSearchKind::Brents; 
        else if(kind=="GoldenSection")
            return LineSearchKind::GoldenSection; 
        else if(kind=="BackTracking")
            return LineSearchKind::BackTracking; 
        else if(kind=="TwoPointA")
            return LineSearchKind::TwoPointA; 
        else if(kind=="TwoPointB")
            return LineSearchKind::TwoPointB; 
        else
            throw;
    }

    // Checks whether or not a string is valid
    bool LineSearchKind::is_valid::operator () (std::string const & name)const {
        if( name=="Brents" ||
            name=="GoldenSection" ||
            name=="BackTracking" ||
            name=="TwoPointA" ||
            name=="TwoPointB" 
        )
            return true;
        else
            return false;
    }

    // Determine whether or not the line-search checks the sufficient
    // decrease condition.
    bool LineSearchKind::is_sufficient_decrease(LineSearchKind::t const & kind){
        switch(kind){
        case LineSearchKind::GoldenSection:
        case LineSearchKind::BackTracking:
        case LineSearchKind::Brents:
            return true;
        case LineSearchKind::TwoPointA:
        case LineSearchKind::TwoPointB:
            return false;
        default:
            throw;
        }
    }
    
    // Converts the optimization location to a string 
    std::string OptimizationLocation::to_string(
        OptimizationLocation::t const & loc
    ){
        switch(loc){
        case OptimizationLocation::BeginningOfOptimization:
            return "BeginningOfOptimization";
        case OptimizationLocation::BeforeInitialFuncAndGrad:
            return "BeforeInitialFuncAndGrad";
        case OptimizationLocation::AfterInitialFuncAndGrad:
            return "AfterInitialFuncAndGrad";
        case OptimizationLocation::BeforeOptimizationLoop:
            return "BeforeOptimizationLoop";
        case OptimizationLocation::BeforeSaveOld:
            return "BeforeSaveOld";
        case OptimizationLocation::BeforeStep:
            return "BeforeStep";
        case OptimizationLocation::BeforeGetStep:
            return "BeforeGetStep";
        case OptimizationLocation::GetStep:
            return "GetStep";
        case OptimizationLocation::AfterStepBeforeGradient:
            return "AfterStepBeforeGradient";
        case OptimizationLocation::AfterGradient:
            return "AfterGradient";
        case OptimizationLocation::BeforeQuasi:
            return "BeforeQuasi";
        case OptimizationLocation::AfterQuasi:
            return "AfterQuasi";
        case OptimizationLocation::EndOfOptimizationIteration:
            return "EndOfOptimizationIteration";
        case OptimizationLocation::BeforeLineSearch:
            return "BeforeLineSearch";
        case OptimizationLocation::AfterRejectedTrustRegion:
            return "AfterRejectedTrustRegion";
        case OptimizationLocation::AfterRejectedLineSearch:
            return "AfterRejectedLineSearch";
        case OptimizationLocation::BeforeActualVersusPredicted:
            return "BeforeActualVersusPredicted";
        case OptimizationLocation::EndOfKrylovIteration:
            return "EndOfKrylovIteration";
        case OptimizationLocation::EndOfOptimization:
            return "EndOfOptimization";
        default:
            throw;
        }
    }
    
    // Converts a string to a line-search kind 
    OptimizationLocation::t OptimizationLocation::from_string(
        std::string const & loc
    ){
        if(loc=="BeginningOfOptimization")
            return OptimizationLocation::BeginningOfOptimization;
        else if(loc=="BeforeInitialFuncAndGrad")
            return OptimizationLocation::BeforeInitialFuncAndGrad;
        else if(loc=="AfterInitialFuncAndGrad")
            return OptimizationLocation::AfterInitialFuncAndGrad;
        else if(loc=="BeforeOptimizationLoop")
            return OptimizationLocation::BeforeOptimizationLoop;
        else if(loc=="BeforeSaveOld")
            return OptimizationLocation::BeforeSaveOld; 
        else if(loc=="BeforeStep")
            return OptimizationLocation::BeforeStep; 
        else if(loc=="BeforeGetStep")
            return OptimizationLocation::BeforeGetStep; 
        else if(loc=="GetStep")
            return OptimizationLocation::GetStep; 
        else if(loc=="AfterStepBeforeGradient")
            return OptimizationLocation::AfterStepBeforeGradient; 
        else if(loc=="AfterGradient")
            return OptimizationLocation::AfterGradient;
        else if(loc=="BeforeQuasi")
            return OptimizationLocation::BeforeQuasi;
        else if(loc=="AfterQuasi")
            return OptimizationLocation::AfterQuasi;
        else if(loc=="EndOfOptimizationIteration")
            return OptimizationLocation::EndOfOptimizationIteration; 
        else if(loc=="BeforeLineSearch")
            return OptimizationLocation::BeforeLineSearch;
        else if(loc=="AfterRejectedTrustRegion")
            return OptimizationLocation::AfterRejectedTrustRegion;
        else if(loc=="AfterRejectedLineSearch")
            return OptimizationLocation::AfterRejectedLineSearch;
        else if(loc=="BeforeActualVersusPredicted")
            return OptimizationLocation::BeforeActualVersusPredicted;
        else if(loc=="EndOfKrylovIteration")
            return OptimizationLocation::EndOfKrylovIteration;
        else if(loc=="EndOfOptimization")
            return OptimizationLocation::EndOfOptimization;
        else
            throw;
    }

    // Checks whether or not a string is valid
    bool OptimizationLocation::is_valid::operator () (
        std::string const & name
    ) const {
        if( name=="BeginningOfOptimization" ||
            name=="BeforeInitialFuncAndGrad" ||
            name=="AfterInitialFuncAndGrad" ||
            name=="BeforeOptimizationLoop" ||
            name=="BeforeSaveOld" || 
            name=="BeforeStep" || 
            name=="BeforeGetStep" || 
            name=="GetStep" || 
            name=="AfterStepBeforeGradient" ||
            name=="AfterGradient" ||
            name=="BeforeQuasi" ||
            name=="AfterQuasi" ||
            name=="EndOfOptimizationIteration" ||
            name=="BeforeLineSearch" ||
            name=="AfterRejectedTrustRegion" ||
            name=="AfterRejectedLineSearch" ||
            name=="BeforeActualVersusPredicted" ||
            name=="EndOfKrylovIteration" ||
            name=="EndOfOptimization"
        )
            return true;
        else
            return false;
    }
    
    // Converts the problem class to a string
    std::string ProblemClass::to_string(ProblemClass::t const & problem_class){
        switch(problem_class){
        case ProblemClass::Unconstrained:
            return "Unconstrained";
        case ProblemClass::EqualityConstrained:
            return "EqualityConstrained";
        case ProblemClass::InequalityConstrained:
            return "InequalityConstrained";
        case ProblemClass::Constrained:
            return "Constrained";
        default:
                throw;
        }
    }

    // Converts a string to a problem class 
    ProblemClass::t ProblemClass::from_string(std::string const &problem_class){
        if(problem_class=="Unconstrained")
            return ProblemClass::Unconstrained;
        else if(problem_class=="EqualityConstrained")
            return ProblemClass::EqualityConstrained;
        else if(problem_class=="InequalityConstrained")
            return ProblemClass::InequalityConstrained;
        else if(problem_class=="Constrained")
            return ProblemClass::Constrained;
        else
            throw;
    }

    // Checks whether or not a string is valid
    bool ProblemClass::is_valid::operator () (std::string const & name) const {
        if( name=="Unconstrained" ||
            name=="EqualityConstrained" ||
            name=="InequalityConstrained" ||
            name=="Constrained" 
        )
            return true;
        else
            return false;
    }
    
    // Converts the Krylov solver to a string
    std::string KrylovSolverTruncated::to_string(
        const KrylovSolverTruncated::t& truncated_krylov
    ){
        switch(truncated_krylov){
        case KrylovSolverTruncated::ConjugateDirection:
            return "ConjugateDirection";
        case KrylovSolverTruncated::MINRES:
            return "MINRES";
        default:
                throw;
        }
    }

    // Converts a string to a problem class 
    KrylovSolverTruncated::t KrylovSolverTruncated::from_string(
        std::string const & truncated_krylov
    ){
        if(truncated_krylov=="ConjugateDirection")
            return KrylovSolverTruncated::ConjugateDirection;
        else if(truncated_krylov=="MINRES")
            return KrylovSolverTruncated::MINRES;
        else
            throw;
    }

    // Checks whether or not a string is valid
    bool KrylovSolverTruncated::is_valid::operator () (
        std::string const & name
    ) const {
        if( name=="ConjugateDirection" ||
            name=="MINRES" 
        )
            return true;
        else
            return false;
    }
    
    // Converts the interior point method to a string 
    std::string InteriorPointMethod::to_string(
        InteriorPointMethod::t const & ipm
    ){
        switch(ipm){
        case InteriorPointMethod::PrimalDual:
            return "PrimalDual";
        case InteriorPointMethod::PrimalDualLinked:
            return "PrimalDualLinked";
        case InteriorPointMethod::LogBarrier:
            return "LogBarrier";
        default:
            throw;
        }
    }
    
    // Converts a string to an interior point method 
    InteriorPointMethod::t InteriorPointMethod::from_string(
        std::string const & ipm
    ){
        if(ipm=="PrimalDual")
            return InteriorPointMethod::PrimalDual; 
        else if(ipm=="PrimalDualLinked")
            return InteriorPointMethod::PrimalDualLinked; 
        else if(ipm=="LogBarrier")
            return InteriorPointMethod::LogBarrier; 
        else
            throw;
    }

    // Checks whether or not a string is valid
    bool InteriorPointMethod::is_valid::operator () (
        std::string const & name
    ) const {
        if( name=="PrimalDual" ||
            name=="PrimalDualLinked" ||
            name=="LogBarrier"
        )
            return true;
        else
            return false;
    }
    
    // Converts the centrality strategy to a string
    std::string CentralityStrategy::to_string(
        CentralityStrategy::t const & cstrat
    ){
        switch(cstrat){
        case CentralityStrategy::Constant:
            return "Constant";
        case CentralityStrategy::StairStep:
            return "StairStep";
        case CentralityStrategy::PredictorCorrector:
            return "PredictorCorrector";
        default:
            throw;
        }
    }
    
    // Converts a string to the cstrat
    CentralityStrategy::t CentralityStrategy::from_string(
        std::string const & cstrat
    ){
        if(cstrat=="Constant")
            return CentralityStrategy::Constant; 
        else if(cstrat=="StairStep")
            return CentralityStrategy::StairStep; 
        else if(cstrat=="PredictorCorrector")
            return CentralityStrategy::PredictorCorrector; 
        else
            throw;
    }

    // Checks whether or not a string is valid
    bool CentralityStrategy::is_valid::operator () (
        std::string const & name
    ) const {
        if( name=="Constant" ||
            name=="StairStep" ||
            name=="PredictorCorrector" 
        )
            return true;
        else
            return false;
    }

    // Converts a variety of basic datatypes to strings
    std::ostream& Utility::formatReal(std::ostream & out) {
        return out<<std::setprecision(2) << std::scientific << std::setw(10)
            << std::left;
    }
    std::ostream& Utility::formatInt(std::ostream & out) {
        return out << std::setw(10) << std::left;
    }
    std::ostream& Utility::formatString(std::ostream & out) {
        return out << std::setw(10) << std::left;
    }
    
    // Converts anything to a string.
    std::string Utility::atos(double const & x){
        std::stringstream ss;
        ss << formatReal << x;
        return ss.str();
    }
    std::string Utility::atos(Natural const & x){
        std::stringstream ss;
        ss << formatInt << x;
        return ss.str();
    }
    std::string Utility::atos(std::string const & x){
        std::stringstream ss;
        ss << formatString << x;
        return ss.str();
    }
    std::string Utility::atos(KrylovStop::t const & x){
        std::stringstream ss;
        // Converts the Krylov stopping condition to a shorter string 
        switch(x){
        case KrylovStop::NegativeCurvature:
            return atos("NegCurv");
        case KrylovStop::RelativeErrorSmall:
            return atos("RelErrSml");
        case KrylovStop::MaxItersExceeded:
            return atos("IterExcd");
        case KrylovStop::TrustRegionViolated:
            return atos("TrstReg");
        case KrylovStop::Instability:
            return atos("Unstable");
        case KrylovStop::InvalidTrustRegionCenter:
            return atos("InvldCnt");
        default:
            throw;
        }
    }
}
