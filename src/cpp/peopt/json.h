#ifndef JSON_H
#define JSON_H

#include <fstream>
#include "peopt.h"
#include "json/json.h"

namespace peopt {
    using namespace peopt;
    struct json {
        // Parses a JSON file and returns the root
        static Json::Value parse(
            const peopt::Messaging& msg,
            const std::string fname
        ) {
            // Read in the input file
            Json::Value root;
            Json::Reader reader;
            std::ifstream file(fname.c_str(),std::ifstream::in);
            bool parsingSuccessful = reader.parse( file, root, true );
            if ( !parsingSuccessful ) 
                msg.error("Failed to parse the optimization parameter "
                    "file:  " + reader.getFormattedErrorMessages());

            // Close everything out and return the root
            file.close();
            return root;
        }

        template <typename Real,template <typename> class XX> 
        struct Unconstrained {
            // Read parameters from file
            static void read_(
                const peopt::Messaging& msg,
                const std::string fname,
                typename peopt::Unconstrained <Real,XX>::State::t& state
            ) {
                // Create a type shortcut
                typedef typename std::vector<Real>::size_type Natural;

                // Base error message
                const std::string base = "Invalid JSON parameter: ";

                // Read in the input file
                Json::Value root=parse(msg,fname);

                // Read in the parameters
                state.eps_g=Real(root["peopt"]
                    .get("eps_g",state.eps_g).asDouble());
                state.eps_dx=Real(root["peopt"]
                    .get("eps_dx",state.eps_dx).asDouble());
                state.stored_history=Natural(root["peopt"]
                    .get("stored_history",
                        Json::Value::UInt64(state.stored_history)).asUInt64());
                state.history_reset=Natural(root["peopt"]
                    .get("history_reset",
                        Json::Value::UInt64(state.history_reset)).asUInt64());
                state.iter_max=Natural(root["peopt"]
                    .get("iter_max",
                        Json::Value::UInt64(state.iter_max)).asUInt64());
                state.krylov_iter_max=Natural(root["peopt"]
                    .get("krylov_iter_max",
                        Json::Value::UInt64(state.krylov_iter_max)).asUInt64());
                state.krylov_orthog_max=Natural(root["peopt"]
                    .get("krylov_orthog_max",
                        Json::Value::UInt64(state.krylov_orthog_max))
                    .asUInt64());
                state.eps_krylov=Real(root["peopt"]
                    .get("eps_krylov",state.eps_krylov).asDouble());

                std::string krylov_solver=root["peopt"]
                    .get("krylov_solver",
                        KrylovSolverTruncated::to_string(state.krylov_solver))
                    .asString();
                if(KrylovSolverTruncated::is_valid()(krylov_solver))
                    state.krylov_solver=
                        KrylovSolverTruncated::from_string(krylov_solver); 
                else
                    msg.error(base + krylov_solver
                        + " is not a valid krylov_solver");

                std::string algorithm_class=root["peopt"]
                    .get("algorithm_class",
                        AlgorithmClass::to_string(state.algorithm_class))
                    .asString();
                if(AlgorithmClass::is_valid()(algorithm_class))
                    state.algorithm_class=
                        AlgorithmClass::from_string(algorithm_class); 
                else
                    msg.error(base + algorithm_class
                        + " is not a valid algorithm_class");

                std::string Minv_type=root["peopt"]
                    .get("Minv_type",
                        Operators::to_string(state.Minv_type))
                    .asString();
                if(Operators::is_valid()(Minv_type))
                    state.Minv_type=Operators::from_string(Minv_type); 
                else
                    msg.error(base + Minv_type + " is not a valid Minv_type");

                std::string H_type=root["peopt"]
                    .get("H_type",
                        Operators::to_string(state.H_type))
                    .asString();
                if(Operators::is_valid()(H_type))
                    state.H_type=Operators::from_string(H_type); 
                else
                    msg.error(base + H_type + " is not a valid H_type");


                state.msg_level=Natural(root["peopt"]
                    .get("msg_level",
                        Json::Value::UInt64(state.msg_level)).asUInt64());
                state.delta=Real(root["peopt"]
                    .get("delta",state.delta).asDouble());
                state.delta_max=Real(root["peopt"]
                    .get("delta_max",state.delta_max).asDouble());

                state.eta1=Real(root["peopt"]
                    .get("eta1",state.eta1).asDouble());
                state.eta2=Real(root["peopt"]
                    .get("eta2",state.eta2).asDouble());
                state.alpha=Real(root["peopt"]
                    .get("alpha",state.alpha).asDouble());
                state.linesearch_iter_max=Natural(root["peopt"]
                    .get("linesearch_iter_max",
                        Json::Value::UInt64(state.linesearch_iter_max))
                    .asUInt64());
                state.eps_ls=Real(root["peopt"]
                    .get("eps_ls",state.eps_ls).asDouble());
                
                std::string dir=root["peopt"]
                    .get("dir",LineSearchDirection::to_string(state.dir))
                    .asString();
                if(LineSearchDirection::is_valid()(dir))
                    state.dir=LineSearchDirection::from_string(dir); 
                else
                    msg.error(base + dir + " is not a valid dir.");
                
                std::string kind=root["peopt"]
                    .get("kind",LineSearchKind::to_string(state.kind))
                    .asString();
                if(LineSearchKind::is_valid()(kind))
                    state.kind=LineSearchKind::from_string(kind); 
                else
                    msg.error(base + kind + " is not a valid kind.");
            }
            static void read(
                const peopt::Messaging& msg,
                const std::string fname,
                typename peopt::Unconstrained <Real,XX>::State::t& state
            ) {
                Unconstrained <Real,XX>::read_(msg,fname,state);
            }

            // Convert parameters to a string 
            static std::string to_string_(
                typename peopt::Unconstrained <Real,XX>::State::t& state
            ) {
                // Create a new root for writing
                Json::Value root;

                // Write the optimization parameters
                root["peopt"]["eps_g"]=state.eps_g;
                root["peopt"]["eps_dx"]=state.eps_dx;
                root["peopt"]["stored_history"]
                    =Json::Value::UInt64(state.stored_history);
                root["peopt"]["history_reset"]
                    =Json::Value::UInt64(state.history_reset);
                root["peopt"]["iter_max"]
                    =Json::Value::UInt64(state.iter_max);
                root["peopt"]["krylov_iter_max"]
                    =Json::Value::UInt64(state.krylov_iter_max);
                root["peopt"]["krylov_orthog_max"]
                    =Json::Value::UInt64(state.krylov_orthog_max);
                root["peopt"]["eps_krylov"]=state.eps_krylov;
                root["peopt"]["krylov_solver"]
                    =KrylovSolverTruncated::to_string(state.krylov_solver);
                root["peopt"]["algorithm_class"]
                    =AlgorithmClass::to_string(state.algorithm_class);
                root["peopt"]["Minv_type"]
                    =Operators::to_string(state.Minv_type);
                root["peopt"]["H_type"]=Operators::to_string(state.H_type);
                root["peopt"]["msg_level"]=Json::Value::UInt64(state.msg_level);
                root["peopt"]["delta"]=state.delta;
                root["peopt"]["delta_max"]=state.delta_max;
                root["peopt"]["eta1"]=state.eta1;
                root["peopt"]["eta2"]=state.eta2;
                root["peopt"]["alpha"]=state.alpha;
                root["peopt"]["linesearch_iter_max"]
                    =Json::Value::UInt64(state.linesearch_iter_max);
                root["peopt"]["eps_ls"]=state.eps_ls;
                root["peopt"]["dir"]=LineSearchDirection::to_string(state.dir);
                root["peopt"]["kind"]=LineSearchKind::to_string(state.kind);

                // Create a string with the above output
                Json::StyledWriter writer;

                return writer.write(root);
            }
            static std::string to_string(
                typename peopt::Unconstrained <Real,XX>::State::t& state
            ) {
                return Unconstrained <Real,XX>::to_string_(state);
            }
        };

        template <
            typename Real,
            template <typename> class XX,
            template <typename> class YY
        > 
        struct EqualityConstrained {
            // Read parameters from file
            static void read_(
                const peopt::Messaging& msg,
                const std::string fname,
                typename peopt::EqualityConstrained <Real,XX,YY>::State::t&
                    state
            ) {
                // Read in the input file
                Json::Value root=parse(msg,fname);

                // Read in the parameters
                state.zeta=Real(root["peopt"]
                    .get("zeta",state.zeta).asDouble());
                state.eta0=Real(root["peopt"]
                    .get("eta0",state.eta0).asDouble());
                state.rho=Real(root["peopt"]
                    .get("rho",state.rho).asDouble());
                state.rho_bar=Real(root["peopt"]
                    .get("rho_bar",state.rho_bar).asDouble());
                state.eps_constr=Real(root["peopt"]
                    .get("eps_constr",state.eps_constr).asDouble());
                state.xi_all(Real(root["peopt"]
                    .get("xi_all",state.xi_qn).asDouble()));
                state.xi_qn=Real(root["peopt"]
                    .get("xi_qn",state.xi_qn).asDouble());
                state.xi_pg=Real(root["peopt"]
                    .get("xi_pg",state.xi_pg).asDouble());
                state.xi_proj=Real(root["peopt"]
                    .get("xi_proj",state.xi_proj).asDouble());
                state.xi_tang=Real(root["peopt"]
                    .get("xi_tang",state.xi_tang).asDouble());
                state.xi_lmh=Real(root["peopt"]
                    .get("xi_lmh",state.xi_lmh).asDouble());
                state.xi_lmg=Real(root["peopt"]
                    .get("xi_lmg",state.xi_lmg).asDouble());
                state.xi_4=Real(root["peopt"]
                    .get("xi_4",state.xi_4).asDouble());
                state.augsys_iter_max=Natural(root["peopt"]
                    .get("augsys_iter_max",
                        Json::Value::UInt64(state.augsys_iter_max)).asUInt64());
                state.augsys_rst_freq=Natural(root["peopt"]
                    .get("augsys_rst_freq",
                        Json::Value::UInt64(state.augsys_rst_freq)).asUInt64());
            }
            static void read(
                const peopt::Messaging& msg,
                const std::string fname,
                typename peopt::EqualityConstrained <Real,XX,YY>::State::t&
                    state
            ) {
                Unconstrained <Real,XX>::read_(msg,fname,state);
                EqualityConstrained <Real,XX,YY>::read_(msg,fname,state);
            }

            // Convert parameters to a string 
            static std::string to_string_(
                typename peopt::EqualityConstrained <Real,XX,YY>::State::t&
                    state
            ) {
                // Create a new root for writing
                Json::Value root;

                // Create a string with the above output
                Json::StyledWriter writer;

                // Write the optimization parameters
                root["peopt"]["zeta"]=state.zeta;
                root["peopt"]["eta0"]=state.eta0;
                root["peopt"]["rho"]=state.rho;
                root["peopt"]["rho_bar"]=state.rho_bar;
                root["peopt"]["eps_constr"]=state.eps_constr;
                root["peopt"]["xi_all"]=state.xi_qn;
                root["peopt"]["xi_qn"]=state.xi_qn;
                root["peopt"]["xi_pg"]=state.xi_pg;
                root["peopt"]["xi_proj"]=state.xi_proj;
                root["peopt"]["xi_tang"]=state.xi_tang;
                root["peopt"]["xi_lmh"]=state.xi_lmh;
                root["peopt"]["xi_lmg"]=state.xi_lmg;
                root["peopt"]["xi_4"]=state.xi_4;
                root["peopt"]["augsys_iter_max"]
                    =Json::Value::UInt64(state.augsys_iter_max);
                root["peopt"]["augsys_rst_freq"]
                    =Json::Value::UInt64(state.augsys_rst_freq);

                return writer.write(root);
            }
            static std::string to_string(
                typename peopt::EqualityConstrained <Real,XX,YY>::State::t&
                    state
            ) {
                std::string ucon
                    = Unconstrained <Real,XX>::to_string_(state);
                std::string econ
                    = EqualityConstrained <Real,XX,YY>::to_string_(state);
                return ucon.substr(0,ucon.size()-8)+",\n"+
                       econ.substr(17,econ.size());
            }
        };

        template < typename Real,
            template <typename> class XX,
            template <typename> class ZZ 
        > 
        struct InequalityConstrained {
            // Read parameters from file
            static void read_(
                const peopt::Messaging& msg,
                const std::string fname,
                typename peopt::InequalityConstrained <Real,XX,ZZ>::State::t&
                    state
            ) {
                // Base error message
                const std::string base = "Invalid JSON parameter: ";

                // Read in the input file
                Json::Value root=parse(msg,fname);

                // Read in the parameters
                state.mu=Real(root["peopt"]
                    .get("mu",state.mu).asDouble());
                state.eps_mu=Real(root["peopt"]
                    .get("eps_mu",state.eps_mu).asDouble());
                state.sigma=Real(root["peopt"]
                    .get("sigma",state.sigma).asDouble());
                state.gamma=Real(root["peopt"]
                    .get("gamma",state.gamma).asDouble());
                
                std::string ipm=root["peopt"]
                    .get("ipm",InteriorPointMethod::to_string(state.ipm))
                    .asString();
                if(InteriorPointMethod::is_valid()(ipm))
                    state.ipm=InteriorPointMethod::from_string(ipm); 
                else
                    msg.error(base + ipm + " is not a valid ipm.");

                std::string cstrat=root["peopt"]
                    .get("cstrat",CentralityStrategy::to_string(state.cstrat))
                    .asString();
                if(CentralityStrategy::is_valid()(cstrat))
                    state.cstrat=CentralityStrategy::from_string(cstrat); 
                else
                    msg.error(base + cstrat + " is not a valid cstrat.");
            }
            static void read(
                const peopt::Messaging& msg,
                const std::string fname,
                typename peopt::InequalityConstrained <Real,XX,ZZ>::State::t&
                    state
            ) {
                Unconstrained <Real,XX>::read_(msg,fname,state);
                InequalityConstrained <Real,XX,ZZ>::read_(msg,fname,state);
            }

            // Convert parameters to a string 
            static std::string to_string_(
                typename peopt::InequalityConstrained <Real,XX,ZZ>::State::t&
                    state
            ) {
                // Create a new root for writing
                Json::Value root;

                // Create a string with the above output
                Json::StyledWriter writer;
                
                // Write the optimization parameters
                root["peopt"]["mu"]=state.mu;
                root["peopt"]["eps_mu"]=state.eps_mu;
                root["peopt"]["sigma"]=state.sigma;
                root["peopt"]["gamma"]=state.gamma;
                root["peopt"]["ipm"]=InteriorPointMethod::to_string(state.ipm);
                root["peopt"]["cstrat"]
                    =CentralityStrategy::to_string(state.cstrat);

                return writer.write(root);
            }
            static std::string to_string(
                typename peopt::InequalityConstrained <Real,XX,ZZ>::State::t&
                    state
            ) {
                std::string ucon
                    = Unconstrained <Real,XX>::to_string_(state);
                std::string icon
                    = InequalityConstrained <Real,XX,ZZ>::to_string_(state);
                return ucon.substr(0,ucon.size()-8)+",\n"+
                       icon.substr(17,icon.size());
            }
        };

        template < typename Real,
            template <typename> class XX,
            template <typename> class YY, 
            template <typename> class ZZ 
        > 
        struct Constrained {
            // Read parameters from file
            static void read_(
                const peopt::Messaging& msg,
                const std::string fname,
                typename peopt::Constrained <Real,XX,YY,ZZ>::State::t&
                    state
            ) {
                // Read in the input file
                Json::Value root=parse(msg,fname);

                // Read in the parameters
            }
            static void read(
                const peopt::Messaging& msg,
                const std::string fname,
                typename peopt::Constrained <Real,XX,YY,ZZ>::State::t&
                    state
            ) {
                Unconstrained <Real,XX>::read_(msg,fname,state);
                EqualityConstrained <Real,XX,YY>::read_(msg,fname,state);
                InequalityConstrained <Real,XX,ZZ>::read_(msg,fname,state);
            }

            // Convert parameters to a string 
            static std::string to_string_(
                typename peopt::Constrained <Real,XX,YY,ZZ>::State::t&
                    state
            ) {
                // Create a new root for writing
                Json::Value root;

                // Create a string with the above output
                Json::StyledWriter writer;
                
                // Write the optimization parameters

                return writer.write(root);
            }
            static std::string to_string(
                typename peopt::Constrained <Real,XX,YY,ZZ>::State::t&
                    state
            ) {
                std::string ucon
                    = Unconstrained <Real,XX>::to_string_(state);
                std::string econ
                    = EqualityConstrained <Real,XX,YY>::to_string_(state);
                std::string icon
                    = InequalityConstrained <Real,XX,ZZ>::to_string_(state);
                return ucon.substr(0,ucon.size()-8)+",\n"+
                       econ.substr(17,econ.size()-8)+",\n";
                       icon.substr(17,icon.size());
            }
        };
    };
}

#endif
