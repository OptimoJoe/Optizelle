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
                // Base error message
                const std::string base = "Invalid JSON parameter: ";

                // Read in the input file
                Json::Value root=parse(msg,fname);

                // Read in the parameters
                state.eps_g=Real(root["peopt"]
                    .get("eps_g",state.eps_g).asDouble());
                state.eps_dx=Real(root["peopt"]
                    .get("eps_dx",state.eps_dx).asDouble());
                state.stored_history=root["peopt"]
                    .get("stored_history",state.stored_history).asInt();
                state.history_reset=root["peopt"]
                    .get("history_reset",state.history_reset).asInt();
                state.iter_max=root["peopt"]
                    .get("iter_max",state.iter_max).asInt();
                state.krylov_iter_max=root["peopt"]
                    .get("krylov_iter_max",state.krylov_iter_max).asInt();
                state.eps_krylov=Real(root["peopt"]
                    .get("eps_krylov",state.eps_krylov).asDouble());

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


                state.msg_level=root["peopt"]
                    .get("msg_level",state.msg_level).asInt();
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
                state.linesearch_iter_max=root["peopt"]
                    .get("linesearch_iter_max",state.linesearch_iter_max)
                    .asInt();
                state.eps_ls=Real(root["peopt"]
                    .get("eps_ls",state.eps_ls).asDouble());
                
                std::string dir=root["peopt"]
                    .get("dir",LineSearchDirection::to_string(state.dir))
                    .asString();
                if(LineSearchDirection::is_valid()(dir))
                    state.dir=LineSearchDirection::from_string(dir); 
                else
                    msg.error(base + dir + " is not a valid dir");
                
                std::string kind=root["peopt"]
                    .get("kind",LineSearchKind::to_string(state.kind))
                    .asString();
                if(LineSearchKind::is_valid()(kind))
                    state.kind=LineSearchKind::from_string(kind); 
                else
                    msg.error(base + kind + " is not a valid kind");
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
                root["peopt"]["stored_history"]=state.stored_history;
                root["peopt"]["history_reset"]=state.history_reset;
                root["peopt"]["iter_max"]=state.iter_max;
                root["peopt"]["krylov_iter_max"]=state.krylov_iter_max;
                root["peopt"]["eps_krylov"]=state.eps_krylov;
                root["peopt"]["algorithm_class"]=
                    AlgorithmClass::to_string(state.algorithm_class);
                root["peopt"]["Minv_type"]=
                    Operators::to_string(state.Minv_type);
                root["peopt"]["H_type"]=
                    Operators::to_string(state.H_type);
                root["peopt"]["msg_level"]=state.msg_level;
                root["peopt"]["delta"]=state.delta;
                root["peopt"]["delta_max"]=state.delta_max;
                root["peopt"]["eta1"]=state.eta1;
                root["peopt"]["eta2"]=state.eta2;
                root["peopt"]["alpha"]=state.alpha;
                root["peopt"]["linesearch_iter_max"]=state.linesearch_iter_max;
                root["peopt"]["eps_ls"]=state.eps_ls;
                root["peopt"]["dir"]=
                    LineSearchDirection::to_string(state.dir);
                root["peopt"]["kind"]=
                    LineSearchKind::to_string(state.kind);

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
                // Base error message
                const std::string base = "Invalid JSON parameter: ";

                // Read in the input file
                Json::Value root=parse(msg,fname);

                // Read in the parameters
                state.zeta=Real(root["peopt"]
                    .get("zeta",state.zeta).asDouble());
                state.eta0=Real(root["peopt"]
                    .get("eta0",state.zeta).asDouble());
                state.rho=Real(root["peopt"]
                    .get("rho",state.zeta).asDouble());
                state.rho_bar=Real(root["peopt"]
                    .get("rho_bar",state.zeta).asDouble());
                state.eps_constr=Real(root["peopt"]
                    .get("eps_constr",state.zeta).asDouble());
                state.xi_all=Real(root["peopt"]
                    .get("xi_all",state.zeta).asDouble());
                state.xi_qn=Real(root["peopt"]
                    .get("xi_qn",state.zeta).asDouble());
                state.xi_pg=Real(root["peopt"]
                    .get("xi_pg",state.zeta).asDouble());
                state.xi_proj=Real(root["peopt"]
                    .get("xi_proj",state.zeta).asDouble());
                state.xi_tang=Real(root["peopt"]
                    .get("xi_tang",state.zeta).asDouble());
                state.xi_lmh=Real(root["peopt"]
                    .get("xi_lmh",state.zeta).asDouble());
                state.xi_lmg=Real(root["peopt"]
                    .get("xi_lmg",state.zeta).asDouble());
                state.xi_4=Real(root["peopt"]
                    .get("xi_4",state.zeta).asDouble());
                
                std::string Minv_right_type=root["peopt"]
                    .get("Minv_right_type",
                        Operators::to_string(state.Minv_right_type))
                    .asString();
                if(Operators::is_valid()(Minv_right_type))
                    state.Minv_right_type
                        = Operators::from_string(Minv_right_type);
                else
                    msg.error(base + Minv_right_type
                        + " is not a valid Minv_right_type");
                
                std::string Minv_left_type=root["peopt"]
                    .get("Minv_left_type",
                        Operators::to_string(state.Minv_left_type))
                    .asString();
                if(Operators::is_valid()(Minv_left_type))
                    state.Minv_left_type
                        = Operators::from_string(Minv_left_type);
                else
                    msg.error(base + Minv_left_type
                        + " is not a valid Minv_left_type");
                
                state.augsys_iter_max=Real(root["peopt"]
                    .get("augsys_iter_max",state.augsys_iter_max).asInt());
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
                root["peopt"]["xi_all"]=state.xi_all;
                root["peopt"]["xi_qn"]=state.xi_qn;
                root["peopt"]["xi_pg"]=state.xi_pg;
                root["peopt"]["xi_proj"]=state.xi_proj;
                root["peopt"]["xi_tang"]=state.xi_tang;
                root["peopt"]["xi_lmh"]=state.xi_lmh;
                root["peopt"]["xi_lmg"]=state.xi_lmg;
                root["peopt"]["xi_4"]=state.xi_4;
                root["peopt"]["Minv_right_type"]=
                    Operators::to_string(state.Minv_right_type);
                root["peopt"]["Minv_left_type"]=
                    Operators::to_string(state.Minv_left_type);
                root["peopt"]["augsys_iter_max"]=state.augsys_iter_max;

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
                // Read in the input file
                Json::Value root=parse(msg,fname);

                // Read in the parameters
                state.mu_tol=Real(root["peopt"]
                    .get("mu_tol",state.mu_tol).asDouble());
                state.sigma=Real(root["peopt"]
                    .get("sigma",state.sigma).asDouble());
                state.gamma=Real(root["peopt"]
                    .get("gamma",state.gamma).asDouble());
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
                root["peopt"]["mu_tol"]=state.mu_tol;
                root["peopt"]["sigma"]=state.sigma;
                root["peopt"]["gamma"]=state.gamma;

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
