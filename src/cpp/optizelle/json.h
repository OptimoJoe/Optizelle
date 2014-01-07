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

#ifndef JSON_H
#define JSON_H

#include <fstream>
#include "optizelle/optizelle.h"
#include "json/json.h"

namespace Optizelle {
    using namespace Optizelle;
    namespace json {
        // Parses a JSON file and returns the root
        Json::Value parse(
            Optizelle::Messaging const & msg,
            std::string const & fname
        ); 
       
        // Writes a JSON spec to file
        void write(
            Optizelle::Messaging const & msg,
            std::string const & fname,
            Json::Value const & root
        ); 

        // A helper class to help with serialization of vectors into json
        // objects.
        template <typename Real,template <typename> class XX>
        struct Serialization {
            static std::string serialize(
                typename XX <Real>::Vector const & x
            ) { }
            static typename XX <Real>::Vector deserialize(
                std::string const & x_json
            ) { }
        };

        // Routines to serialize lists of elements for restarting
        namespace Serialize{
            // Vectors 
            template <typename Real,template <typename> class XX>
            void vectors(
                std::pair <
                    std::list <std::string>,
                    std::list <typename XX <Real>::Vector>
                > const & xs,
                std::string const & vs,
                Json::Value & root
            ) {
                // Create some type shortcuts
                typedef XX <Real> X;
                typedef typename X::Vector X_Vector;

                // Create a reader object to parse a json tree
                Json::Reader reader;

                // Loop over all the vectors and serialize things
                typename std::list <X_Vector>::const_iterator x
                    =xs.second.begin();
                for(typename std::list <std::string>::const_iterator
                        name=xs.first.begin();
                    name!=xs.first.end();
                    name++, x++
                ) {
                    // Grab the json string of the vector
                    std::string x_json_(Serialization <Real,XX>::serialize(*x));
                   
                    // Parse the string
                    Json::Value x_json;
                    reader.parse(x_json_,x_json,true);

                    // Insert the information into the correct place
                    root[vs][*name]=x_json;
                }
            }

            // Reals 
            template <typename Real>
            void reals(
                std::pair <
                    std::list <std::string>,
                    std::list <Real>
                > const & reals,
                std::string const & vs,
                Json::Value & root
            ) {
                // Loop over all the reals and serialize things
                typename std::list <Real>::const_iterator real 
                    =reals.second.begin();
                for(typename std::list <std::string>::const_iterator
                        name=reals.first.begin();
                    name!=reals.first.end();
                    name++, real++
                )
                    root[vs][*name]=*real;
            }

            // Naturals 
            void naturals(
                std::pair <
                    std::list <std::string>,
                    std::list <Natural>
                > const & nats,
                std::string const & vs,
                Json::Value & root
            );

            // Parameters 
            void parameters(
                std::pair <
                    std::list <std::string>,
                    std::list <std::string>
                > const & params,
                std::string const & vs,
                Json::Value & root
            );
        }
        
        // Routines to deserialize lists of elements for restarting
        namespace Deserialize{

            // Vectors
            template <typename Real,template <typename> class XX>
            static void vectors(
                typename XX <Real>::Vector const & x,
                Json::Value const & root,
                std::string const & vs,
                std::pair <
                    std::list <std::string>,
                    std::list <typename XX <Real>::Vector>
                > & xs
            ) {
                // Create some type shortcuts
                typedef XX <Real> X;
                typedef typename X::Vector X_Vector;
                
                // Create a writer so that we can tranlate json objects into
                // strings
                Json::StyledWriter writer;

                // Loop over all the names in the root
                for(Json::ValueIterator itr=root[vs].begin();
                    itr!=root[vs].end();
                    itr++
                ){
                    // Grab the vector
                    std::string name(itr.key().asString());
                    xs.first.emplace_back(name);
                    xs.second.emplace_back(std::move(
                        Serialization <Real,XX>::deserialize(
                            writer.write(root[vs][name]))));
                }
            }
            
            // Reals 
            template <typename Real>
            void reals(
                Json::Value const & root,
                std::string const & vs,
                std::pair <
                    std::list <std::string>,
                    std::list <Real>
                > & reals
            ) {
                // Loop over all the names in the root
                for(Json::ValueIterator itr=root[vs].begin();
                    itr!=root[vs].end();
                    itr++
                ){
                    // Grab the real 
                    reals.first.emplace_back(itr.key().asString());
                    reals.second.emplace_back(
                        Real(root[vs][reals.first.back()].asDouble()));
                }
            }
            
            // Naturals 
            void naturals(
                Json::Value const & root,
                std::string const & vs,
                std::pair <
                    std::list <std::string>,
                    std::list <Natural>
                > & nats
            );
            
            // Parameters 
            void parameters(
                Json::Value const & root,
                std::string const & vs,
                std::pair <
                    std::list <std::string>,
                    std::list <std::string>
                > & params
            );
        }

        template <typename Real,template <typename> class XX> 
        struct Unconstrained {
            // Create some type shortcuts
            typedef typename Optizelle::Unconstrained <Real,XX>::X_Vectors
                X_Vectors; 
            typedef typename Optizelle::Unconstrained <Real,XX>::Reals Reals;
            typedef typename Optizelle::Unconstrained <Real,XX>::Nats Nats;
            typedef typename Optizelle::Unconstrained <Real,XX>::Params Params; 

            // Read parameters from file
            static void read_(
                Optizelle::Messaging const & msg,
                std::string const & fname,
                typename Optizelle::Unconstrained <Real,XX>::State::t & state
            ) {
                // Base error message
                std::string const base = "Invalid JSON parameter: ";

                // Read in the input file
                Json::Value root=parse(msg,fname);

                // Read in the parameters
                state.eps_grad=Real(root["Optizelle"]
                    .get("eps_grad",state.eps_grad).asDouble());
                state.eps_dx=Real(root["Optizelle"]
                    .get("eps_dx",state.eps_dx).asDouble());
                state.stored_history=Natural(root["Optizelle"]
                    .get("stored_history",
                        Json::Value::UInt64(state.stored_history)).asUInt64());
                state.history_reset=Natural(root["Optizelle"]
                    .get("history_reset",
                        Json::Value::UInt64(state.history_reset)).asUInt64());
                state.iter_max=Natural(root["Optizelle"]
                    .get("iter_max",
                        Json::Value::UInt64(state.iter_max)).asUInt64());
                state.krylov_iter_max=Natural(root["Optizelle"]
                    .get("krylov_iter_max",
                        Json::Value::UInt64(state.krylov_iter_max)).asUInt64());
                state.krylov_orthog_max=Natural(root["Optizelle"]
                    .get("krylov_orthog_max",
                        Json::Value::UInt64(state.krylov_orthog_max))
                    .asUInt64());
                state.eps_krylov=Real(root["Optizelle"]
                    .get("eps_krylov",state.eps_krylov).asDouble());

                std::string krylov_solver=root["Optizelle"]
                    .get("krylov_solver",
                        KrylovSolverTruncated::to_string(state.krylov_solver))
                    .asString();
                if(KrylovSolverTruncated::is_valid()(krylov_solver))
                    state.krylov_solver=
                        KrylovSolverTruncated::from_string(krylov_solver); 
                else
                    msg.error(base + krylov_solver
                        + " is not a valid krylov_solver");

                std::string algorithm_class=root["Optizelle"]
                    .get("algorithm_class",
                        AlgorithmClass::to_string(state.algorithm_class))
                    .asString();
                if(AlgorithmClass::is_valid()(algorithm_class))
                    state.algorithm_class=
                        AlgorithmClass::from_string(algorithm_class); 
                else
                    msg.error(base + algorithm_class
                        + " is not a valid algorithm_class");

                std::string PH_type=root["Optizelle"]
                    .get("PH_type",Operators::to_string(state.PH_type))
                    .asString();
                if(Operators::is_valid()(PH_type))
                    state.PH_type=Operators::from_string(PH_type); 
                else
                    msg.error(base + PH_type + " is not a valid PH_type");

                std::string H_type=root["Optizelle"]
                    .get("H_type",Operators::to_string(state.H_type))
                    .asString();
                if(Operators::is_valid()(H_type))
                    state.H_type=Operators::from_string(H_type); 
                else
                    msg.error(base + H_type + " is not a valid H_type");

                state.msg_level=Natural(root["Optizelle"]
                    .get("msg_level",
                        Json::Value::UInt64(state.msg_level)).asUInt64());
                state.delta=Real(root["Optizelle"]
                    .get("delta",state.delta).asDouble());

                state.eta1=Real(root["Optizelle"]
                    .get("eta1",state.eta1).asDouble());
                state.eta2=Real(root["Optizelle"]
                    .get("eta2",state.eta2).asDouble());
                state.alpha0=Real(root["Optizelle"]
                    .get("alpha0",state.alpha0).asDouble());
                state.c1=Real(root["Optizelle"].get("c1",state.c1).asDouble());
                state.linesearch_iter_max=Natural(root["Optizelle"]
                    .get("linesearch_iter_max",
                        Json::Value::UInt64(state.linesearch_iter_max))
                    .asUInt64());
                state.eps_ls=Real(root["Optizelle"]
                    .get("eps_ls",state.eps_ls).asDouble());
                
                std::string dir=root["Optizelle"]
                    .get("dir",LineSearchDirection::to_string(state.dir))
                    .asString();
                if(LineSearchDirection::is_valid()(dir))
                    state.dir=LineSearchDirection::from_string(dir); 
                else
                    msg.error(base + dir + " is not a valid dir.");
                
                std::string kind=root["Optizelle"]
                    .get("kind",LineSearchKind::to_string(state.kind))
                    .asString();
                if(LineSearchKind::is_valid()(kind))
                    state.kind=LineSearchKind::from_string(kind); 
                else
                    msg.error(base + kind + " is not a valid kind.");
            }
            static void read(
                Optizelle::Messaging const & msg,
                std::string const & fname,
                typename Optizelle::Unconstrained <Real,XX>::State::t& state
            ) {
                Unconstrained <Real,XX>::read_(msg,fname,state);
            }

            // Convert parameters to a string 
            static std::string to_string_(
                typename Optizelle::Unconstrained <Real,XX>::State::t & state
            ) {
                // Create a new root for writing
                Json::Value root;

                // Write the optimization parameters
                root["Optizelle"]["eps_grad"]=state.eps_grad;
                root["Optizelle"]["eps_dx"]=state.eps_dx;
                root["Optizelle"]["stored_history"]
                    =Json::Value::UInt64(state.stored_history);
                root["Optizelle"]["history_reset"]
                    =Json::Value::UInt64(state.history_reset);
                root["Optizelle"]["iter_max"]
                    =Json::Value::UInt64(state.iter_max);
                root["Optizelle"]["krylov_iter_max"]
                    =Json::Value::UInt64(state.krylov_iter_max);
                root["Optizelle"]["krylov_orthog_max"]
                    =Json::Value::UInt64(state.krylov_orthog_max);
                root["Optizelle"]["eps_krylov"]=state.eps_krylov;
                root["Optizelle"]["krylov_solver"]
                    =KrylovSolverTruncated::to_string(state.krylov_solver);
                root["Optizelle"]["algorithm_class"]
                    =AlgorithmClass::to_string(state.algorithm_class);
                root["Optizelle"]["PH_type"]
                    =Operators::to_string(state.PH_type);
                root["Optizelle"]["H_type"]=Operators::to_string(state.H_type);
                root["Optizelle"]["msg_level"]
                    =Json::Value::UInt64(state.msg_level);
                root["Optizelle"]["delta"]=state.delta;
                root["Optizelle"]["eta1"]=state.eta1;
                root["Optizelle"]["eta2"]=state.eta2;
                root["Optizelle"]["alpha0"]=state.alpha0;
                root["Optizelle"]["c1"]=state.c1;
                root["Optizelle"]["linesearch_iter_max"]
                    =Json::Value::UInt64(state.linesearch_iter_max);
                root["Optizelle"]["eps_ls"]=state.eps_ls;
                root["Optizelle"]["dir"]
                    =LineSearchDirection::to_string(state.dir);
                root["Optizelle"]["kind"]=LineSearchKind::to_string(state.kind);

                // Create a string with the above output
                Json::StyledWriter writer;

                return writer.write(root);
            }
            static std::string to_string(
                typename Optizelle::Unconstrained <Real,XX>::State::t& state
            ) {
                return Unconstrained <Real,XX>::to_string_(state);
            }

            // Write all parameters to file
            static void write_restart(
                Optizelle::Messaging const & msg,
                std::string const & fname,
                typename Optizelle::Unconstrained <Real,XX>::State::t & state
            ) {
                // Do a release 
                X_Vectors xs;
                Reals reals;
                Nats nats;
                Params params;
                Optizelle::Unconstrained <Real,XX>::Restart::release(
                    state,xs,reals,nats,params);

                // Serialize everything
                Json::Value root;
                Serialize::vectors <Real,XX>(xs,"X_Vectors",root);
                Serialize::reals <Real> (reals,"Reals",root);
                Serialize::naturals(nats,"Naturals",root);
                Serialize::parameters(params,"Parameters",root);
                
                // Create a string with the above output
                write(msg,fname,root);

                // Recapture the state
                Optizelle::Unconstrained <Real,XX>::Restart::capture(
                    msg,state,xs,reals,nats,params);
            }

            // Read all the parameters from file
            static void read_restart(
                Optizelle::Messaging const & msg,
                std::string const & fname,
                typename Optizelle::Unconstrained <Real,XX>::State::t & state
            ) {
                // Read in the input file
                Json::Value root=parse(msg,fname);

                // Extract everything from the parsed json file 
                X_Vectors xs;
                Reals reals;
                Nats nats;
                Params params;
                Deserialize::vectors <Real,XX>(root,"X_Vectors",xs);
                Deserialize::reals <Real> (root,"Reals",reals);
                Deserialize::naturals(root,"Naturals",nats);
                Deserialize::parameters(root,"Parameters",params);
               
                // Move this information into the state
                Optizelle::Unconstrained <Real,XX>::Restart::capture(
                    state,xs,reals,nats,params);
            }
        };

        template <
            typename Real,
            template <typename> class XX,
            template <typename> class YY
        > 
        struct EqualityConstrained {
            // Create some type shortcuts
            typedef
                typename Optizelle::EqualityConstrained<Real,XX,YY>::X_Vectors
                X_Vectors; 
            typedef
                typename Optizelle::EqualityConstrained<Real,XX,YY>::Y_Vectors
                Y_Vectors; 
            typedef typename Optizelle::EqualityConstrained <Real,XX,YY>::Reals
                Reals;
            typedef typename Optizelle::EqualityConstrained <Real,XX,YY>::Nats
                Nats;
            typedef typename Optizelle::EqualityConstrained <Real,XX,YY>::Params
                Params; 

            // Read parameters from file
            static void read_(
                Optizelle::Messaging const & msg,
                std::string const & fname,
                typename Optizelle::EqualityConstrained <Real,XX,YY>::State::t &
                    state
            ) {
                // Base error message
                std::string const base = "Invalid JSON parameter: ";

                // Read in the input file
                Json::Value root=parse(msg,fname);

                // Read in the parameters
                state.zeta=Real(root["Optizelle"]
                    .get("zeta",state.zeta).asDouble());
                state.eta0=Real(root["Optizelle"]
                    .get("eta0",state.eta0).asDouble());
                state.rho=Real(root["Optizelle"]
                    .get("rho",state.rho).asDouble());
                state.rho_bar=Real(root["Optizelle"]
                    .get("rho_bar",state.rho_bar).asDouble());
                state.eps_constr=Real(root["Optizelle"]
                    .get("eps_constr",state.eps_constr).asDouble());
                state.xi_all(Real(root["Optizelle"]
                    .get("xi_all",state.xi_qn).asDouble()));
                state.xi_qn=Real(root["Optizelle"]
                    .get("xi_qn",state.xi_qn).asDouble());
                state.xi_pg=Real(root["Optizelle"]
                    .get("xi_pg",state.xi_pg).asDouble());
                state.xi_proj=Real(root["Optizelle"]
                    .get("xi_proj",state.xi_proj).asDouble());
                state.xi_tang=Real(root["Optizelle"]
                    .get("xi_tang",state.xi_tang).asDouble());
                state.xi_lmh=Real(root["Optizelle"]
                    .get("xi_lmh",state.xi_lmh).asDouble());
                state.xi_lmg=Real(root["Optizelle"]
                    .get("xi_lmg",state.xi_lmg).asDouble());
                state.xi_4=Real(root["Optizelle"]
                    .get("xi_4",state.xi_4).asDouble());
                state.augsys_iter_max=Natural(root["Optizelle"]
                    .get("augsys_iter_max",
                        Json::Value::UInt64(state.augsys_iter_max)).asUInt64());
                state.augsys_rst_freq=Natural(root["Optizelle"]
                    .get("augsys_rst_freq",
                        Json::Value::UInt64(state.augsys_rst_freq)).asUInt64());
                std::string PSchur_left_type=root["Optizelle"]
                    .get("PSchur_left_type",
                        Operators::to_string(state.PSchur_left_type))
                    .asString();
                if(Operators::is_valid()(PSchur_left_type))
                    state.PSchur_left_type
                        = Operators::from_string(PSchur_left_type); 
                else
                    msg.error(base + PSchur_left_type
                        + " is not a valid PSchur_left_type");

                std::string PSchur_right_type=root["Optizelle"]
                    .get("PSchur_right_type",
                        Operators::to_string(state.PSchur_right_type))
                    .asString();
                if(Operators::is_valid()(PSchur_right_type))
                    state.PSchur_right_type
                        = Operators::from_string(PSchur_right_type); 
                else
                    msg.error(base + PSchur_right_type
                        + " is not a valid PSchur_right_type");
            }
            static void read(
                Optizelle::Messaging const & msg,
                std::string const & fname,
                typename Optizelle::EqualityConstrained <Real,XX,YY>::State::t &
                    state
            ) {
                Unconstrained <Real,XX>::read_(msg,fname,state);
                EqualityConstrained <Real,XX,YY>::read_(msg,fname,state);
            }

            // Convert parameters to a string 
            static std::string to_string_(
                typename Optizelle::EqualityConstrained <Real,XX,YY>::State::t &
                    state
            ) {
                // Create a new root for writing
                Json::Value root;

                // Create a string with the above output
                Json::StyledWriter writer;

                // Write the optimization parameters
                root["Optizelle"]["zeta"]=state.zeta;
                root["Optizelle"]["eta0"]=state.eta0;
                root["Optizelle"]["rho"]=state.rho;
                root["Optizelle"]["rho_bar"]=state.rho_bar;
                root["Optizelle"]["eps_constr"]=state.eps_constr;
                root["Optizelle"]["xi_all"]=state.xi_qn;
                root["Optizelle"]["xi_qn"]=state.xi_qn;
                root["Optizelle"]["xi_pg"]=state.xi_pg;
                root["Optizelle"]["xi_proj"]=state.xi_proj;
                root["Optizelle"]["xi_tang"]=state.xi_tang;
                root["Optizelle"]["xi_lmh"]=state.xi_lmh;
                root["Optizelle"]["xi_lmg"]=state.xi_lmg;
                root["Optizelle"]["xi_4"]=state.xi_4;
                root["Optizelle"]["augsys_iter_max"]
                    =Json::Value::UInt64(state.augsys_iter_max);
                root["Optizelle"]["augsys_rst_freq"]
                    =Json::Value::UInt64(state.augsys_rst_freq);
                root["Optizelle"]["PSchur_left_type"]
                    =Operators::to_string(state.PSchur_left_type);
                root["Optizelle"]["PSchur_right_type"]
                    =Operators::to_string(state.PSchur_right_type);

                return writer.write(root);
            }
            static std::string to_string(
                typename Optizelle::EqualityConstrained <Real,XX,YY>::State::t &
                    state
            ) {
                std::string ucon
                    = Unconstrained <Real,XX>::to_string_(state);
                std::string econ
                    = EqualityConstrained <Real,XX,YY>::to_string_(state);
                return ucon.substr(0,ucon.size()-8)+",\n"+
                       econ.substr(17,econ.size());
            }

            // Write all parameters to file
            static void write_restart(
                Optizelle::Messaging const & msg,
                std::string const & fname,
                typename Optizelle::EqualityConstrained <Real,XX,YY>::State::t &
                    state
            ) {
                // Do a release 
                X_Vectors xs;
                Y_Vectors ys;
                Reals reals;
                Nats nats;
                Params params;
                Optizelle::EqualityConstrained <Real,XX,YY>::Restart::release(
                    state,xs,ys,reals,nats,params);

                // Serialize everything
                Json::Value root;
                Serialize::vectors <Real,XX>(xs,"X_Vectors",root);
                Serialize::vectors <Real,YY>(ys,"Y_Vectors",root);
                Serialize::reals <Real> (reals,"Reals",root);
                Serialize::naturals(nats,"Naturals",root);
                Serialize::parameters(params,"Parameters",root);
                
                // Create a string with the above output
                write(msg,fname,root);

                // Recapture the state
                Optizelle::EqualityConstrained<Real,XX,YY>::Restart::capture(
                    msg,state,xs,ys,reals,nats,params);
            }

            // Read all the parameters from file
            static void read_restart(
                Optizelle::Messaging const & msg,
                std::string const & fname,
                typename Optizelle::EqualityConstrained <Real,XX,YY>::State::t &
                    state
            ) {
                // Read in the input file
                Json::Value root=parse(msg,fname);

                // Extract everything from the parsed json file 
                X_Vectors xs;
                Y_Vectors ys;
                Reals reals;
                Nats nats;
                Params params;
                Deserialize::vectors <Real,XX>(root,"X_Vectors",xs);
                Deserialize::vectors <Real,YY>(root,"Y_Vectors",ys);
                Deserialize::reals <Real> (root,"Reals",reals);
                Deserialize::naturals(root,"Naturals",nats);
                Deserialize::parameters(root,"Parameters",params);
               
                // Move this information into the state
                Optizelle::EqualityConstrained <Real,XX,YY>::Restart::capture(
                    state,xs,ys,reals,nats,params);
            }
        };

        template < typename Real,
            template <typename> class XX,
            template <typename> class ZZ 
        > 
        struct InequalityConstrained {
            // Create some type shortcuts
            typedef 
                typename Optizelle::InequalityConstrained<Real,XX,ZZ>::X_Vectors
                X_Vectors; 
            typedef
                typename Optizelle::InequalityConstrained<Real,XX,ZZ>::Z_Vectors
                Z_Vectors; 
            typedef 
                typename Optizelle::InequalityConstrained <Real,XX,ZZ>::Reals
                Reals;
            typedef typename Optizelle::InequalityConstrained <Real,XX,ZZ>::Nats
                Nats;
            typedef 
                typename Optizelle::InequalityConstrained <Real,XX,ZZ>::Params
                Params; 

            // Read parameters from file
            static void read_(
                Optizelle::Messaging const & msg,
                std::string const & fname,
                typename Optizelle::InequalityConstrained<Real,XX,ZZ>::State::t&
                    state
            ) {
                // Base error message
                std::string const base = "Invalid JSON parameter: ";

                // Read in the input file
                Json::Value root=parse(msg,fname);

                // Read in the parameters
                state.eps_mu=Real(root["Optizelle"]
                    .get("eps_mu",state.eps_mu).asDouble());
                state.sigma=Real(root["Optizelle"]
                    .get("sigma",state.sigma).asDouble());
                state.gamma=Real(root["Optizelle"]
                    .get("gamma",state.gamma).asDouble());
                
                std::string ipm=root["Optizelle"]
                    .get("ipm",InteriorPointMethod::to_string(state.ipm))
                    .asString();
                if(InteriorPointMethod::is_valid()(ipm))
                    state.ipm=InteriorPointMethod::from_string(ipm); 
                else
                    msg.error(base + ipm + " is not a valid ipm.");

                std::string cstrat=root["Optizelle"]
                    .get("cstrat",CentralityStrategy::to_string(state.cstrat))
                    .asString();
                if(CentralityStrategy::is_valid()(cstrat))
                    state.cstrat=CentralityStrategy::from_string(cstrat); 
                else
                    msg.error(base + cstrat + " is not a valid cstrat.");
            }
            static void read(
                Optizelle::Messaging const & msg,
                std::string const & fname,
                typename Optizelle::InequalityConstrained<Real,XX,ZZ>::State::t&
                    state
            ) {
                Unconstrained <Real,XX>::read_(msg,fname,state);
                InequalityConstrained <Real,XX,ZZ>::read_(msg,fname,state);
            }

            // Convert parameters to a string 
            static std::string to_string_(
                typename Optizelle::InequalityConstrained<Real,XX,ZZ>::State::t&
                    state
            ) {
                // Create a new root for writing
                Json::Value root;

                // Create a string with the above output
                Json::StyledWriter writer;
                
                // Write the optimization parameters
                root["Optizelle"]["eps_mu"]=state.eps_mu;
                root["Optizelle"]["sigma"]=state.sigma;
                root["Optizelle"]["gamma"]=state.gamma;
                root["Optizelle"]["ipm"]
                    =InteriorPointMethod::to_string(state.ipm);
                root["Optizelle"]["cstrat"]
                    =CentralityStrategy::to_string(state.cstrat);

                return writer.write(root);
            }
            static std::string to_string(
                typename Optizelle::InequalityConstrained<Real,XX,ZZ>::State::t&
                    state
            ) {
                std::string ucon
                    = Unconstrained <Real,XX>::to_string_(state);
                std::string icon
                    = InequalityConstrained <Real,XX,ZZ>::to_string_(state);
                return ucon.substr(0,ucon.size()-8)+",\n"+
                       icon.substr(17,icon.size());
            }

            // Write all parameters to file
            static void write_restart(
                Optizelle::Messaging const & msg,
                std::string const & fname,
                typename Optizelle::InequalityConstrained<Real,XX,ZZ>::State::t&
                    state
            ) {
                // Do a release 
                X_Vectors xs;
                Z_Vectors zs;
                Reals reals;
                Nats nats;
                Params params;
                Optizelle::InequalityConstrained <Real,XX,ZZ>::Restart::release(
                    state,xs,zs,reals,nats,params);

                // Serialize everything
                Json::Value root;
                Serialize::vectors <Real,XX>(xs,"X_Vectors",root);
                Serialize::vectors <Real,ZZ>(zs,"Z_Vectors",root);
                Serialize::reals <Real> (reals,"Reals",root);
                Serialize::naturals(nats,"Naturals",root);
                Serialize::parameters(params,"Parameters",root);
                
                // Create a string with the above output
                write(msg,fname,root);

                // Recapture the state
                Optizelle::InequalityConstrained<Real,XX,ZZ>::Restart::capture(
                    msg,state,xs,zs,reals,nats,params);
            }

            // Read all the parameters from file
            static void read_restart(
                Optizelle::Messaging const & msg,
                std::string const & fname,
                typename Optizelle::InequalityConstrained<Real,XX,ZZ>::State::t&
                    state
            ) {
                // Read in the input file
                Json::Value root=parse(msg,fname);

                // Extract everything from the parsed json file 
                X_Vectors xs;
                Z_Vectors zs;
                Reals reals;
                Nats nats;
                Params params;
                Deserialize::vectors <Real,XX>(root,"X_Vectors",xs);
                Deserialize::vectors <Real,ZZ>(root,"Z_Vectors",zs);
                Deserialize::reals <Real> (root,"Reals",reals);
                Deserialize::naturals(root,"Naturals",nats);
                Deserialize::parameters(root,"Parameters",params);
               
                // Move this information into the state
                Optizelle::InequalityConstrained <Real,XX,ZZ>::Restart::capture(
                    state,xs,zs,reals,nats,params);
            }
        };

        template < typename Real,
            template <typename> class XX,
            template <typename> class YY, 
            template <typename> class ZZ 
        > 
        struct Constrained {
            // Create some type shortcuts
            typedef typename Optizelle::Constrained<Real,XX,YY,ZZ>::X_Vectors
                X_Vectors; 
            typedef typename Optizelle::Constrained<Real,XX,YY,ZZ>::Y_Vectors
                Y_Vectors; 
            typedef typename Optizelle::Constrained<Real,XX,YY,ZZ>::Z_Vectors
                Z_Vectors; 
            typedef typename Optizelle::Constrained <Real,XX,YY,ZZ>::Reals
                Reals;
            typedef typename Optizelle::Constrained <Real,XX,YY,ZZ>::Nats
                Nats;
            typedef typename Optizelle::Constrained <Real,XX,YY,ZZ>::Params
                Params; 

            // Read parameters from file
            static void read(
                Optizelle::Messaging const & msg,
                std::string const & fname,
                typename Optizelle::Constrained <Real,XX,YY,ZZ>::State::t&
                    state
            ) {
                Unconstrained <Real,XX>::read_(msg,fname,state);
                EqualityConstrained <Real,XX,YY>::read_(msg,fname,state);
                InequalityConstrained <Real,XX,ZZ>::read_(msg,fname,state);
            }

            // Convert parameters to a string 
            static std::string to_string(
                typename Optizelle::Constrained <Real,XX,YY,ZZ>::State::t&
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
            
            // Write all parameters to file
            static void write_restart(
                Optizelle::Messaging const & msg,
                std::string const & fname,
                typename Optizelle::Constrained <Real,XX,YY,ZZ>::State::t &
                    state
            ) {
                // Do a release 
                X_Vectors xs;
                Y_Vectors ys;
                Z_Vectors zs;
                Reals reals;
                Nats nats;
                Params params;
                Optizelle::Constrained <Real,XX,YY,ZZ>::Restart::release(
                    state,xs,ys,zs,reals,nats,params);

                // Serialize everything
                Json::Value root;
                Serialize::vectors <Real,XX>(xs,"X_Vectors",root);
                Serialize::vectors <Real,YY>(ys,"Y_Vectors",root);
                Serialize::vectors <Real,ZZ>(zs,"Z_Vectors",root);
                Serialize::reals <Real> (reals,"Reals",root);
                Serialize::naturals(nats,"Naturals",root);
                Serialize::parameters(params,"Parameters",root);
                
                // Create a string with the above output
                write(msg,fname,root);

                // Recapture the state
                Optizelle::Constrained<Real,XX,YY,ZZ>::Restart::capture(
                    msg,state,xs,ys,zs,reals,nats,params);
            }
            
            // Read all the parameters from file
            static void read_restart(
                Optizelle::Messaging const & msg,
                std::string const & fname,
                typename Optizelle::Constrained <Real,XX,YY,ZZ>::State::t& state
            ) {
                // Read in the input file
                Json::Value root=parse(msg,fname);

                // Extract everything from the parsed json file 
                X_Vectors xs;
                Y_Vectors ys;
                Z_Vectors zs;
                Reals reals;
                Nats nats;
                Params params;
                Deserialize::vectors <Real,XX>(root,"X_Vectors",xs);
                Deserialize::vectors <Real,YY>(root,"Y_Vectors",ys);
                Deserialize::vectors <Real,ZZ>(root,"Z_Vectors",zs);
                Deserialize::reals <Real> (root,"Reals",reals);
                Deserialize::naturals(root,"Naturals",nats);
                Deserialize::parameters(root,"Parameters",params);
               
                // Move this information into the state
                Optizelle::Constrained <Real,XX,YY,ZZ>::Restart::capture(
                    state,xs,ys,zs,reals,nats,params);
            }
        };
    }
}

#endif
