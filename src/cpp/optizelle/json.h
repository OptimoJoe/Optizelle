#pragma once

#include <fstream>
#include <typeinfo>
#include "optizelle/optizelle.h"
#include "json/json.h"

namespace Optizelle {
    using namespace Optizelle;
    namespace json {
        // Parses a JSON file and returns the root
        Json::Value parse(
            std::string const & fname
        ); 
       
        // Writes a JSON spec to file
        void write_to_file(
            std::string const & fname,
            Json::Value const & root
        ); 

        // Safely reads from a json tree 
        namespace read {
            // Read a real
            template <typename Real>
            Real real(
                Json::Value const & json,
                std::string const & name
            ) {
                // Set the error message
                std::string const err_msg = "Invalid JSON parameter: "
                    + name + " contains an invalid real.";

                // Check if we have a NaN or Inf 
                if(json.isString()) {
                    // Extract the string
                    std::string val=json.asString();

                    // Figure out if we have NaN, Inf, or -Inf
                    if(val=="NaN")
                        return std::numeric_limits <Real>::quiet_NaN();
                    else if(val=="Inf")
                        return std::numeric_limits <Real>::infinity();
                    else if(val=="-Inf")
                        return -std::numeric_limits <Real>::infinity();
                    else
                        throw Exception::t(__LOC__ + ", " + err_msg);

                // As long as we have a number, grab it
                } else if (json.isNumeric())
                    return Real(json.asDouble());

                // Anything else is an error
                else
                    throw Exception::t(__LOC__ + ", " + err_msg);
                
                // We should not hit this point
                throw;
            }

            // Read a natural 
            Natural natural(
                Json::Value const & json,
                std::string const & name
            );

            // Read a paramter
            template <typename enum_t>
            enum_t param(
                Json::Value const & json,
                std::function<bool(std::string const &)> const & is_valid,
                std::function<enum_t(std::string const&)>const& from_string,
                std::string const & name
            ) {

                // Set the error message
                std::string const err_msg = "Invalid JSON parameter: "
                    + name + " contains an invalid parameter.";

                // All parameters start off as strings
                if(json.isString()) {

                    // Grab the string
                    std::string val=json.asString();

                    // If we have a valid parameter, return it 
                    if(is_valid(val))
                        return from_string(val);

                    // Otherwise, raise an error
                    throw Exception::t(__LOC__ + ", " + err_msg);

                // If we don't start with a string, raise an error
                } else
                    throw Exception::t(__LOC__ + ", " + err_msg);
                
                // We should not hit this point
                throw;
            }

            // Read a string 
            std::string string(
                Json::Value const & json,
                std::string const & name
            );
        }
        
        // Writes into a json tree 
        namespace write {
            // Write a real
            template <typename Real>
            Json::Value real(Real const & val) {

                // Write out a NaN
                if(val!=val)
                    return Json::Value("NaN");

                // Positive infinity
                else if(val > std::numeric_limits <Real>::max())
                    return Json::Value("Inf");
                
                // Negative infinity
                else if(val < std::numeric_limits <Real>::lowest())
                    return Json::Value("-Inf");

                // A plain old number
                else
                    return Json::Value(val);
            }

            // Write a natural 
            Json::Value natural(Natural const & val);

            // Write a paramter
            template <typename enum_t>
            enum_t param(
                std::function<std::string(enum_t const &)>const& to_string,
                enum_t const & val
            ) {
                return Json::Value(to_string(val));
            }
        }

        // A helper class to help with serialization of vectors into json
        // objects.
        template <typename Real,template <typename> class XX>
        struct Serialization {
            static std::string serialize(
                typename XX <Real>::Vector const & x,
                std::string const & name,
                Natural const & iter
            ) { 
                throw Exception::t(__LOC__
                    + "Optizelle::json::Serialization <>::serialize "
                    + "undefined for the type: "
                    + typeid(XX <Real>).name());
            }
            static typename XX <Real>::Vector deserialize(
                typename XX <Real>::Vector const & x,
                std::string const & x_json
            ) { 
                throw Exception::t(__LOC__
                    + "Optizelle::json::Serialization <>::deserialize "
                    + "undefined for the type: "
                    + typeid(XX <Real>).name());
            }
        };

        // Routines to serialize lists of elements for restarting
        namespace Serialize{
            // Vectors 
            template <typename Real,template <typename> class XX>
            void vectors(
                typename RestartPackage<typename XX<Real>::Vector>::t const& xs,
                std::string const & vs,
                Natural const & iter,
                Json::Value & root
            ) {
                // Create some type shortcuts
                typedef XX <Real> X;
                typedef typename X::Vector X_Vector;
                typedef typename RestartPackage <X_Vector>::t X_Vectors;

                // Create a reader object to parse a json tree
                Json::Reader reader;

                // Loop over all the vectors and serialize things
                for(typename X_Vectors::const_iterator item = xs.cbegin();
                    item!=xs.cend();
                    item++
                ){
                    // Grab the json string of the vector
                    std::string x_json_(Serialization <Real,XX>::serialize(
                        item->second,item->first,iter));
                   
                    // Parse the string
                    Json::Value x_json;
                    reader.parse(x_json_,x_json,true);

                    // Insert the information into the correct place
                    root[vs][item->first]=x_json;
                }
            }

            // Reals 
            template <typename Real>
            void reals(
                typename RestartPackage <Real>::t const & reals,
                std::string const & vs,
                Json::Value & root
            ) {
                // Create some type shortcuts
                typedef typename RestartPackage <Real>::t Reals; 

                // Loop over all the reals and serialize things
                for(typename Reals::const_iterator item = reals.cbegin();
                    item!=reals.cend();
                    item++
                )
                    root[vs][item->first]=write::real(item->second);
            }

            // Naturals 
            void naturals(
                typename RestartPackage <Natural>::t const & nats,
                std::string const & vs,
                Json::Value & root
            );

            // Parameters 
            void parameters(
                typename RestartPackage <std::string>::t const & params,
                std::string const & vs,
                Json::Value & root
            );
        }
        
        // Routines to deserialize lists of elements for restarting
        namespace Deserialize{

            // Vectors
            template <typename Real,template <typename> class XX>
            static void vectors(
                Json::Value const & root,
                std::string const & vs,
                typename XX <Real>::Vector const & x,
                typename RestartPackage<typename XX<Real>::Vector>::t & xs
            ) {
                // Create a writer so that we can tranlate json objects into
                // strings
                Json::StyledWriter writer;

                // Loop over all the names in the root
                for(Json::ValueConstIterator itr=root[vs].begin();
                    itr!=root[vs].end();
                    itr++
                ){
                    // Grab the vector
                    std::string name(itr.key().asString());
                    xs.emplace_back(name,std::move(
                        Serialization <Real,XX>::deserialize(
                            x,writer.write(root[vs][name]))));
                }
            }
            
            // Reals 
            template <typename Real>
            void reals(
                Json::Value const & root,
                std::string const & vs,
                typename RestartPackage <Real>::t & reals
            ) {
                // Loop over all the names in the root
                for(Json::ValueConstIterator itr=root[vs].begin();
                    itr!=root[vs].end();
                    itr++
                ){
                    // Grab the real 
                    std::string name(itr.key().asString());
                    reals.emplace_back(name,std::move(
                        read::real <Real> (root[vs][name],name)));
                }
            }
            
            // Naturals 
            void naturals(
                Json::Value const & root,
                std::string const & vs,
                typename RestartPackage <Natural>::t & nats
            );
            
            // Parameters 
            void parameters(
                Json::Value const & root,
                std::string const & vs,
                typename RestartPackage <std::string>::t & params 
            );
        }

        template <typename Real,template <typename> class XX> 
        struct Unconstrained {
            // Create some type shortcuts
            typedef typename Optizelle::Unconstrained <Real,XX>
                ::X_Vector X_Vector;

            typedef typename Optizelle::Unconstrained <Real,XX>::Restart
                ::X_Vectors X_Vectors; 
            typedef typename Optizelle::Unconstrained <Real,XX>::Restart
                ::Reals Reals;
            typedef typename Optizelle::Unconstrained <Real,XX>::Restart
                ::Naturals Naturals;
            typedef typename Optizelle::Unconstrained <Real,XX>::Restart
                ::Params Params; 

            // Read parameters from file
            static void read_(
                std::string const & fname,
                typename Optizelle::Unconstrained <Real,XX>::State::t & state
            ) {
                // Read in the input file
                Json::Value root=parse(fname);

                // Read in the parameters
                state.eps_grad=read::real <Real> (
                    root["Optizelle"].get("eps_grad",state.eps_grad),
                    "eps_grad");
                state.eps_dx=read::real <Real> (
                    root["Optizelle"].get("eps_dx",state.eps_dx),
                    "eps_dx");
                state.stored_history=read::natural(
                    root["Optizelle"].get(
                        "stored_history",
                        Json::Value::UInt64(state.stored_history)),
                    "stored_history");
                state.iter_max=read::natural(
                    root["Optizelle"].get(
                        "iter_max",
                        Json::Value::UInt64(state.iter_max)),
                    "iter_max");
                state.glob_iter_max=read::natural(
                    root["Optizelle"].get(
                        "glob_iter_max",
                        Json::Value::UInt64(state.glob_iter_max)),
                    "glob_iter_max");
                state.trunc_iter_max=read::natural(
                    root["Optizelle"].get(
                        "trunc_iter_max",
                        Json::Value::UInt64(state.trunc_iter_max)),
                    "trunc_iter_max");
                state.trunc_orthog_storage_max=read::natural(
                    root["Optizelle"].get(
                        "trunc_orthog_storage_max",
                        Json::Value::UInt64(state.trunc_orthog_storage_max)),
                    "trunc_orthog_storage_max");
                state.trunc_orthog_iter_max=read::natural(
                    root["Optizelle"].get(
                        "trunc_orthog_iter_max",
                        Json::Value::UInt64(state.trunc_orthog_iter_max)),
                    "trunc_orthog_iter_max");
                state.eps_trunc=read::real <Real> (
                    root["Optizelle"].get("eps_trunc",state.eps_trunc),
                    "eps_trunc");
                state.algorithm_class=read::param <AlgorithmClass::t> (
                    root["Optizelle"].get("algorithm_class",
                        AlgorithmClass::to_string(state.algorithm_class)),
                    AlgorithmClass::is_valid,
                    AlgorithmClass::from_string,
                    "algorithm_class");
                state.PH_type=read::param <Operators::t> (
                    root["Optizelle"].get("PH_type",
                        Operators::to_string(state.PH_type)),
                    Operators::is_valid,
                    Operators::from_string,
                    "PH_type");
                state.H_type=read::param <Operators::t> (
                    root["Optizelle"].get("H_type",
                        Operators::to_string(state.H_type)),
                    Operators::is_valid,
                    Operators::from_string,
                    "H_type");
                state.msg_level=read::natural(
                    root["Optizelle"].get(
                        "msg_level",
                        Json::Value::UInt64(state.msg_level)),
                    "msg_level");
                state.safeguard_failed_max=read::natural(
                    root["Optizelle"].get(
                        "safeguard_failed_max",
                        Json::Value::UInt64(state.safeguard_failed_max)),
                    "safeguard_failed_max");
                state.delta=read::real <Real> (
                    root["Optizelle"].get("delta",state.delta),
                    "delta");
                state.eta1=read::real <Real> (
                    root["Optizelle"].get("eta1",state.eta1),
                    "eta1");
                state.eta2=read::real <Real> (
                    root["Optizelle"].get("eta2",state.eta2),
                    "eta2");
                state.alpha0=read::real <Real> (
                    root["Optizelle"].get("alpha0",state.alpha0),
                    "alpha0");
                state.c1=read::real <Real> (
                    root["Optizelle"].get("c1",state.c1),
                    "c1");
                state.ls_iter_max=read::natural(
                    root["Optizelle"].get(
                        "ls_iter_max",
                        Json::Value::UInt64(state.ls_iter_max)),
                    "ls_iter_max");
                state.eps_ls=read::real <Real> (
                    root["Optizelle"].get("eps_ls",state.eps_ls),
                    "eps_ls");
                state.dir=read::param <LineSearchDirection::t> (
                    root["Optizelle"].get("dir",
                        LineSearchDirection::to_string(state.dir)),
                    LineSearchDirection::is_valid,
                    LineSearchDirection::from_string,
                    "dir");
                state.kind=read::param <LineSearchKind::t> (
                    root["Optizelle"].get("kind",
                        LineSearchKind::to_string(state.kind)),
                    LineSearchKind::is_valid,
                    LineSearchKind::from_string,
                    "kind");
                state.f_diag=read::param <FunctionDiagnostics::t> (
                    root["Optizelle"].get("f_diag",
                        FunctionDiagnostics::to_string(state.f_diag)),
                    FunctionDiagnostics::is_valid,
                    FunctionDiagnostics::from_string,
                    "f_diag");
                state.L_diag=read::param <FunctionDiagnostics::t> (
                    root["Optizelle"].get("L_diag",
                        FunctionDiagnostics::to_string(state.L_diag)),
                    FunctionDiagnostics::is_valid,
                    FunctionDiagnostics::from_string,
                    "L_diag");
                state.x_diag=read::param <VectorSpaceDiagnostics::t> (
                    root["Optizelle"].get("x_diag",
                        VectorSpaceDiagnostics::to_string(state.x_diag)),
                    VectorSpaceDiagnostics::is_valid,
                    VectorSpaceDiagnostics::from_string,
                    "x_diag");
                state.dscheme=read::param <DiagnosticScheme::t> (
                    root["Optizelle"].get("dscheme",
                        DiagnosticScheme::to_string(state.dscheme)),
                    DiagnosticScheme::is_valid,
                    DiagnosticScheme::from_string,
                    "dscheme");
                state.eps_kind=read::param
                    <ToleranceKind::t> (
                    root["Optizelle"].get("eps_kind",
                        ToleranceKind::to_string(
                            state.eps_kind)),
                    ToleranceKind::is_valid,
                    ToleranceKind::from_string,
                    "eps_kind");
            }
            static void read(
                std::string const & fname,
                typename Optizelle::Unconstrained <Real,XX>::State::t& state
            ) {
                Unconstrained <Real,XX>::read_(fname,state);
                Optizelle::Unconstrained<Real,XX>::State::check(state);
            }

            // Convert parameters to a string 
            static std::string to_string_(
                typename Optizelle::Unconstrained <Real,XX>::State::t & state
            ) {
                // Create a new root for writing
                Json::Value root;

                // Write the optimization parameters
                root["Optizelle"]["eps_grad"]=write::real(state.eps_grad);
                root["Optizelle"]["eps_dx"]=write::real(state.eps_dx);
                root["Optizelle"]["stored_history"]=write::natural(
                    state.stored_history);
                root["Optizelle"]["iter_max"]=write::natural(state.iter_max);
                root["Optizelle"]["glob_iter_max"]=write::natural(
                    state.glob_iter_max);
                root["Optizelle"]["trunc_iter_max"]=write::natural(
                    state.trunc_iter_max);
                root["Optizelle"]["trunc_orthog_storage_max"]=write::natural(
                    state.trunc_orthog_storage_max);
                root["Optizelle"]["trunc_orthog_iter_max"]=write::natural(
                    state.trunc_orthog_iter_max);
                root["Optizelle"]["eps_trunc"]=write::real(state.eps_trunc);
                root["Optizelle"]["algorithm_class"]=write_param(
                    AlgorithmClass::to_string,state.algorithm_class);
                root["Optizelle"]["PH_type"]=write_param(
                    Operators::to_string,state.PH_type);
                root["Optizelle"]["H_type"]=write_param(
                    Operators::to_string,state.H_type);
                root["Optizelle"]["msg_level"]=write::natural(state.msg_level);
                root["Optizelle"]["safeguard_failed_max"]=write::natural(
                    state.safeguard_failed_max);
                root["Optizelle"]["delta"]=write::real(state.delta);
                root["Optizelle"]["eta1"]=write::real(state.eta1);
                root["Optizelle"]["eta2"]=write::real(state.eta2);
                root["Optizelle"]["alpha0"]=write::real(state.alpha0);
                root["Optizelle"]["c1"]=write::real(state.c1);
                root["Optizelle"]["ls_iter_max"]=write::natural(
                    state.ls_iter_max);
                root["Optizelle"]["eps_ls"]=write::real(state.eps_ls);
                root["Optizelle"]["dir"]=write_param(
                    LineSearchDirection::to_string,state.dir);
                root["Optizelle"]["kind"]=write_param(
                    LineSearchKind::to_string,state.kind);
                root["Optizelle"]["f_diag"]=write_param(
                    FunctionDiagnostics::to_string,state.f_diag);
                root["Optizelle"]["L_diag"]=write_param(
                    FunctionDiagnostics::to_string,state.L_diag);
                root["Optizelle"]["x_diag"]=write_param(
                    VectorSpaceDiagnostics::to_string,state.x_diag);
                root["Optizelle"]["dscheme"]=write_param(
                    DiagnosticScheme::to_string,state.dscheme);
                root["Optizelle"]["eps_kind"]=write_param(
                    ToleranceKind::to_string,state.eps_kind);

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
                std::string const & fname,
                typename Optizelle::Unconstrained <Real,XX>::State::t & state
            ) {
                // Grab the iteration number
                Natural iter = state.iter;

                // Do a release 
                X_Vectors xs;
                Reals reals;
                Naturals nats;
                Params params;
                Optizelle::Unconstrained <Real,XX>::Restart::release(
                    state,xs,reals,nats,params);

                // Serialize everything
                Json::Value root;
                Serialize::vectors <Real,XX>(xs,"X_Vectors",iter,root);
                Serialize::reals <Real> (reals,"Reals",root);
                Serialize::naturals(nats,"Naturals",root);
                Serialize::parameters(params,"Parameters",root);
                
                // Write everything to file 
                write_to_file(fname,root);

                // Recapture the state
                Optizelle::Unconstrained <Real,XX>::Restart::capture(
                    state,xs,reals,nats,params);
            }

            // Read all the parameters from file
            static void read_restart(
                std::string const & fname,
                X_Vector const & x,
                typename Optizelle::Unconstrained <Real,XX>::State::t & state
            ) {
                // Read in the input file
                Json::Value root=parse(fname);

                // Extract everything from the parsed json file 
                X_Vectors xs;
                Reals reals;
                Naturals nats;
                Params params;
                Deserialize::vectors <Real,XX>(root,"X_Vectors",x,xs);
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
            typedef typename Optizelle::EqualityConstrained<Real,XX,YY>
                ::X_Vector X_Vector;
            typedef typename Optizelle::EqualityConstrained<Real,XX,YY>
                ::Y_Vector Y_Vector;

            typedef typename Optizelle::EqualityConstrained<Real,XX,YY>::Restart
                ::X_Vectors X_Vectors; 
            typedef typename Optizelle::EqualityConstrained<Real,XX,YY>::Restart
                ::Y_Vectors Y_Vectors; 
            typedef typename Optizelle::EqualityConstrained<Real,XX,YY>::Restart
                ::Reals Reals;
            typedef typename Optizelle::EqualityConstrained<Real,XX,YY>::Restart
                ::Naturals Naturals;
            typedef typename Optizelle::EqualityConstrained<Real,XX,YY>::Restart
                ::Params Params; 

            // Read parameters from file
            static void read_(
                std::string const & fname,
                typename Optizelle::EqualityConstrained <Real,XX,YY>::State::t &
                    state
            ) {
                // Read in the input file
                Json::Value root=parse(fname);

                // Read in the parameters
                state.zeta=read::real <Real> (
                    root["Optizelle"].get("zeta",state.zeta),
                    "zeta");
                state.eta0=read::real <Real> (
                    root["Optizelle"].get("eta0",state.eta0),
                    "eta0");
                state.rho=read::real <Real> (
                    root["Optizelle"].get("rho",state.rho),
                    "rho");
                state.rho_bar=read::real <Real> (
                    root["Optizelle"].get("rho_bar",state.rho_bar),
                    "rho_bar");
                state.eps_constr=read::real <Real> (
                    root["Optizelle"].get("eps_constr",state.eps_constr),
                    "eps_constr");
                state.xi_all(read::real <Real> (
                    root["Optizelle"].get("xi_all",state.xi_qn),
                    "xi_all"));
                state.eps_constr=read::real <Real> (
                    root["Optizelle"].get("eps_constr",state.eps_constr),
                    "eps_constr");
                state.xi_qn=read::real <Real> (
                    root["Optizelle"].get("xi_qn",state.xi_qn),
                    "xi_qn");
                state.xi_pg=read::real <Real> (
                    root["Optizelle"].get("xi_pg",state.xi_pg),
                    "xi_pg");
                state.xi_proj=read::real <Real> (
                    root["Optizelle"].get("xi_proj",state.xi_proj),
                    "xi_proj");
                state.xi_tang=read::real <Real> (
                    root["Optizelle"].get("xi_tang",state.xi_tang),
                    "xi_tang");
                state.xi_lmh=read::real <Real> (
                    root["Optizelle"].get("xi_lmh",state.xi_lmh),
                    "xi_lmh");
                state.xi_lmg=read::real <Real> (
                    root["Optizelle"].get("xi_lmg",state.xi_lmg),
                    "xi_lmg");
                state.xi_4=read::real <Real> (
                    root["Optizelle"].get("xi_4",state.xi_4),
                    "xi_4");
                state.augsys_iter_max=read::natural(
                    root["Optizelle"].get(
                        "augsys_iter_max",
                        Json::Value::UInt64(state.augsys_iter_max)),
                    "augsys_iter_max");
                state.augsys_rst_freq=read::natural(
                    root["Optizelle"].get(
                        "augsys_rst_freq",
                        Json::Value::UInt64(state.augsys_rst_freq)),
                    "augsys_rst_freq");
                state.PSchur_left_type=read::param <Operators::t> (
                    root["Optizelle"].get("PSchur_left_type",
                        Operators::to_string(state.PSchur_left_type)),
                    Operators::is_valid,
                    Operators::from_string,
                    "PSchur_left_type");
                state.PSchur_right_type=read::param <Operators::t> (
                    root["Optizelle"].get("PSchur_right_type",
                        Operators::to_string(state.PSchur_right_type)),
                    Operators::is_valid,
                    Operators::from_string,
                    "PSchur_right_type");
                state.g_diag=read::param <FunctionDiagnostics::t> (
                    root["Optizelle"].get("g_diag",
                        FunctionDiagnostics::to_string(state.g_diag)),
                    FunctionDiagnostics::is_valid,
                    FunctionDiagnostics::from_string,
                    "g_diag");
                state.y_diag=read::param <VectorSpaceDiagnostics::t> (
                    root["Optizelle"].get("y_diag",
                        VectorSpaceDiagnostics::to_string(state.y_diag)),
                    VectorSpaceDiagnostics::is_valid,
                    VectorSpaceDiagnostics::from_string,
                    "y_diag");
            }
            static void read(
                std::string const & fname,
                typename Optizelle::EqualityConstrained <Real,XX,YY>::State::t &
                    state
            ) {
                Unconstrained <Real,XX>::read_(fname,state);
                EqualityConstrained <Real,XX,YY>::read_(fname,state);
                Optizelle::EqualityConstrained<Real,XX,YY>::State::check(state);
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
                root["Optizelle"]["zeta"]=write::real(state.zeta);
                root["Optizelle"]["eta0"]=write::real(state.eta0);
                root["Optizelle"]["rho"]=write::real(state.rho);
                root["Optizelle"]["rho_bar"]=write::real(state.rho_bar);
                root["Optizelle"]["eps_constr"]=write::real(state.eps_constr);
                root["Optizelle"]["xi_qn"]=write::real(state.xi_qn);
                root["Optizelle"]["xi_pg"]=write::real(state.xi_pg);
                root["Optizelle"]["xi_proj"]=write::real(state.xi_proj);
                root["Optizelle"]["xi_tang"]=write::real(state.xi_tang);
                root["Optizelle"]["xi_lmh"]=write::real(state.xi_lmh);
                root["Optizelle"]["xi_lmg"]=write::real(state.xi_lmg);
                root["Optizelle"]["xi_4"]=write::real(state.xi_4);
                root["Optizelle"]["augsys_iter_max"]=write::natural(
                    state.augsys_iter_max);
                root["Optizelle"]["augsys_rst_freq"]=write::natural(
                    state.augsys_rst_freq);
                root["Optizelle"]["PSchur_left_type"]=write_param(
                    Operators::to_string,state.PSchur_left_type);
                root["Optizelle"]["PSchur_right_type"]=write_param(
                    Operators::to_string,state.PSchur_right_type);
                root["Optizelle"]["g_diag"]=write_param(
                    FunctionDiagnostics::to_string,state.g_diag);
                root["Optizelle"]["y_diag"]=write_param(
                    VectorSpaceDiagnostics::to_string,state.y_diag);

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
                std::string const & fname,
                typename Optizelle::EqualityConstrained <Real,XX,YY>::State::t &
                    state
            ) {
                // Grab the iteration number
                Natural iter = state.iter;

                // Do a release 
                X_Vectors xs;
                Y_Vectors ys;
                Reals reals;
                Naturals nats;
                Params params;
                Optizelle::EqualityConstrained <Real,XX,YY>::Restart::release(
                    state,xs,ys,reals,nats,params);

                // Serialize everything
                Json::Value root;
                Serialize::vectors <Real,XX>(xs,"X_Vectors",iter,root);
                Serialize::vectors <Real,YY>(ys,"Y_Vectors",iter,root);
                Serialize::reals <Real> (reals,"Reals",root);
                Serialize::naturals(nats,"Naturals",root);
                Serialize::parameters(params,"Parameters",root);
                
                // Write everything to file 
                write_to_file(fname,root);

                // Recapture the state
                Optizelle::EqualityConstrained<Real,XX,YY>::Restart::capture(
                    state,xs,ys,reals,nats,params);
            }

            // Read all the parameters from file
            static void read_restart(
                std::string const & fname,
                X_Vector const & x,
                Y_Vector const & y,
                typename Optizelle::EqualityConstrained <Real,XX,YY>::State::t &
                    state
            ) {
                // Read in the input file
                Json::Value root=parse(fname);

                // Extract everything from the parsed json file 
                X_Vectors xs;
                Y_Vectors ys;
                Reals reals;
                Naturals nats;
                Params params;
                Deserialize::vectors <Real,XX>(root,"X_Vectors",x,xs);
                Deserialize::vectors <Real,YY>(root,"Y_Vectors",y,ys);
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
            typedef typename Optizelle::InequalityConstrained<Real,XX,ZZ>
                ::X_Vector X_Vector;
            typedef typename Optizelle::InequalityConstrained<Real,XX,ZZ>
                ::Z_Vector Z_Vector;

            typedef typename Optizelle::InequalityConstrained<Real,XX,ZZ>
                ::Restart::X_Vectors X_Vectors; 
            typedef typename Optizelle::InequalityConstrained<Real,XX,ZZ>
                ::Restart::Z_Vectors Z_Vectors; 
            typedef typename Optizelle::InequalityConstrained<Real,XX,ZZ>
                ::Restart::Reals Reals;
            typedef typename Optizelle::InequalityConstrained<Real,XX,ZZ>
                ::Restart::Naturals Naturals;
            typedef typename Optizelle::InequalityConstrained<Real,XX,ZZ>
                ::Restart::Params Params; 

            // Read parameters from file
            static void read_(
                std::string const & fname,
                typename Optizelle::InequalityConstrained<Real,XX,ZZ>::State::t&
                    state
            ) {
                // Read in the input file
                Json::Value root=parse(fname);

                // Read in the parameters
                state.eps_mu=read::real <Real> (
                    root["Optizelle"].get("eps_mu",state.eps_mu),
                    "eps_mu");
                state.mu=read::real <Real> (
                    root["Optizelle"].get("mu",state.mu),
                    "mu");
                state.sigma=read::real <Real> (
                    root["Optizelle"].get("sigma",state.sigma),
                    "sigma");
                state.gamma=read::real <Real> (
                    root["Optizelle"].get("gamma",state.gamma),
                    "gamma");
                state.h_diag=read::param <FunctionDiagnostics::t> (
                    root["Optizelle"].get("h_diag",
                        FunctionDiagnostics::to_string(state.h_diag)),
                    FunctionDiagnostics::is_valid,
                    FunctionDiagnostics::from_string,
                    "h_diag");
                state.z_diag=read::param <VectorSpaceDiagnostics::t> (
                    root["Optizelle"].get("z_diag",
                        VectorSpaceDiagnostics::to_string(state.z_diag)),
                    VectorSpaceDiagnostics::is_valid,
                    VectorSpaceDiagnostics::from_string,
                    "z_diag");
            }
            static void read(
                std::string const & fname,
                typename Optizelle::InequalityConstrained<Real,XX,ZZ>::State::t&
                    state
            ) {
                Unconstrained <Real,XX>::read_(fname,state);
                InequalityConstrained <Real,XX,ZZ>::read_(fname,state);
                Optizelle::InequalityConstrained<Real,XX,ZZ>::State::check(
                    state);
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
                root["Optizelle"]["eps_mu"]=write::real(state.eps_mu);
                root["Optizelle"]["mu"]=write::real(state.mu);
                root["Optizelle"]["sigma"]=write::real(state.sigma);
                root["Optizelle"]["gamma"]=write::real(state.gamma);
                root["Optizelle"]["h_diag"]=write_param(
                    FunctionDiagnostics::to_string,state.h_diag);
                root["Optizelle"]["z_diag"]=write_param(
                    VectorSpaceDiagnostics::to_string,state.z_diag);

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
                std::string const & fname,
                typename Optizelle::InequalityConstrained<Real,XX,ZZ>::State::t&
                    state
            ) {
                // Grab the iteration number
                Natural iter = state.iter;

                // Do a release 
                X_Vectors xs;
                Z_Vectors zs;
                Reals reals;
                Naturals nats;
                Params params;
                Optizelle::InequalityConstrained <Real,XX,ZZ>::Restart::release(
                    state,xs,zs,reals,nats,params);

                // Serialize everything
                Json::Value root;
                Serialize::vectors <Real,XX>(xs,"X_Vectors",iter,root);
                Serialize::vectors <Real,ZZ>(zs,"Z_Vectors",iter,root);
                Serialize::reals <Real> (reals,"Reals",root);
                Serialize::naturals(nats,"Naturals",root);
                Serialize::parameters(params,"Parameters",root);
                
                // Write everything to file 
                write_to_file(fname,root);

                // Recapture the state
                Optizelle::InequalityConstrained<Real,XX,ZZ>::Restart::capture(
                    state,xs,zs,reals,nats,params);
            }

            // Read all the parameters from file
            static void read_restart(
                std::string const & fname,
                X_Vector const & x,
                Z_Vector const & z,
                typename Optizelle::InequalityConstrained<Real,XX,ZZ>::State::t&
                    state
            ) {
                // Read in the input file
                Json::Value root=parse(fname);

                // Extract everything from the parsed json file 
                X_Vectors xs;
                Z_Vectors zs;
                Reals reals;
                Naturals nats;
                Params params;
                Deserialize::vectors <Real,XX>(root,"X_Vectors",x,xs);
                Deserialize::vectors <Real,ZZ>(root,"Z_Vectors",z,zs);
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
            typedef typename Optizelle::Constrained<Real,XX,YY,ZZ>
                ::X_Vector X_Vector;
            typedef typename Optizelle::Constrained<Real,XX,YY,ZZ>
                ::Y_Vector Y_Vector;
            typedef typename Optizelle::Constrained<Real,XX,YY,ZZ>
                ::Z_Vector Z_Vector;

            typedef typename Optizelle::Constrained<Real,XX,YY,ZZ>::Restart
                ::X_Vectors X_Vectors; 
            typedef typename Optizelle::Constrained<Real,XX,YY,ZZ>::Restart
                ::Y_Vectors Y_Vectors; 
            typedef typename Optizelle::Constrained<Real,XX,YY,ZZ>::Restart
                ::Z_Vectors Z_Vectors; 
            typedef typename Optizelle::Constrained<Real,XX,YY,ZZ>::Restart
                ::Reals Reals;
            typedef typename Optizelle::Constrained<Real,XX,YY,ZZ>::Restart
                ::Naturals Naturals;
            typedef typename Optizelle::Constrained<Real,XX,YY,ZZ>::Restart
                ::Params Params; 

            // Read parameters from file
            static void read(
                std::string const & fname,
                typename Optizelle::Constrained <Real,XX,YY,ZZ>::State::t&
                    state
            ) {
                Unconstrained <Real,XX>::read_(fname,state);
                EqualityConstrained <Real,XX,YY>::read_(fname,state);
                InequalityConstrained <Real,XX,ZZ>::read_(fname,state);
                Optizelle::Constrained<Real,XX,YY,ZZ>::State::check(state);
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
                std::string const & fname,
                typename Optizelle::Constrained <Real,XX,YY,ZZ>::State::t &
                    state
            ) {
                // Grab the iteration number
                Natural iter = state.iter;

                // Do a release 
                X_Vectors xs;
                Y_Vectors ys;
                Z_Vectors zs;
                Reals reals;
                Naturals nats;
                Params params;
                Optizelle::Constrained <Real,XX,YY,ZZ>::Restart::release(
                    state,xs,ys,zs,reals,nats,params);

                // Serialize everything
                Json::Value root;
                Serialize::vectors <Real,XX>(xs,"X_Vectors",iter,root);
                Serialize::vectors <Real,YY>(ys,"Y_Vectors",iter,root);
                Serialize::vectors <Real,ZZ>(zs,"Z_Vectors",iter,root);
                Serialize::reals <Real> (reals,"Reals",root);
                Serialize::naturals(nats,"Naturals",root);
                Serialize::parameters(params,"Parameters",root);
                
                // Write everything to file 
                write_to_file(fname,root);

                // Recapture the state
                Optizelle::Constrained<Real,XX,YY,ZZ>::Restart::capture(
                    state,xs,ys,zs,reals,nats,params);
            }
            
            // Read all the parameters from file
            static void read_restart(
                std::string const & fname,
                X_Vector const & x,
                Y_Vector const & y,
                Z_Vector const & z,
                typename Optizelle::Constrained <Real,XX,YY,ZZ>::State::t& state
            ) {
                // Read in the input file
                Json::Value root=parse(fname);

                // Extract everything from the parsed json file 
                X_Vectors xs;
                Y_Vectors ys;
                Z_Vectors zs;
                Reals reals;
                Naturals nats;
                Params params;
                Deserialize::vectors <Real,XX>(root,"X_Vectors",x,xs);
                Deserialize::vectors <Real,YY>(root,"Y_Vectors",y,ys);
                Deserialize::vectors <Real,ZZ>(root,"Z_Vectors",z,zs);
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
