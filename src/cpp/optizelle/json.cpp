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

#include "optizelle/json.h"

namespace Optizelle {
    namespace json {
        // Parses a JSON file and returns the root
        Json::Value parse(
            Optizelle::Messaging const & msg,
            std::string const & fname
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
       
        // Writes a JSON spec to file
        void write_to_file(
            Optizelle::Messaging const & msg,
            std::string const & fname,
            Json::Value const & root
        ) {
            // Create a string with the above output
            Json::StyledWriter writer;
            std::string output = writer.write(root);

            // Open a file for writing
            std::ofstream fout(fname.c_str());
            if(fout.fail())
                msg.error("While writing the restart file, unable to open "
                    "the file: " + fname + ".");

            // Write out the json tree
            fout << output;
            if(fout.fail())
                msg.error("While writing the restart file, unable to write "
                    "the json tree.");

            // Close the file
            fout.close();
        }
        
        // Safely reads from a json tree 
        namespace read {
            // Read a natural 
            Natural natural(
                Optizelle::Messaging const & msg,
                Json::Value const & json,
                std::string const & name
            ) {
                // Set the error message
                std::string const err_msg = "Invalid JSON parameter: "
                    + name + " contains an invalid natural.";

                // As long as we have an unsigned integer, grab it
                if(json.isUInt())
                    return Natural(Json::Value::UInt64(json.asUInt64()));
                
                // If we have an integer, grab it if it's positive
                else if(json.isInt()) {
                    Integer val(json.asInt64());
                    if(val>=0)
                        return Natural(val);
                    else
                        msg.error(err_msg);

                // Anything else is an error
                } else
                    msg.error(err_msg);

                // We should not hit this point
                throw;
            }
            
            // Read a string 
            std::string string(
                Optizelle::Messaging const & msg,
                Json::Value const & json,
                std::string const & name
            ) {
                // Set the error message
                std::string const err_msg = "Invalid JSON parameter: "
                    + name + " is not a valid parameter.";

                // As long as we have an unsigned integer, grab it
                if(json.isString())
                    return json.asString();

                // Anything else is an error
                else
                    msg.error(err_msg);

                // We should not hit this point
                throw;
            }
        }

        // Writes into a json tree 
        namespace write {
            // Write a natural 
            Json::Value natural(Natural const & val) {
                return Json::Value::UInt64(val);
            }
        }

        // Routines to serialize lists of elements for restarting
        namespace Serialize{
       
            // Naturals 
            void naturals(
                typename RestartPackage <Natural>::t const & nats,
                std::string const & vs,
                Json::Value & root
            ) {
                // Create some type shortcuts
                typedef typename RestartPackage <Natural>::t Naturals; 

                // Loop over all the naturals and serialize things
                for(typename Naturals::const_iterator item = nats.cbegin();
                    item!=nats.cend();
                    item++
                )
                    root[vs][item->first]=write::natural(item->second);
            }

            // Parameters 
            void parameters(
                typename RestartPackage <std::string>::t const & params,
                std::string const & vs,
                Json::Value & root
            ) {
                // Create some type shortcuts
                typedef typename RestartPackage <std::string>::t Params; 

                // Loop over all the params and serialize things
                for(typename Params::const_iterator item = params.cbegin();
                    item!=params.cend();
                    item++
                )
                    root[vs][item->first]=item->second;
            }
        }
        
        // Routines to deserialize lists of elements for restarting
        namespace Deserialize{

            // Naturals 
            void naturals(
                Optizelle::Messaging const & msg,
                Json::Value const & root,
                std::string const & vs,
                typename RestartPackage <Natural>::t & nats
            ) {
                // Loop over all the names in the root
                for(Json::ValueConstIterator itr=root[vs].begin();
                    itr!=root[vs].end();
                    itr++
                ){
                    // Grab the natural 
                    std::string name(itr.key().asString());
                    nats.emplace_back(name,std::move(
                        read::natural(msg,root[vs][name],name))); 
                }
            }
            
            // Parameters 
            void parameters(
                Optizelle::Messaging const & msg,
                Json::Value const & root,
                std::string const & vs,
                typename RestartPackage <std::string>::t & params 
            ) {
                // Loop over all the names in the root
                for(Json::ValueConstIterator itr=root[vs].begin();
                    itr!=root[vs].end();
                    itr++
                ){
                    // Grab the parameter 
                    std::string name(itr.key().asString());
                    params.emplace_back(name,std::move(
                        read::string(msg,root[vs][name],name)));
                }
            }
        }
    }
}
