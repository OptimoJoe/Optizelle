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
        void write(
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

        // Routines to serialize lists of elements for restarting
        namespace Serialize{
       
            // Naturals 
            void naturals(
                std::pair <
                    std::list <std::string>,
                    std::list <Natural>
                > const & nats,
                std::string const & vs,
                Json::Value & root
            ) {
                // Loop over all the naturals and serialize things
                typename std::list <Natural>::const_iterator nat 
                    =nats.second.begin();
                for(typename std::list <std::string>::const_iterator
                        name=nats.first.begin();
                    name!=nats.first.end();
                    name++, nat++
                )
                    root[vs][*name]=Json::Value::UInt64(*nat);
            }

            // Parameters 
            void parameters(
                const std::pair <
                    std::list <std::string>,
                    std::list <std::string>
                >& params,
                std::string const & vs,
                Json::Value & root
            ) {
                // Loop over all the parameters and serialize things
                typename std::list <std::string>::const_iterator param 
                    =params.second.begin();
                for(typename std::list <std::string>::const_iterator
                        name=params.first.begin();
                    name!=params.first.end();
                    name++, param++
                )
                    root[vs][*name]=*param;
            }
        }
        
        // Routines to deserialize lists of elements for restarting
        namespace Deserialize{

            // Naturals 
            void naturals(
                Json::Value const & root,
                std::string const & vs,
                std::pair <
                    std::list <std::string>,
                    std::list <Natural>
                >& nats
            ) {
                // Loop over all the names in the root
                for(Json::ValueIterator itr=root[vs].begin();
                    itr!=root[vs].end();
                    itr++
                ){
                    // Grab the natural 
                    nats.first.emplace_back(itr.key().asString());
                    nats.second.emplace_back(
                        root[vs][nats.first.back()].asUInt64());
                }
            }
            
            // Parameters 
            void parameters(
                Json::Value const & root,
                std::string const & vs,
                std::pair <
                    std::list <std::string>,
                    std::list <std::string>
                >& params
            ) {
                // Loop over all the names in the root
                for(Json::ValueIterator itr=root[vs].begin();
                    itr!=root[vs].end();
                    itr++
                ){
                    // Grab the parameter 
                    params.first.emplace_back(itr.key().asString());
                    params.second.emplace_back(
                        root[vs][params.first.back()].asString());
                }
            }
        }
    }
}
