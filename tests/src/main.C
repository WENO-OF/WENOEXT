/*---------------------------------------------------------------------------*\
       ██╗    ██╗███████╗███╗   ██╗ ██████╗     ███████╗██╗  ██╗████████╗
       ██║    ██║██╔════╝████╗  ██║██╔═══██╗    ██╔════╝╚██╗██╔╝╚══██╔══╝
       ██║ █╗ ██║█████╗  ██╔██╗ ██║██║   ██║    █████╗   ╚███╔╝    ██║   
       ██║███╗██║██╔══╝  ██║╚██╗██║██║   ██║    ██╔══╝   ██╔██╗    ██║   
       ╚███╔███╔╝███████╗██║ ╚████║╚██████╔╝    ███████╗██╔╝ ██╗   ██║   
        ╚══╝╚══╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝     ╚══════╝╚═╝  ╚═╝   ╚═╝   
-------------------------------------------------------------------------------                                                                                                                                                         
License
    This file is part of WENO Ext.

    WENO Ext is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    WENO Ext is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with  WENO Ext.  If not, see <http://www.gnu.org/licenses/>.

Description
    Custom main function for Catch2 to provide serial and parallel runs 
    Switch between serial and parallel using the -p or --parallel flag

    e.g:
      # Run in serial
      ./unitTest

      # Run in parallel 
      mpirun -np 8 ./unitTest -p 

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/


#include <catch2/catch_session.hpp>
#include <iostream>
#include "fvCFD.H"
#include "globalFoamArgs.H"

Foam::argList* argsPtr = nullptr;

int main(int argc, char* argv[])
{
    Catch::Session session;          // There must be exactly one instance

    // Add the command line argument to switch between parallel mode
    // Build a new parser on top of Catch2's
    bool parallelRun = false;

    using namespace Catch::Clara;
    auto cli
    = session.cli()                  // Get Catch2's command line parser
    | Opt( parallelRun ) // bind variable to a new option, with a hint string
        ["-p"]["--parallel"]         // the option names it will respond to
        ("parallel run");            // description string for the help output

    // Now pass the new composite back to Catch2 so it uses that
    session.cli( cli );

    // Let Catch2 (using Clara) parse the command line
    int returnCode = session.applyCommandLine( argc, argv );
    if( returnCode != 0 ) // Indicates a command line error
        return returnCode;


    // Create OpenFOAM arguments 
    // Has to be done here for MPI support
    int argcOF = 1;
    if (parallelRun)
        argcOF++;
    
    char **argvOF = static_cast<char**>(malloc(sizeof(char*)));
    char executable[] = {'p','a','r','a','l','l','e','l','T','e','s','t'};
    char flags[] = "-parallel";
    argvOF[0] = executable;
    argvOF[1] = flags;
    setFoamArgs(argcOF, argvOF);

    // Start the session
    const int result = session.run();

    // Clean up MPI 
    clearFoamArgs();

    return result;
}
