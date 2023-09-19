/** @file
 * This file is left mostly empty so that catch doesn't need 
 * slow (~16s) recompilation each time unit tests are edited
 *
 * @author Tyson Jones
 */

/** Use our modified Catch in custom-main mode (main defined below).
 * catch.hpp was modified to, in distributed mode, output only once.
 */
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include "QuEST.h"

// #ifdef _OPENMP
#include <omp.h>
// #endif

/** The global QuESTEnv instance, to be created and destroyed once in this 
 * main(), so that the MPI environment is correctly created once when running 
 * distributed unit tests 
 */
QuESTEnv QUEST_ENV;

/** Redefinition of QuEST_validation's invalidQuESTInputError function, called when a 
 * user passes an incorrect parameter (e.g. an negative qubit index). This is 
 * redefined here to, in lieu of printing and exiting, throw a C++ exception
 * which can be caught (and hence unit tested for) by Catch2
 */
 extern "C" void invalidQuESTInputError(const char* errMsg, const char* errFunc) {
     throw errMsg;
 }
 
/** Explicit declaration of main to create (destroy) the QuESTEnv before (after)
 * invoking the Catch unit tests 
 */
int main( int argc, char* argv[] ) {
  Catch::Session session; // There must be exactly one instance

  int numThreads = 0; // Some user variable you want to be able to set

  // Build a new parser on top of Catch2's
  using namespace Catch::clara;
  auto cli
    = session.cli()           // Get Catch2's command line parser
    | Opt( numThreads, "numThreads" ) // bind variable to a new option, with a hint string
        ["-j"]["--numThreads"]    // the option names it will respond to
        ("Number of threads");        // description string for the help output

  // Now pass the new composite back to Catch2 so it uses that
  session.cli( cli );

  // Let Catch2 (using Clara) parse the command line
  int returnCode = session.applyCommandLine( argc, argv );
  if( returnCode != 0 ) // Indicates a command line error
      return returnCode;

  if (numThreads != 0) {
    omp_set_num_threads(numThreads);
  }

  QUEST_ENV = createQuESTEnv();

  int result = session.run();
  destroyQuESTEnv(QUEST_ENV);
  return result;
}