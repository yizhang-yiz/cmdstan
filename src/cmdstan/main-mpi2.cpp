#include <cmdstan/command.hpp>
#include <stan/services/error_codes.hpp>
#include <boost/exception/diagnostic_information.hpp> 
#include <boost/exception_ptr.hpp>

#include <stan/math/prim/arr/functor/mpi_cluster.hpp>

int main(int argc, const char* argv[]) {
  boost::mpi::environment env;
  int exitcode = 0;
  //boost::mpi::communicator world;
  
  // on non-root processes this makes the workers listen to commands
  // send from the root
  stan::math::mpi_cluster cluster;

  //const char** argv_const = argv;

  try {
    exitcode = cmdstan::command<stan_model>(argc,argv);
  } catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
    exitcode = stan::services::error_codes::SOFTWARE;
  }

  return(exitcode);
}
