#include <cmdstan/command.hpp>
#include <stan/services/error_codes.hpp>
#include <boost/exception/diagnostic_information.hpp> 
#include <boost/exception_ptr.hpp>

#include <stan/math/prim/arr/functor/mpi_cluster.hpp>

int main(int argc, const char* argv[]) {
  boost::mpi::environment env;
  int exitcode = 0;
  
  // on non-root processes this makes the workers listen to commands
  // send from the root
  stan::math::mpi_cluster cluster;

  boost::mpi::communicator world;

  const std::size_t rank = world.rank();

  cluster.listen();
  
  //const char** argv_const = argv;

  if (rank == 0) {
    try {
      exitcode = cmdstan::command<stan_model>(argc,argv);
    } catch (const std::exception& e) {
      std::cout << e.what() << std::endl;
      exitcode = stan::services::error_codes::SOFTWARE;
    }
  }

  cluster.stop_listen();

  return(exitcode);
}
