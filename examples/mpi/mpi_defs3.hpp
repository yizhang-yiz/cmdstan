//#include <boost/smart_ptr/shared_ptr.hpp>
#include <stan/math/prim/arr/functor/mpi_command.hpp>
#include <stan/math/prim/arr/functor/mpi_cluster.hpp>
#include <stan/math/prim/mat/functor/map_rect.hpp>
#include <stan/math/prim/mat/functor/map_rect_mpi.hpp>
#include <stan/math/rev/mat/functor/map_rect_reduce.hpp>
#include <stan/math/prim/mat/functor/mpi_parallel_call.hpp>

namespace mpi_model_namespace {

  template <typename T0__, typename T1__, typename T2__>
  Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic, 1>
  map_rect_mpi(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& eta,
               const std::vector<Eigen::Matrix<T1__, Eigen::Dynamic, 1> >& Theta,
               const std::vector<std::vector<T2__> >& X_r,
               const std::vector<std::vector<int> >& X_i, std::ostream* pstream__) {
    return stan::math::map_rect_mpi<0,mpi_function_functor__>(eta, Theta, X_r, X_i, pstream__);
  }
  
  template <typename T0__, typename T1__, typename T2__>
  Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic, 1>
  map_rect_serial(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& eta,
               const std::vector<Eigen::Matrix<T1__, Eigen::Dynamic, 1> >& Theta,
               const std::vector<std::vector<T2__> >& X_r,
               const std::vector<std::vector<int> >& X_i, std::ostream* pstream__) {
    return stan::math::map_rect_serial<0,mpi_function_functor__>(eta, Theta, X_r, X_i, pstream__);
  }

}

STAN_REGISTER_MAP_RECT(0, mpi_model_namespace::mpi_function_functor__)

