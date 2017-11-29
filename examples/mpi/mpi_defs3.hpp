//#include <boost/smart_ptr/shared_ptr.hpp>
#include <stan/math/prim/arr/functor/mpi_command.hpp>
#include <stan/math/prim/arr/functor/mpi_cluster.hpp>
#include <stan/math/prim/mat/functor/map_rect.hpp>
#include <stan/math/prim/mat/functor/map_rect_mpi.hpp>
#include <stan/math/rev/mat/functor/map_rect_reduce.hpp>

struct mpi_call {
    template <typename T0__, typename T1__, typename T2__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic,1>
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic,1>& eta,
                 const Eigen::Matrix<T1__, Eigen::Dynamic,1>& theta,
                 const std::vector<T2__>& x_r,
                 const std::vector<int>& x_i) const {
      return mpi_model_namespace::mpi_function(eta, theta, x_r, x_i, 0);
    }
};

namespace mpi_model_namespace {

  template <typename T0__, typename T1__, typename T2__>
  Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic, 1>
  map_rect(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& eta,
           const std::vector<Eigen::Matrix<T1__, Eigen::Dynamic, 1> >& Theta,
           const std::vector<std::vector<T2__> >& X_r,
           const std::vector<std::vector<int> >& X_i, std::ostream* pstream__) {
    return stan::math::map_rect<mpi_call>(eta, Theta, X_r, X_i, 0);
  }
}

/*
BOOST_CLASS_EXPORT(stan::math::mpi_distributed_apply<stan::math::internal::distributed_map_rect_data>)
BOOST_CLASS_TRACKING(stan::math::mpi_distributed_apply<stan::math::internal::distributed_map_rect_data>,track_never)
BOOST_SERIALIZATION_FACTORY_0(stan::math::mpi_distributed_apply<stan::math::internal::distributed_map_rect_data>)
*/

// argh: class names must be shorter than 128 characters. Otherwise
// boost::serialization throws an exception

typedef stan::math::map_rect_reduce<mpi_call, stan::math::var, stan::math::var> mpi_function_functor_reducer_vv;
typedef stan::math::map_rect_combine<mpi_call, stan::math::var, stan::math::var> mpi_function_functor_combiner_vv;
typedef stan::math::mpi_parallel_call<mpi_function_functor_reducer_vv,mpi_function_functor_combiner_vv> mpi_function_functor_parallel_call_vv;

BOOST_CLASS_EXPORT(stan::math::mpi_distributed_apply<mpi_function_functor_parallel_call_vv>);
BOOST_CLASS_TRACKING(stan::math::mpi_distributed_apply<mpi_function_functor_parallel_call_vv>,track_never);
BOOST_SERIALIZATION_FACTORY_0(stan::math::mpi_distributed_apply<mpi_function_functor_parallel_call_vv>);


typedef stan::math::map_rect_reduce<mpi_call, double, double> mpi_function_functor_reducer_dd;
typedef stan::math::map_rect_combine<mpi_call, double, double> mpi_function_functor_combiner_dd;
typedef stan::math::mpi_parallel_call<mpi_function_functor_reducer_dd,mpi_function_functor_combiner_dd> mpi_function_functor_parallel_call_dd;

BOOST_CLASS_EXPORT(stan::math::mpi_distributed_apply<mpi_function_functor_parallel_call_dd>);
BOOST_CLASS_TRACKING(stan::math::mpi_distributed_apply<mpi_function_functor_parallel_call_dd>,track_never);
BOOST_SERIALIZATION_FACTORY_0(stan::math::mpi_distributed_apply<mpi_function_functor_parallel_call_dd>);

