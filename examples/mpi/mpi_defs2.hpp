#include <stan/math/prim/arr/functor/mpi_command.hpp>
#include <stan/math/prim/arr/functor/mpi_cluster.hpp>
#include <stan/math/prim/arr/functor/map_rect_mpi.hpp>
#include <stan/math/rev/arr/functor/map_rect_mpi.hpp>

struct mpi_apply {
  template<typename T>
  static
  std::vector<T> apply(std::vector<T> theta, std::vector<double> x_r, std::vector<int> x_i) {
    return oral_2cmt_mpi2_model_namespace::mpi_function(theta, x_r, x_i, 0);
  }
};

namespace oral_2cmt_mpi2_model_namespace {


  template <typename T0__, typename T1__>
  std::vector<typename boost::math::tools::promote_args<T0__, T1__>::type>
  run_mpi_function(const std::vector<std::vector<T0__> >& Theta,
                   const std::vector<std::vector<T1__> >& X_r,
                   const std::vector<std::vector<int> >& X_i, std::ostream* pstream__) {
    return stan::math::map_rect_mpi<mpi_apply>(Theta, X_r, X_i, 0);
  }

}


BOOST_CLASS_EXPORT(stan::math::distributed_apply<stan::math::internal::distributed_map_rect_data>)
BOOST_CLASS_TRACKING(stan::math::distributed_apply<stan::math::internal::distributed_map_rect_data>,track_never)
BOOST_SERIALIZATION_FACTORY_0(stan::math::distributed_apply<stan::math::internal::distributed_map_rect_data>)

// argh: class names must be shorter than 128 characters. Otherwise
// boost::serialization throws an exception
BOOST_CLASS_EXPORT(stan::math::distributed_apply<stan::math::internal::distributed_map_rect<mpi_apply> >)
BOOST_CLASS_TRACKING(stan::math::distributed_apply<stan::math::internal::distributed_map_rect<mpi_apply> >,track_never)
BOOST_SERIALIZATION_FACTORY_0(stan::math::distributed_apply<stan::math::internal::distributed_map_rect<mpi_apply> >)


/* does not help for an attempt to get static linking working
#include <boost/mpi/packed_iarchive.hpp>
#include <boost/mpi/packed_oarchive.hpp>

template void stan::math::stop_worker::serialize<boost::mpi::packed_iarchive>(boost::mpi::packed_iarchive & ar, const unsigned int version);
template void stan::math::stop_worker::serialize<boost::mpi::packed_oarchive>(boost::mpi::packed_oarchive & ar, const unsigned int version);

template void stan::math::distributed_apply<stan::math::internal::distributed_map_rect<mpi_apply> >::serialize<boost::mpi::packed_iarchive>(boost::mpi::packed_iarchive & ar, const unsigned int version);
template void stan::math::distributed_apply<stan::math::internal::distributed_map_rect<mpi_apply> >::serialize<boost::mpi::packed_oarchive>(boost::mpi::packed_oarchive & ar, const unsigned int version);

template void stan::math::distributed_apply<stan::math::internal::distributed_map_rect<mpi_apply> >::serialize<boost::archive::detail::basic_iarchive>(boost::archive::detail::basic_iarchive & ar, const unsigned int version);
template void stan::math::distributed_apply<stan::math::internal::distributed_map_rect<mpi_apply> >::serialize<boost::archive::detail::basic_oarchive>(boost::archive::detail::basic_oarchive & ar, const unsigned int version);

template void stan::math::distributed_apply<stan::math::internal::distributed_map_rect_data>::serialize<boost::mpi::packed_iarchive>(boost::mpi::packed_iarchive & ar, const unsigned int version);
template void stan::math::distributed_apply<stan::math::internal::distributed_map_rect_data>::serialize<boost::mpi::packed_oarchive>(boost::mpi::packed_oarchive & ar, const unsigned int version);
*/
