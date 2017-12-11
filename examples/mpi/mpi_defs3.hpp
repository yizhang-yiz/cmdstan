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

  // run a map_rect_serial which uses just like the MPI code a
  // precomputed gradient approach
  /*
  template <>
  Eigen::Matrix<var, Eigen::Dynamic, 1>
  map_rect_serial(const Eigen::Matrix<var, Eigen::Dynamic, 1>& shared_params,
                  const std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1> >& job_params,
                  const std::vector<std::vector<double> >& X_r,
                  const std::vector<std::vector<int> >& X_i, std::ostream* pstream__) {
    typedef map_rect_reduce<mpi_function_functor__, var, var> reducer_t;
    const size_type num_shared_params = shared_params.rows();
    const size_type num_job_specific_params = job_params[0].rows();
    const size_type num_jobs = job_params.size();
    const size_type num_params = num_shared_params  + num_job_specific_params;
    std::vector<int> world_f_out(num_jobs, 1);
    // WARNING: HARD CODE THAT WE HAVE 1 OUTPUT PER JOB!!!
    Eigen::Matrix<var, Eigen::Dynamic, 1> out(num_jobs);

    for(size_type i = 0, ij=0; i < num_jobs; i++) {
      matrix_d world_result = reducer_t::apply(value_of(shared_params), value_of(job_params[i]), X_r[i], X_i[i]);
      const std::size_t offset_job_params = stan::is_constant_struct<var>::value ? 1 : 1+num_shared_params ;

      operands_and_partials<Eigen::Matrix<var, Eigen::Dynamic, 1>,
                            Eigen::Matrix<var, Eigen::Dynamic, 1> >
        ops_partials(shared_params, job_params[i]);

      for(std::size_t j=0; j != world_f_out[i]; ++j, ++ij) {
        // check if the outputs flags a failure
        if(unlikely(world_result(0,j) == std::numeric_limits<double>::max())) {
          //std::cout << "THROWING ON RANK " << rank_ << std::endl;
          throw std::runtime_error("MPI error.");
        }

        if (!stan::is_constant_struct<var>::value) {
          ops_partials.edge1_.partials_ = world_result.block(1,j,num_shared_params,1);
        }
              
        if (!stan::is_constant_struct<var>::value) {
          ops_partials.edge2_.partials_ = world_result.block(offset_job_params,j,num_job_specific_params,1);
        }
            
        out(ij) = ops_partials.build(world_result(0,j));
      }
    }
    return(out);
  }
  */
}

STAN_REGISTER_MAP_RECT(0, mpi_model_namespace::mpi_function_functor__)

