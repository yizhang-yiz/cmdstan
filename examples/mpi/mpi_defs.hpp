#include <boost/mpi.hpp>

namespace mpi_util {
  // can be used to ship static data around; obsolete as MPI jobs
  // start such that the data is distributed with NFS. Left here to
  // show how the boost serialization library can be used to share
  // data with MPI. However, its much faster to copy raw memory
  // buffers around which avoids a lot of expensive abstraction.
  struct function_data {
    std::vector<double> x_r_;
    std::vector<int> x_i_;
    
    function_data()
      : x_r_(), x_i_() {}
  
    function_data(const std::vector<double>& x_r, const std::vector<int>& x_i)
      : x_r_(x_r), x_i_(x_i) {}
    
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & x_r_;
      ar & x_i_;
    }
  };

  template <typename F>
  void jacobian(const F& f,
                const Eigen::VectorXd& theta,
                const std::vector<double>& x_r,
                const std::vector<int>& x_i,
                Eigen::Block<Eigen::MatrixXd> FJ) {
      const size_t M = theta.size();

      stan::math::start_nested();
      try {
        std::vector<stan::math::var> theta_v(theta.data(), theta.data() + M);
        std::vector<stan::math::var> fx_v = f(theta_v, x_r, x_i, 0);

        const std::size_t N = fx_v.size();

        // FJ is filled columns wise with the first entry being the
        // value f and in the same column all partials and so on
      
        for (std::size_t i = 0; i != N; ++i) {
          FJ(0,i) = fx_v[i].val();
          stan::math::set_zero_all_adjoints_nested();
          stan::math::grad(fx_v[i].vi_);
          for (std::size_t k = 0; k != M; ++k)
            FJ(k+1,i) = theta_v[k].adj();
        }
      } catch (const std::exception& e) {
        stan::math::recover_memory_nested();
        throw;
      }
      stan::math::recover_memory_nested();
   }
};

namespace oral_2cmt_mpi_model_namespace {
  
  template <typename T0__, typename T2__>
  std::vector<typename boost::math::tools::promote_args<T0__, T2__>::type>
  run_mpi_function(const std::vector<std::vector<T0__> >& Theta,
                   const std::vector<int >& chunks,
                   const std::vector<std::vector<T2__> >& X_r,
                   const std::vector<std::vector<int> >& X_i,
                   std::ostream* pstream__) {
    typedef Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T2__>::type, Eigen::Dynamic,1> return_t;

    return_t res(10);

    typedef stan::is_var<T0__> T0_var;
    typedef stan::is_var<T2__> T2_var;

    if(T0_var::value) {
      std::cout << "T0 = var" << std::endl;
    } else {
      std::cout << "T0 != var" << std::endl;
    }

    if(T2_var::value) {
      std::cout << "T2 = var" << std::endl;
    } else {
      std::cout << "T2 != var" << std::endl;
    }

    std::cout << "WRONG FUNCTION!!!" << std::endl;

    return(res);
  }

  // TODO: abstract out MPI common code !!!
  template <>
  std::vector<double>
  run_mpi_function(const std::vector<std::vector<double> >& Theta,
                   const std::vector<std::vector<int> >& work_map,
                   const std::vector<std::vector<double> >& X_r,
                   const std::vector<std::vector<int> >& X_i, std::ostream* pstream__) {
    boost::mpi::communicator world;

    const std::size_t M = Theta[0].size();
    const std::size_t J = Theta.size();
    const std::size_t W = world.size();

    std::vector<int> chunks(W, 0);
    std::vector<int> chunks_Theta(W, 0);
    std::vector<int> chunks_result(W, 0);

    std::vector<int > N_node;
    int N_total = 0;
    int N_world = 0;

    for(std::size_t j = 0; j != J; j++) {
      const int worker = work_map[j][0] - 1;
      const int N_job = work_map[j][1];
      
      chunks[worker] += 1;
      chunks_Theta[worker] += M;
      chunks_result[worker] += (1 + M) * N_job;
      
      if(worker == 0) {
        N_total += N_job;
        N_node.push_back(N_job);
      }

      N_world += N_job;
    }
    
    const std::size_t C = chunks[0];

    // scatter parameters
    Eigen::MatrixXd world_Theta(M, J);

    for(std::size_t j = 0; j != J; ++j)
      for(std::size_t m = 0; m != M; ++m)
        world_Theta(m,j) = stan::math::value_of(Theta[j][m]);

    //std::cout << "Shipping parameters to childs:" << std::endl << world_Theta << std::endl;

    Eigen::MatrixXd local_Theta(M, C);
    boost::mpi::scatterv(world, world_Theta.data(), chunks_Theta, local_Theta.data(), 0);
    
    // now work on the chunk for the root
    // calculate the chunk of the root

    // allocate storage for local function_result
    Eigen::MatrixXd local_result(M + 1, N_total);

    // allocate memory on AD stack for final result. Note that the
    // gradients and the function results will land there
    double* final_result
      = ChainableStack::memalloc_.alloc_array<double>( (M + 1) * N_world );

    try {
      std::size_t nc = 0;
      for(std::size_t i = 0; i != C; i++) {
        mpi_util::jacobian(mpi_function_functor__(),
                           local_Theta.col(i),
                           X_r[i], X_i[i],
                           local_result.block(0, nc, M+1, N_node[i]));
        nc += N_node[i];
      }
    } catch(const std::exception& e) {
      // we have to gather results from childs, otherwise they go into
      // an undefined state. Then we throw.
      boost::mpi::gatherv(world, local_result.data(), (M + 1) * N_total, final_result, chunks_result, 0);
      throw;
    }

    // collect results from all workers
    boost::mpi::gatherv(world, local_result.data(), (M + 1) * N_total, final_result, chunks_result, 0);

    std::vector<double> res(N_world);

    Eigen::Map<Eigen::MatrixXd> final_result_eig(final_result, M+1, N_world);
    
    Eigen::Map<Eigen::RowVectorXd>(&res[0], N_world) = final_result_eig.row(0);

    for(std::size_t i=0; i != N_world; i++)
      if(unlikely(final_result_eig(0,i) == std::numeric_limits<double>::quiet_NaN()))
        throw std::runtime_error("MPI error");

    return(res);
  }

  template <>
  std::vector<boost::math::tools::promote_args<stan::math::var, double>::type>
  run_mpi_function(const std::vector<std::vector<stan::math::var> >& Theta,
                   const std::vector<std::vector<int> >& work_map,
                   const std::vector<std::vector<double> >& X_r,
                   const std::vector<std::vector<int> >& X_i,
                   std::ostream* pstream__) {
    boost::mpi::communicator world;

    const std::size_t M = Theta[0].size();
    const std::size_t J = Theta.size();
    const std::size_t W = world.size();


    std::vector<int> chunks(W, 0);
    std::vector<int> chunks_Theta(W, 0);
    std::vector<int> chunks_result(W, 0);

    std::vector<int > N_node;
    int N_total = 0;
    int N_world = 0;

    for(std::size_t j = 0; j != J; j++) {
      const int worker = work_map[j][0] - 1;
      const int N_job = work_map[j][1];
      
      chunks[worker] += 1;
      chunks_Theta[worker] += M;
      chunks_result[worker] += (1 + M) * N_job;
      
      if(worker == 0) {
        N_total += N_job;
        N_node.push_back(N_job);
      }

      N_world += N_job;
    }
    
    const std::size_t C = chunks[0];

    // scatter parameters
    Eigen::MatrixXd world_Theta(M, J);

    for(std::size_t j = 0; j != J; ++j)
      for(std::size_t m = 0; m != M; ++m)
        world_Theta(m,j) = stan::math::value_of(Theta[j][m]);

    //std::cout << "Shipping parameters to childs:" << std::endl << world_Theta << std::endl;

    Eigen::MatrixXd local_Theta(M, C);
    boost::mpi::scatterv(world, world_Theta.data(), chunks_Theta, local_Theta.data(), 0);
    
    // now work on the chunk for the root
    // calculate the chunk of the root

    // allocate storage for local function_result
    Eigen::MatrixXd local_result(M + 1, N_total);

    // allocate memory on AD stack for final result. Note that the
    // gradients and the function results will land there
    double* final_result
      = ChainableStack::memalloc_.alloc_array<double>( (M + 1) * N_world );

    // build up pointers in advance
    stan::math::vari** varis
      = ChainableStack::memalloc_.alloc_array<stan::math::vari*>(J*M);
    for(std::size_t j=0; j != J; j++) {
      for (int m = 0; m != M; ++m)
        varis[j * M + m] = Theta[j][m].vi_;
    }
    
    try {
      std::size_t nc = 0;
      for(std::size_t i = 0; i != C; i++) {
        mpi_util::jacobian(mpi_function_functor__(),
                           local_Theta.col(i),
                           X_r[i], X_i[i],
                           local_result.block(0, nc, M+1, N_node[i]));
        nc += N_node[i];
      }
    } catch(const std::exception& e) {
      // we have to gather results from childs, otherwise they go into
      // an undefined state. Then we throw.
      boost::mpi::gatherv(world, local_result.data(), (M + 1) * N_total, final_result, chunks_result, 0);
      throw;
    }

    // collect results from all workers
    boost::mpi::gatherv(world, local_result.data(), (M + 1) * N_total, final_result, chunks_result, 0);

    // finally build AD tree
    std::vector<stan::math::var> res;
    res.reserve(N_world);
    
    std::size_t nc = 0;
    for(std::size_t j=0; j != J; j++) {
      for(std::size_t k=0; k != work_map[j][1]; k++) {
        const double val = *(final_result + (M+1) * nc);
        if(unlikely(val == std::numeric_limits<double>::quiet_NaN()))
          throw std::runtime_error("MPI error");
        
        res.push_back(stan::math::var(new precomputed_gradients_vari(val, M, varis + j * M, final_result + (M+1) * nc + 1)));
        ++nc;
      }
    }

    return(res);
  }

  template <typename T0__, typename T1__>
  std::vector<std::vector<int> >
  setup_mpi_function(const std::vector<std::vector<T0__> >& Theta,
                     const std::vector<std::vector<T1__> >& X_r,
                     const std::vector<std::vector<int> >& X_i, std::ostream* pstream__) {

    boost::mpi::communicator world;

    const std::size_t J = Theta.size();
    const std::size_t M = Theta[0].size();
    const std::size_t W = world.size();
    const std::size_t C = J / W;
    const std::size_t CR = J % W;
    const std::size_t R = world.rank();

    // perform one double only function evaluation to learn the
    // resulting output sizes per job.
    std::vector<std::vector<int> > work_map(J, std::vector<int>(2, 0));

    for(std::size_t j=0; j != J; j++) {
      const std::vector<double> f = mpi_function(stan::math::value_of(Theta[j]), X_r[j], X_i[j], pstream__);
      work_map[j][0] = 1;
      work_map[j][1] = f.size();
    }
    
    if(W <= 1) {
      return(work_map);
    }
    
    // map the jobs to nodes
    std::vector<int> chunks(W, C);
    std::vector<int> chunk_start(W, 0);

    for(std::size_t j=0; j != CR; j++)
      ++chunks[j+1];

    std::size_t jc = 0;
    for(std::size_t w = 0; w != W; w++) {
      for(std::size_t i = 0; i != chunks[w]; i++) {
        work_map[jc + i][0] = w + 1;
      }
      chunk_start[w] = jc;
      jc += chunks[w];
    }

    std::cout << "Worker " << R+1 << " of " << W << " will handle " << chunks[R] << " of " << J << " chunks." << std::endl;

    if(R != 0) {
      std::cout << "Worker " << R+1 << " waits for work to come..." << std::endl;
      const std::size_t C_node = chunks[R];

      std::vector<std::vector<T0__> > X_r_node;
      std::vector<std::vector<int> > X_i_node;
      std::vector<int > N_node;
      int N_total = 0;

      for(std::size_t j=chunk_start[R]; j != chunk_start[R]+chunks[R]; ++j) {
        X_r_node.push_back(X_r[j]);
        X_i_node.push_back(X_i[j]);
        N_node.push_back(work_map[j][1]);
        N_total += work_map[j][1];
      }
      
      Eigen::MatrixXd local_Theta(M, C_node);
      
      // allocate storage for local function_result
      Eigen::MatrixXd local_result(M + 1, C_node * N_total);

      while(1) {
        // wait for some parameters 
        boost::mpi::scatterv(world, local_Theta.data(), M*C_node, 0);

        // now work on the chunk of parameters
        try {
          int nc = 0;
          for(std::size_t i = 0; i != C_node; i++) {
            mpi_util::jacobian(mpi_function_functor__(),
                               local_Theta.col(i),
                               X_r[i], X_i[i],
                               local_result.block(0, nc, M+1, N_node[i]));
            nc += N_node[i];
          }
        } catch(const std::exception& e) {
          // we only abort further processing. The root node will
          // detect unfinished runs and throw. Flag failure with a
          // NaN.
          local_result(0,0) = std::numeric_limits<double>::quiet_NaN();
        }
        
        // sent results to root
        boost::mpi::gatherv(world, local_result.data(), (M+1) * N_total, 0);
      }
    }
    return(work_map);
  }
  
};
