#ifndef CMDSTAN_NUTS_WITH_MASS_MATRIX_HPP
#define CMDSTAN_NUTS_WITH_MASS_MATRIX_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/callbacks/interrupt.hpp>
#include <stan/callbacks/writer.hpp>
#include <stan/mcmc/fixed_param_sampler.hpp>
#include <stan/services/error_codes.hpp>
#include <stan/mcmc/hmc/nuts/adapt_dense_e_nuts.hpp>
#include <stan/mcmc/hmc/nuts/dense_e_nuts.hpp>
#include <stan/services/util/run_sampler.hpp>
#include <stan/services/util/rng.hpp>
#include <stan/services/util/initialize.hpp>
#include <vector>

namespace cmdstan {

  /**
   * Runs HMC with NUTS with dense Euclidean
   * metric without adapatation.
   *
   * @tparam Model Model class
   *
   * @param model Input model to test (with data already instantiated)
   * @param init var context for initialization
   * @param random_seed random seed for the pseudo random number generator
   * @param chain chain id to advance the pseudo random number generator
   * @param init_radius radius to initialize
   * @param num_warmup Number of warmup samples
   * @param num_samples Number of samples
   * @param num_thin Number to thin the samples
   * @param save_warmup Indicates whether to save the warmup iterations
   * @param refresh Controls the output
   * @param stepsize initial stepsize for discrete evolution
   * @param stepsize_jitter uniform random jitter of stepsize
   * @param max_depth Maximum tree depth
   * @param interrupt Callback for interrupts
   * @param message_writer Writer for messages
   * @param error_writer Writer for messages
   * @param init_writer Writer callback for unconstrained inits
   * @param sample_writer Writer for draws
   * @param diagnostic_writer Writer for diagnostic information
   * @return error code; 0 if no error
   */
  template <class Model>
  int nuts_with_mass_matrix(Model& model,
                            stan::io::var_context& init,
                            stan::io::var_context& mass_matrix,
                            unsigned int random_seed,
                            unsigned int chain,
                            double init_radius,
                            int num_samples,
                            int num_thin,
                            int refresh,
                            double stepsize,
                            double stepsize_jitter,
                            int max_depth,
                            stan::callbacks::interrupt&
                            interrupt,
                            stan::callbacks::writer&
                            message_writer,
                            stan::callbacks::writer&
                            error_writer,
                            stan::callbacks::writer&
                            init_writer,
                            stan::callbacks::writer&
                            sample_writer,
                            stan::callbacks::writer&
                            diagnostic_writer) {
    stan::services::sample::mcmc_writer
      writer(sample_writer, diagnostic_writer, message_writer);

    boost::ecuyer1988 rng = stan::services::util::rng(random_seed, chain);

    std::vector<int> disc_vector;
    std::vector<double> cont_vector
      = stan::services::util::initialize(model, init, rng, init_radius,
                                         true,
                                         message_writer, init_writer);
    
    Eigen::Map<Eigen::VectorXd> cont_params(cont_vector.data(),
                                            cont_vector.size());
    stan::mcmc::sample s(cont_params, 0, 0);

    
    //stan::mcmc::adapt_dense_e_nuts<Model, boost::ecuyer1988> sampler(model, rng);
    stan::mcmc::dense_e_nuts<Model, boost::ecuyer1988> sampler(model, rng);
    sampler.set_nominal_stepsize(stepsize);
    sampler.set_stepsize_jitter(stepsize_jitter);
    sampler.set_max_depth(max_depth);


    std::vector<double> step_size = mass_matrix.vals_r("step_size");
    std::vector<double> inverse_mass_matrix = mass_matrix.vals_r("inverse_mass_matrix");

    sampler.set_nominal_stepsize(step_size[0]);
    sampler.z().mInv = Eigen::Map<Eigen::MatrixXd>(&inverse_mass_matrix[0],
                                                   model.num_params_r(),
                                                   model.num_params_r());
     // Headers
    writer.write_sample_names(s, sampler, model);
    writer.write_diagnostic_names(s, sampler, model);
    
    writer.write_adapt_finish(sampler);

    clock_t start = clock();
    stan::services::util::generate_transitions
      (sampler, num_samples, 0, num_samples, num_thin,
       refresh, true, false,
       writer,
       s, model, rng,
       interrupt, message_writer, error_writer);
    clock_t end = clock();
    double sample_delta_t
      = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    
    writer.write_timing(0, sample_delta_t);

    return stan::services::error_codes::OK;
  }

}
#endif
