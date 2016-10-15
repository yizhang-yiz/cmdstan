#ifndef CMDSTAN_COMMAND_HPP
#define CMDSTAN_COMMAND_HPP

#include <cmdstan/arguments/argument_parser.hpp>
#include <cmdstan/arguments/arg_data.hpp>
#include <cmdstan/arguments/arg_id.hpp>
#include <cmdstan/arguments/arg_init.hpp>
#include <cmdstan/arguments/arg_output.hpp>
#include <cmdstan/arguments/arg_random.hpp>
#include <cmdstan/write_model.hpp>
#include <cmdstan/write_stan.hpp>
#include <stan/callbacks/interrupt.hpp>
#include <stan/callbacks/noop_interrupt.hpp>
#include <stan/callbacks/writer.hpp>
#include <stan/callbacks/noop_writer.hpp>
#include <stan/callbacks/stream_writer.hpp>
#include <stan/io/dump.hpp>
#include <stan/services/diagnose/diagnose.hpp>
#include <stan/services/optimize/bfgs.hpp>
#include <stan/services/optimize/lbfgs.hpp>
#include <stan/services/optimize/newton.hpp>
#include <stan/services/sample/fixed_param.hpp>
#include <stan/services/sample/hmc_nuts_dense_e.hpp>
#include <stan/services/sample/hmc_nuts_dense_e_adapt.hpp>
#include <stan/services/sample/hmc_nuts_diag_e.hpp>
#include <stan/services/sample/hmc_nuts_diag_e_adapt.hpp>
#include <stan/services/sample/hmc_nuts_unit_e.hpp>
#include <stan/services/sample/hmc_nuts_unit_e_adapt.hpp>
#include <stan/services/sample/hmc_static_dense_e.hpp>
#include <stan/services/sample/hmc_static_dense_e_adapt.hpp>
#include <stan/services/sample/hmc_static_diag_e.hpp>
#include <stan/services/sample/hmc_static_diag_e_adapt.hpp>
#include <stan/services/sample/hmc_static_unit_e.hpp>
#include <stan/services/sample/hmc_static_unit_e_adapt.hpp>
#include <stan/services/experimental/advi/fullrank.hpp>
#include <stan/services/experimental/advi/meanfield.hpp>
#include <cmdstan/arguments/arg_adapt.hpp>
#include <cmdstan/arguments/arg_adapt_delta.hpp>
#include <cmdstan/arguments/arg_adapt_engaged.hpp>
#include <cmdstan/arguments/arg_adapt_gamma.hpp>
#include <cmdstan/arguments/arg_adapt_init_buffer.hpp>
#include <cmdstan/arguments/arg_adapt_kappa.hpp>
#include <cmdstan/arguments/arg_adapt_t0.hpp>
#include <cmdstan/arguments/arg_adapt_term_buffer.hpp>
#include <cmdstan/arguments/arg_adapt_window.hpp>
#include <cmdstan/arguments/arg_bfgs.hpp>
#include <cmdstan/arguments/arg_data.hpp>
#include <cmdstan/arguments/arg_data_file.hpp>
#include <cmdstan/arguments/arg_dense_e.hpp>
#include <cmdstan/arguments/arg_diag_e.hpp>
#include <cmdstan/arguments/arg_diagnose.hpp>
#include <cmdstan/arguments/arg_diagnostic_file.hpp>
#include <cmdstan/arguments/arg_engine.hpp>
#include <cmdstan/arguments/arg_fail.hpp>
#include <cmdstan/arguments/arg_fixed_param.hpp>
#include <cmdstan/arguments/arg_history_size.hpp>
#include <cmdstan/arguments/arg_hmc.hpp>
#include <cmdstan/arguments/arg_id.hpp>
#include <cmdstan/arguments/arg_init.hpp>
#include <cmdstan/arguments/arg_init_alpha.hpp>
#include <cmdstan/arguments/arg_int_time.hpp>
#include <cmdstan/arguments/arg_iter.hpp>
#include <cmdstan/arguments/arg_lbfgs.hpp>
#include <cmdstan/arguments/arg_max_depth.hpp>
#include <cmdstan/arguments/arg_method.hpp>
#include <cmdstan/arguments/arg_metric.hpp>
#include <cmdstan/arguments/arg_newton.hpp>
#include <cmdstan/arguments/arg_num_samples.hpp>
#include <cmdstan/arguments/arg_num_warmup.hpp>
#include <cmdstan/arguments/arg_nuts.hpp>
#include <cmdstan/arguments/arg_optimize.hpp>
#include <cmdstan/arguments/arg_optimize_algo.hpp>
#include <cmdstan/arguments/arg_output.hpp>
#include <cmdstan/arguments/arg_output_file.hpp>
#include <cmdstan/arguments/arg_random.hpp>
#include <cmdstan/arguments/arg_refresh.hpp>
#include <cmdstan/arguments/arg_rwm.hpp>
#include <cmdstan/arguments/arg_sample.hpp>
#include <cmdstan/arguments/arg_sample_algo.hpp>
#include <cmdstan/arguments/arg_save_iterations.hpp>
#include <cmdstan/arguments/arg_save_warmup.hpp>
#include <cmdstan/arguments/arg_seed.hpp>
#include <cmdstan/arguments/arg_static.hpp>
#include <cmdstan/arguments/arg_stepsize.hpp>
#include <cmdstan/arguments/arg_stepsize_jitter.hpp>
#include <cmdstan/arguments/arg_test.hpp>
#include <cmdstan/arguments/arg_test_grad_eps.hpp>
#include <cmdstan/arguments/arg_test_grad_err.hpp>
#include <cmdstan/arguments/arg_test_gradient.hpp>
#include <cmdstan/arguments/arg_thin.hpp>
#include <cmdstan/arguments/arg_tolerance.hpp>
#include <cmdstan/arguments/arg_unit_e.hpp>
#include <cmdstan/arguments/argument.hpp>
#include <cmdstan/arguments/argument_parser.hpp>
#include <cmdstan/arguments/argument_probe.hpp>
#include <cmdstan/arguments/categorical_argument.hpp>
#include <cmdstan/arguments/list_argument.hpp>
#include <cmdstan/arguments/singleton_argument.hpp>
#include <cmdstan/arguments/unvalued_argument.hpp>
#include <cmdstan/arguments/valued_argument.hpp>
#include <cmdstan/nuts_with_mass_matrix.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <fstream>
#include <string>
#include <vector>

namespace cmdstan {

  stan::io::dump get_var_context(const std::string file) {
    std::fstream stream(file.c_str(), std::fstream::in);
    stan::io::dump var_context(stream);
    stream.close();
    return var_context;
  }

  template <class Model>
  int command(int argc, const char* argv[]) {
    stan::callbacks::stream_writer info(std::cout);
    stan::callbacks::stream_writer err(std::cout);


    // Read arguments
    std::vector<argument*> valid_arguments;
    valid_arguments.push_back(new arg_id());
    valid_arguments.push_back(new arg_data());
    valid_arguments.push_back(new arg_init());
    valid_arguments.push_back(new arg_random());
    valid_arguments.push_back(new arg_output());
    argument_parser parser(valid_arguments);
    int err_code = parser.parse_args(argc, argv, info, err);
    if (err_code != 0) {
      std::cout << "Failed to parse arguments, terminating Stan" << std::endl;
      return err_code;
    }
    if (parser.help_printed())
      return err_code;
    u_int_argument* random_arg = dynamic_cast<u_int_argument*>(parser.arg("random")->arg("seed"));
    if (random_arg->is_default()) {
      random_arg->set_value((boost::posix_time::microsec_clock::universal_time() - boost::posix_time::ptime(boost::posix_time::min_date_time)).total_milliseconds());
    }
    parser.print(info);
    info();


    stan::callbacks::noop_writer init_writer;
    stan::callbacks::noop_interrupt interrupt;

    stan::io::dump data_var_context(get_var_context(dynamic_cast<string_argument*>(parser.arg("data")->arg("file"))->value()));

    std::fstream output_stream(dynamic_cast<string_argument*>(parser.arg("output")->arg("file"))->value().c_str(),
                               std::fstream::out);
    stan::callbacks::stream_writer sample_writer(output_stream, "# ");

    std::fstream diagnostic_stream(dynamic_cast<string_argument*>(parser.arg("output")->arg("diagnostic_file"))->value().c_str(),
                                   std::fstream::out);
    stan::callbacks::stream_writer diagnostic_writer(diagnostic_stream, "# ");


    //////////////////////////////////////////////////
    //                Initialize Model              //
    //////////////////////////////////////////////////
    Model model(data_var_context, &std::cout);
    write_stan(sample_writer);
    write_model(sample_writer, model.model_name());
    parser.print(sample_writer);

    write_stan(diagnostic_writer);
    write_model(diagnostic_writer, model.model_name());
    parser.print(diagnostic_writer);


    int refresh = dynamic_cast<int_argument*>(parser.arg("output")->arg("refresh"))->value();
    unsigned int id = dynamic_cast<int_argument*>(parser.arg("id"))->value();
    unsigned int random_seed = dynamic_cast<u_int_argument*>(parser.arg("random")->arg("seed"))->value();

    std::string init = dynamic_cast<string_argument*>(parser.arg("init"))->value();
    double init_radius = 2.0;
    try {
      init_radius = boost::lexical_cast<double>(init);
      init = "";
    } catch (const boost::bad_lexical_cast& e) {
    }
    stan::io::dump init_context(get_var_context(init));
stan::io::dump mass_matrix_context(get_var_context(dynamic_cast<string_argument*>(parser.arg("data")->arg("mass_matrix_file"))->value()));    

    int return_code = stan::services::error_codes::CONFIG;
    int num_samples = dynamic_cast<int_argument*>(parser.arg("method")->arg("sample")->arg("num_samples"))->value();
    list_argument* algo
      = dynamic_cast<list_argument*>
      (parser.arg("method")->arg("sample")->arg("algorithm"));

    categorical_argument* hmc
      = dynamic_cast<categorical_argument*>
      (algo->arg("hmc"));
    double stepsize
      = dynamic_cast<real_argument*>
      (hmc->arg("stepsize"))->value();
    double stepsize_jitter = dynamic_cast<real_argument*>
      (hmc->arg("stepsize_jitter"))->value();
    int num_thin = dynamic_cast<int_argument*>(parser.arg("method")->arg("sample")->arg("thin"))->value();
    categorical_argument* base
      = dynamic_cast<categorical_argument*>
      (algo->arg("hmc")->arg("engine")->arg("nuts"));
    int max_depth
      = dynamic_cast<int_argument*>(base->arg("max_depth"))
      ->value();
    
    return_code = nuts_with_mass_matrix(model, init_context, mass_matrix_context,
                                        random_seed, id, init_radius,
                                        num_samples, num_thin,
                                        refresh, stepsize, stepsize_jitter,
                                        max_depth, interrupt,
                                        info, err, init_writer,
                                        sample_writer, diagnostic_writer);
    
    output_stream.close();
    diagnostic_stream.close();
    for (size_t i = 0; i < valid_arguments.size(); ++i)
      delete valid_arguments.at(i);
    return return_code;
  }

}
#endif
