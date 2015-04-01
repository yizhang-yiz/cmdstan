#include <algorithm>
#include <iostream>
#include <iomanip>
#include <ios>
#include <stan/mcmc/chains.hpp>
#include <stan/command/print.hpp>

/**
 * The Stan print function.
 *
 * @param argc Number of arguments
 * @param argv Arguments
 * 
 * @return 0 for success, 
 *         non-zero otherwise
 */
int main(int argc, const char* argv[]) {
  
  if (argc == 1) {
    return 0;
  }
  
  // Parse any arguments specifying filenames
  std::ifstream ifstream;
  std::vector<std::string> filenames;
  
  for (int i = 1; i < argc; i++) {
    ifstream.open(argv[i]);
    
    if (ifstream.good()) {
      filenames.push_back(argv[i]);
      ifstream.close();
    } else {
      std::cout << "File " << argv[i] << " not found" << std::endl;
    }
  }
  
  if (!filenames.size()) {
    std::cout << "No valid input files, exiting." << std::endl;
    return 0;
  }
  
  // Parse specified files
  ifstream.open(filenames[0].c_str());
  
  stan::io::stan_csv stan_csv = stan::io::stan_csv_reader::parse(ifstream);
  
  stan::mcmc::chains<> chains(stan_csv);
  ifstream.close();
  
  for (std::vector<std::string>::size_type chain = 1; 
       chain < filenames.size(); chain++) {
    ifstream.open(filenames[chain].c_str());
    stan_csv = stan::io::stan_csv_reader::parse(ifstream);
    chains.add(stan_csv);
    ifstream.close();
  }
  
  // Diagnostics header
  if (   stan_csv.metadata.algorithm != "hmc"
      || stan_csv.metadata.engine != "nuts") {
    std::cout << "Diagnostics available for only samples generated with NUTS" << std::endl;
    return 0;
  }
  
  std::cout << "WARNING: The following diagnostics are necessary for well-behaved"
            << std::endl;
  std::cout << "Markov Chain Monte Carlo estimators but they are not necessarily sufficient!"
            << std::endl << std::endl;
  
  // Divergences
  Eigen::VectorXd n_divergent = chains.samples("n_divergent__");

  int n_divergent_iterations = 0;
  for (int n = 0; n < n_divergent.size(); ++n)
    if (n_divergent(n) != 0) ++n_divergent_iterations;
  
  if (n_divergent_iterations == 0) {
    std::cout << "No divergences prevented NUTS from exploring the posterior."
              << std::endl;
  } else {
    std::cout << n_divergent_iterations << " of " << n_divergent.size()
              << " Markov transitions ended when the sampler encountered"
              << " a divergence.  This indicates that the Markov chain was"
              << " not able to explore all of the posterior and Markov Chain"
              << " Monte Carlo estimators may be biased.  Try running again"
              << " with a higher target acceptance probability."
              << std::endl;
  }
  std::cout << std::endl;
  
  // Treedepth
  int max_observed_depth = chains.samples("treedepth__").maxCoeff();
  int max_treedepth = 10; // FIXME: Need to grab this from file which requires updating stan_csv_reader
  
  if (max_observed_depth < max_treedepth) {
    std::cout << "No treedepth limits prevents NUTS from exploring the posterior."
              << std::endl;
  } else {
    std::cout << "At least one Markov transition was terminated prematurely "
                 "because the NUTS tree hit the max tree depth limit.  This "
                 "can lead to inefficient exploration -- try rerunning with "
                 "a higher max_depth setting."
              << std::endl;
  }
  std::cout << std::endl;
  
  // R_hat
  std::vector<std::string> slow_params;
  for (int n = 0; n < chains.num_params(); n++) {
    if (chains.split_potential_scale_reduction(n) > 1.1)
      slow_params.push_back(chains.param_name(n));
  }
  
  if (slow_params.size() == 0) {
    std::cout << "All parameters had potential scale reduction factors, or R_hats, "
              << "below 1.1, hence there is no indication of poor mixing."
              << std::endl;
  } else {
    if (slow_params.size() == 1) {
      std::cout << "The parameter " << slow_params.at(0) << " had a potential scale "
                << "reduction factor, or R_hat, above 1.1 indicating that the "
                << "Markov chain is not mixing well and Markov Chain Monte Carlo "
                << "estimators may be biased.  Try running again with more samples."
                << std::endl;
    } else if (slow_params.size() == 2) {
      std::cout << "The parameters " << slow_params.at(0) << " and " << slow_params.at(1)
                << "had potential scale "
                << "reduction factors, or R_hats, above 1.1 indicating that the "
                << "Markov chain is not mixing well and Markov Chain Monte Carlo "
                << "estimators may be biased.  Try running again with more samples."
                << std::endl;
    } else {
      std::cout << "The parameters ";
      for (int n = 0; n < slow_params.size() - 1; ++n) {
        std::cout << slow_params.at(n) << ", ";
      }
      std::cout << " and " << slow_params.back() << " all had potential scale "
                << "reduction factors, or R_hats, above 1.1 indicating that the "
                << "Markov chain is not mixing well and Markov Chain Monte Carlo "
                << "estimators may be biased.  Try running again with more samples."
                << std::endl;
    }
  }
  
  return 0;
        
}



