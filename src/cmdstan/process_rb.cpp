#include <boost/lexical_cast.hpp>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

// clang++ -O3 -Istan/lib/stan_math/lib/eigen_3.2.8/ -Istan/lib/stan_math/lib/boost_1.60.0/ -o process_rb src/cmdstan/process_rb.cpp

void print_usage() {
  std::cout << "USAGE: process_rb input=<inputname> [output=<outputname>]"
            << std::endl;
  std::cout << "OPTIONS:" << std::endl;
  std::cout << "  input=<inputname>\tInput file created by a "
            << "Rao-Blackwellized sampler"
            << std::endl;
  std::cout << "  output=<outputname>\tReduced output file "
            << "(Defaults to \"red_<inputname>\")"
            << std::endl;
}

/**
 * Reduce Rao-Blackwellized sampler output
 *
 * @param argc Number of arguments
 * @param argv Arguments
 *
 * @return 0 for success,
 *         non-zero otherwise
 */
int main(int argc, const char* argv[]) {

  if (!(argc == 2 || argc == 3)) {
    print_usage();
    return 1;
  }

  // Parse any arguments
  std::string input_name("");
  std::string input_arg(argv[1]);

  if (input_arg.find("input=") == std::string::npos) {
    print_usage();
    return 1;
  } else {
    size_t equal = input_arg.find("=");
    input_name = input_arg.substr(equal + 1, input_arg.size());
  }

  std::string output_name("");

  if (argc == 3) {
    std::string output_arg(argv[2]);
    if (output_arg.find("output=") == std::string::npos) {
      print_usage();
      return 1;
    } else {
      size_t equal = output_arg.find("=");
      output_name = output_arg.substr(equal + 1, input_arg.size());
    }
  } else {
    output_name = "red_" + input_name;
  }

  std::ifstream input_file;
  input_file.open(input_name);
  if (!input_file.good()) {
    std::cout << "Can't open " << input_name << " for reading" << std::endl;
    return 1;
  }

  std::ofstream output_file;
  output_file.open(output_name);
  if (!output_file.good()) {
    std::cout << "Can't open " << output_name << " for writing" << std::endl;
    return 1;
  }

  bool samples = false;
  bool first_sample = false;
  size_t n_params = 0;
  size_t n_sampler_params = 0;

  double n_samples = 0;
  double sum_weights = 0;
  double running_ave_accept = 0;
  Eigen::VectorXd running_ave_params;

  std::vector<std::string> sampler_params;

  while (input_file.good()) {
    std::string line;
    std::getline(input_file, line);

    // Grab sizes from parameter header
    if (line.substr(0, 4) == "lp__") {
      size_t n_total = std::count(line.begin(), line.end(), ',') + 1;
      size_t n_under = (std::count(line.begin(), line.end(), '_') - 1) / 2;
      n_params = n_total - n_under + 1;
      n_sampler_params = n_under - 2;

      running_ave_params.resize(n_params);

      output_file << line << std::endl;
      continue;
    }

    if (line == "# BEGIN SAMPLES") {
      samples = true;
      first_sample = true;
      continue;
    }

    if (line == "# END SAMPLES") {
      samples = false;
      continue;
    }

    if (!samples) {
      output_file << line << std::endl;
      continue;
    }

    if (line == "# ") {
      output_file << running_ave_params(0) << ",";
      output_file << running_ave_accept;

      // output sampler params
      for (size_t n = 0; n < n_sampler_params; ++n)
        output_file << "," << sampler_params.at(n);

      for (size_t n = 0; n < n_params - 2; ++n)
        output_file << "," << running_ave_params(n + 1);
      output_file << "," << running_ave_params(n_params - 1) << std::endl;

      // Reset averages
      n_samples = 0;
      running_ave_accept = 0;

      sum_weights = 0;
      running_ave_params.setZero();

      sampler_params.clear();
      first_sample = true;

      continue;
    }

    std::string value;
    std::stringstream lstream(line);
    Eigen::VectorXd params(n_params);

    // log_prob
    std::getline(lstream, value, ',');
    params(0) = boost::lexical_cast<double>(value);

    // weight
    std::getline(lstream, value, ',');
    double weight = boost::lexical_cast<double>(value);
    double accept = weight > 1 ? 1 : weight;

    ++n_samples;
    double delta = accept - running_ave_accept;
    running_ave_accept += delta / n_samples;

    // sampler parameters
    if (first_sample) {
      for (size_t n = 0; n < n_sampler_params; ++n) {
        std::getline(lstream, value, ',');
        sampler_params.push_back(value);
      }
      first_sample = false;
    } else {
      for (size_t n = 0; n < n_sampler_params; ++n)
        std::getline(lstream, value, ',');
    }

    // model parameters
    for (size_t n = 0; n < n_params - 2; ++n) {
      std::getline(lstream, value, ',');
      params(n + 1) = boost::lexical_cast<double>(value);
    }
    std::getline(lstream, value, '\n');
    params(n_params - 1) = boost::lexical_cast<double>(value);

    // weighted Welford
    sum_weights += weight;
    params -= running_ave_params;
    running_ave_params += params / sum_weights;
  }

  input_file.close();
  output_file.close();

  return 0;

}
