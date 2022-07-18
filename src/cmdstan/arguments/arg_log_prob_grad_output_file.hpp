#ifndef CMDSTAN_ARGUMENTS_ARG_LOG_PROB_GRAD_OUTPUT_FILE_HPP
#define CMDSTAN_ARGUMENTS_ARG_LOG_PROB_GRAD_OUTPUT_FILE_HPP

#include <cmdstan/arguments/singleton_argument.hpp>

namespace cmdstan {

class arg_log_prob_grad_output_file : public string_argument {
 public:
  arg_log_prob_grad_output_file() : string_argument() {
    _name = "grad_output_file";
    _description = "Gradients output file";
    _validity = "Path to existing file";
    _default = "grad_output.csv";
    _default_value = "grad_output.csv";
    _constrained = false;
    _good_value = "good";
    _value = _default_value;
  }
};

}  // namespace cmdstan

#endif
