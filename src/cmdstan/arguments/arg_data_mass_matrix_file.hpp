#ifndef CMDSTAN_ARGUMENTS_ARG_DATA_MASS_MATRIX_FILE_HPP
#define CMDSTAN_ARGUMENTS_ARG_DATA_MASS_MATRIX_FILE_HPP

#include <cmdstan/arguments/singleton_argument.hpp>

namespace cmdstan {

  class arg_data_mass_matrix_file : public string_argument {
  public:
    arg_data_mass_matrix_file(): string_argument() {
      _name = "mass_matrix_file";
      _description = "Input mass matrix file";
      _validity = "Path to existing file";
      _default = "\"mass.R\"";
      _default_value = "";
      _constrained = false;
      _good_value = "good";
      _value = _default_value;
    }
  };

}
#endif
