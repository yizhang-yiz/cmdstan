file(GLOB_RECURSE MODELS
     RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
     *.stan)

set(COMPILED_MODELS "")
foreach(stanfile ${MODELS})
  string(REPLACE ".stan" "${empty_string}" executable ${stanfile})
  string(REPLACE "/"  "_" executable ${executable})
  string(REPLACE "\\"  "_" executable ${executable})
  set(executable cmdstan_unitmodels_${executable})

  stan_model(${executable} ${stanfile} CMDSTAN)
  list(APPEND COMPILED_MODELS ${executable})
endforeach()
add_custom_target(CmdStan_Test_Models DEPENDS ${COMPILED_MODELS})

