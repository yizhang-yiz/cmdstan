#!/bin/sh

## new argument: mass_matrix_file
## mass_matrix_file must have step_size and inverse_mass_matrix


make examples/bernoulli/bernoulli examples/blocker/blocker

./examples/bernoulli/bernoulli sample data file=examples/bernoulli/bernoulli.data.R mass_matrix_file=examples/bernoulli/bernoulli.mass.R output file=bernoulli.csv

./examples/blocker/blocker sample data file=examples/blocker/blocker.data.R mass_matrix_file=examples/blocker/blocker.mass.R output file=blocker.csv
