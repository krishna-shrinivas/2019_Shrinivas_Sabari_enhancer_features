# Analyses post HOOMD simulations
## Dependencies
* HOOMD associated libraries - Fresnel, GSD, Freud, HOOMD
* Python Dependencies - Numpy, PIL, JSON, matplotlib

## Functions
### Parsing output data organization
Simulation runs are stored under an *out_folder* (refer simulation package).

parse_MD_output.py contains a function **return_data** which takes in the path to *out_folder* as input and returns the following:
*  param_name & param_range
* input_parameters and unique paths for each different simulation condition within out_folder

If only one condition exists (i.e. no parameter range), then param_name ***'default'*** and param_value = **1.0** are returned.
