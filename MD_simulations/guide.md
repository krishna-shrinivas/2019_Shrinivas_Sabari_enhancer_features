# Guide to running simulations

## Broad overview of code

This code runs coarse-grained Langevin dynamics simulations of a three-component system of DNA (species A), transcription factors (TFs or species B), and coactivators (species C). Specific and monovalent binding between TFs and "active" DNA binding sites is achieved through rigid body + embedded patches on the core of each particle. Further, any two particles can potentially interact through a LJ-type potential, and this mimics IDR-like interactions. Refer guide.md for more details.


## How to run code

The file md_poly.py is the main function that intializes and runs the molecular dynamics simulations. Various *required* and optional parameters are passed as flags to md_poly.py. Detailed listing of the various files are provided below. Examples of each type of file is provided in the repository along with the code.

## Setting simulation parameters

### Input parameters
Many overall simulation parameters are set by passing the flag "--i
input_file". The input_file contains a list of all input parameters, and each line contains the parameter_name, parameter_value, i.e.

* **parameter_name**, **parameter_value**

The description of various parameters are provided below:

| Parameter     | Description   | Requirements |
| ------------- |:-------------:|:-------------:|
| N_A     | Number of DNA molecules | non-negative |
| N_B     | Number of TF molecules  |non-negative |
| N_C     | Number of coactivator molecules|non-negative |
| seq_A   | Sequence of DNA molecule | must include 'A' type monomer i.e. ones constructed with binding sites for TFs (refer 'ABS' type)|
| seq_B     | Sequence of TF molecule | must include 'B' type - structured part of TF molecule with DNA binding domain  ( refer 'BDNA' type), IDR can incorporated as polymeric tail attached to B, for e.g. BDDDDDD or implicitly through the B type monomer    |
| seq_C     | Sequence of coactivator | must include 'C' type - structured part of coactivator with potential DNA binding domain  ( refer 'CDNA' type), IDR can incorporated as polymeric tail attached to C, for e.g. CDDDDDD or implicitly through the C type monomer       |
| T_A     | Relative temperature of A-type monomers      | 1 |
| dt     | Time step     | |
| L     | Length of periodic cube      | |
| L_A_box     | Length of sub-cube centered on origin for initializing DNA molecule | L_A < L|
| T     | Temperature of system      | Typically 1 (i.e kT=1) - set scales relative to that |
| N_bs_AB     | Number of binding sites on each A monomer for B    ||
| N_bs_AC     | Number of binding sites on each A monomer for C    | 0 in all reported simulations|
| log_step     | Observables logged every log_step    | |
| log_trajectory_step     | Trajectory data logged every log_trajectory_step    | |
| t_end  |  N_steps to run simulation |   |
| seed_position   | seed for randomly initializing positions of molecules  | seeds numpy.random |
|  drug_test |  Flag for disrupting TF-DNA interactions | Set 1 for True   |
|  t_post_drug |  Number of simulation steps after disrupting TF-DNA interactions | Only used if drug_test is 1   |
| p_drug  | Fraction of TF-DNA interactions to disrupt  |  1 disrupts all interactions b/w TF-DNA |
| movie_flag  | Creates and saves movies of trajectories  |  Stored under Trajectory/Movies/ |
| gamma_A  | Damping coefficient for DNA molecules  |  Typically higher than 1.0 that is set for TFs, coactivators |
| warm_up_flag  | If 1, warm up simulations with smaller dt steps are run  |  |
|  plot_flag |   Generates standard output plots of observables if set to 1|   |
| lb     | Bjeurrum length | Deprecated feature - no support|
| kappa     | Screening parameter for salt-mediated electrostatics   |Deprecated feature - no support|
| seed   | seed for Langevin dynamics  | Deprecated- instead randomly generated now  |


### Input particle information
Parameters for each monomer type are set by passing "--f input_particle_information". input_particle_information contains information for each monomer and patch-type (for binding sites) about type, mass, diameter, charge, and is_active_DNA. All lines starting with % are not read. Example line is:

* A,1.0,1.0,0,1

This sets the particle type to 'A', mass to 1.0 units, diameter to 1.0 units (same unit as length of simulation box), charge (variable is currently not supported and always zero) to netural('0'), and is_active_DNA set to 1.

Setting is_active_DNA to 1 ensures that A type monomers are constructed to potentially contain binding sites for TFs and coactivators. See above for setting number of binding sites/monomer for TF, and method for details.
The following particle types are required:
* A - DNA monomer with binding sites, is_active_DNA FLAG SET TO 1.
* B - DNA-binding structured domain of TFs
* C - Structured domain of coactivators
* ABS,ACS - Binding site patches on 'A' type monomer for TFs/coactivators. Size is one-fourth the length of a monomer. DO NOT SET is_active_DNA flag to 1 for patch, only for parent monomer.
* BDNA,CDNA - cognate binding patches within structured domains of B,C monomers. Size same as ABS/ACS patches in reported simulations.
* AW - Scaffold particle for constructing excluded volume around DNA monomer, size should be 1/3 diameter of A unit
* GS - Mutated binding patch. If drug_test is set to 1, then BDNA,CDNA patches are converted to GS type patches. Energetic interactions are set to hard-core between GS patches and cognate binding sites.

### Energetics parameters
Interaction parameters (LJ) between monomers are passed through "--e energetics_parameters". energetics_parameters contains the 'epsilon' LJ parameters between patches. All LJ parameters not-specified are assumed to be purely hard sphere, always cut off at ~2.5 sigma. Example line reads
*  **type_1,type_2,E12**
* ABS,BDNA,16

This line sets the energy between the ABS type binding site and BDNA binding patch in the structured domain of TFs to have a well-depth of ~ 16kT.

### Param list
If a parameter sweep is desired, than an optional parameter can be passed with  "--p param_list". param_list has the first line with parameter name (refer above table for modifiable parameters), and each subsequent line represents a distinct value. See example file provided with code.

## Simulation output organization

### Output folder + file name
The destination of the folder is set by passing ""--outfolder out_folder". Output will be written to /Output/out_folder/.
The prefix output_file passed in "--o output_file" will be added to all output files.

### Output data organization
All outputs are stored under the Output folder. A brief overview of the output structure is presented below:

Output/
* out_folder/
  - Each simulation condition with a different meta_parameter has its own sub-folder under out_folder. Automatically generated when param_list parameter is a meta_parameter.
  - parameter_condition/ -  standard format is 'NA_NB_NC_L_seed_position_seqA_seqB_seqC'
    - Data/
      - Observables/
        - parameter_condition/ (could be meta-parameter or other parameter such as N_bs_AB)
          - **output_file_seed_cluster_size.dat** - This file contains 4 columns and stores the timestep, size of largest cluster, radius of largest cluster, and density of largest cluster.
          - **output_file_seed_PE_T_KE.dat** - This file contains 4 columns and stores the timestep,potential energy, temperature (in kT units), and kinetic energy.
          - **output_file_seed_input_params** - This file contains a consolidated list of all system parameters used to run MD simulations.
          - **output_file_seed_warm_up_T.log** - This file contains the list initial time steps and system temperature during warm up.
      - Trajectory/
        - parameter_condition/
          - **output_file_seed.gsd** - Complete trajectory data
      - Figure/
        - parameter_condition/
          - 6 pngs for each simulation trajectory plotting the 6 measured observables.

### Analyzing simulation output
The attached python codes and files in Simulation_tools have several parsers and processing files that automatically navigate the above file  naming conventions and data organization.
