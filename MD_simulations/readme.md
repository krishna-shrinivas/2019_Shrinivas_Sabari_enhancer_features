# Coarse grained molecular dynamics of chromatin, transcription factors, and coactivators

## Dependencies

* Python3
* Libraries and versions (since many have developed very fast)
  - NumPy - 1.11
  - HOOMD - 2.2.1
  - Fresnel - 0.5.0
  - Freud - 0.6.2
  - GSD - 1.2.0
  - Pillow - 4.2.1

## How to run code

The file md_poly.py is the main function that intializes and runs the molecular dynamics simulations. Various *required* and optional parameters are passed as flags to md_poly.py.

### Required parameters

* --i input_file             

input_file is a list of all input parameters including the number of chains, sequence of different types, size of simulation of box,
size of simulation box for chain A (DNA) initialization, length of simulation, dt for simulation, bond strength, bjeurrum length, screening parameter, drug characterization, and p_drug. Please refer guide.MD for more comprehensive documentation or refer example files.

* --o output_file

output_file is the name of the output_file to which three types of output are written: observables (cluster_size, PE, input_params),images of observables, and trajectory_file.

* --e energy_file

energy_file contains the LJ parameters between patches. All LJ parameters not-specified are assumed to be purely hard sphere.

* --f monomer_information_file

monomer_information_file contains information for each monomer and patch-type (for binding sites) about mass, diameter, and charge.

### Optional parameters

* --p params_file

File with parameter name to sweep on first line, and each subsequent line is the value to sweep

* --outfolder

Destination folder for output files

## Broad overview of code

This code runs coarse-grained Langevin dynamics simulations of a three-component system of DNA (species A), transcription factors (species B), and coactivators (species C). Specific and monovalent binding between TFs and "active" DNA binding sites is achieved through rigid body + embedded patches on the core of each particle. Further, any two particles can potentially interact through a LJ-type potential, and this mimics IDR-like interactions. Refer guide.md for more details.

## Code output
Typically, successful simulations create 3 types of subfolders under outfolder for any given simulation condition. The Observables sub-folder contains text files of key observables for each trajectory, Figures/ contains default figures of observables, and Trajectory/ contains gsd files of complete trajectory data. To analyze and process output data, please use python provided tools in Simulation_tools.
