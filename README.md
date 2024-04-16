# WestAfricaModel

A model to simulate mosquito populations in West Africa

The scripts in this repository were written to study genetic (gene-drive) methods of controlling mosquito vectors of malaria, as desribed in Hancock, P.A, North, A et al.

Installation

CM.cpp is a c++ code for a simulation model of mosquito population dynamics, accompanied by the header file head.h. To run the model, download these source files and compile them locally (though we cannot guarentee they are compatible with all compilers). Alternatively, you may download RunModel, a compiled executable version of the c++ code, compiled with The Intel(R) C++ Compiler Classic (ICC) (version 2021.9.0). RunModel has been compiled and tested on a linux desktop computer (Ubuntu 20.04.6 LTS). The sub-directories ‘InputFiles’, ‘MuFiles’, ‘RainFiles’ and ‘ParameterFiles’ are needed to run the model in sixteen spatial settings, each an area of ~12,000km^2 in West Africa. 

Usage

To run the compiled simulation with default parameters, copy this repository, including sub-directories, onto your local machine. 

The simulation begins by reading numerous parameters from an input file, including the names of input files that will be required. 48 input files are provided in the folder ‘ParameterFiles’, which differ by the mosquito species and spatial setting. With respect to species, the files named with “gam” simulate the species pair Anopheles gambiae and An. Coluzzii, those named with “arab” simulate An. Arabiensis, and those named with “fun” simulate An. funestis. The spatial settings are (1 degree latitude by 1 degree longitude) squares in the following countries/regions:

    1. Western Mali
    2. Gambia/ Senegal
    3. Southern Mali
    4. Southern Niger
    5. Guinea Bissau/Senegal
    6. Northern Nigeria
    7. Western Burkina Faso
    8. Benin/ Burkina Faso
    9. Guinea
    10. Liberia/ Sierra Leone
    11. Côte d'Ivoire
    12. Benin/Togo
    13. Nigeria (Lagos)
    14. Cameroon
    15. Liberia
    16. Ghana

The model can be simulated with these parameter input files from the directory containing the input file directories (‘InputFiles’, ‘MuFiles’, ‘RainFiles’).

Example

For example, with this file architecture, a simulation of Anopheles arabiensis in Southern Mali can be executed in a bash script with the command:

./RunModel<./ParameterFiles/Parameters_arab_3.csv>Output_arab_3.txt &


Output

A simulation run will write the following files:

    • “Totals_[species]_[setting]_index_X.txt” where [species] is either ‘gam’, ‘arab’, or ‘fun’ as above; setting is a number from 1 to 16 to represent the spatial setting as above, and the index ‘X’ is an integer index specified in the parameter file which may be useful if there are multiple runs with the same species and setting. This file outputs the total numbers of adult female mosquitoes (for the species being simulated), with each of six genotypes considered in the simulation model, for each day of the simulation run. The columns of the txt file are: day number (column 1), W/W females (c2), W/D females (c3), D/D females (c4), W/R females (c5), R/R females (c6) and D/R females (c7). Here W refers to the ‘wildtype’ allele, ‘D’ refers to the gene-drive allele, and ‘R’ refers to the non-functional resistant allele.
    • “Emergence_[species]_[setting]_index_X.txt”. This file outputs a list of the daily number of biting female mosquitoes (W/W plus W/D plus W/R; for the species being simulated) that emerge in the study area (one value per day).
    • “Local_[species]_[setting]_index_X.txt”. This file writes the numbers of males (in each of the six genotype categories) and females (also for the six genotypes) in each local population, at time points specified in the input file. The default parameter files are set to not output into these files.



Parameter specifics

The default parameters in the ‘ParameterFiles’ directory have the following settings.

    • Duration. The maximum possible duration is 46 years. The parameter files are set to start at the start of year 30 (day 10951) and run until the end of year 45 (day 16425).
    • Gene drive introduction. This occurs on day 150 of the third year after the simulation begins (day 11830). The introduction is the addition of 1000 W/D type male mosquitoes into 50 local sites.
    • Gene drive characteristics. The default assumed fitness cost of W/D females is 0.35, and it is assumed that D/D, R/R and D/R females are all sterile (and non-biting). The default homing rate is 0.95.
    • Replication. The default files are set to execute one run of the simulation and then stop.
