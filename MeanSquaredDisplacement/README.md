### This C++ code calculates the following:

1) Averaged mean squared displacement (MSD) of the inner monomers in a chain (g<sub>1</sub>(t)) (or) the center of mass (COM) MSD of the chain (g<sub>3</sub>(t)) in a file named "MSDValues.txt"
2) Non-Gaussian parameter (alpha2) values in a file named "MSDValues.txt"

Note that the "MSDValues.txt" will have 5 columns: The first column is the timestep, the second column is timestep * time unit conversion, the third column is the MSD values, the fourth column is the square of the MSD values, and the fifth column is the non-Gaussian parameter (alpha2) values.

### Requirements

* This code requires dump files named "DumpFile.TimeStep" sorted based on atom-ids and has to be in the following format: "id mol type q mass xu yu zu ix iy iz"  
* xu, yu, zu refers to unwrapped coordinates
* For the center of mass MSD, each chain must have a unique molecule ID.  That is, if the system contains 500 chains, molecule ID should go from 1 to 500.  

### Required files

* This code comes with 3 files: (1) one .cpp file, (2) one .h header file, (3) one parameter text file

* In the parameter file, one needs to input the information about the system and the simulation details like trajectory file name, number of chains, beads per chain, number of inner monomers, start time, end time, interval time, number of atom types of interest, list of atom types, time unit conversion, and the COM flag.  If COMFlag == 'no', averaged mean squared displacement (MSD) of the individual beads in a chain is computed and if COMFlag == 'yes', COM MSD of the chain is computed

* The parameter 'NumberofInnerMonomers' in the parameter text file gives the user the flexibility to choose the number of inner monomers in a polymer chain to compute g<sub>1</sub>(t).  For example, if the number of monomers in a polymer chain is 50, 30 inner monomers will be selected as 15 each on either side of the central monomer in each polymer chain when specified.  For COM MSD, all the monomers in a polymer chain will be selected.  

### Compile and run the C++ code as below

* If using a GCC compiler:  
    gcc -lstdc++ -lm  MSD.cpp   
    ./a.out
    
* If using a ICC compiler:  
    icc MSD.cpp  
    ./a.out
    
### Other details

This code does moving average over all the frames and also averages over all the atoms of interest selected based on the atom types given by the user in the parameter file


