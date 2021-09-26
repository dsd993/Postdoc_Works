### This C++ code calculates the following:

1) Averaged mean squared displacement (MSD) of the individual beads in a chain (or) the COM MSD of the chain in a file named "MSDValues.txt"
2) Non-Gaussian parameter (alpha2) values in a file named "MSDValues.txt"

Note that the "MSDValues.txt" will have 5 columns: The first column is the timestep, the second column is timestep * time unit conversion, the third column is the MSD values, the fourth column is the square of the MSD values, and the fifth column is the non-Gaussian parameter (alpha2) values.

### Requirements

* This code requires dump files named "DumpFile.TimeStep" sorted based on atom-ids and has to be in the following format: "id mol type q mass xu yu zu ix iy iz"  
* xu, yu, zu refers to unwrapped coordinates

### Required Files

* This code comes with 3 files: (1) one .cpp file, (2) one .h header file, (3) parameter text file

* In the parameter file, one needs to input the information about the system and the simulation details like trajectory file name, start time, end time, interval time, number of atom types of interest, list of atom types, time unit conversion, and the COM flag.  If COMFlag == 'no', averaged mean squared displacement (MSD) of the individual beads in a chain is computed and if COMFlag == 'yes', COM MSD of the chain is computed

### Compile and run the C++ code as below:

* If using a GCC compiler:  
    gcc -lstdc++ -lm  BondAutoCorr_P2.cpp   
    ./a.out
    
* If using a ICC compiler:  
    icc BondAutoCorr_P2.cpp  
    ./a.out
    
### Other details:

This code does moving average over all the frames and also averages over all the atoms of interest selected based on the atom types given by the user in the parameter file


