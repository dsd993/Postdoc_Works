### This C++ code calculates the following:

1) Bond autocorrelation function for all the bonds of a chain in a file named "BondAutoCorr_ForBond.BOND-ID"
2) P2 (second order Legendre polynomial) to quantify the orientational dynamics in a file named "P2_WholeChain.txt"

Note that the "BondAutoCorr_ForBond.BONDID" will have 3 columns: The first column is the timestep, the second column is just cos(theta), and the third column is the dot product of the bond vectors divided by the initial dot product value at time t = 0.

### Requirements

* This code requires dump files named "DumpFile.TimeStep" sorted based on atom-ids and has to be in the following format: "id mol type q mass xu yu zu ix iy iz"  
* xu, yu, zu refers to unwrapped coordinates

### Required files

* This code comes with 3 files: (1) one .cpp file, (2) one .h header file, (3) one parameter text file

* In the parameter file, one needs to input the information about the system and the simulation details like trajectory file name, number of chains, number of beads per chain, start time, end time, interval time, P2 end time.  P2 end time is nothing but let's say one has trajectory information for 10<sup>8</sup> steps but knows that the correlation goes to zero in 10<sup>6</sup> steps, one can then specify P2 end time as 10<sup>6</sup> and the code will terminate once it has calculated the bond autocorrelation/P2 value until 10<sup>6</sup> steps

### Compile and run the C++ code as below

* If using a GCC compiler:  
    gcc -lstdc++ -lm  BondAutoCorr_P2.cpp   
    ./a.out
    
* If using a ICC compiler:  
    icc BondAutoCorr_P2.cpp  
    ./a.out
    
### Other details

This code does moving average over all the frames and also averages over chains if there are more than 1 chain in the system

