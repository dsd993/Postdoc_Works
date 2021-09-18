### This C++ code calculates the following:

1) End-to-end distance as a function of time
2) End-to-end autocorrelation for inter-chains as a function of time (only when there are multiple chains in the system)
3) End-to-end autocorrelation for intra-chains as a function of time

### Required Files

* This code comes with 3 files: (1)one .cpp file, (2) one .h header file, (3) parameter text file.  

* In the parameter file, one needs to input the information about the system and the simulation like trajectory file name, number of chains, number of beads per chain, start time, end time, interval time, autocorrelation end time.  Autocorrelation end time is nothing but lets say you have trajectory information for 10^8 steps but you know that the correlation goes to zero in 10^6 steps, you can specify autocorrelation end time as 10^6 and the code will terminate once it has calculated the autocorrelation value until 10^6.

### Compile and run the C++ code as below:

* If using a GCC compuler:
    gcc -lstdc++ -lm  EndtoEndOrientation.cpp
    ./a.out
    
* If using a ICC compiler:
    icc EndtoEndOrientation.cpp
    ./a.out
