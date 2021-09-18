### This C++ code calculates the following:

1) End-to-end distance as a function of time
2) End-to-end autocorrelation for inter-chains as a function of time (only when there are multiple chains in the system)
3) End-to-end autocorrelation for intra-chains as a function of time

Compile and run the C++ code as below:

* If using a GCC compuler:
    gcc -lstdc++ -lm  EndtoEndOrientation.cpp
    ./a.out
    
* If using a ICC compiler:
    icc EndtoEndOrientation.cpp
    ./a.out
