### This C++ code calculates and outputs the following:

1) Radius of gyration, Rg
2) Hydrodynamic radius based on Kirkwood approximation, Rh
3) Shape factors: Asphericity (b), Acyndricity (c), and relative anisotropy (\kappa squared)
4) Eigen values
5) Principal axis vector both from gyration tensor and MOI (momnet of inertia) tensor

The eigen vector corresponding to the maximum eigen value in the case of gyration tensor or the eigen vector corresponding to the 
minimum eigen value in the case of the moment of inertia tensor was chosen to be the principle axis of a given polymer chain.

### Output files:

* In the 'Rg2Shapes_VS_Time.txt' file that the code outputs, the first column is timestep, second column is Rg squared, third column is Rg, fourth column is b, fifth column is c, sixth column is \kappa squared, seventh column is (b/Rg squared), and eight column is (c/Rg squared)

* In the 'Rh_VS_Time.txt' file that the code outputs, the first column is timestep and the second column is the Rh value

* In the 'EigenValues_GyrationTensor.txt' file that the code outputs, the first column is timestep, second column is maximum Eigen value from Gyration tensor, third column is intermediate Eigen value from Gyration tensor, and fourth column is minimum Eigen value from Gyration tensor

* In the 'EigenValues_MOITensor.txt' file that the code outputs, the first column is timestep, second column is maximum Eigen value from MOI tensor, third column is intermediate Eigen value from MOI tensor, and fourth column is minimum Eigen value from MOI tensor

* In the 'EigenVectors_GyrationTensor.txt' file that the code outputs, the first column is timestep and the second, third, and fourth column are the principal axis vector components corresponding to the maximum Eigen value obtained from the Gyration tensor.

* In the 'EigenVectors_MOITensor.txt' file that the code outputs, the first column is timestep and the second, third, and fourth column are the principal axis vector components corresponding to the minimum Eigen value obtained from the MOI tensor.

* In the 'Rg2Shapes_Averages.txt' file that the code outputs, the first and second columns are averaged Rg squared and standard deviation of Rg squared, the third and fourth columns are averaged b and standard deviation of b, the fifth and sixth columns are averaged c and standard deviation of c, and finally, the seventh and eigth columns are averaged \kappa squared and standard deviation of \kappa squared.

* In the 'AverageRh.txt' file that the code outputs, the first and second columns are averaged Rh and standard deviation of Rh.


### Requirements

* This code requires dump files named "DumpFile.TimeStep" sorted based on atom-ids and has to be in the following format: "id mol type q mass xu yu zu ix iy iz"  
* xu, yu, zu refers to unwrapped coordinates
* If there are multiple chains, each chain must have a unique molecule ID.  That is, if the system contains 500 chains, molecule ID should go from 1 to 500.  

### Required files

* This code comes with 7 files: (1) one .cpp file, (2) five .h header files, (3) one parameter text file

* In the parameter file, one needs to input the information about the system and the simulation details like trajectory file name, start time, end time, interval time, number of chains, and beads per chain.  

### Compile and run the C++ code as below

* If using a GCC compiler:  
    gcc -lstdc++ -lm  GyrationTensor_SizeShape.cpp   
    ./a.out
    
* If using a ICC compiler:  
    icc GyrationTensor_SizeShape.cpp   
    ./a.out
    
### Other details

This code averages over all the chains



