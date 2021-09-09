### This FORTRAN code generates sequences of E-K residues based on SCD optimization. 

* Added getSCD() function to Young's FORTRAN code to generate sequences of E-K residues based on SCD optimization.
* This code has the flexibility to generate sequences between a user preferred minimum and maximum SCD, if needed. I have added two input parameters 'minSCD' and 'maxSCD' in the 'input.in' file. 
 --> If the user does not have a preferred choice of SCD range within which the sequences has to be generated, specifiy a number greater than 10000 for 'minSCD' and 'maxSCD' in the 'input.in' file. 

#### Required files:

1) input.in: Contains input parameters to generate SV series
2) q.dat: containes the charge information for E-K residues.

#### Compile and run the FORTRAN code as below:

gfortran mtfort90.f90 getSequence.f90 

./a.out input.in > seq_scd.dat

* Note that seq_scd.dat is the output file to which the generated sequences and their corresponding SCD values are written. The file name is user's choice.
* One can make appropriate changes in the input.in and q.dat files to generates sequences of other different residues as well.
