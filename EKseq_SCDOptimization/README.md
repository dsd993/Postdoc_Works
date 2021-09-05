#This FORTRAN code generates sequences of E-K residues based on SCD optimization. 

#Required files:

1) input.in: contains input paramters to generate SV series
2) q.dat: containes the charge information for E-K residues.

#Compile and execute the FORTRAN code as below

gfortran mtfort90.f90 getSequence.f90 

./a.out input.in > seq_scd.dat

#Note that seq_scd.dat is the output file to which the generated sequences and their corresponding SCD values are written. The file name is user's choice.
