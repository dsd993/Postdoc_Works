// Written by Dinesh S. Devarajan, dated 10/16/2021

/*
Code outputs the following:

(1) Rg**2 and Rg
(2) Rh
(3) Asphericity, b and (b/Rg**2)
(4) Acyndricity, c and (c/Rg**2)
(5) Relative anisotropy, k**2
(6) Eigen values
(7) Principal axis vector both from gyration tensor and MOI tensor

The eigen vector corresponding to the maximum eigen value in the case of gyration tensor or the eigen vector corresponding to the 
minimum eigen value in the case of the moment of inertia tensor was chosen to be the principle axis of a given polymer chain.
*/

// For definitions of the radius of gyration tensor components (S) and moment of inertia tensor components(I), refer to the folloing paper: https://doi.org/10.1021/jp2065612

//Header files
#include "ClassLAMMPSDumpFile.h"
#include "nr3.h"
#include "eigen_sym.h"
#include "eigen_unsym.h"
#include "Statistic.h"

int main()
{
	//***************Reading the parameter file*****************
	ifstream ReadConfig("ParameterFile.txt");
	double aveRg2 = 0; double aveb = 0; double avec = 0; double aveK2 = 0;
	double sdevRg2, sdevb, sdevc, sdevK2, sdevRh;
	int NumberofMolecules, NumberofAtomPerMolecule;
	int StartTime, EndTime, IntervalTime;
	string TempString, DumpFilename, COMFilename, AverageResults, InstantaneousResults, GyrationTensor;
	
	ReadConfig >> TempString>>DumpFilename>> TempString>>COMFilename>>TempString>> StartTime>>TempString>>EndTime>>TempString>>IntervalTime>>TempString>> NumberofMolecules >> TempString>> NumberofAtomPerMolecule;
	ReadConfig.close();
	// Finished reading the parameter file

	//*********************Creating all necessary Output files**************************
	ofstream WriteShapeTime ("Rg2Shapes_VS_Time.txt");
	ofstream OutAverageRG("Rg2Shapes_Averages.txt");
	ofstream OutInertiaEigenvaluesI("EigenValues_MOITensor.txt");
	ofstream OutInertiaEigenvaluesS("EigenValues_GyrationTensor.txt");
	ofstream OutInertiaEigenvectorsI("EigenVectors_MOITensor.txt");
	ofstream OutInertiaEigenvectorsS("EigenVectors_GyrationTensor.txt");
	ofstream WriteRhTime("Rh_VS_Time.txt");
	ofstream AverageRh("AverageRh.txt");

	int NumberofIntervals=int(floor((EndTime-StartTime)/IntervalTime)) + 1;
	
	//***********************************Master loop over Files*****************************************
	int IntervalCounter = 0;
	VecDoub Rg2(NumberofIntervals); VecDoub Rg2xx(NumberofIntervals); VecDoub Rg2yy(NumberofIntervals); VecDoub Rg2zz(NumberofIntervals);
	VecDoub b(NumberofIntervals); VecDoub c(NumberofIntervals); VecDoub K2(NumberofIntervals);
	VecDoub MaxEigI(NumberofIntervals); VecDoub MinEigI(NumberofIntervals); VecDoub MidEigI(NumberofIntervals);
	VecDoub MaxEigS(NumberofIntervals); VecDoub MinEigS(NumberofIntervals); VecDoub MidEigS(NumberofIntervals);
	VecDoub MaxEigVecI_x(NumberofIntervals); VecDoub MaxEigVecI_y(NumberofIntervals); VecDoub MaxEigVecI_z(NumberofIntervals);
	VecDoub MinEigVecI_x(NumberofIntervals); VecDoub MinEigVecI_y(NumberofIntervals); VecDoub MinEigVecI_z(NumberofIntervals);
	VecDoub MidEigVecI_x(NumberofIntervals); VecDoub MidEigVecI_y(NumberofIntervals); VecDoub MidEigVecI_z(NumberofIntervals);
	VecDoub MaxEigVecS_x(NumberofIntervals); VecDoub MaxEigVecS_y(NumberofIntervals); VecDoub MaxEigVecS_z(NumberofIntervals); 
	VecDoub MinEigVecS_x(NumberofIntervals); VecDoub MinEigVecS_y(NumberofIntervals); VecDoub MinEigVecS_z(NumberofIntervals); 
	VecDoub MidEigVecS_x(NumberofIntervals); VecDoub MidEigVecS_y(NumberofIntervals); VecDoub MidEigVecS_z(NumberofIntervals);
	
	double **Rg2Mol; double **bMol; double **cMol; double **K2Mol; double **RhMol;
	Rg2Mol = new double *[NumberofIntervals]; bMol = new double *[NumberofIntervals]; cMol = new double *[NumberofIntervals];
	K2Mol = new double *[NumberofIntervals]; RhMol = new double *[NumberofIntervals];
	
	for (int FileIndex =0; FileIndex < NumberofIntervals; FileIndex++)
	{
		Rg2Mol[FileIndex] = new double [NumberofMolecules]; bMol[FileIndex] = new double [NumberofMolecules];
		cMol[FileIndex] = new double [NumberofMolecules]; K2Mol[FileIndex] = new double [NumberofMolecules];
		RhMol[FileIndex] = new double [NumberofMolecules];
	}
	
	for (int FileIndex =0; FileIndex < NumberofIntervals; FileIndex++)
	{
			for (int MolIndex =0; MolIndex < NumberofMolecules; MolIndex++)
			{
				Rg2Mol[FileIndex][MolIndex] = 0.0; bMol[FileIndex][MolIndex] = 0.0; cMol[FileIndex][MolIndex] = 0.0;
				K2Mol[FileIndex][MolIndex] = 0.0; RhMol[FileIndex][MolIndex] = 0.0;
			}
	}

	MatDoub S(3,3);  // Gyration tensor
	MatDoub I(3,3);  // Moment of Inertia Tensor
	DumpFileData D;
	DumpFileData COM;
	
	float mass;
	for (int FileIndex = 0; FileIndex < NumberofIntervals;FileIndex ++)
	{
		D.ReadfromFile(FilenameDotTimestep(DumpFilename,StartTime + FileIndex*IntervalTime).c_str());
		InitializeDumpFile(D,COM); // Just initialize box and timestep info
		COM.NumberofAtoms = NumberofMolecules;
		COM.Atoms = new AtomsData[NumberofMolecules];
 		for (int MolIndex = 0; MolIndex < NumberofMolecules; MolIndex ++)
		{
			COM.Atoms[MolIndex].AtomID = MolIndex + 1;
			COM.Atoms[MolIndex].MoleculeID = MolIndex + 1;
			COM.Atoms[MolIndex].AtomTypeID = 1;
			COM.Atoms[MolIndex].Charge = 0.0;
			COM.Atoms[MolIndex].AtomCoord.x = 0.0; // Initialize
			COM.Atoms[MolIndex].AtomCoord.y = 0.0; // Initialize
			COM.Atoms[MolIndex].AtomCoord.z = 0.0; // Initialize
			mass = 0.0;
			
			for (int AtomIndex = 0; AtomIndex < NumberofAtomPerMolecule; AtomIndex++)
			{
				int Index1 = MolIndex*NumberofAtomPerMolecule + AtomIndex; // Atom Index
				COM.Atoms[MolIndex].AtomCoord.x += D.Atoms[Index1].AtomCoord.x * D.Atoms[Index1].mass;
				COM.Atoms[MolIndex].AtomCoord.y += D.Atoms[Index1].AtomCoord.y * D.Atoms[Index1].mass;
				COM.Atoms[MolIndex].AtomCoord.z += D.Atoms[Index1].AtomCoord.z * D.Atoms[Index1].mass;
				mass += D.Atoms[Index1].mass; // Accumulate mass of the atoms to calculate mass of the molecule
			}
			// Calculate the center of mass of the molecule
			COM.Atoms[MolIndex].mass = mass;
			COM.Atoms[MolIndex].AtomCoord.x /= mass;
			COM.Atoms[MolIndex].AtomCoord.y /= mass;
			COM.Atoms[MolIndex].AtomCoord.z /= mass;
		}
		// Perform Rg2 Calculation...
		
		for (int MolIndex = 0; MolIndex < NumberofMolecules; MolIndex ++)
		{
			S[0][0] = 0.0;
			S[1][0] = 0.0;
			S[2][0] = 0.0;
			S[0][1] = 0.0;
			S[1][1] = 0.0;
			S[2][1] = 0.0;
			S[0][2] = 0.0;
			S[1][2] = 0.0;
			S[2][2] = 0.0;

			I[0][0] = 0.0;
			I[1][0] = 0.0;
			I[2][0] = 0.0;
			I[0][1] = 0.0;
			I[1][1] = 0.0;
			I[2][1] = 0.0;
			I[0][2] = 0.0;
			I[1][2] = 0.0;
			I[2][2] = 0.0;
			
			//***********************Computing gyration tensor (S) and moment of inertia tensor (I) components******************************
			for (int AtomIndex = 0; AtomIndex < NumberofAtomPerMolecule; AtomIndex++)
			{
				int Index = MolIndex*NumberofAtomPerMolecule + AtomIndex; // Atom Index
				S[0][0] += ((D.Atoms[Index].AtomCoord.x-COM.Atoms[MolIndex].AtomCoord.x)*(D.Atoms[Index].AtomCoord.x-COM.Atoms[MolIndex].AtomCoord.x));
				S[1][1] += ((D.Atoms[Index].AtomCoord.y-COM.Atoms[MolIndex].AtomCoord.y)*(D.Atoms[Index].AtomCoord.y-COM.Atoms[MolIndex].AtomCoord.y));
				S[2][2] += ((D.Atoms[Index].AtomCoord.z-COM.Atoms[MolIndex].AtomCoord.z)*(D.Atoms[Index].AtomCoord.z-COM.Atoms[MolIndex].AtomCoord.z));
				S[0][1] += ((D.Atoms[Index].AtomCoord.x-COM.Atoms[MolIndex].AtomCoord.x)*(D.Atoms[Index].AtomCoord.y-COM.Atoms[MolIndex].AtomCoord.y));
				S[0][2] += ((D.Atoms[Index].AtomCoord.x-COM.Atoms[MolIndex].AtomCoord.x)*(D.Atoms[Index].AtomCoord.z-COM.Atoms[MolIndex].AtomCoord.z));
				S[1][2] += ((D.Atoms[Index].AtomCoord.y-COM.Atoms[MolIndex].AtomCoord.y)*(D.Atoms[Index].AtomCoord.z-COM.Atoms[MolIndex].AtomCoord.z));

				I[0][0] += ((D.Atoms[Index].AtomCoord.y * D.Atoms[Index].AtomCoord.y + D.Atoms[Index].AtomCoord.z * D.Atoms[Index].AtomCoord.z)*D.Atoms[Index].mass); 
				I[1][1] += ((D.Atoms[Index].AtomCoord.x * D.Atoms[Index].AtomCoord.x + D.Atoms[Index].AtomCoord.z * D.Atoms[Index].AtomCoord.z)*D.Atoms[Index].mass); 
				I[2][2] += ((D.Atoms[Index].AtomCoord.x * D.Atoms[Index].AtomCoord.x + D.Atoms[Index].AtomCoord.y * D.Atoms[Index].AtomCoord.y)*D.Atoms[Index].mass); 
				I[0][1] += -((D.Atoms[Index].AtomCoord.x * D.Atoms[Index].AtomCoord.y)*D.Atoms[Index].mass); 
				I[0][2] += -((D.Atoms[Index].AtomCoord.x * D.Atoms[Index].AtomCoord.z)*D.Atoms[Index].mass); 
				I[1][2] += -((D.Atoms[Index].AtomCoord.y * D.Atoms[Index].AtomCoord.z)*D.Atoms[Index].mass); 
			}
			
			// Normalizing the gyration tensor components wtr number of atoms per molecule
			S[0][0] = S[0][0]/NumberofAtomPerMolecule;
			S[1][1] = S[1][1]/NumberofAtomPerMolecule;
			S[2][2] = S[2][2]/NumberofAtomPerMolecule;
			S[0][1] = S[0][1]/NumberofAtomPerMolecule;
			S[0][2] = S[0][2]/NumberofAtomPerMolecule;
			S[1][2] = S[1][2]/NumberofAtomPerMolecule;
			S[1][0] = S[0][1];
			S[2][0] = S[0][2];
			S[2][1] = S[1][2];

			I[1][0] = I[0][1];
			I[2][0] = I[0][2];
			I[2][1] = I[1][2];

			Symmeig EigenProcessor(S);// Eigenvalue calculations
			Symmeig InertiaEigenProcessor(I); // Eigenvalue calculations for I
			
			// Eigen value of gyration tensor and moment of inertia tensor
			MaxEigI[FileIndex] = InertiaEigenProcessor.d[0]; MinEigI[FileIndex] = InertiaEigenProcessor.d[2]; MidEigI[FileIndex] = InertiaEigenProcessor.d[1];
			MaxEigS[FileIndex] = EigenProcessor.d[0]; MinEigS[FileIndex] = EigenProcessor.d[2]; MidEigS[FileIndex] = EigenProcessor.d[1];
			
			// Eigen vectors of gyration tensor and moment of inertia tensor
			MaxEigVecI_x[FileIndex] = InertiaEigenProcessor.z[0][0]; MaxEigVecI_y[FileIndex] = InertiaEigenProcessor.z[1][0]; MaxEigVecI_z[FileIndex] = InertiaEigenProcessor.z[2][0]; 
			MinEigVecI_x[FileIndex] = InertiaEigenProcessor.z[0][2]; MinEigVecI_y[FileIndex] = InertiaEigenProcessor.z[1][2]; MinEigVecI_z[FileIndex] = InertiaEigenProcessor.z[2][2];
			MidEigVecI_x[FileIndex] = InertiaEigenProcessor.z[0][1]; MidEigVecI_y[FileIndex] = InertiaEigenProcessor.z[1][1]; MidEigVecI_z[FileIndex] = InertiaEigenProcessor.z[2][1];
			MaxEigVecS_x[FileIndex] = EigenProcessor.z[0][0]; MaxEigVecS_y[FileIndex] = EigenProcessor.z[1][0]; MaxEigVecS_z[FileIndex] = EigenProcessor.z[2][0];
			MinEigVecS_x[FileIndex] = EigenProcessor.z[0][2]; MinEigVecS_y[FileIndex] = EigenProcessor.z[1][2]; MinEigVecS_z[FileIndex] = EigenProcessor.z[2][2];
			MidEigVecS_x[FileIndex] = EigenProcessor.z[0][1]; MidEigVecS_y[FileIndex] = EigenProcessor.z[1][1]; MidEigVecS_z[FileIndex] = EigenProcessor.z[2][1];
			
			Rg2[FileIndex] = EigenProcessor.d[0]+EigenProcessor.d[1]+EigenProcessor.d[2]; // Rg2 calculation
			b[FileIndex] = Calculate_b(EigenProcessor.d[0],EigenProcessor.d[1],EigenProcessor.d[2]); // Asphericity calculation
			c[FileIndex] = Calculate_c(EigenProcessor.d[1],EigenProcessor.d[2]); // Acylindricity calculation
			K2[FileIndex] = Calculate_K2(b[FileIndex],c[FileIndex],Rg2[FileIndex]); // Relative anisotropy calulation
			Rg2Mol[FileIndex][MolIndex] = Rg2[FileIndex];
			bMol[FileIndex][MolIndex] = b[FileIndex];
			cMol[FileIndex][MolIndex] = c[FileIndex];
			K2Mol[FileIndex][MolIndex] = K2[FileIndex];
		}
	}
	COM.DestroyDumpFileData();
    D.DestroyDumpFileData();
	for (int FileIndex = 0; FileIndex < NumberofIntervals; FileIndex++)
	{
		Rg2[FileIndex] = 0; b[FileIndex] = 0.0; c[FileIndex] = 0; K2[FileIndex] = 0;
		for (int MolIndex = 0; MolIndex < NumberofMolecules; MolIndex++)
		{
			Rg2[FileIndex] += Rg2Mol[FileIndex][MolIndex];
			b[FileIndex] += bMol[FileIndex][MolIndex];
			c[FileIndex] += cMol[FileIndex][MolIndex];
			K2[FileIndex] += K2Mol[FileIndex][MolIndex];
		}
		
		// Normalizing the quantities by number of molecules in the system
		Rg2[FileIndex]/=NumberofMolecules;
		b[FileIndex]/=NumberofMolecules;
		c[FileIndex]/=NumberofMolecules;
		K2[FileIndex]/=NumberofMolecules;
		
		//Outputting all the quantities as a function of time
		WriteShapeTime << StartTime + FileIndex * IntervalTime <<" "<<Rg2[FileIndex]<<" "<<sqrt(Rg2[FileIndex])<<" "<<b[FileIndex]<< " " <<c[FileIndex]<<" "<<K2[FileIndex]<<" "<<b[FileIndex]/Rg2[FileIndex]<<" "<<c[FileIndex]/Rg2[FileIndex]<<endl;
		OutInertiaEigenvaluesI << StartTime + FileIndex * IntervalTime << " " <<	MaxEigI[FileIndex] << " " << MidEigI[FileIndex] << " " << MinEigI[FileIndex] << endl;
		OutInertiaEigenvaluesS << StartTime + FileIndex * IntervalTime << " " << MaxEigS[FileIndex] << " " << MidEigS[FileIndex] << " " << MinEigS[FileIndex] << endl;
		OutInertiaEigenvectorsI << StartTime + FileIndex * IntervalTime << " " << MinEigVecI_x[FileIndex] << " " << MinEigVecI_y[FileIndex] << " " << MinEigVecI_z[FileIndex] << endl;
		OutInertiaEigenvectorsS << StartTime + FileIndex * IntervalTime << " " << MaxEigVecS_x[FileIndex] << " " << MaxEigVecS_y[FileIndex] << " " << MaxEigVecS_z[FileIndex] << endl;
	}
	
	// Computing the averages and the standard deviation of the quantities
	moment(Rg2,aveRg2,sdevRg2);
	moment(b,aveb,sdevb);
	moment(c,avec,sdevc);
	moment(K2,aveK2,sdevK2);
	
	OutAverageRG << aveRg2 <<" "<<sdevRg2<<" " <<aveb<<" "<<sdevb<<" "<<avec<<" "<<sdevc<<" "<<aveK2<<" "<<sdevK2;
	// Rg calculations has finished ...All OK!

	//**************************Begining Rh calculations ...based on the Kirkwood approximation ...**********************************

	// MasterLoop over files		
	double aveRH = 0;
	VecDoub RH(NumberofIntervals);
	for (int FileIndex = 0; FileIndex < NumberofIntervals;FileIndex ++)
	{
		DumpFileData D; // Dump file
		D.ReadfromFile(FilenameDotTimestep(DumpFilename,StartTime + FileIndex*IntervalTime).c_str());
		RH[IntervalCounter] = 0.0; // Initialize
		for (int MolIndex = 0; MolIndex < NumberofMolecules; MolIndex ++)
		{
			int Index1 = 0;int Index2 = 0;
			for ( int AtomIndexI = 0; AtomIndexI < NumberofAtomPerMolecule; ++AtomIndexI)
				{
					int Index1 = MolIndex*NumberofAtomPerMolecule + AtomIndexI; // Atom Index
					for (int AtomIndexJ = 0; AtomIndexJ < NumberofAtomPerMolecule; ++AtomIndexJ)
					{
						int Index2 = MolIndex*NumberofAtomPerMolecule + AtomIndexJ; // Atom Index
						if (Index1 != Index2)
						{
							RhMol[FileIndex][MolIndex] += 1/(CalculateDistance(D.Atoms[Index1].AtomCoord,D.Atoms[Index2].AtomCoord));
						}
					}
				}
			RhMol[FileIndex][MolIndex] = (NumberofAtomPerMolecule*NumberofAtomPerMolecule)/RhMol[FileIndex][MolIndex];		
		}
		 D.DestroyDumpFileData();
	}
	
	for (int FileIndex = 0; FileIndex < NumberofIntervals; FileIndex++)
	{
		RH[FileIndex] = 0.0; // Initialize
		for (int MolIndex = 0; MolIndex < NumberofMolecules; MolIndex++)
		{
			RH[FileIndex] += RhMol[FileIndex][MolIndex];
		}
		RH[FileIndex]/=NumberofMolecules;
        
		// Outputting Rh values as a function of time
		WriteRhTime << StartTime + FileIndex * IntervalTime <<" "<<RH[FileIndex]<<endl;
	}
	moment(RH,aveRH,sdevRh);
	AverageRh << aveRH<<" "<<sdevRh;
	return 9111989;
}
