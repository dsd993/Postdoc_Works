//Written by Dinesh S. Devarajan, dated 10/12/2021

//Code to calculate the averaged mean squared displacement (MSD) of the individual beads in chains and the COM MSD of the chains

//Code also outputs the non-Gaussian parameter (alpha2) values 

#include "ClassLAMMPSDumpFile.h"

struct COMInfo
{
	double COMx, COMy, COMz, mass;
};

int main()
{
	string TempStr, TrajFileName, COMFlag;
	int NumberofChains, BeadsPerChain, NumberofInnerMonomers, StartTime, EndTime, IntervalTime, NumberofAtomTypes, *ListofAtomTypes;
	double TimeUnitConversion, MSD, MSD2;
	
	//****************Reading parameter file********************
	ifstream ReadConfig("ParameterFile.txt");
	cout<<"Reading the Parameter File... \n \n";

	if(!ReadConfig) {cout<<"Error..Parameter file does not exist \n";exit(-1);}

	//Reading the parameter file 
	while (!ReadConfig.eof())
	{
		ReadConfig >> TempStr;
		
		if (TempStr == "TrajFileName:")
			ReadConfig >> TrajFileName;
		
		if (TempStr == "NumberofChains:")
		    ReadConfig >> NumberofChains;
		
		if (TempStr == "BeadsPerChain:")
		    ReadConfig >> BeadsPerChain;
		
		if (TempStr == "NumberofInnerMonomers:")
		    ReadConfig >> NumberofInnerMonomers;
			
		if (TempStr == "Start:")
		    ReadConfig >> StartTime;
		
		if (TempStr == "Stop:")
		    ReadConfig >> EndTime;
		
		if (TempStr == "Interval:")
		    ReadConfig >> IntervalTime;
			
		if (TempStr == "NoofAtomTypes:")
		{
		    ReadConfig >> NumberofAtomTypes;
		    ListofAtomTypes = new int[NumberofAtomTypes];
		}
			
		if (TempStr == "AtomTypeList:")
		{
		    for (int AtomTypeCounter = 0; AtomTypeCounter < NumberofAtomTypes; ++AtomTypeCounter)
		        ReadConfig >> ListofAtomTypes[AtomTypeCounter];
		}
			
		if (TempStr == "TimeUnitConversion:")
			ReadConfig >> TimeUnitConversion;
		
		if (TempStr == "COMFlag(yes/no):")
			ReadConfig >> COMFlag;
	}
	ReadConfig.close();
	
	int NumberofIntervals = ((EndTime - StartTime) / IntervalTime) + 1;
	
	//****************Counting the number of atoms of interest based on the atom types given by the user********************
	DumpFileData AtomsofInterest;
	AtomsofInterest.ReadfromFile(FilenameDotTimestep(TrajFileName, StartTime));
	
	int ChainID = -1;
	int NumberofAtomsofInterest = 0;
	if (COMFlag == "no") //Picking inner monomers based on atom type id and atom id
	{
		for (int AtomCounter = 0; AtomCounter < AtomsofInterest.NumberofAtoms; ++AtomCounter)
		{
			for (int AtomTypeCounter = 0; AtomTypeCounter < NumberofAtomTypes; ++AtomTypeCounter)
			{
				if ((ListofAtomTypes[AtomTypeCounter] == AtomsofInterest.Atoms[AtomCounter].AtomTypeID) && (AtomsofInterest.Atoms[AtomCounter].AtomID >= (((ChainID*BeadsPerChain)+floor(BeadsPerChain/2))-(NumberofInnerMonomers/2))) && (AtomsofInterest.Atoms[AtomCounter].AtomID <= (((ChainID*BeadsPerChain)+floor(BeadsPerChain/2))+(NumberofInnerMonomers/2))))
				{
					NumberofAtomsofInterest++; //Counting the atoms of interest
				}
			}
			
			if (AtomCounter%BeadsPerChain == 0) {ChainID += 1;}
		}
	}
	else //Picking all the monomers based on atom type id for COM MSD calculation
	{
		for (int AtomCounter = 0; AtomCounter < AtomsofInterest.NumberofAtoms; ++AtomCounter)
		{
			for (int AtomTypeCounter = 0; AtomTypeCounter < NumberofAtomTypes; ++AtomTypeCounter)
			{
				if (ListofAtomTypes[AtomTypeCounter] == AtomsofInterest.Atoms[AtomCounter].AtomTypeID)
				{
					NumberofAtomsofInterest++; //Counting the atoms of interest
				}
			}

		}
	}

	//****************Selecting the atoms of interest based on the atom types given by the user********************
	
	//Selected atom ids are printed to a file named 'AtomsSelected.txt'
        int *IndexAtomsofInterest; IndexAtomsofInterest = new int[NumberofAtomsofInterest];
	
	ofstream Print("AtomsSelected.txt");
	Print << "Following Atoms Selected for MSD Calculations" << endl;
	
	ChainID = -1;
	int AtomsofInterestCounter = 0;
        if (COMFlag == "no") //Picking inner monomers based on atom type id and atom id
	{
		for (int AtomCounter = 0; AtomCounter < AtomsofInterest.NumberofAtoms; ++AtomCounter)
		{
			for (int AtomTypeCounter = 0; AtomTypeCounter < NumberofAtomTypes; ++AtomTypeCounter)
			{

				if ((ListofAtomTypes[AtomTypeCounter] == AtomsofInterest.Atoms[AtomCounter].AtomTypeID) && (AtomsofInterest.Atoms[AtomCounter].AtomID >= (((ChainID*BeadsPerChain)+floor(BeadsPerChain/2))-(NumberofInnerMonomers/2))) && (AtomsofInterest.Atoms[AtomCounter].AtomID <= (((ChainID*BeadsPerChain)+floor(BeadsPerChain/2))+(NumberofInnerMonomers/2))))
				{
				    IndexAtomsofInterest[AtomsofInterestCounter] = AtomCounter;
					Print << AtomCounter <<' '<< AtomsofInterest.Atoms[AtomCounter].AtomID << endl; //Outputting the selected atom ids
					++AtomsofInterestCounter;
				}
			}
			
			if (AtomCounter%BeadsPerChain == 0) {ChainID += 1;}
		}
	}
	else //Picking all the monomers based on atom type id for COM MSD calculation
	{
		for (int AtomCounter = 0; AtomCounter < AtomsofInterest.NumberofAtoms; ++AtomCounter)
		{
			for (int AtomTypeCounter = 0; AtomTypeCounter < NumberofAtomTypes; ++AtomTypeCounter)
			{
				if (ListofAtomTypes[AtomTypeCounter] == AtomsofInterest.Atoms[AtomCounter].AtomTypeID)
				{
				    IndexAtomsofInterest[AtomsofInterestCounter] = AtomCounter;
					Print << AtomCounter <<' '<< AtomsofInterest.Atoms[AtomCounter].AtomID << endl; //Outputting the selected atom ids
					++AtomsofInterestCounter;
				}
			}
		}
	}
				
	NumberofAtomsofInterest = AtomsofInterestCounter;
	AtomsofInterest.DestroyDumpFileData();
	
	//******************Reading all the trajectory files first to speed up the calculations**********************
	DumpFileData *D1;
	D1 = new DumpFileData[NumberofIntervals];
	cout << "Reading Dump Files" << endl;
	for (int FileCounter = 0; FileCounter < NumberofIntervals; ++FileCounter)
	{
		D1[FileCounter].ReadfromFile(FilenameDotTimestep(TrajFileName, (StartTime + FileCounter * IntervalTime)));
		cout << "Reading dump file at time: " << (StartTime + (FileCounter * IntervalTime)) << endl;
	}
	
	//******************Calculating the center of mass of atoms of interest in each frame**********************
	
	//COM values at each time frame is printed to a file named 'COM.txt'
	ofstream COMPrint("COM.txt");

	COMInfo *COM;
	COM = new COMInfo[NumberofIntervals];
	
	for (int FileCounter = 0; FileCounter < NumberofIntervals; ++FileCounter)
	{
		double d1 = 0, d2 = 0, d3 = 0, m = 0;
		for (int AI = 0; AI < NumberofAtomsofInterest; ++AI)
		{
			d1 += D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].AtomCoord.x * D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].Mass;
			d2 += D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].AtomCoord.y * D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].Mass;
			d3 += D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].AtomCoord.z * D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].Mass;
			m += D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].Mass;
		}
		COM[FileCounter].COMx = d1/m; //COM in X
		COM[FileCounter].COMy = d2/m; //COM in Y
		COM[FileCounter].COMz = d3/m; //COM in Z
		COMPrint << FileCounter*IntervalTime << " " << COM[FileCounter].COMx << " " << COM[FileCounter].COMy << " " << COM[FileCounter].COMz << endl; //Outputting the COM values
	}
	
	//******************Calculating the center of mass of all the chains in each frame**********************
	
	//COM values of all the chains at each time frame is printed to a file named 'COM_EachChain.txt'
	ofstream COMEachChainPrint("COM_EachChain.txt");
	
	COMInfo **COM_EachChain;
	COM_EachChain = new COMInfo *[NumberofIntervals];
	
	for (int FileCounter = 0; FileCounter < NumberofIntervals; ++FileCounter)
	{
		COM_EachChain[FileCounter] = new COMInfo[NumberofChains];
	}
	
	for (int FileCounter = 0; FileCounter < NumberofIntervals; ++FileCounter)
	{
		for (int ChainCounter = 0; ChainCounter < NumberofChains; ++ChainCounter)
		{
			double d1 = 0, d2 = 0, d3 = 0, m = 0;
			for (int AI = 0; AI < NumberofAtomsofInterest; ++AI)
			{
				if (D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].MoleculeID == ChainCounter + 1)
				{
					d1 += D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].AtomCoord.x * D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].Mass;
					d2 += D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].AtomCoord.y * D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].Mass;
					d3 += D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].AtomCoord.z * D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].Mass;
					m += D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].Mass;
				}
			}
			COM_EachChain[FileCounter][ChainCounter].COMx = d1/m; //COM of each chain in X
			COM_EachChain[FileCounter][ChainCounter].COMy = d2/m; //COM of each chain in Y
			COM_EachChain[FileCounter][ChainCounter].COMz = d3/m; //COM of each chain in Z
			COMEachChainPrint << FileCounter*IntervalTime << " " << COM_EachChain[FileCounter][ChainCounter].COMx << " " << COM_EachChain[FileCounter][ChainCounter].COMy << " " <<  COM_EachChain[FileCounter][ChainCounter].COMz << endl; //Outputting the COM values of each chain at each timestep
		}
	}
	
	ofstream FinalOutput("MSDValues.txt"); //Final MSD values are outputted to a file named 'MSDValues.txt'
	
	if (COMFlag == "no") //Checking for the COM flag
	{
		//If COMFlag == 'no', averaged mean squared displacement (MSD) of the individual beads in chains is computed
		int n = 1;
		NumberofIntervals -= 1;
		while (1)
		{
			if (NumberofIntervals < 1)
				break;
			
			MSD = 0;
			MSD2 = 0;
			for (int FileCounter = 0; FileCounter < NumberofIntervals; FileCounter++)
			{
				for (int AI = 0; AI < NumberofAtomsofInterest; AI++)
				{
					MSD += CalculateSquaredDisplacement(D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].AtomCoord, D1[FileCounter+n].Atoms[IndexAtomsofInterest[AI]].AtomCoord);
					MSD2 += pow(CalculateSquaredDisplacement(D1[FileCounter].Atoms[IndexAtomsofInterest[AI]].AtomCoord, D1[FileCounter+n].Atoms[IndexAtomsofInterest[AI]].AtomCoord),2);
				}
				
			}
			
			MSD /= (NumberofIntervals*NumberofAtomsofInterest);
			MSD2 /= (NumberofIntervals*NumberofAtomsofInterest);
			
			FinalOutput << n*IntervalTime << " " << n*IntervalTime*TimeUnitConversion << " " << MSD << " " << MSD2 << " " << ((3*MSD2)/(5*MSD*MSD))-1 << endl;
			n++;
			NumberofIntervals -= 1;
		}
	}
	else
	{
		//If COMFlag == 'yes', COM MSD of the chains is computed
		int n = 1;
		NumberofIntervals -= 1;
		while (1)
		{
			if (NumberofIntervals < 1)
				break;
			
			MSD = 0;
			MSD2 = 0;
			
			for (int FileCounter = 0; FileCounter < NumberofIntervals; FileCounter++)
			{
				for (int ChainCounter = 0; ChainCounter < NumberofChains; ++ChainCounter)
				{
					MSD += pow((COM_EachChain[FileCounter][ChainCounter].COMx - COM_EachChain[FileCounter+n][ChainCounter].COMx),2) + pow((COM_EachChain[FileCounter][ChainCounter].COMy - COM_EachChain[FileCounter+n][ChainCounter].COMy),2) + pow((COM_EachChain[FileCounter][ChainCounter].COMz - COM_EachChain[FileCounter+n][ChainCounter].COMz),2);
					MSD2 += pow((pow((COM_EachChain[FileCounter][ChainCounter].COMx - COM_EachChain[FileCounter+n][ChainCounter].COMx),2) + pow((COM_EachChain[FileCounter][ChainCounter].COMy - COM_EachChain[FileCounter+n][ChainCounter].COMy),2) + pow((COM_EachChain[FileCounter][ChainCounter].COMz - COM_EachChain[FileCounter+n][ChainCounter].COMz),2)),2);
			    }
			}
			
			MSD /= (NumberofIntervals*NumberofChains);
			MSD2 /= (NumberofIntervals*NumberofChains);
			
			FinalOutput << n*IntervalTime << " " << n*IntervalTime*TimeUnitConversion << " " << MSD << " " << MSD2 << " " << ((3*MSD2)/(5*MSD*MSD))-1 << endl;
			n++;
			NumberofIntervals -= 1;
		}
	}
return(420);
}
