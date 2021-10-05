//Written by Dinesh S. Devarajan, dated 09/19/2021

//Code to calculate the local bond autocorrelation of each bond in a chain and P2 - second order Legendre polynomial for the overall chain

#include "ClassLAMMPSDump.h"

//To calculate cos(theta)
#define cos(x1, y1, z1, x2, y2, z2) ((x1*x2)+(y1*y2)+(z1*z2))/(sqrt((x1*x1)+(y1*y1)+(z1*z1))*sqrt((x2*x2)+(y2*y2)+(z2*z2)))

//To calculate only the dot product between inter chains and intra chains
#define dot_prod(x1, y1, z1, x2, y2, z2) ((x1*x2)+(y1*y2)+(z1*z2))

struct SOPInfo
{
	double Value;
};

int main()
{
	//Declaring variables for reading the parameter file
	string TempStr, TrajFileName;
	int NumberofChains, NumberofBeadsPerChain, StartTime, EndTime, IntervalTime, P2_EndTime;
	
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

		if (TempStr == "NumberofBeadsPerChain:")
			ReadConfig >> NumberofBeadsPerChain;

		if (TempStr == "StartTime:")
			ReadConfig >> StartTime;

		if (TempStr == "EndTime:")
			ReadConfig >> EndTime;

		if (TempStr == "IntervalTime:")
			ReadConfig >> IntervalTime;
		
		if (TempStr == "P2_EndTime:")
			ReadConfig >> P2_EndTime;
	}
	ReadConfig.close();
	
	//Variables declaration
	int NumberofFiles, id1, id2;
	double xv1, yv1, zv1, xv2, yv2, zv2, E2E_Dist, OPV, TotalValue, TotalValue_DotProd;
	
	int NumberofAtoms = NumberofChains * NumberofBeadsPerChain;  //Gettin the number of atoms in the system
	int NumberofBonds = NumberofAtoms - 1;  //Gettin the number of bonds in the system
	NumberofFiles = ((EndTime - StartTime) / IntervalTime) + 1; //Getting the number of files
	
	int n = 0;
	SOPInfo *SOPV, *SOPV_DotProd, *InitValue_DotProd;
	SOPV = new SOPInfo[NumberofChains];
	SOPV_DotProd = new SOPInfo[NumberofChains];
	InitValue_DotProd = new SOPInfo[NumberofChains];

	for (int i = 0; i < NumberofChains; i++)
	{
		InitValue_DotProd[i].Value = 0;
	}
	

	//******************Reading all the trajectory files first to speed up the calculations**********************
	DumpFileData *D1;
	D1 = new DumpFileData[NumberofFiles];
	cout << "Reading Dump Files" << endl;
	for (int FileCounter = 0; FileCounter < NumberofFiles; ++FileCounter)
	{
		D1[FileCounter].ReadfromFile(FilenameDotTimestep(TrajFileName, (StartTime + FileCounter * IntervalTime)));
		cout << "Reading dump file at time: " << (StartTime + (FileCounter * IntervalTime)) << endl;
	}

	//***************************Loop to calculate bond autocorrelation for each bonds in a chain*******************************
	
	//Moving time average is done and also averaged over the number of chains in the system!

	//Output files
	for (int id1 = 0; id1 < NumberofBonds; id1++)
	{
	id2 = id1 + 1;
	cout << "Calculating bond autocorrelation for bond between atom id " << " " << id1+1 << " " << "and atom id" << " " << id2+1 << endl;
	ofstream OutputResults1(FilenameDotTimestep("BondAutoCorr_ForBond", id2));

	n = 0;

	NumberofFiles = ((EndTime - StartTime) / IntervalTime) + 1;

	while(1)
	{
		for (int i = 0; i < NumberofChains; i++)
		{
		SOPV[i].Value = 0;
		SOPV_DotProd[i].Value = 0;
		}
	
		TotalValue = 0;	
		TotalValue_DotProd = 0;
	
		if (NumberofFiles < 1)
			break;
		
		if (n*IntervalTime > P2_EndTime) //Breaking out the loop once the bond autocorrelation is calculated until the user defined P2_EndTime
			break;

		for (int i = 0; i < NumberofFiles; i++)
		{	
			for (int j = id1; j < NumberofAtoms; j = j + NumberofBeadsPerChain)
			{
				xv1 = D1[i].AtomsCoord[j].AtomPosition.x - D1[i].AtomsCoord[j+1].AtomPosition.x;
				yv1 = D1[i].AtomsCoord[j].AtomPosition.y - D1[i].AtomsCoord[j+1].AtomPosition.y;
				zv1 = D1[i].AtomsCoord[j].AtomPosition.z - D1[i].AtomsCoord[j+1].AtomPosition.z;
				xv2 = D1[i+n].AtomsCoord[j].AtomPosition.x - D1[i+n].AtomsCoord[j+1].AtomPosition.x;
				yv2 = D1[i+n].AtomsCoord[j].AtomPosition.y - D1[i+n].AtomsCoord[j+1].AtomPosition.y;
				zv2 = D1[i+n].AtomsCoord[j].AtomPosition.z - D1[i+n].AtomsCoord[j+1].AtomPosition.z;
				SOPV[int(j/NumberofBeadsPerChain)].Value += cos(xv1, yv1, zv1, xv2, yv2, zv2);
				SOPV_DotProd[int(j/NumberofBeadsPerChain)].Value += dot_prod(xv1, yv1, zv1, xv2, yv2, zv2);
			}
		}
		for (int i = 0; i < NumberofChains; i++)
		{
			TotalValue += SOPV[i].Value;
		}
		TotalValue /=  (NumberofFiles*NumberofChains); //Dividing the accumulated value by number of files and total number of chains in the system
		
		for (int i = 0; i < NumberofChains; i++)
		{
			SOPV_DotProd[i].Value /= NumberofFiles; //Dividing each chain dot product value by number of files in the system
		}

		if (n == 0)
		{
			for (int i = 0; i < NumberofChains; i++)
			{
				InitValue_DotProd[i].Value = SOPV_DotProd[i].Value; //Storing the initial dot product value of each chain in an array
			}
		}

		for (int i = 0; i < NumberofChains; i++)
		{
			SOPV_DotProd[i].Value /= InitValue_DotProd[i].Value; //Dividing the dot product value of each chain at any given time by their respective initial values at t = 0
		}

		for (int i = 0; i < NumberofChains; i++)
		{
		        TotalValue_DotProd += SOPV_DotProd[i].Value;
		}

		TotalValue_DotProd /= NumberofChains; //Dividing the accumulated value by number of chains in the system

		OutputResults1 << n * IntervalTime << ' ' << TotalValue << ' ' << TotalValue_DotProd << endl;
		n++;
		NumberofFiles -= 1;
	}
	}
	
	//*************************Loop to calculate the P2 for whole chain******************************
	
	//Moving time average is done and also averaged over the number of chains in the system if there are multiple chains!
	//Output files
	cout << "Calculating P2 for whole chain" << endl;
	ofstream OutputResults2("P2_WholeChain.txt");
	
	n = 0;
	NumberofFiles = ((EndTime - StartTime) / IntervalTime) + 1;
	
	while(1)
	{
		if (NumberofFiles < 1)
			break;
		
		if (n*IntervalTime > P2_EndTime) //Breaking out the loop once the P2 is calculated until the user defined P2_EndTime
			break;
		
		for (int i = 0; i < NumberofChains; i++)
		{
			SOPV[i].Value = 0;
		}
		TotalValue = 0;
		for (int i = 0; i < NumberofFiles; i++)
		{
		for (int j = 0; j < NumberofAtoms; j = j + NumberofBeadsPerChain)
		{
			xv1 = D1[i].AtomsCoord[j].AtomPosition.x - D1[i].AtomsCoord[j + (NumberofBeadsPerChain - 1)].AtomPosition.x;
			yv1 = D1[i].AtomsCoord[j].AtomPosition.y - D1[i].AtomsCoord[j + (NumberofBeadsPerChain - 1)].AtomPosition.y;
			zv1 = D1[i].AtomsCoord[j].AtomPosition.z - D1[i].AtomsCoord[j + (NumberofBeadsPerChain - 1)].AtomPosition.z;
			xv2 = D1[i+n].AtomsCoord[j].AtomPosition.x - D1[i+n].AtomsCoord[j + (NumberofBeadsPerChain - 1)].AtomPosition.x;
			yv2 = D1[i+n].AtomsCoord[j].AtomPosition.y - D1[i+n].AtomsCoord[j + (NumberofBeadsPerChain - 1)].AtomPosition.y;
			zv2 = D1[i+n].AtomsCoord[j].AtomPosition.z - D1[i+n].AtomsCoord[j + (NumberofBeadsPerChain - 1)].AtomPosition.z;
			SOPV[int(j/NumberofBeadsPerChain)].Value += pow(cos(xv1, yv1, zv1, xv2, yv2, zv2),2);
		}
		}
		for (int i = 0; i < NumberofChains; i++)
		{
			TotalValue += SOPV[i].Value;
		}
		TotalValue /=  (NumberofFiles*NumberofChains); //Dividing the accumulated value by number of files and total number of chains in the system
		TotalValue = (3 * TotalValue  - 1) * 0.5; //P2 second order Legendre polynomial functional form
	
		OutputResults2 << n * IntervalTime << ' ' << TotalValue << endl; //Outputs
		n++;
		NumberofFiles -= 1;
	}        
return(420);
}
