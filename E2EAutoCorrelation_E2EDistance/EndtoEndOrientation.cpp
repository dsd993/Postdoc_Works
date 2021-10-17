//Written by Dinesh S. Devarajan, dated 09/15/2021

//Code to calculate the end-to-end autocorrelation between both inter chains and intra chains
//Code also calculates the end-to-end distance

#include "ClassLAMMPSDump.h"

//To calculate cos(theta) between inter chains and intra chains
#define cos(x1, y1, z1, x2, y2, z2) ((x1*x2)+(y1*y2)+(z1*z2))/(sqrt((x1*x1)+(y1*y1)+(z1*z1))*sqrt((x2*x2)+(y2*y2)+(z2*z2)))

//To calculate only the dot product between inter chains and intra chains
#define dot_prod(x1, y1, z1, x2, y2, z2) ((x1*x2)+(y1*y2)+(z1*z2))

//To calculate the end-to-end distance
#define Distance(a,b,c) sqrt((a*a)+(b*b)+(c*c))

struct SOPInfo
{
	double Value;
};

int main()
{
	//Declaring variables for reading the parameter file
	string TempStr, TrajFileName;
	int NumberofChains, NumberofBeadsPerChain, StartTime, EndTime, IntervalTime, AutoCorr_EndTime;
	
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
		
		if (TempStr == "AutoCorr_EndTime:")
			ReadConfig >> AutoCorr_EndTime;
	}
	ReadConfig.close();
	
	//Variables declaration
	int NumberofFiles;
	double xv1, yv1, zv1, xv2, yv2, zv2, E2E_Dist, OPV, TotalValue, TotalValue_DotProd;
	
	int NumberofPairs = (NumberofChains * (NumberofChains - 1)) / 2;
	int NumberofAtoms = NumberofChains * NumberofBeadsPerChain;  //Gettin the number of atoms in the system
	NumberofFiles = ((EndTime - StartTime) / IntervalTime) + 1; //Getting the number of files
	
	int m, n = 0;
	OPV = 0;
	SOPInfo *SOPV, *SOPV_DotProd, *InitValue_DotProd;
	SOPV = new SOPInfo[NumberofChains];
	SOPV_DotProd = new SOPInfo[NumberofChains];
	InitValue_DotProd = new SOPInfo[NumberofChains];
	
	for (int i = 0; i < NumberofChains; i++)
	{
		InitValue_DotProd[i].Value = 0;
	}
	
	//Reading all the trajectory files first to speed up the calculations
	DumpFileData *D1;
	D1 = new DumpFileData[NumberofFiles];
	cout << "Reading Dump Files" << endl;
	for (int FileCounter = 0; FileCounter < NumberofFiles; ++FileCounter)
	{
		D1[FileCounter].ReadfromFile(FilenameDotTimestep(TrajFileName, (StartTime + FileCounter * IntervalTime)));
		cout << "Reading dump file at time: " << (StartTime + (FileCounter * IntervalTime)) << endl;
	}

	//****************Loop to calculate the end-to-end distance************************
	
	cout << "Calculating end-to-end distance in nm as a function of time" << endl;
	
	//Output file named 'E2EDistance.txt'
	ofstream OutputResults_E2EDistance("E2EDistance.txt");
	OutputResults_E2EDistance << "Time" << ' ' << "E2EDistance in nm" << endl;
	
	for (int i = 0; i < NumberofFiles; ++i)
	{
		E2E_Dist = 0;
		for (int j = 0; j < NumberofChains; ++j)
		{
			xv1 = D1[i].AtomsCoord[NumberofBeadsPerChain * j].AtomPosition.x - D1[i].AtomsCoord[NumberofBeadsPerChain * j + (NumberofBeadsPerChain - 1)].AtomPosition.x;
			yv1 = D1[i].AtomsCoord[NumberofBeadsPerChain * j].AtomPosition.y - D1[i].AtomsCoord[NumberofBeadsPerChain * j + (NumberofBeadsPerChain - 1)].AtomPosition.y;
			zv1 = D1[i].AtomsCoord[NumberofBeadsPerChain * j].AtomPosition.z - D1[i].AtomsCoord[NumberofBeadsPerChain * j + (NumberofBeadsPerChain - 1)].AtomPosition.z;
			E2E_Dist += Distance(xv1,yv1,zv1);
		}
		//E2E_Dist /= (NumberofChains*10);
		E2E_Dist /= (NumberofChains);
		OutputResults_E2EDistance << (StartTime + (i * IntervalTime)) << ' ' << E2E_Dist << endl; //Outputs
	}

	//************************Loop to calculate the end-to-end autocorrelation for inter chains (only if there are multiple chains)**************************
	
	//Within chains at a particular time and moves along time
	if (NumberofChains != 1)
	{	
	
		cout << "Calculating end-to-end autocorrelation for inter chains" << endl;
		//Output files
		ofstream OutputResults_Inter_E2EAutoCorr("Inter_E2EAutoCorr.txt");
		OutputResults_Inter_E2EAutoCorr << "Time" << ' ' << "Inter_E2EAutoCorr" << ' ' << "Pairs" << endl;
		for (int i = 0; i < NumberofFiles; ++i)
		{
			OPV = 0;
			m = 0;
			for (int j = 0; j < NumberofChains; ++j)
			{
				xv1 = D1[i].AtomsCoord[NumberofBeadsPerChain * j].AtomPosition.x - D1[i].AtomsCoord[NumberofBeadsPerChain * j + (NumberofBeadsPerChain - 1)].AtomPosition.x;
				yv1 = D1[i].AtomsCoord[NumberofBeadsPerChain * j].AtomPosition.y - D1[i].AtomsCoord[NumberofBeadsPerChain * j + (NumberofBeadsPerChain - 1)].AtomPosition.y;
				zv1 = D1[i].AtomsCoord[NumberofBeadsPerChain * j].AtomPosition.z - D1[i].AtomsCoord[NumberofBeadsPerChain * j + (NumberofBeadsPerChain - 1)].AtomPosition.z;
					for (int k = j + 1; k < NumberofChains; ++k)
					{
						xv2 = D1[i].AtomsCoord[NumberofBeadsPerChain * k].AtomPosition.x - D1[i].AtomsCoord[NumberofBeadsPerChain * k + (NumberofBeadsPerChain - 1)].AtomPosition.x;
						yv2 = D1[i].AtomsCoord[NumberofBeadsPerChain * k].AtomPosition.y - D1[i].AtomsCoord[NumberofBeadsPerChain * k + (NumberofBeadsPerChain - 1)].AtomPosition.y;
						zv2 = D1[i].AtomsCoord[NumberofBeadsPerChain * k].AtomPosition.z - D1[i].AtomsCoord[NumberofBeadsPerChain * k + (NumberofBeadsPerChain - 1)].AtomPosition.z;
						++m;
						OPV += cos(xv1, yv1, zv1, xv2, yv2, zv2);
					}
			}
		OPV /= NumberofPairs;
		OutputResults_Inter_E2EAutoCorr << (StartTime + (i * IntervalTime)) << ' ' << OPV << ' ' << m << endl;  //Outputs
		}
	}

	//*************************Loop to calculate the end-to-end autocorrelation for intra chains******************************
	
	//Moving time average is done and also averaged over the number of chains in the system if there are multiple chains!

	//Output files
	cout << "Calculating end-to-end autocorrelation for intra chains" << endl;
	ofstream OutputResults_Intra_E2EAutoCorr("Intra_E2EAutoCorr.txt");
	OutputResults_Intra_E2EAutoCorr << "Time" << ' ' << "Intra_E2EAutoCorr_Cos" << ' ' << "Intra_E2EAutoCorr_Init" << endl;
	while(1)
	{
		if (NumberofFiles < 1)
			break;
		
		if (n*IntervalTime > AutoCorr_EndTime)
			break;
		
		for (int i = 0; i < NumberofChains; i++)
		{
			SOPV[i].Value = 0;
			SOPV_DotProd[i].Value = 0;
		}
		TotalValue = 0;
		TotalValue_DotProd = 0;
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
	
		OutputResults_Intra_E2EAutoCorr << n * IntervalTime << ' ' << TotalValue << ' ' << TotalValue_DotProd << endl; //Outputs
		n++;
		cout << ((NumberofFiles)/(double)(((EndTime - StartTime) / IntervalTime) + 1))*100 << ' ' << "% left for completion" << endl;
		NumberofFiles -= 1;
	}    
return(420);
}
