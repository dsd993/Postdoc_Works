// This header file reads dump files at different time steps

#include<fstream>
#include<iostream>
#include<cstdio>
#include<string>
#include<cstdlib>
#include<iomanip>
#include<limits>
#include<math.h>
#include<sstream>

using namespace std;

class Coord
{
public:
	double x,y,z,ix,iy,iz;
};

inline double CalculateDistance(Coord &A,Coord &B)
{
	double distance = (A.x - B.x)*(A.x - B.x) + (A.y - B.y)*(A.y - B.y) + (A.z - B.z)*(A.z - B.z);
	return sqrt(distance);
}

class AtomsData
{
public:	
	int AtomID,MoleculeID,AtomTypeID;
	float Charge,mass;
	Coord AtomCoord;
	Coord AtomImageFlag;
		
	void rd(ifstream &in)
	{
		in>>AtomID>>MoleculeID>>AtomTypeID>>Charge>>mass>>AtomCoord.x>>AtomCoord.y>>AtomCoord.z>>AtomImageFlag.ix>>AtomImageFlag.iy>>AtomImageFlag.iz;
	}
};

// Class to read the dump files
class DumpFileData
{
public:
	int TimeStep, NumberofAtoms;
	Coord BoxMin,BoxMax;
	AtomsData *Atoms;

	void ReadfromFile(string Filename)
	{
		ifstream in(Filename.c_str());
		if(!in)
		{
			cout<<"Error reading "<<Filename;
			cin>>Filename;
		}

		string tempstring;
		in>>tempstring>>tempstring;
		in>>TimeStep;
		in>>tempstring>>tempstring>>tempstring>>tempstring;
		in>>NumberofAtoms;
		in>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring;
		in>>BoxMin.x>>BoxMax.x>>BoxMin.y>>BoxMax.y>>BoxMin.z>>BoxMax.z;
		in>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring;
		Atoms = new AtomsData[NumberofAtoms];
		for(int AtomsCounter=0;AtomsCounter<NumberofAtoms;++AtomsCounter)
		{
			Atoms[AtomsCounter].rd(in);
		}

		in.close();
	}
			
	void DestroyDumpFileData()
	{
		delete [] Atoms; Atoms=NULL;
	}
};

// This function is useful in combining the string part of the dump file name with the timestep value
string FilenameDotTimestep(string InitialFilename, int TimeStep)
{
	/* Convert Timestep to append to filename */
	stringstream TimeStream;
	TimeStream<<TimeStep;
	return(InitialFilename+'.'+TimeStream.str());
}

void InitializeDumpFile(DumpFileData &File1, DumpFileData &File2)
{
	File2.BoxMax.x = File1.BoxMax.x;
	File2.BoxMax.y = File1.BoxMax.y;
	File2.BoxMax.z = File1.BoxMax.z;
	File2.BoxMin.x = File1.BoxMin.x;
	File2.BoxMin.y = File1.BoxMin.y;
	File2.BoxMin.z = File1.BoxMin.z;
	File2.TimeStep = File1.TimeStep;
}

// Function to compute asphericity based on the Eigen values
double Calculate_b(double &a,double &b,double &c)
{
	return a-0.5*(b+c);
}

// Function to compute acylindricity based on the Eigen values
double Calculate_c(double &a,double &b)
{
	return (a-b);
}

// Function to compute relative anisotropy based on the Eigen values and Rg**2
double Calculate_K2 (double &b,double &c, double &Rg2)
{
	return (b*b + (3/4)*c*c)/(Rg2*Rg2);
}