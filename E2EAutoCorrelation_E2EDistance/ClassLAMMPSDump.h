//This header file reads dump files at different time steps

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
	double x, y, z;
	int ix, iy, iz;
};

class DumpCoord
{
public:
	int AtomID, MoleculeID, AtomTypeID;
	float Charge, Mass;
	Coord AtomPosition;
	Coord AtomImageFlags;

	void ReadCoord(ifstream &in)
	{
		in >> AtomID >> MoleculeID >> AtomTypeID >> Charge >> Mass >> AtomPosition.x >> AtomPosition.y >> AtomPosition.z >> AtomImageFlags.ix >> AtomImageFlags.iy >> AtomImageFlags.iz;
	}
};

//Class to read the dump files
class DumpFileData
{
public:
	int TimeStep, NumberofAtoms;
	Coord BoxMin, BoxMax;
	DumpCoord *AtomsCoord;

	void ReadfromFile(string Filename)
	{
		ifstream in(Filename.c_str());
		if (!in)
		{
			cout << "Error reading " << Filename;
			cin >> Filename;
		}

		string tempstring;
		in >> tempstring >> tempstring;
		in >> TimeStep;
		in >> tempstring >> tempstring >> tempstring >> tempstring;
		in >> NumberofAtoms;
		in >> tempstring >> tempstring >> tempstring >> tempstring >> tempstring >> tempstring;
		in >> BoxMin.x >> BoxMax.x >> BoxMin.y >> BoxMax.y >> BoxMin.z >> BoxMax.z;

		in >> tempstring >> tempstring >> tempstring >> tempstring >> tempstring >> tempstring >> tempstring >> tempstring >> tempstring >> tempstring >> tempstring >> tempstring >> tempstring;
		AtomsCoord = new DumpCoord[NumberofAtoms];
		for (int AtomsCounter = 0; AtomsCounter<NumberofAtoms; ++AtomsCounter)
		{
			AtomsCoord[AtomsCounter].ReadCoord(in);
		}
			
		in.close();
	}

	void DestroyDumpFileData()
	{
		delete[] AtomsCoord; AtomsCoord = NULL;
	}
};

//This function is useful in combining the string part of the dump file name with the timestep value
string FilenameDotTimestep(string InitialFilename, int TimeStep)
{
	/* Convert Timestep to append to filename */
	stringstream TimeStream;
	TimeStream << TimeStep;
	return(InitialFilename + '.' + TimeStream.str());
}

