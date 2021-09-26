//This header file reads dump files at different time steps

#include<fstream>
#include<iostream>
#include<cstdio>
#include<cstring>
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

class AtomsData
{
public:
	int AtomID, MoleculeID, AtomTypeID;
	float Charge, Mass;
	Coord AtomCoord;
	Coord AtomImageFlags;

	void rd(ifstream &in)
	{
		in >> AtomID >> MoleculeID >> AtomTypeID >> Charge >> Mass >> AtomCoord.x >> AtomCoord.y >> AtomCoord.z >> AtomImageFlags.ix >> AtomImageFlags.iy >> AtomImageFlags.iz;
	}
};

//Class to read the dump files
class DumpFileData
{
public:
	int TimeStep, NumberofAtoms;
	Coord BoxMin, BoxMax;
	AtomsData *Atoms;

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
		Atoms = new AtomsData[NumberofAtoms];
		for (int AtomsCounter = 0; AtomsCounter<NumberofAtoms; ++AtomsCounter)
		{
			Atoms[AtomsCounter].rd(in);
		}

		in.close();
	}

	void DestroyDumpFileData()
	{
		delete[] Atoms; Atoms = NULL;
	}
};

string FilenameDotTimestep(string InitialFilename, int TimeStep)
{
	/* Convert Timestep to append to filename */
	stringstream TimeStream;
	TimeStream << TimeStep;
	return(InitialFilename + '.' + TimeStream.str());
}

inline double CalculateSquaredDisplacement(Coord A, Coord B)
{
	Coord SquaredDisplacement;
	SquaredDisplacement.x = (A.x - B.x)*(A.x - B.x);
	SquaredDisplacement.y = (A.y - B.y)*(A.y - B.y);
	SquaredDisplacement.z = (A.z - B.z)*(A.z - B.z);

	return((SquaredDisplacement.x + SquaredDisplacement.y + SquaredDisplacement.z));
}