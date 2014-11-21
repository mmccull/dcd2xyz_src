/*
 * 
 * 
 * 
 */
#include "xyzlib.h"


void write_xyz(FILE *xyzFile, double **coord, char** atomNames, int nAtoms)
{
	int   atom;
	int   k;
	int   junk;
	char buffer[1024];
	char posChar[10];
	
	// need to print the per step xyz header info (consists of two lines one with nAtoms the other is a title.  I forget which is which so I print nAtoms on both)
	fprintf(xyzFile,"%12d\n", nAtoms);
	fprintf(xyzFile,"%12d\n", nAtoms);
	for(atom=0;atom<nAtoms;atom++) {
		fprintf(xyzFile,"%5s %12.6f%12.6f%12.6f\n",atomNames[atom],coord[atom][0],coord[atom][1],coord[atom][2]);
	}

}

