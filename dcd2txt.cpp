
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include "dcdlib.h"
#include "psflib.h"
#include "xyzlib.h"
#include "stringlib.h"

using namespace std;

void shift_pos_positive(double **, int, int, double *); 
void dump_pos(double **, int, int, char *);
void matmul(double **, int, int, double **, int, int, double **);
double det_3_3(double **);
void read_cfg_data(char *,char *, char *, char *,char *,int *); 
void read_par_file(char **, int, char *, double *, double *,int *, int);

int main (char* argv[]) {
	
	FILE *psfFile;
	FILE *dcdCoorFile;
	FILE *dcdForceFile;
	FILE *coorOutFile;
	FILE *forceOutFile;
	char psfFileName[1024];
	char dcdCoorFileName[1024];
	char dcdForceFileName[1024];
	char forceOutFileName[1024];
	char coorOutFileName[1024];
	char buffer[1024];

	int i,j,k;

	double **coord;
	double **force;
	double *axis;
	double *charge;
	char **atomNames;
	int box;

	int nAtoms;
	int nAtomTypes;
	int nSteps;
	int step;
	long int dcdCoorFilePosition;
	long int dcdForceFilePosition;


	time_t currentTime;
	time_t previousTime;

	dcdCoorFileName[0] = '\0';
	dcdForceFileName[0] = '\0';

	// read cfg file
	read_cfg_data(psfFileName,dcdCoorFileName,dcdForceFileName,coorOutFileName,forceOutFileName,&box);

	// read in number of atoms
	psfFile = fopen(psfFileName,"r");
	read_psf_header(psfFile,&nAtoms);
	printf("Number of atoms in psf file= %d\n",nAtoms);
	//now can allocate arrays
	charge = new double [nAtoms];
	atomNames = new char* [nAtoms];
	for (i=0;i<nAtoms;i++) {
		atomNames[i] = new char [5];
	}
	//read atom type data from psf file
	read_psf_data(psfFile, nAtoms, atomNames, charge);
	//close psf file
	fclose(psfFile);

	//allocate GB arrays
	coord = new double* [nAtoms];
	force = new double* [nAtoms];
	for(i=0;i<nAtoms;i++) {
		coord[i] = new double [3];
		force[i] = new double [3];
	}


	//open output files
	forceOutFile = fopen(forceOutFileName,"w");
	coorOutFile = fopen(coorOutFileName,"w");

	// read heard for dcd coor file
	dcdCoorFile = fopen(dcdCoorFileName,"r");
	read_dcd_header(dcdCoorFile,&nAtoms,&nSteps);	
	printf("Number of atoms in coor dcd file= %d\n",nAtoms);
	printf("Number of steps in coor dcd file= %d\n",nSteps);
	// read heard for dcd force file
	dcdForceFile = fopen(dcdForceFileName,"r");
	read_dcd_header(dcdForceFile,&nAtoms,&nSteps);	
	printf("Number of atoms in force dcd file= %d\n",nAtoms);
	printf("Number of steps in force dcd file= %d\n",nSteps);
	
	// allocate arrays
	axis = new double [3];
	
	previousTime = clock();
	for (step=0;step<nSteps;step++) {
	
		//read in positions
		dcdCoorFilePosition = ftell(dcdCoorFile);
		read_dcd_step(dcdCoorFile,&dcdCoorFilePosition,axis,coord,nAtoms,box);
		//read in force
		dcdForceFilePosition = ftell(dcdForceFile);
		read_dcd_step(dcdForceFile,&dcdForceFilePosition,axis,force,nAtoms,0);

		//write coordinate xyz file
		write_xyz(coorOutFile,coord,atomNames,nAtoms);
		//write force xyz file
		write_xyz(forceOutFile,force,atomNames,nAtoms);


		/*	
		if(step%1==0) {
			printf("Reading step %d from trajectory file\n",step);
			// save time	
			currentTime = clock();
			printf("Step %4d to %4d took %f\n",step+1-1000, step+1,( (double) (currentTime-previousTime)) / CLOCKS_PER_SEC);
			previousTime = currentTime;
		}
		*/

	
	}
	fclose(dcdCoorFile);
	fclose(dcdForceFile);

	fclose(coorOutFile);
	fclose(forceOutFile);


}

/*
 *                SUBROUTINES
 */


void read_cfg_data(char *psfFileName,char *dcdCoorFileName, char *dcdForceFileName, char *coorOutFileName,char *forceOutFileName, int *box) {

	double C = 0.000395394; // Angstroms^2*mole/(K*L) combined constants in Debye screening length term

	char buffer[1024];
	char tempBuffer[1024];
	char check[15];
	char *firstWord;


	while (fgets(buffer,1024,stdin) != NULL) {

		strncpy(tempBuffer,buffer,1024);
		firstWord=string_firstword(tempBuffer);
//		printf ("First word = %s\n",firstWord);
		if (strncmp(firstWord,"psfFile",7)==0) {
			strcpy(psfFileName,string_secondword(buffer));
		} else if (strncmp(firstWord,"box",3)==0) {
			*box = atoi(string_secondword(buffer));
		} else if (strncmp(firstWord,"dcdCoorFile",11)==0) {
			strcpy(dcdCoorFileName,string_secondword(buffer));
		} else if (strncmp(firstWord,"dcdForceFile",12)==0) {
			strcpy(dcdForceFileName,string_secondword(buffer));
		} else if (strncmp(firstWord,"coorOutFile",11)==0) {
			strcpy(coorOutFileName,string_secondword(buffer));
		} else if (strncmp(firstWord,"forceOutFile",12)==0) {
			strcpy(forceOutFileName,string_secondword(buffer));
		}
	

	}	

	printf("psf file: %s\n",psfFileName);
	printf("dcd coor file: %s\n",dcdCoorFileName);
	printf("dcd force file: %s\n",dcdForceFileName);
	printf("coor out: %s\n",coorOutFileName);
	printf("force out: %s\n",forceOutFileName);


}


void dump_pos(double **coord, int nAtoms, int nSteps, char *fileName) {

	int i, j, k;
	int atom;
	int step;
	FILE *xyzOut;

	xyzOut = fopen(fileName,"w");

	for (step=0;step<nSteps;step++) {

		for (atom=0;atom<nAtoms;atom++) {

			fprintf(xyzOut,"%12.6f%12.6f%12.6f\n",coord[atom+step*nAtoms][0],coord[atom+step*nAtoms][1],coord[atom+step*nAtoms][2]);

		}

	}

	fclose(xyzOut);

}

double det_3_3(double **mat) {

	double det;

	det = mat[0][0]*mat[1][1]*mat[2][2] - mat[0][0]*mat[1][2]*mat[2][1] - mat[0][1]*mat[1][0]*mat[2][2] + mat[0][1]*mat[1][2]*mat[2][0];
	det+= mat[0][2]*mat[1][0]*mat[2][1] - mat[0][2]*mat[1][1]*mat[2][0];

	return det;
}

void matmul(double **mat1, int n, int m, double **mat2, int p, int q, double **mat3) {

	int i,j,k;


	if (m!=p) {
		printf ("Error in matmul: dimensions do not match (m=%d,p=%d)\n",m,p);
		exit(1);
	}

	for(i=0;i<n;i++) {
		for(j=0;j<q;j++) {
			mat3[i][j]=0;
			for(k=0;k<m;k++) {
				mat3[i][j]+=mat1[i][k]*mat2[k][j];
			}
		}
	}

}


	// shift reference frame to be in quadrant 1 for hilbert curve
void shift_pos_positive(double **pos,int nRows, int nCols, double *max) {

	int i, j;
	double min[nCols];

	for (i=0;i<nRows;i++) {
		for (j=0;j<nCols;j++) {
			if (i==0) {
				min[j] = pos[i][j];
				max[j] = pos[i][j];
			} else if (pos[i][j] < min[j]) {
				min[j] = pos[i][j];
			} else if (pos[i][j] > max[j]) {
				max[j] = pos[i][j];
			}
		}
	}
	for (i=0;i<nRows;i++) {
		for (j=0;j<nCols;j++) {
			pos[i][j] -= min[j];
		}
	}

	for(j=0;j<nCols;j++) {
		max[j] -= min[j];
	}
}




