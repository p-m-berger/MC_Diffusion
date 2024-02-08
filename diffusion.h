
#ifdef DIFF_MAIN
#define DIFF_EXTERN 
#else
#define DIFF_EXTERN extern
#endif

#define BUFLEN 512
#pragma warning(disable:4996)


#include "species.h"
#include "node.h"

//number of species and nodes
DIFF_EXTERN int Nspe, Nnode;

//mineral masses
DIFF_EXTERN double **Min_mol;
//mineral names
DIFF_EXTERN char **Min_name;
//mineral rates and surface areas
DIFF_EXTERN double *Min_rate, *Min_surf;


//constant node spacing in meters
DIFF_EXTERN double Dx;

//info about each species
DIFF_EXTERN species *Spe;

//holdes the diffusion coefficients for each node
DIFF_EXTERN node *Dnode;

//index of Cl - should be at the end of all of the arrays
DIFF_EXTERN int Cl;

DIFF_EXTERN int Nin, Ndom;

//is the left boundary constant, if not it is closed
DIFF_EXTERN bool Fixedbound;
//is the left boundary constant with regard to CO2, if not it is closed
DIFF_EXTERN bool FixedCO2;

//main nodes for modeling

DIFF_EXTERN char Inlet[256][512];
DIFF_EXTERN char Domain[256][512];
