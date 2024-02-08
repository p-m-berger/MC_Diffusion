#include "diffusion.h"
#include "diffusion_proto.h"

#include <fstream>
#include <stdio.h>
#include <math.h>
#include <string.h>
using namespace std;

#define NTEMP 7
#define NCLASS 6
double TempC[NTEMP] = {25, 60, 100, 150, 200, 250, 300};

double AT[NCLASS][NTEMP] = {{0, 35.43, 93.76, 153.6, 217.71, 278.00, 340.12},
{0, 7.47, 19.17, 44.34, 55.31, 68.91, 83.14},
{0, -48.83, -142.25, -244.77, -334.43, -429.2, -522.76},
{0, 21.74, 44.79, 67.27, 93.5, 121.23, 154.10},
{0, 36.63, 77.88, 119.01, 164.41, 210.3, 250.65},
{0, -15.00, -33.7, -50.97, -72.85, -93.89, -115.00}};
double BT[NCLASS][NTEMP] = {{1, 1.15, 1.38, 1.45, 1.75, 1.98, 2.21},
{1, 1.66, 2.76, 3.57, 4.80, 5.91, 7.03},
{1, 2.49, 5.16, 7.93, 10.65, 13.37, 16.10},
{1, 1.78, 2.17, 3.03, 3.91, 4.71, 5.71},
{1, 1.24, 1.65, 2.19, 2.57, 3.03, 3.51},
{1, 2.12, 3.48, 4.77, 6.39, 7.95, 9.52}};

//gas constant in J/mol/K
#define GCONST 8.314472
//faraday constant in j/equiv
#define FARAD 96485.3399

void node::init(int nsp, int ind, bool act) {
	int i;
	nspe = nsp;
	index = ind;
	d1 = new double *[nspe];
	h	= new double [nspe];
	dt = new double [nspe];
	spe = new double [nspe];
   bas = new double[nspe];
	delts = new double[nspe];
	wmass = 1;
	for (i=0; i<nspe; i++) {
		d1[i] = new double [nspe];
		delts[i] = 0;
      bas[i] = -1;
	}
	calc_act = act;
	//initialize the calculation node
	cp = new ChemPlugin("", "-d \"C:/Program Files/Gwb/Gtdata/thermo_ymp.R2.tdat\"");
}

node::~node() {
	int i;
	for (i=0; i<nspe; i++) {
		delete [] d1[i];
	}
	delete [] d1;
	delete [] h;
	delete [] spe;
   delete [] bas;
	delete [] delts;
	delete cp;
}

//calcualtes the salt fraction then calls the functions to calculate
//fluid density and viscosity
void node::calc_props() {
	int i;
	double s, na=0.0;
	for (i=0; i<nspe; i++)
		if (Spe[i].name == "Na+")
			na = spe[i];
	s = (35.453*spe[Cl] + 22.99*na)/1000000000;

	calc_density(s);
	calc_viscos(s);
}


//calculates the density of brine in a node from the equations of
//Batzle and Wang, 1992, Seismic properties of pore fluids, Geophysics
//v. 57, n. 11, p 1396-1408.
void node::calc_density(double s) {
	density = 1 + 1e-6*(-80.0*temp - 3.3*temp*temp + 0.00175*temp*temp*temp 
		+ 489.0*press - 2.0*temp*press + 0.016*temp*temp*press - 
		1.3e-5*temp*temp*temp*press - 0.333*press*press - 0.002*temp*press*press);
	density += s*(0.668 + 0.44*s + 1e-6*(300*press - 2400.0*press*s + temp*(80.0 +
		3.0*temp - 3300.0*s - 13.0*press + 47.0*press*s)));
}

//calculates the viscosity of brine in a node from the equations of
//Batzle and Wang, 1992, Seismic properties of pore fluids, Geophysics
//v. 57, n. 11, p 1396-1408.
void node::calc_viscos(double s) {
	viscos = .1 + .333*s + (1.65 + 91.9*s*s*s)*exp(-(.42*(pow(s, .8) - .17)*(pow(s, .8) - .17)
		+.045)*pow(temp, .8));
}


//calculate the potential derivatives for a given node and puts them in the given array
//the derivative is unitless
//the formula is from 
//Oelkers, 1996, Physical and Chemical Properties of Rocks
//and Fluids for Chemical Mass Transport Calculations, Reviews in Mineralogy, v 34.
void node::calc_potd() {
	int i, k;
	static double *act1, *act2, *bas1, *c1;
   static bool *pertb;
	static bool init = false;

	if (!calc_act) {
		for (k=0; k<nspe; k++) {
			h[k] = 1;
		}
		return;
	}
	if (!init) {
		init = true;
		act1 = new double[nspe];
		act2 = new double[nspe];
		bas1 = new double[nspe];
		c1 = new double[nspe];
      pertb = new bool[nspe];
	}
	
   for (i = 0; i < nspe; i++) {
      pertb[i] = false;
      this->delts[i] = spe[i] * 1e-2;
   }

	//get the unperturbed activity coefficients
	run_model(nspe, this, c1, act1, bas1, pertb, false);

	//for each species, perturb the concentration and get the new
	//activity the potential derivative is 1 + dln[a]/dc
	for (k=0; k<nspe; k++) {

      //perturb one value and see how everything changes
      for (i = 0; i < nspe; i++)
         pertb[i] = false;
      pertb[k] = true;

      //get the perturbed activity 
		run_model(nspe, this, c1, act2, bas1, pertb, false);
		
		h[k] = 1 + (log(act2[k]) - log(act1[k]))/log(this->delts[k]);
	}
   for (int j = 0; j < nspe; j++)
      this->delts[j] = 0;
}




//calculate the diffusion matrix for a given node and puts them in the given array
//the diffusion coefficients have units of cm^2/sec
//the formula is from 
//Oelkers, 1996, Physical and Chemical Properties of Rocks
//and Fluids for Chemical Mass Transport Calculations, Reviews in Mineralogy, v 34.
void node::calc_dc() {
	int i, j, k;
	double denom;

	calc_potd();

	denom = 0;
	for (k=0; k<nspe; k++)
		denom += Spe[k].z*Spe[k].z*dt[k]*spe[k];
	denom = 1/denom;
	for (i=0; i<nspe; i++) {
		for (j=0; j<nspe; j++) {
			if (i == j)
				d1[i][j] = dt[i]*h[i];
			else
				d1[i][j] = 0;
			d1[i][j] -= Spe[i].z*dt[i]*spe[i]*denom*Spe[j].z*dt[j]*h[j];
 		}
	}
   if (tort!=1.0) {
      for (i = 0; i < nspe; i++) {
         for (j = 0; j < nspe; j++) {
            d1[i][j] *= tort;
         }
      }
   }

}



//calculates tracer diffusion coefficients from the nerst einstein equation
//using conductivities v temperature relations and puts them in global array Dt
//the diffusion coefficients have units of cm^2/sec
//the formula is from 
//Nigrini, 1970, Diffusion in Rock Alteration Systems: I. Prediction of Limiting Equivalent 
//Ionic Conductances at Elevated Temperatures, American Journal of Science, v 269, p 65-91.
void node::calc_dt() {
	int j, close, ind1, ind2;
	double intrpol, at, bt, cond;

	//figure out which indices to use for the interpolation
	close = NTEMP - 1;
	for (j=0; j<NTEMP; j++) {
		if (temp > TempC[j]) {
			close = j - 1;
			break;
		}
	}
	if (close == -1) {
		ind1 = 0;
		ind2 = 1;
	} 
	else if (close == NTEMP - 1) {
		ind1 = close - 2;
		ind2 = close - 1;
	}
	else {
		ind1 = close;
		ind2 = close + 1;
	}
	//calculate the correct factors to do the interpolation
	intrpol = (temp - 273.16 - TempC[ind1])/(TempC[ind2] - TempC[ind1]);

	//build the new tracer diffusion matrix
	for (j=0; j<nspe; j++) {
		if (Spe[j].cl == 6)//handle hydrogen
			cond = 4.1354*(temp - 273.16) + 246.265;
		else if (Spe[j].cl == 7) //handle hydroxyl
			cond = 3.3795*(temp - 273.16) + 113.04;
		else {
			//calculate the temperature corrected conductivity
			at = AT[Spe[j].cl][ind1] + intrpol*(AT[Spe[j].cl][ind2] - AT[Spe[j].cl][ind1]);
			bt = BT[Spe[j].cl][ind1] + intrpol*(BT[Spe[j].cl][ind2] - BT[Spe[j].cl][ind1]);
			cond = at + bt*Spe[j].cond;
		}
		if (Spe[j].cl != 8)
			dt[j] = cond*GCONST*temp/abs(Spe[j].z)/FARAD/FARAD;
		else // diffusion coefficient for CO2 taken from Lu et al. Geochem. et. Cosmo. Acta 115 (2013) 183-204.
			dt[j] = 13.942e-5 * pow(temp/227.0 - 1, 1.7094);
	}

	return;
}

//compares the composition of one node to another, but not temperature
bool node::same_as(node* other) {
	bool match = true;
		
	for (int i = 0; i < nspe; i++) {
		if (spe[i] != other->spe[i])
			match = false;
	}

	return match;
}
