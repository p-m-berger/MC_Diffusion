#include "Chemplugin.h"

class node {
public:
	void init(int, int, bool);
	~node();
	double **d1, *h; //diffusion coefficient, onsager coefficients, hessian of chemical potential
	double *spe, *bas, *dt, *delts; //species concentrations, basis concentrations, tracer diffusion coefficients, 
									  //and changes in chemical composition
	double press, temp, tort; //pressure, temperature, tortuosity
	double wmass, por; //mass of water, porosity
	ChemPlugin *cp; //for doing the calculations
	void calc_props();
	double getvis() {return viscos;};
	double getk() {return k;};
	void calc_dc();
	void calc_dt();
	bool same_as(node *);
private:
	double density, viscos, poros, k;
	void calc_density(double);
	void calc_viscos(double);
	void calc_potd();
	int nspe, index;
	bool calc_act;
};