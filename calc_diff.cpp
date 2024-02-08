#include "diffusion.h"
#include "diffusion_proto.h"
#include <math.h>
#include <algorithm>
#include <string.h>

#include <omp.h>



//this function performs the crank-nicholson finite difference algorithm
//to caculate changes in concentration due to diffusion. However, the 
//diffusion coefficients are dependent on concentration. The function uses
//the fully explicit finite difference algorithm calculate the concentration
void calc_diff(double *dt, double tend, double time, double maxstep) {
	int i, j, k, count=0;
	bool converge = true, same;
	double dtinv_dx, maxd=0, d1, d2;

	//get the diffusion coefficients for each node
	//if it has the same concentrations as the one beside it, then just copy
#pragma omp parallel for
	for (i = 0; i < Nnode; i++)
		Dnode[i].calc_dc();

	//find the absolute value of the largest diffusion coefficient
	//they can be negative
	for (i=1; i<Nnode-1; i++)
		for (j=0; j<Nspe; j++)
			for (k=0; k<Nspe; k++)
				if (fabs(Dnode[i].d1[j][k]) > maxd)
					maxd = fabs(Dnode[i].d1[j][k]);

	//reset all the delta concentrations to 0
	for (i = 0; i < Nnode; i++)
		for (j = 0; j < Nspe; j++)
			Dnode[i].delts[j] = 0;

	//probably need to update this
	//calculate the time step based on the largest diffusion coefficient
	*dt = std::min(Dx*Dx/maxd/4.0,*dt);
	if (time + *dt > tend)
		*dt = tend - time;
	if (*dt > maxstep)
		*dt = maxstep;
	dtinv_dx = *dt/Dx/Dx;
	out_d(Dnode[1].d1, Nspe, Dx * Dx / maxd);

	//(c+1 - c) = d1*dc
	//solve for v+1
   //skipping index because it may be constant and the far boundary because it is constant
	for (k = 1; k < Nnode - 2; k++) {
		for (i = 0; i < Nspe; i++) {
			for (j = 0; j < Nspe; j++) {
            d1 = 0.5 * (Dnode[k-1].d1[i][j] + Dnode[k].d1[i][j]);
            d2 = 0.5 * (Dnode[k].d1[i][j] + Dnode[k+1].d1[i][j]);
				Dnode[k].delts[i] += dtinv_dx*(d1*Dnode[k-1].spe[j] + d2*Dnode[k+1].spe[j] - Dnode[k].spe[j]*(d1 + d2));
			}
		}
	}

	//Are we having a fixed or a closed boundary
	if (!Fixedbound) {
		//calculating the 0 flux near boundary
		for (i = 0; i < Nspe; i++) {
			if ((Spe[i].name != "HCO3-" && Spe[i].name != "CO2(aq)") || !FixedCO2) {
				for (j = 0; j < Nspe; j++) {
					d1 = 0.5 * (Dnode[0].d1[i][j] + Dnode[1].d1[i][j]);
					Dnode[0].delts[i] += dtinv_dx * (-d1 * Dnode[0].spe[j] + d1 * Dnode[1].spe[j]);
				}
			}
			else
				Dnode[0].delts[i] = 0;
		}
	}

   //calculating the 0 flux far boundary
   for (i = 0; i < Nspe; i++) {
      for (j = 0; j < Nspe; j++) {
			d1 = 0.5 * (Dnode[Nnode-2].d1[i][j] + Dnode[Nnode-1].d1[i][j]);
         Dnode[Nnode - 1].delts[i] += dtinv_dx * (-d1* Dnode[Nnode-1].spe[j] + d1*Dnode[Nnode-2].spe[j]);
      }
   }

	return;
}