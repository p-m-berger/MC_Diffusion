#define DIFF_MAIN 1
#include "diffusion.h"
#include "diffusion_proto.h"


#include <omp.h>


int main(int argc, char **argv) {
	//time is the time at the current time step, tend is the maximum time, dt is the length of the time step
	double time, tend, tmax, dt, *holdspe, *holdbas, *holdact, tot;
	int step = 0;
	bool *update;

	if (!load_input(argv[1], &time, &tend, &tmax))
		return 0;
	
	//all time variables in seconds
   dt = tend;
	write_step(time, 0);
	holdspe = new double[Nspe];
   holdbas = new double[Nspe];
	holdact = new double[Nspe];
	update = new bool[Nspe];

	for (int i = 0; i < Nspe; i++)
		update[i] = true;

	while (time < tend) {

		calc_diff(&dt, tend, time, tmax);
		time+=dt;

		//add the fluxes into each node and recalc the concentrations

#pragma omp parallel for
		for (int i = Nnode - 1; i >= 0; i--) {
			run_model(Nspe, &Dnode[i], holdspe, holdact, holdbas, update, false);
         for (int j = 0; j < Nspe; j++) {
            Dnode[i].spe[j] = holdspe[j];
            Dnode[i].bas[j] = holdbas[j];
         }
      }

      //handle the near boundary which is 0 flux but constant CO2 fugacity
    /*  run_CO2(Nspe, &Dnode[0], holdspe, holdact, holdbas, update, false);
      for (int j = 0; j < Nspe; j++) {
         Dnode[0].spe[j] = holdspe[j];
         Dnode[0].bas[j] = holdbas[j];
      }*/
      //left hand boundary does nothing, constant 

		if (step % 1000 == 0 || step < 10) {
         if (step == 0)
            write_delts(0);
			write_step(time, step);
		}

		step++;
	}
	

	write_step(time, step);
	write_step(-1, 0);

	delete[] holdspe;
   delete[] holdbas;
	delete[] holdact;
	delete[] update;
	return 0;
}