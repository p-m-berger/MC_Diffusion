
#include "diffusion.h"
#include "diffusion_proto.h"

#include <fstream>
#include <stdio.h>
#include <string>

#include "Chemplugin.h"

//utility to help with precision problems in going to strings
std::string double2hexstr(double x) {
   uint64_t i;
   std::memcpy(&i, &x, 8);

   char buf[17];
   snprintf(buf, sizeof(buf), "%016llx", i);
   buf[16] = 0; //make sure it is null terminated.

   return std::string(buf);
}

void run_model(int nspe, node *anode, double *concs, double *act, double *basis, bool *change, bool save) {
   using namespace std;
   string cmd;
	int cnt = 0;
	ChemPlugin cp("", "-d \"C:/Program Files/Gwb/Gtdata/thermo_ymp.R2.tdat\"");

   cmd = "temperature " + to_string(anode->temp - 273.15) + "; suppress all; balance off; plot off; printout = off";
   cp.Config(cmd.c_str());
   cp.Config("suppress CO3-- CaCO3(aq) CaCl+ CaOH+; reset reactants; cpu_max 360");

   for (int i = 0; i < Nspe; i++) {
      //change the value for the new node
      if (basis[i] != -1) {
         cmd = Spe[i].name + " mol 0x" + double2hexstr(anode->bas[i]);
         cp.Config(cmd.c_str());
      }
      if (change[i]) {
         cmd = "react mol 0x" + double2hexstr(anode->delts[i]) + " " + Spe[i].name;
         cp.Config(cmd.c_str());
      }
   }

   cp.Initialize();
	while (true) {
		double deltat = cp.ReportTimeStep();
		if (cp.AdvanceTimeStep(deltat)) break;
		if (cp.AdvanceChemical()) break;
		cnt++;
		if (cnt > 100000) {
			cp.Config("save broken.rea");
			panic("did not converge in run_model");
		}
	}

   //get the concentrations
   for (int i = 0; i < nspe; i++) {
      //get results
      cmd = "mass aqueous " + Spe[i].name;
      concs[i] = cp.Report1(cmd.c_str(), "mol");
      if (concs[i] < 0)
         panic("crash in run_model");
   }

   //get the activities
   for (int i = 0; i < nspe; i++) {
      //get results
      cmd = "activity " + Spe[i].name;
      act[i] = cp.Report1(cmd.c_str());
      if (concs[i] < 0)
         panic("crash in run_model");
   }

   //get the basis concentrations
   for (int i = 0; i < nspe; i++) {
      //get results
      if (Spe[i].name != "OH-" && Spe[i].name != "CO2(aq)") {
         cmd = "concentration original " + Spe[i].name;
         basis[i] = cp.Report1(cmd.c_str(), "mol");
         if (basis[i] < 0 && Spe[i].name != "H+")
            panic("crash in run_model");
      }
      else {
         basis[i] = -1;
      }
   }
}



//runs the initial model to determine the initial concentrations
void initc(int nspe, node *anode, int ncmd, double temp, bool isinlet, double *concs, double *basis) {
   using namespace std;
   string cmd;

   cmd = "temperature " + to_string(temp - 273.15) + "; plot off; printout off";
   anode->cp->Config(cmd.c_str());
   anode->cp->Config("suppress CO3-- CaCO3(aq) CaCl+ CaOH+");

   //issue the commands
   for (int i = 0; i < ncmd; i++) {
      if (isinlet)
         anode->cp->Config(Inlet[i]);
      else
         anode->cp->Config(Domain[i]);
   }

	anode->cp->Config("save test.cp");
   if (anode->cp->Initialize())
		panic("failure to init");

   //get the concentrations
   for (int i = 0; i < nspe; i++) {
      //get results
      cmd = "mass aqueous " + Spe[i].name;
		concs[i] = anode->cp->Report1(cmd.c_str(), "mol/kg");
      if (concs[i] < 0) {
         panic("crash getting the concs");
		}
   }

   //get the basis concentrations
   for (int i = 0; i < nspe; i++) {
      //get results
      if (Spe[i].name != "OH-" && Spe[i].name != "CO2(aq)") {
         cmd = "concentration original " + Spe[i].name;
			basis[i] = anode->cp->Report1(cmd.c_str(), "mol/kg");
         if (basis[i] < 0 && Spe[i].name != "H+")
            panic("crash getting the basis concs");
      }
      else {
         basis[i] = -1;
      }
   }

}
