
#include "diffusion.h"
#include "diffusion_proto.h"

#include <fstream>
#include <stdio.h>
#include <iomanip>

using namespace std;


//loads in the input from a script with a semi-fixed format
//arguments: path - character array with the file path, tend - end time
//    of the model run (everything else is stored in globals allocated in this function)
int load_input(char *path, double *tstart, double *tend, double *tmax) {
	int i;
	char junk[512];
	ifstream input(path, ios::in), prv;
	double hold, tor, *holdc, *holdb;
	bool indata = true, prev, act;

	//see if its there
	if (input.fail())
		return 0;

	//load in the number of species and initialize the arrays
	input>>junk>>Nspe;

	Spe = new species [Nspe];
   holdc = new double[Nspe];
   holdb = new double[Nspe];

	//load number of nodes
	input>>junk>>Nnode;

	//load whether or not to calcuate the acitivity coefficient derivative
	input>>junk>>i;
	act = (bool)i;

	Dnode = new node [Nnode];
	for (i=0; i<Nnode; i++) {
		Dnode[i].init(Nspe, i, act);
	}

	//node spacing (m) and end time in seconds
	input>>junk>>Dx;
	input>>junk>>*tend;
	//convert the node spacing to cm
	Dx *= 100;

	//input temp is in C, convert it to K
	input>>junk>>i;
	if (i) {//is the system isothermal
		input>>junk>>hold;
		hold += 273.15;
		for (i=0; i<Nnode; i++)
			Dnode[i].temp = hold;
	}
	else {//or not
		input>>junk>>hold;
		hold += 273.15;
		Dnode[0].temp = hold;
		//Get the increment
		input >> hold;
		//add the increment
		for (i=1; i<Nnode; i++) {
			Dnode[i].temp = Dnode[0].temp + hold*(double)i;
		}
	}

	input>>junk>>*tmax;
	input>>junk>>i;
	Fixedbound = (bool)i;
   FixedCO2 = false;
   if (!Fixedbound) {
      input>>junk>>i;
      FixedCO2 = (bool)i;
   }

	input >> junk >> tor;
	for (i = 0; i<Nnode; i++)
		Dnode[i].tort = tor;

	//load in the species information
	Cl = -1;
	for (i=0; i<Nspe; i++) {
		input>>Spe[i].name>>Spe[i].z>>Spe[i].cond>>Spe[i].cl;
		if (Spe[i].name == "Cl-")
			Cl = i;
	}

	input >> junk >> i;
	prev = (bool)i;

	if (!prev) {
		//load in the commands for each domain
		*tstart = 0;
		Nin = Ndom = 0;
		input.ignore();
		while (!input.eof()) {
			input.getline(junk, 512, '\n');
			if (!strcmp(junk, "scope inlet"))
				indata = true;
			else if (!strcmp(junk, "scope domain"))
				indata = false;
			else if (!strcmp(junk, "end"))
				break;
			else {
				if (indata) {
					strcpy(Inlet[Nin], junk);
					Nin++;
				}
				else {
					strcpy(Domain[Ndom], junk);
					Ndom++;
				}
			}
		}

		//get the initial concentrations at the inlet
		initc(Nspe, &Dnode[0], Nin, Dnode[0].temp, true, holdc, holdb);
      for (i = 0; i < Nspe; i++) {
         Dnode[0].spe[i] = holdc[i];
         Dnode[0].bas[i] = holdb[i];
      }

		//get the initial concentrations for the domain
		initc(Nspe, &Dnode[1], Nin, Dnode[1].temp, false, holdc, holdb);


		//once the species information is loaded, we can calculate the tracer diffusion coefficients
		//run the model to get the initial concentrations
		for (i = 0; i < Nnode; i++) {
			Dnode[i].calc_dt();
			if (i != 0) {
            for (int j = 0; j < Nspe; j++) {
               Dnode[i].spe[j] = holdc[j];
               Dnode[i].bas[j] = holdb[j];
            }
			}
		}

		delete[] holdc;
      delete[] holdb;
	}
	else {
		//just load in the concentrations and calculate the tracer diffusion coef.
		prv.open("previous.txt", ios::in);
		prv >> *tstart;
		for (i = 0; i < Nnode; i++) {
			for (int j = 0; j < Nspe; j++)
				prv >> Dnode[i].spe[j];
			Dnode[i].calc_dt();
		}

		prv.close();
	}

	return 1;
}

//called when a fatal error occurs
//arguements: err - character string with an error description
void panic(char *err) {
	static ofstream output("err.dat", ios::out);
   output << err;
	output.close();
	exit(1);
}

//writes the concentrations at each node at a time
//arguments: time - current time step, -1 is a flag to close the file
void write_step(double time, int iter) {
	int i, j;
	static ofstream output("output.dat", ios::out);
   static ofstream outputb("outputb.dat", ios::out);

   //close the file
	if (time == -1) {
		output.close();
      outputb.close();
		return;
	}

	//output the time
	output <<time<<"\t"<<iter<<endl<<"m\t";
   outputb<<time<<"\t"<<iter<<endl<< "m\t";
	//output the nodal spacing
	for (i=0; i<Nnode; i++) {
		output <<i*Dx/100.0<<"\t";
      outputb<<i*Dx/100.0<<"\t";
	}
	output<<endl;
	//output the species concentrations
	for (j=0; j<Nspe; j++) {
		output<<Spe[j].name<<"\t";
		for (i=0; i<Nnode; i++) {
			output <<setprecision(14)<<Dnode[i].spe[j]<<"\t";
         outputb<<setprecision(14)<<Dnode[i].bas[j]<<"\t";
		}
		output <<endl;
      outputb<<endl;
	}
	output<<"Temp"<<"\t";
	for (i=0; i<Nnode; i++) {
		output <<setprecision(14)<<Dnode[i].temp-273.15<<"\t";
	}
	output <<endl;

	return;
}

//this outputs the mineralogy at each node into a file
//arguments: name - character string with the file name, cp - pointer to the first chemplugin node,
//    nnode - number of chemplugin nodes
/*void write_mins(char *name, int nnode) {
   int i, j, k, maxmin = 0, back, curmin;
   //I concede that this is a bit of a hack
   char minlist[1048], **mins, cmd[256], amin[256];
   bool found;
   ofstream output;

   //first see how many minerals there are
   //create a list and check to see if new minerals need to be added
   minlist[0] = '\0';
   for (i = 0; i < nnode; i++) {
      back = (int)cp[i].Report1("nminerals");
      for (j = 0; j < back; j++) {
         sprintf(cmd, "minerals %i", j);
         cp[i].Report(amin, cmd);
         if (strstr(minlist, amin)!=NULL) {
            strcat(minlist, " ");
            strcat(minlist, amin);
            maxmin++;
         }
      }
   }

   //allocate for the mineral names
   mins = new char *[maxmin];
   for (i = 0; i < maxmin; i++)
      mins[i] = new char[256];

   //copy over the mineral names
   curmin = 0;
   for (i = 0; i < nnode; i++) {
      found = false;
      back = (int)cp[i].Report1("nminerals");
      for (j = 0; j < back; j++) {
         sprintf(cmd, "minerals %i", j);
         cp[i].Report(amin, cmd);
         for (k = 0; k < curmin; k++) {
            if (!strcmp(amin, mins[i])) {
               found = true;
               break;
            }
         }
         if (!found) {
            strcpy(mins[curmin], amin);
            curmin++;
         }
      }
   }

   //print out the values
   output.open(name, ios::out);
   for (i = 0; i < maxmin; i++) {
      output << mins[i] << "\t";
      for (j = 0; j < nnode; j++) {
         sprintf(cmd, "mass minerals %s", mins[i]);
         output << cp[j].Report1(cmd) << "\t";
      }
   }
   if (maxmin == 0) output << "none";
   output.close();

   //clean up
   for (i = 0; i < maxmin; i++)
      delete[] mins[i];
   delete[] mins;

   return;
}*/

//a function for debugging, writes out the calculated changes in concentration
//arguments: time - current time step, -1 is a flag to close the file
void write_delts(double time) {
	int i, j;
	static ofstream output("outdelta.dat", ios::out);

   //close the file
	if (time == -1) {
		output.close();
		return;
	}

	//output the time
	output << time << endl << "m\t";
	//output the nodal spacing
	for (i = 0; i<Nnode; i++) {
		output << i*Dx / 100 << "\t";
	}
	output << endl;
	//output the species concentrations
	for (j = 0; j<Nspe; j++) {
		output << Spe[j].name << "\t";
		for (i = 0; i<Nnode; i++) {
			output << Dnode[i].delts[j] << "\t";
		}
		output << endl;
	}


	return;
}

//another debugging tool, outputs the contents of an array
//arguments: array of doubles to output, size of the array (assumes it is square)
void out_d(double **arr, int size, double val) {
	int i, j;
	ofstream output("dees.dat", ios::out);

	for (i=0; i<size; i++) {
		for (j=0; j<size; j++) {
			output<<arr[i][j]<<"\t";
		}
		output<<endl;
	}
	output << val;

	return;
}