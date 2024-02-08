
#include <string>

//a class containing all of the necessary info about a species
class species {
public:
	//the class- 0: +1 cation, 1: +2 cations, 2: +3 cations, 3: halogens, 4: -1 oxyanions,
	//5: -2 oxyanions, 6: Hydrogen, 7: hydroxyl
	int cl;
	//conductivity at infinite dilution and 25 C units of cm2 Siemans / mole
	//from the handbook of chemistry and physics
	double cond;
	//the charge
	double z;
	//species name in a form SpecE8 will accept
	std::string name;
};