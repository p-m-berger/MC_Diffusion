Coder

Peter M. Berger
GFZ German Research Centre for Geosciences
Telegrafenberg, 14473, Germany 

Software capabilities

The software is built to calculate multi-component diffusion of ions over a user defined domain based on the equations of previous workers which are summarized in Oelkers (1996). These equations account for charge balance and the derivative of the species activity coefficient with respect to changes in concentration. It uses a simple, fully explicit finite difference approximation to calculate diffusion profiles and has been compared to a previous benchmark used by other software packages (Rasouli et al., 2015).
The software is capable of simulating three different boundary conditions at the left boundary: closed, closed with a constant CO2 concentration, and constant. The right boundary is always closed. One can also apply a temperature gradient over the domain, though the software has not been tested against a Soret effect benchmark. 

Software compilation

The software is written in C++ and while it was originally developed using Visual Studio, it should compile with any suitable compiler. The software does use ChemPlugin to do the activity coefficient calculations. This is a third-party software library that must be purchased separately. The use of ChemPlugin by this software is in no way an endorsement, rather it is a question of personal preference. Regardless, the code does require ChemPlugin.h to compile and the associated library to link. No make file is provided.

Software alteration

The user is free to alter this software as they see fit. If they do not like using ChemPlugin, several other geochemical software libraries, such as PHREEQC	, should be easy enough to swap out. It is also possible for the software to be extended to include mineral reactions and different boundary conditions. If you do this, please document your code as well any changes in the input file format.

Software use

The pdf manual contains instructions on how to set up and run the software.

Software license

This software is released under the GNU GPL License 3.0. See the license.txt file included with this software for details or check www.gnu.org/licenses/

References

Oelkers, E.H., 1996, Chapter 3. PHYSICAL AND CHEMICAL PROPERTIES OF ROCKS AND FLUIDS FOR CHEMICAL MASS TRANSPORT CALCULATIONS, in Lichtner, P.C., Steefel, C.I., and Oelkers, E.H. eds., Reactive Transport in Porous Media, Berlin, Boston, De Gruyter, p. 131–192, doi:10.1515/9781501509797-006.

Rasouli, P., Steefel, C.I., Mayer, K.U., and Rolle, M., 2015, Benchmarks for multicomponent diffusion and electrochemical migration: Computational Geosciences, v. 19, p. 523–533, doi:10.1007/s10596-015-9481-z.
