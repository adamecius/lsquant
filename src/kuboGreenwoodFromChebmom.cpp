
// C & C++ libraries
#include <iostream>		/* for std::cout mostly */
#include <string>		/* for std::string class */
#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/
#include <vector>		/* for std::vector mostly class*/
#include <complex>		/* for std::vector mostly class*/

#include "chebyshev_moments.hpp"
#include "observable.hpp"

// Thin wrapper (Phase 5): Kubo-Greenwood is the Fermi-surface (N x 1) degenerate case of
// Kubo-Bastin -- same double-Chebyshev sum, half the prefactor, no Fermi-sea integral. It now
// runs through the unified lsquant::reconstruct_kubo with the KUBO_GREENWOOD observable.
// (The Bug-1 fix -- prefactor -SystemSize*ScaleFactor^2, not *ShiftFactor -- lives in the
// KUBO_GREENWOOD definition; the greenwood_regression test guards it.)

void printHelpMessage();
void printWelcomeMessage();

int main(int argc, char *argv[])
{
	if (argc != 3) { printHelpMessage(); return 0; }
	printWelcomeMessage();

	const double broadening = stod(argv[2]);

	chebyshev::Moments2D mu(argv[1]);
	mu.ApplyJacksonKernel(broadening, broadening);

	std::cout << "Computing the Kubo-Greenwood reconstruction using "
	          << mu.HighestMomentNumber(0) << " x " << mu.HighestMomentNumber(1) << " moments" << std::endl;

	const std::vector<std::pair<double,double> > data = lsquant::reconstruct_kubo(mu, lsquant::KUBO_GREENWOOD);

	const std::string outputName = "KuboGreenWood_" + mu.SystemLabel() + "JACKSON.dat";
	std::cout << "Saving the data in " << outputName << std::endl;
	std::cout << "PARAMETERS: " << mu.SystemSize() << " " << mu.HalfWidth() << std::endl;
	std::ofstream outputfile(outputName.c_str());
	for (size_t i = 0; i < data.size(); ++i)
		outputfile << data[i].first << " " << data[i].second << std::endl;
	outputfile.close();

	std::cout << "The program finished succesfully." << std::endl;
	return 0;
}

void printHelpMessage()
{
	std::cout << "The program should be called with the following options: moments_filename broadening(meV)" << std::endl << std::endl;
	std::cout << "moments_filename will be used to look for .chebmom2D file" << std::endl;
	std::cout << "broadening in (meV) will define the broadening of the delta functions" << std::endl;
}

void printWelcomeMessage()
{
	std::cout << "WELCOME: This program will compute the chebyshev sum of the kubo-greenwood formula" << std::endl;
}
