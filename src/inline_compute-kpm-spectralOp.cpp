// C & C++ libraries
#include <iostream> /* for std::cout mostly */
#include <string>   /* for std::string class */
#include <fstream>  /* for std::ofstream and std::ifstream functions classes*/
#include <stdlib.h>
#include <chrono>


#include "kpm_noneqop.hpp" //Message functions
#include "compute.hpp"

namespace spectral
{
	void printHelpMessage();
	void printWelcomeMessage();
};


int main(int argc, char *argv[])
{
	if ( !(argc == 4 || argc == 5 ) )
	{
		spectral::printHelpMessage();
		return 0;
	}
	else
		spectral::printWelcomeMessage();

	const std::string
		LABEL  = argv[1],
		S_OP   = argv[2],
		S_NMOM = argv[3];

	const int numMoms = atoi(S_NMOM.c_str());

	// The orchestration lives in lsquant::compute_spectral (shared with `lsquant compute --config`
	// mode "spectral"). This driver only parses argv and delegates.
	std::string state_file;
	if ( argc == 5 ) state_file = argv[4];

	lsquant::compute_spectral(LABEL, S_OP, numMoms, state_file);

	std::cout << "End of program" << std::endl;
	return 0;
};


	void spectral::printHelpMessage()
	{
		std::cout << "The program should be called with the following options: Label Op numMom " << std::endl
				  << std::endl;
		std::cout << "Label will be used to look for Label.Ham, Label.Op" << std::endl;
		std::cout << "Op will be used to located the  matrix file of the operators for the spectral " << std::endl;
		std::cout << "numMom will be used to set the number of moments in the chebyshev table" << std::endl;
	};


	void spectral::printWelcomeMessage()
	{
		std::cout << "WELCOME: This program will compute a list needed for expanding the spectral function in Chebyshev polynomials" << std::endl;
	};
