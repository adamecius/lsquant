// C & C++ libraries
#include <iostream> /* for std::cout mostly */
#include <string>   /* for std::string class */
#include <fstream>  /* for std::ofstream and std::ifstream functions classes*/
#include <stdlib.h>
#include <chrono>


#include "kpm_noneqop.hpp" //Message functions
#include "compute.hpp"


namespace time_evolution
{
	void printHelpMessage();
	void printWelcomeMessage();
}


int main(int argc, char *argv[])
{
	if ( !(argc == 6 || argc == 7 ) )
	{
		time_evolution::printHelpMessage();
		return 0;
	}
	else
		time_evolution::printWelcomeMessage();

	const std::string
		LABEL  = argv[1],
		S_OP   = argv[2],
		S_NMOM = argv[3],
		S_NTIME= argv[4],
		S_TMAX = argv[5];

	const int numMoms = atoi(S_NMOM.c_str() );
	const int numTimes= atoi(S_NTIME.c_str() );
	const double tmax = stod(S_TMAX );

	// The orchestration lives in lsquant::compute_msd (shared with `lsquant compute --config` mode
	// "msd"). This driver only parses argv and delegates.
	std::string state_file;
	if ( argc == 7 ) state_file = argv[6];

	lsquant::compute_msd(LABEL, S_OP, numMoms, numTimes, tmax, state_file);

	std::cout << "End of program" << std::endl;
	return 0;
}

void time_evolution::printHelpMessage()
	{
		std::cout << "The program should be called with the following options: Label Op numMom numTimeSteps MaxTime (optional) num_states_file" << std::endl
				  << std::endl;
		std::cout << "Label will be used to look for Label.Ham, Label.Op " << std::endl;
		std::cout << "Op will be the velocity used to compute the evolution" << std::endl;
		std::cout << "numMom will be used to set the number of moments in the chebyshev table" << std::endl;
		std::cout << "numTimeSteps  will be used to set the number of timesteps in the chebyshev table" << std::endl;
		std::cout << "TimeMax  will be set the maximum time where the correlation will be evaluted " << std::endl;
	};

void time_evolution::printWelcomeMessage()
	{
		std::cout << "WELCOME: This program will compute a table needed for expanding the correlation function in Chebyshev polynomialms" << std::endl;
	};
