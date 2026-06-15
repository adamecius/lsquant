// C & C++ libraries
#include <iostream> /* for std::cout mostly */
#include <string>   /* for std::string class */
#include <fstream>  /* for std::ofstream and std::ifstream functions classes*/
#include <stdlib.h>
#include <chrono>


#include "kpm_noneqop.hpp" //Message functions
#include "chebyshev_moments.hpp"
#include "sparse_matrix.hpp"
#include "quantum_states.hpp"
#include "chebyshev_solver.hpp"
#include "run_config.hpp"
#include "compute.hpp"

namespace kpmKubo
{

void printHelpMessage();

void printWelcomeMessage();

}
int main(int argc, char *argv[])
{
	// Two input paths, same computation:
	//   lsquant ... --config run.json          (Phase 4: JSON config)
	//   LABEL OPR OPL numMom [state]           (legacy positional argv -- fallback)
	std::string LABEL, S_OPR, S_OPL, S_NUM_MOM, state_file;
	bool have_state = false;

	if ( argc == 3 && std::string(argv[1]) == "--config" )
	{
		const lsquant::RunConfig c = lsquant::read_run_config(argv[2]);
		if ( !c.valid ) { std::cerr << "config error: " << c.error << std::endl; return 1; }
		kpmKubo::printWelcomeMessage();
		LABEL = c.label; S_OPR = c.operator_right; S_OPL = c.operator_left;
		S_NUM_MOM = std::to_string(c.num_moments);
		if ( c.state != "default" && !c.state.empty() ) { state_file = c.state; have_state = true; }
	}
	else if ( argc == 5 || argc == 6 )
	{
		kpmKubo::printWelcomeMessage();
		LABEL = argv[1]; S_OPR = argv[2]; S_OPL = argv[3]; S_NUM_MOM = argv[4];
		if ( argc == 6 ) { state_file = argv[5]; have_state = true; }
	}
	else
	{
		kpmKubo::printHelpMessage();
		return 0;
	}

	// The orchestration lives in lsquant::compute_noneq (shared with the `lsquant compute`
	// subcommand). This driver only parses inputs (argv or --config) and delegates.
	const int numMoms = atoi(S_NUM_MOM.c_str());
	lsquant::compute_noneq(LABEL, S_OPR, S_OPL, numMoms, have_state ? state_file : std::string());

	std::cout<<"End of program"<<std::endl;
	return 0;
}

void kpmKubo::printHelpMessage()
{
	std::cout << "The program should be called with the following options: Label Op1 Op2 numMom num_states (default 1)" << std::endl
			  << std::endl;
	std::cout << "Label will be used to look for Label.Ham, Label.Op1 and Label.Op2" << std::endl;
	std::cout << "Op1 and Op2 will be used to located the sparse matrix file of two operators for the correlation" << std::endl;
	std::cout << "numMom will be used to set the number of moments in the chebyshev table" << std::endl;
};

void kpmKubo::printWelcomeMessage()
{
	std::cout << "WELCOME: This program will compute a table needed for expanding the correlation function in Chebyshev polynomialms" << std::endl;
};
