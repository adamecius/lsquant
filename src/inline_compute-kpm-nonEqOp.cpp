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

	const int numMoms= atoi(S_NUM_MOM.c_str());
	chebyshev::Moments2D chebMoms(numMoms,numMoms); //load number of moments

	SparseMatrixType OP[3];
	OP[0].SetID("HAM");
	OP[1].SetID(S_OPR);
	OP[2].SetID(S_OPL);

	// Build the operators from Files
	SparseMatrixBuilder builder;
	std::array<double,2> spectral_bounds;	
	for (int i = 0; i < 3; i++)
	{
		std::string input = "operators/" + LABEL + "." + OP[i].ID() + ".CSR";
		builder.setSparseMatrix(&OP[i]);
		builder.BuildOPFromCSRFile(input);
		
		if( i == 0 ) //is hamiltonian
		//Obtain automatically the energy bounds
		 spectral_bounds = chebyshev::utility::SpectralBounds(OP[0], "operators/" + LABEL + ".HAM.desc");
	};
	//CONFIGURE THE CHEBYSHEV MOMENTS
	chebMoms.SystemLabel(LABEL);
	chebMoms.BandWidth ( (spectral_bounds[1]-spectral_bounds[0])*1.0);
	chebMoms.BandCenter( (spectral_bounds[1]+spectral_bounds[0])*0.5);
	chebMoms.SetAndRescaleHamiltonian(OP[0]);
	chebMoms.Print();


	//Define thes states youll use
	//Factory state_factory ;

	//Compute the chebyshev expansion table
	qstates::generator gen;
	if( have_state )
		gen  = qstates::LoadStateFile(state_file);

	chebyshev::CorrelationExpansionMoments( OP[1], OP[2], chebMoms, gen );

	//Save the table in a file
	std::string outputfilename="NonEqOp"+S_OPR+"-"+S_OPL+LABEL+"KPM_M"+S_NUM_MOM+"x"+S_NUM_MOM+"_state"+gen.StateLabel()+".chebmom2D";
	chebMoms.saveIn(outputfilename);

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
