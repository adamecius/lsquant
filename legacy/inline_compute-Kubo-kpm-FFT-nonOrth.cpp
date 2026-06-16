// [LEGACY -- NOT BUILT] previous-working inline_compute-Kubo-kpm-FFT-nonOrth. Reference only.
// Supported path: lsquant <cmd>  |  port + golden-gate before any reuse.
// C& C++ libraries
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

#include "Kubo_solver_FFT.hpp"

namespace kpmKubo
{

void printHelpMessage();

void printWelcomeMessage();

}
int main(int argc, char *argv[])
{
	if ( !(argc == 5 || argc == 6) )
	{
		kpmKubo::printHelpMessage();
		return 0;
	}
	else
		kpmKubo::printWelcomeMessage();
	
	const std::string
		LABEL = argv[1],
		S_OPR = argv[2],
		S_OPL = argv[3],
		S_NUM_MOM = argv[4];

	
	const int numMoms= atoi(argv[4]);
        int num_sections = 1, nump = numMoms;
	chebyshev::formula sym_formula = chebyshev::KUBO_GREENWOOD;
	//chebyshev::Moments Hamiltonian_dummyMoms; //load number of moments
                

	SparseMatrixType OP[5];
	OP[0].SetID("HAM");
	OP[1].SetID( "S" );
	OP[2].SetID( "dSdK" );
	OP[3].SetID( "dHdK" );

	//OP[2].SetID(S_OPR);
	//OP[3].SetID(S_OPL);

	// Build the operators from Files
	SparseMatrixBuilder builder;
	std::array<double,2> spectral_bounds;	
	for (int i = 0; i < 4; i++)
	{
		std::string input = "operators/" + LABEL + "." + OP[i].ID() + ".CSR";
		builder.setSparseMatrix(&OP[i]);
		builder.BuildOPFromCSRFile(input);
		
		if( i == 0 ) //is hamiltonian
		//Obtain automatically the energy bounds
		 spectral_bounds = chebyshev::utility::SpectralBounds(OP[0]);
	};


	/*Testing
	Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType>
	  H = OP[0].eigen_matrix(),
	  S = OP[1].eigen_matrix(),
	  x1 = OP[2].eigen_matrix(),
	  x2 = OP[3].eigen_matrix();

	Eigen::SparseLU<Eigen::SparseMatrix<complex<double>>> lin_solver;
        lin_solver.compute(S);

	
	Eigen::SparseMatrix<complex<double>> I(S.rows(),S.cols());
        I.setIdentity();
        Eigen::SparseMatrix<complex<double>> S_inv = lin_solver.solve(I);

	H=S_inv * H;
	x1=H*x1 - x1*H;
	x2=H*x2 - x2*H;
	*/


  
	chebyshev::Vectors_sliced_nonOrthogonal 
	  chebVec( OP[1],OP[2], OP[3], numMoms, num_sections ),
          chebVec_2( OP[1],OP[2], OP[3], numMoms, num_sections );

	
	//CONFIGURE THE CHEBYSHEV MOMENTS
	chebVec.SystemLabel(LABEL);
	chebVec.BandWidth ( (spectral_bounds[1]-spectral_bounds[0])*1.0);
        chebVec.BandCenter( (spectral_bounds[1]+spectral_bounds[0])*0.5);
	chebVec.SetAndRescaleHamiltonian(OP[0]);

	chebVec_2 = chebVec;
	

	//Define thes states youll use
	//Factory state_factory ;


	
	//Compute the chebyshev expansion table
	qstates::generator gen;
	if( argc == 6)	
		gen  = qstates::LoadStateFile(argv[5]);

	std::string outputfilename="Greenwood_FFT"+S_OPR+"-"+S_OPL+LABEL+"KPM_M"+S_NUM_MOM+"x"+S_NUM_MOM+"_state"+gen.StateLabel()+".conductivity_nonOrth";

	chebyshev::Kubo_solver_FFT solver(numMoms,  num_sections, nump, sym_formula, chebVec,chebVec_2,  outputfilename);
	solver.compute( OP[0], OP[0], gen );

	



	
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
