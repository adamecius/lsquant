// [LEGACY -- NOT BUILT] previous-working inline_kpm-S-inversion. Reference only.
// Supported path: lsquant <cmd>  |  port + golden-gate before any reuse.
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

namespace spectral
{
	void printHelpMessage();
	void printWelcomeMessage();
};

void postProcess(chebyshev::Moments1D& );


int main(int argc, char *argv[])
{
	if ( !(argc == 3 || argc == 4 ) )
	{
		spectral::printHelpMessage();
		return 0;
	}
	else
		spectral::printWelcomeMessage();
	
	const std::string
		LABEL  = argv[1],
	        S_NMOM = argv[2];

	const int numMoms   = atoi(S_NMOM.c_str() );
	chebyshev::Moments1D chebMoms( numMoms ); //load number of moments

	const int NOPS = 1;
	SparseMatrixType OP[NOPS];
	OP[0].SetID("HAM");
	//OP[1].SetID( S_OP );
	
	// Build the operators from Files
	SparseMatrixBuilder builder;
	std::array<double,2> spectral_bounds;	
	
	for (int i = 0; i < NOPS; i++)
	if( OP[i].isIdentity() )
			std::cout<<"The operator "<<OP[i].ID()<<" is treated as the identity"<<std::endl;
	else
	{
		std::string input = "operators/" + LABEL + "." + OP[i].ID() + ".CSR";
		builder.setSparseMatrix(&OP[i]);
		builder.BuildOPFromCSRFile(input);
		
		if( i == 0 ) //is hamiltonian, obtain automatically the energy bounds
		 spectral_bounds = chebyshev::utility::SpectralBounds(OP[0]);
	}
	
	//CONFIGURE THE CHEBYSHEV MOMENTS
	chebMoms.SystemLabel(LABEL);
	chebMoms.BandWidth ( (spectral_bounds[1]-spectral_bounds[0])*1.0);
	chebMoms.BandCenter( (spectral_bounds[1]+spectral_bounds[0])*0.5);
	chebMoms.SetAndRescaleHamiltonian(OP[0]);
	chebMoms.Print();


	//Compute the chebyshev expansion table
	qstates::generator gen;
	if( argc == 4)
		gen  = qstates::LoadStateFile(argv[3]);

		
	//gen.NumberOfStates(std::stoi(argv[4]));



	chebyshev::SpectralMoments_S_inv(OP[0],chebMoms, gen, spectral_bounds);

	/*
	auto prefix="SpectralOp"+OP[1].ID();
	if( OP[1].isIdentity() )
		prefix="SpectralOp"+OP[1].ID();
	std::string outputfilename=prefix+LABEL+"KPM_M"+S_NMOM+"_state"+gen.StateLabel()+".chebmom1D";	

	std::cout<<"Saving the moments in  "<<outputfilename<<std::endl;
	chebMoms.saveIn(outputfilename);
	*/

	
	
	std::cout<<"End of program"<<std::endl;
	return 0;
};

void postProcess(chebyshev::Moments1D& mu ){
	
  /*	if (argc != 3)
	{
		printHelpMessage();
		return ;
	}
	else
		printWelcomeMessage();
  */
	//const double broadening  = stod(argv[2]);


	//Read the chebyshev moments from the file
 	//chebyshev::Moments1D mu(argv[1]); 
	//and apply the appropiated kernel

	const int numMoms = mu.HighestMomentNumber();
	double broadening = 1.0 / double(numMoms);
	
	mu.ApplyJacksonKernel(broadening);

	const int num_div = 30*mu.HighestMomentNumber();
	
	const double
	xbound = chebyshev::safety_factors().recon_cutoff;
		
	std::vector< double >  energies(num_div,0);
	for( int i=0; i < num_div; i++)
		energies[i] = -xbound + i*(2*xbound)/(num_div-1) ;



	std::cout<<"Computing the Spectral function using "<<mu.HighestMomentNumber()<<" moments "<<std::endl;
	std::cout<<"Over a grid normalized grid  ("<<-xbound<<","<<xbound<<") in steps of "<<energies[1] - energies[0]<< std::endl;
	std::cout<<"The first 10 moments are"<<std::endl;
	for(int m0 =0 ; m0< 1; m0++ )
		std::cout<<mu(m0).real()<<" "<<mu(m0).imag()<<std::endl;
	

	std::string
	outputName  ="mean"+mu.SystemLabel()+"JACKSON.dat";

	std::cout<<"Saving the data in "<<outputName<<std::endl;
	std::cout<<"PARAMETERS: "<< mu.SystemSize()<<" "<<mu.HalfWidth()<<std::endl;
	std::ofstream outputfile( outputName.c_str() );

	std::vector<double> output(num_div,0.0);
	#pragma omp parallel for
	for( int i=0; i < num_div; i++)
	{
		output[i] = 0;
		const double energ = energies[i];
		for( int m0 = 0 ; m0 < mu.HighestMomentNumber() ; m0++)
			output[i] += delta_chebF(energ,m0)*mu(m0).real() ;
		output[i] *= mu.SystemSize()/mu.HalfWidth();
	}

	for( int i=0; i < num_div; i++)
		outputfile<<energies[i]*mu.HalfWidth() + mu.BandCenter() <<" "<< output[i] <<std::endl;

	outputfile.close();

	std::cout<<"The program finished succesfully."<<std::endl;
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



