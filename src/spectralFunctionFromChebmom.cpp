		
// C & C++ libraries
#include <iostream>		/* for std::cout mostly */
#include <string>		/* for std::string class */
#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/
#include <vector>		/* for std::vector mostly class*/
#include <complex>		/* for std::vector mostly class*/
#include <omp.h>

#include "chebyshev_moments.hpp"
#include "observable.hpp"

// Thin wrapper (Phase 5): the spectral-function DOS sum is the unified
// lsquant::reconstruct_density_grid; this driver applies its prefactor (N*ScaleFactor) and
// energy axis (x/ScaleFactor + ShiftFactor). Byte-identical to the legacy loop.

void printHelpMessage();
void printWelcomeMessage();

int main(int argc, char *argv[])
{	
	if (argc != 3)
	{
		printHelpMessage();
		return 0;
	}
	else
		printWelcomeMessage();

	const double broadening  = stod(argv[2]);


	//Read the chebyshev moments from the file
 	chebyshev::Moments1D mu(argv[1]); 
	//and apply the appropiated kernel
	mu.ApplyJacksonKernel(broadening);

	const int num_div = 30*mu.HighestMomentNumber();

	std::vector<double> mu_real( mu.HighestMomentNumber() );
	for( int m = 0 ; m < mu.HighestMomentNumber() ; m++ ) mu_real[m] = mu(m).real();
	const std::vector<std::pair<double,double> > grid =
		lsquant::reconstruct_density_grid( mu_real, chebyshev::safety_factors().recon_cutoff, num_div );

	std::string
	outputName  ="mean"+mu.SystemLabel()+"JACKSON.dat";

	std::cout<<"Computing the Spectral function using "<<mu.HighestMomentNumber()<<" moments "<<std::endl;
	std::cout<<"Saving the data in "<<outputName<<std::endl;
	std::cout<<"PARAMETERS: "<< mu.SystemSize()<<" "<<mu.HalfWidth()<<std::endl;
	std::ofstream outputfile( outputName.c_str() );

	for( size_t i = 0; i < grid.size(); ++i )
		outputfile<< grid[i].first/mu.ScaleFactor() + mu.ShiftFactor()
		          <<" "<< grid[i].second*( mu.SystemSize()*mu.ScaleFactor() ) <<std::endl;

	outputfile.close();

	std::cout<<"The program finished succesfully."<<std::endl;
return 0;
}
	

void printHelpMessage()
{
	std::cout << "The program should be called with the following options: moments_filename broadening(meV)" << std::endl
			  << std::endl;
	std::cout << "moments_filename will be used to look for .chebmom1D file" << std::endl;
	std::cout << "broadening in (meV) will define the broadening of the delta functions" << std::endl;
};

void printWelcomeMessage()
{
	std::cout << "WELCOME: This program will compute the chebyshev sum of the spectral function" << std::endl;
};
