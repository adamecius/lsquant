
// C & C++ libraries
#include <iostream>		/* for std::cout mostly */
#include <string>		/* for std::string class */
#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/
#include <vector>		/* for std::vector mostly class*/
#include <complex>		/* for std::vector mostly class*/
#include <omp.h>

#include "chebyshev_solver.hpp"
#include "chebyshev_coefficients.hpp"
#include "chebyshev_moments.hpp"
#include "observable.hpp"   // lsquant::reconstruct_kubo_bastin_sea (shared grid/integral/write)

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
 	chebyshev::Moments2D mu(argv[1]);
	//and apply the appropiated kernel
	mu.ApplyJacksonKernel(broadening,broadening);

	std::cout<<"Computing the KuboBastin (Sea) kernel using "<<mu.HighestMomentNumber(0)<<" x "<<mu.HighestMomentNumber(1)<<" moments "<<std::endl;

	std::string
	outputName  ="KuboBastinSea_"+mu.SystemLabel()+"JACKSON.dat";

	std::cout<<"Saving the data in "<<outputName<<std::endl;
	std::cout<<"PARAMETERS: "<< mu.SystemSize()<<" "<<mu.HalfWidth()<<std::endl;
	std::ofstream outputfile( outputName.c_str() );

	// Green-difference, cumulative Fermi-sea integral -- the shared reconstruction routine.
	const auto data = lsquant::reconstruct_kubo_bastin_sea(mu);
	for( const auto& point : data )
		outputfile<<point.first<<" "<<point.second<<std::endl;

	outputfile.close();

	std::cout<<"The program finished succesfully."<<std::endl;
return 0;
}


void printHelpMessage()
{
	std::cout << "The program should be called with the following options: moments_filename broadening(meV)" << std::endl
			  << std::endl;
	std::cout << "moments_filename will be used to look for .chebmom2D file" << std::endl;
	std::cout << "broadening in (meV) will define the broadening of the delta functions" << std::endl;
};

void printWelcomeMessage()
{
	std::cout << "WELCOME: This program will compute the chebyshev sum of the kubo-bastin formula for non equlibrium properties" << std::endl;
};
