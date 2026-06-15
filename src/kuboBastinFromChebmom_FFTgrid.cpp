
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

void printHelpMessage();
void printWelcomeMessage();


void rearrange_crescent_order( std::vector<double>& rearranged){//The point choice of the FFT has unconvenient ordering for file saving and integrating; This fixes that.
  int nump = rearranged.size();
  double original[nump];

  for( int k=0; k < nump; k++ )
    original[k] = rearranged[k];

  
  for( int k=0; k < nump / 2; k++ ){
    rearranged[ 2 * k ]   = original[ k ];
    rearranged[ 2 * k + 1 ] = original[ nump - k - 1 ]; 
  }
    
  for( int k=0; k < nump / 2; k++ ){
    double tmp = rearranged[ k ]; 
    rearranged[ k ]   = rearranged[ nump-k-1 ];
    rearranged[ nump-k-1 ] = tmp;    
    }  
}


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

	const int num_div = mu.HighestMomentNumber();
	
	const double
	  xbound = chebyshev::KPM_ALPHA,
	  disp = 0.5;
		
	std::vector< double >  energies(num_div,0);
	for( int i=0; i < num_div; i++)
	  energies[i] = cos( M_PI * ( 2 * i + disp ) /num_div);//-xbound + i*(2*xbound)/(num_div-1) ;


	rearrange_crescent_order(energies);
	



	std::cout<<"Computing the KuboBastin kernel using "<<mu.HighestMomentNumber(0)<<" x "<<mu.HighestMomentNumber(1)<<" moments "<<std::endl;
	std::cout<<"Integrated over a grid normalized grid  ("<<-xbound<<","<<xbound<<") in steps of "<<energies[1] - energies[0]<< std::endl;
	std::cout<<"The first 10 moments are"<<std::endl;
	for(int m0 =0 ; m0< 1; m0++ )
	for(int m1 =0 ; m1< 10; m1++ )
		std::cout<<mu(m0,m1).real()<<" "<<mu(m0,m1).imag()<<std::endl;
	

	std::string
	outputName  ="KuboBastin_"+mu.SystemLabel()+"JACKSON_FFTgrid.dat";

	std::cout<<"Saving the data in "<<outputName<<std::endl;
	std::cout<<"PARAMETERS: "<< mu.SystemSize()<<" "<<mu.HalfWidth()<<std::endl;
	std::ofstream outputfile( outputName.c_str() );

	std::vector<double> kernel(num_div,0.0);
	#pragma omp parallel for
	for( int i=0; i < num_div; i++)
	{
		const double energ = energies[i];
		for( int m0 = 0 ; m0 < mu.HighestMomentNumber(0) ; m0++)
		for( int m1 = 0 ; m1 < mu.HighestMomentNumber(1) ; m1++)
			kernel[i] += delta_chebF(energ,m0)*( DgreenR_chebF(energ,m1)*mu(m0,m1) ).imag() ;
		kernel[i] *= -2.0* mu.SystemSize()/mu.HalfWidth()/mu.HalfWidth();
	}

	double acc = 0;
	for( int i=0; i < num_div-1; i++)
	{
		const double energ  = energies[i];
		const double denerg = energies[i+1]-energies[i];
		acc +=kernel[i]*denerg;
		outputfile<<energ*mu.HalfWidth() + mu.BandCenter() <<" "<<acc <<std::endl;
	}

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
