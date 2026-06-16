// [LEGACY -- NOT BUILT] previous-working opticalConductivityFromChebmom. Reference only.
// Supported path: lsquant <cmd>  |  port + golden-gate before any reuse.

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
#include "optical_conductivity_tools.hpp"


void printHelpMessage();
void printWelcomeMessage();




/*------------------------------

This program computes the current-current correlation function $S_{\alpha \beta}$:

\begin{equation}
S_{\alpha \beta}(\omega) = \int_{-\infty}^{\infty}\frac{d\epsilon}{2\pi}th\frac{\beta \epsilon}{2}Tr[j_\alpha Im G(\epsilon)j_{\beta}G^+(\epsilon+\hbar \omega) + j_\alpha G^-(\epsilon -\hbar \omega)j_\beta Im G(\epsilon)] \label{eq1}
\end{equation}

Which, according to Bastin/1971, is related to the electrical conductivity by:

\begin{equation}
\sigma_{\alpha \beta}(\omega) = \sigma^g_{\alpha \beta}(\omega)+\frac{e^2}{iw}S_{\alpha \beta}(\omega) \label{eq2}
\end{equation}
Where the first term is the gauge current:

\begin{equation}
\sigma^g_{\alpha \beta}(\omega) = -\frac{iNe^2}{m^*\omega}\delta_{\alpha \beta}
\end{equation}


Here, the current-current correlation function is computed through Chebyshev expansions of the involved Greeen's functions. The expansions are reconstructed from the 2D expansion moments in the chebMoms2D file that this code takes as input. \\
 The explicit summation that this code does to compute $S_{\alpha \beta}$ by means of the Chebyshev moments is as follows:

\begin{equation}
    S_{\alpha \beta}(\mu,T, \omega )\frac{4e^2\hbar}{\pi \Omega}\frac{4}{\Delta E^2}\int_{-1}^{1} d\tilde{\epsilon}\frac{f(\tilde{\epsilon})}{(1-\epsilon^2)^2}\sum_{m,n}\Gamma_{nm}(\tilde{\epsilon}, \omega)\mu^{\alpha \beta}_{mn}
\end{equation}

where $\Gamma_{m,n}(\epsilon,\omega)$ are the expansion coefficients:

\begin{equation}
    \Gamma_{m,n}(\epsilon,\omega) = [(\tilde{\epsilon}+\hbar\tilde{\omega} - in \sqrt{1-(\tilde{\epsilon}+\hbar\tilde{\omega}^2})]e^{in\arccos(\tilde{\epsilon}+\hbar\tilde{\omega})}T_m(\tilde{\epsilon}) +  \\
                                    [(\tilde{\epsilon}-\hbar\tilde{\omega} + im \sqrt{1-(\tilde{\epsilon}-\hbar\tilde{\omega}^2})]e^{im\arccos(\tilde{\epsilon}-\hbar\tilde{\omega})}T_n(\tilde{\epsilon})
\end{equation}

The gauge term is not taken into account, which may be innacurate outside gapped regions of the spectrum but yields correct result inside gaps.   



----------------------------- */


int main(int argc, char *argv[])
{	
    // Checking correctness of provided arguments
	if (argc != 4) {
		printHelpMessage();
		return 0;
	} else {
		printWelcomeMessage();
	}

	// Reading broadening parameter for Jackson kernel
	const double broadening  = stod(argv[3]);

	// Check if the system is spinless or not
	const int spinless = true;

	// Reading energy bounds
	std::array<double,2> matterBounds = AC_utils::MatterBounds();

	// Loading the Chebyshev moments from chebMom2D file after checking the file is in place
	chebyshev::MomentsAC mu(argv[1]); 
	if (!mu.Freq(argv[2])) {
		return 0;
	}

	// Rescale moments to [-1,1] (with a safety factor)
	mu.RescaleFreq();

	// Jackson kernel for supressing Gibbs oscilations
	mu.ApplyJacksonKernel(broadening, broadening);

	// Setting the Spin degeneracy factor
	int prefactorSpin = 1;
	if (spinless) {
		prefactorSpin *= 2;
	}

	// Number of energy points in the integration grid
	const int num_div = 2 * mu.HighestMomentNumber() + 1;

	// Physical energy bounds
	const double xboundm = matterBounds[0] * mu.ScaleFactor() + mu.ShiftFactor();
	const double xboundp = matterBounds[1] * mu.ScaleFactor() + mu.ShiftFactor();

	// List of energies energy points for the integration
	std::vector<double> energies(num_div, 0);
	for (int i = 0; i < num_div; i++) {
		energies[i] = xboundm + i * (xboundp - xboundm) / (num_div - 1);
	}

	// Info print
	std::cout << "Computing the Optical Conductivity kernel using "
	          << mu.HighestMomentNumber(0) << " x " << mu.HighestMomentNumber(1)
	          << " moments " << std::endl;

	std::cout << "Integrated over grid (" << xboundm << ", " << xboundp
	          << ") with step size " << energies[1] - energies[0] << std::endl;

	std::cout << "The first 10 moments are:" << std::endl;
	for (int m0 = 0; m0 < 1; m0++)
	for (int m1 = 0; m1 < 10; m1++)
		std::cout << mu(m0, m1).real() << " " << mu(m0, m1).imag() << std::endl;

	// Output file to store conductivity results
	std::string outputName = "OpticalConductivity_" + mu.SystemLabel() + "JACKSON.dat";
	std::cout << "Saving the data in " << outputName << std::endl;

	std::cout << "PARAMETERS: " << mu.SystemSize() << " " << mu.HalfWidth() << std::endl;

	std::ofstream outputfile(outputName.c_str());

	// On the file header: number of frequencies and number of energy points
	outputfile << mu.NumFreq() << ' ' << num_div - 1 << std::endl;

	// Temporary buffer to hold the kernel of the integration
	std::vector<std::complex<double>> kernel(num_div, 0.0);

	// Loop over each frequency
	for (int nomega = 0; nomega < mu.NumFreq(); nomega++) {
		double current_omega = mu.CurrentFreq();

		std::cout << "Computation for Freq " << nomega << std::endl;
		std::cout << "Freq = " << current_omega / mu.ScaleFactor() << std::endl;

		// Compute S_{ab}(\omega) at each energy using the moments and GammaAC kernel
		#pragma omp parallel for
		for (int i = 0; i < num_div; i++) {
			kernel[i] = 0.0;
			const double energ = energies[i];

			for (int m0 = 0; m0 < mu.HighestMomentNumber(0); m0++)
			for (int m1 = 0; m1 < mu.HighestMomentNumber(1); m1++) {
				kernel[i] += GammaAC_chebF(m0, m1, current_omega, energ) * mu(m0, m1);
			}

			// Apply prefactor (system size and rescaling)
			kernel[i] *= mu.SystemSize() * mu.ScaleFactor() * mu.ScaleFactor();
		}

		// Integrate over energy to obtain frequency-dependent conductivity
		std::complex<double> acc = 0;
		for (int i = 0; i < num_div - 1; i++) {
			const double energ  = energies[i];
			const double denerg = energies[i + 1] - energies[i];

			acc += kernel[i] * denerg;

			// Save: omega, energy, Re[S], Im[S]
			outputfile << current_omega / mu.ScaleFactor() << ' '
			           << energ / mu.ScaleFactor() - mu.ShiftFactor() << " "
			           << acc.real() << ' ' << acc.imag() << std::endl;
		}

		// Move to next frequency
		mu.IterateFreq();
	}

	outputfile.close();
	std::cout << "The program finished successfully." << std::endl;

	return 0;
}



	

void printHelpMessage()
{
	std::cout << "The program should be called with the following options: moments_filename Freq_file broadening(meV)" << std::endl
			  << std::endl;
	std::cout << "moments_filename will be used to look for .chebmom2D file" << std::endl;
	std::cout << "Freq_file is the file that stores the Frequencies. (See documentation for format)" << std::endl;
	std::cout << "broadening in (meV) will define the broadening of the delta functions" << std::endl;
};

void printWelcomeMessage()
{
	std::cout << "WELCOME: This program will compute the chebyshev sum of the kubo-bastin formula for non equlibrium properties" << std::endl;
};
