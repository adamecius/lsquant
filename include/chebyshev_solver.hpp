/**
 * @file chebyshev_solver.hpp
 * @brief KPM solver entry points — the Chebyshev recursions that fill the moment tables.
 *
 * Declares the `chebyshev` routines that drive the Kernel Polynomial Method over a
 * @ref SparseMatrixType Hamiltonian and a `qstates::generator` of source states, writing
 * the result into the moment containers of `chebyshev_moments.hpp`:
 *  - spectral / DOS moments (`SpectralMoments`),
 *  - Kubo–Greenwood / non-equilibrium correlations (`CorrelationExpansionMoments`),
 *  - time evolution and time-dependent / projected correlations (`TimeEvolved*`, MSD).
 *
 * `chebyshev::utility::SpectralBounds` supplies the rescaling window: it reads a `BOUNDS`
 * file when present, else estimates a safe Gershgorin enclosure with a 10% pad.
 */
#ifndef CHEBYSHEV_SOLVER
#define CHEBYSHEV_SOLVER

// C & C++ libraries
#include <cassert>   //for assert
#include <array>

#include <vector>    //for std::vector mostly class
#include <numeric>   //for std::accumulate *
#include <algorithm> //for std::max_elem
#include <complex>   ///for std::complex
#include <fstream>   //For ofstream
#include <limits>    //For getting machine precision limit
#include "sparse_matrix.hpp"
#include "chebyshev_moments.hpp"
#include "chebyshev_coefficients.hpp"
#include "linear_algebra.hpp"
#include <omp.h>
#include <chrono>
#include "quantum_states.hpp"
#include "kpm_noneqop.hpp" //Get Batch function
#include "special_functions.hpp"


namespace chebyshev
{
	typedef std::complex<double> value_t;	
	typedef std::vector<value_t> vector_t ;


	namespace utility
	{
		/// @brief Spectral rescaling window [lo,hi]: from a `BOUNDS` file if present,
		///        else a Gershgorin enclosure of the Hamiltonian padded by 10% (safe default).
		std::array<double,2> SpectralBounds( SparseMatrixType& HAM, const std::string& descriptor_path = "");
	};


	int MeanSquareDisplacement_DeviceSim(chebyshev::MomentsTD &, int , int , double , qstates::generator&   );
  
        std::vector<value_t> calculate_cn_bessel(double , double , double );

        void evolution(chebyshev::MomentsTD &, std::vector<value_t>&, std::vector<value_t>&, std::vector<value_t>&, std::vector<value_t>& );      

        void get_dXY2(chebyshev::MomentsTD&, std::vector<value_t>& , std::vector<value_t>&, int , int );

        std::vector<value_t> MomentosDelta(chebyshev::MomentsTD&, const std::vector<value_t>&, int );

  

  
  

	int ComputeMomTable( chebyshev::Vectors &chevVecL, chebyshev::Vectors& chevVecR , chebyshev::Moments2D& sub);

	int CorrelationExpansionMoments( 	const vector_t& PhiR, const vector_t& PhiL,
										SparseMatrixType &OPL,
										SparseMatrixType &OPR,  
										chebyshev::Vectors &chevVecL,
										chebyshev::Vectors &chevVecR,
										chebyshev::Moments2D &chebMoms
										);

	
	/// @brief Fill the 2D Kubo–Greenwood moment table mu_{m,n} = Tr[T_m(H) OPL T_n(H) OPR]
	///        (stochastic trace over @p gen). @p OPL,@p OPR are the operators (e.g. velocities).
	int CorrelationExpansionMoments( SparseMatrixType &OPL, SparseMatrixType &OPR,  chebyshev::Moments2D &chebMoms, qstates::generator& gen );

	/// @brief Time-dependent two-operator correlation moments (velocity autocorrelation, ...).
	int TimeDependentCorrelations( SparseMatrixType &OPL, SparseMatrixType &OPR,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen);

	/// @brief Fill the 1D spectral moment table mu_m = Tr[T_m(H) OP] (DOS when OP is the identity).
	int SpectralMoments(SparseMatrixType &OP,  chebyshev::Moments1D &chebMoms, qstates::generator& gen);

        int SpectralMoments_S_inv(SparseMatrixType &OP,  chebyshev::Moments1D &chebMoms, qstates::generator& gen, std::array<double,2> spectral_bounds);

  
        int SpectralMoments_nonOrth(SparseMatrixType &OP, chebyshev::Moments1D_nonOrth &chebMoms, qstates::generator& gen);
	
        int TimeEvolvedProjectedOperator(SparseMatrixType &OP, SparseMatrixType &OPPRJ,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen  );

        int TimeEvolvedOperator(SparseMatrixType &OP,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen  );

	int MeanSquareDisplacement(chebyshev::MomentsTD &chebMoms, qstates::generator& gen  );

	int TimeEvolvedProjectedOperatorWF(SparseMatrixType &OP, SparseMatrixType &OPRJ,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen  );
    
      	int TimeEvolvedOperatorWithWF(chebyshev::MomentsTD &chebMoms, qstates::generator& gen  );

	int ComputeDeltaPhi( SparseMatrixType &OP,  chebyshev::MomentsLocal &chebMoms, qstates::generator& gen);

	int ComputeDeltaPhi2( SparseMatrixType &OP,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen ,double Energy);

        int SpectralMoments_nonOrth_test(SparseMatrixType &OP, SparseMatrixType &orth_Ham, chebyshev::Moments1D_nonOrth &chebMoms, qstates::generator& gen);

}; // namespace chebyshev

#endif
