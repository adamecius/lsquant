// Used for OPENMP functions
#include "chebyshev_solver.hpp"
#include "units.hpp"
#include "recon_grid.hpp"
#include "operator_descriptor.hpp"
#include <eigen3/Eigen/Core>
#include <boost/math/special_functions/bessel.hpp>


using namespace chebyshev;

// ─────────────────────────────────────────────────────────────────────────────
// calculate_cn_bessel
//
// Direct translation of the Fortran subroutine of the same name.
//
// Computes the Chebyshev expansion coefficients for U(T) = exp(-i*H_bar*T)
// where H_bar = (H - bc) / ac  is the adimensionalised Hamiltonian.
//
// Convention (note: argument of Bessel functions is 2*ac*T, NOT ac*T):
//
//   c[0]   =          J_0(2*ac*T) * exp(-i*bc*T)
//   c[n>0] = sqrt(2) * (-i)^n * J_n(2*ac*T) * exp(-i*bc*T)
//
// Arguments:
//   T   : dimensionless time  =  t_phys * ac / hbar
//           (in Fortran: tstep_fs * 1e-15 * a * teV / hbar_SI,
//            or equivalently tstep_fs * H.a / hbar_eVfs)
//   ac  : half-bandwidth (= H.a, energy scale used in adimensionalisation)
//   bc  : band centre  (= H.b, energy shift used in adimensionalisation)
//
// Returns:
//   c   : coefficient vector; size = npol determined automatically so that
//         |c[npol-1]| < 1e-20  (convergence criterion from the Fortran)
//
// The factor sqrt(2) for n>0 (vs. the common factor of 2) comes from this
// code's normalisation choice for the Chebyshev expansion on [-1,1].
// ─────────────────────────────────────────────────────────────────────────────
	typedef std::complex<double> cdouble;	

std::vector<cdouble> chebyshev::calculate_cn_bessel(double T, double ac, double bc)
{
    std::cout << "\nCalculating the Chebyshev expansion coefficients...\n";


    const double hbar    =  chebyshev::HBAR;// eV*fs (units.hpp)
    const cdouble minus_i  = cdouble(0.0, -1.0);
    const cdouble phase     = std::exp(minus_i * (bc/ac) * T / hbar);   // exp(-i*ac*T)
    const double  bessel_arg =  T * ac / hbar;
    const double  tol       = 1.0e-20;
    

    // ── Determine npol dynamically ────────────────────────────────────────────
    // Find the first n where sqrt(2) * J_n(2*bc*T) < tol.
    // (|phase| = 1 and |(-i)^n| = 1, so only the Bessel value matters.)
    // Start at n=1 (Fortran loop starts at i=1, checking c for i-1=0):
    //   Fortran checks: |(-i)^(i-1) * J_{i-1}(2bc*T) * exp(-iac*T) * sqrt(2)| < 1e-20
    //   For i=1 that's J_0 without the sqrt(2), but the Fortran still
    //   multiplies by sqrt(2) in the convergence check — we match that exactly.
    int npol = 1;
    while (true) {
        cdouble Jn  = 2.0 * phase * boost::math::cyl_bessel_j(npol - 1, bessel_arg);
        double mag = std::abs(Jn);   // |(-i)^n| = 1, |phase| = 1
        

	if (mag < tol) break;
        ++npol;
    }

    std::cout << "  npol = " << npol << "\n";

    // ── Fill coefficient vector ───────────────────────────────────────────────
    std::vector<cdouble> c(npol);

    //c(1) = CDEXP((0.d0,-1.d0)*bc/ac*T) * DBESJ0(T)
    // n=0: c[0] = J_0(2*ac*T) * exp(-i*ac*T)   (no 2.0 factor)
    c[0] =  boost::math::cyl_bessel_j(0, bessel_arg) * phase;

    // n>0: c[n] = 2.0 * (-i)^n * J_n(2*ac*T) * exp(-i*bc/ac*T)
    cdouble im_pow = minus_i;   // (-i)^1
    for (int n = 1; n < npol; ++n) {
        double Jn = boost::math::cyl_bessel_j(n, bessel_arg);
        c[n]      = 2.0 * im_pow * Jn * phase;
        im_pow   *= minus_i;
    }

    return c;
}




int chebyshev::MeanSquareDisplacement_DeviceSim(chebyshev::MomentsTD & H, int nMom_DOS, int nT, double tstep, qstates::generator& gen  )
{
  // ── 1. Get Hamiltonian and simulation parameters ──────────────────────────

    const int    N       = H.SystemSize();
    std::vector<cdouble> c = calculate_cn_bessel(tstep, H.HalfWidth(), H.BandCenter());

    
    int    nMom_tEvol    = c.size();     // nrecurs in Fortran
   

    // ── 2. Compute Chebyshev coefficients for U(T) ───────────────────────────
    // Fortran: CALL calculate_cn_bessel(tstep*1e-15*a*teV/hbar, a, b, npol, c)
    // The dimensionless time argument is  T_bar = tstep / hbar = 1/eV, where tstep is in femtoseconds (hbar is in eV . fs)
    // (H.a already absorbs teV since Adimensionalize was called with a*teV)
    

    // npol is determined by when J_n(T_bar) becomes negligibly small.
    // A safe upper bound: npol ~ T_bar + 10*cbrt(T_bar) (standard KPM rule).
    // Here we let nMom serve as npol, consistent with the Fortran usage.
    
    
    std::cout << "\nComputed " << nMom_tEvol
              << " Chebyshev coefficients for timestep = " << tstep << "\n";
    
    
    // ── 3. Allocate wave-function vectors ─────────────────────────────────────
    std::vector<cdouble> psi(N),   xupsi(N, cdouble{0,0}),
                                   yupsi(N, cdouble{0,0});

    // ── 4. Initialise random state ────────────────────────────────────────────
    // gen.random_state fills psi with a normalised random vector.

    //Initialize the Random Phase vector used for the Trace approximation
    gen.SystemSize(N);	
    gen.getQuantumState();
    psi = gen.State();


    auto mu_0 = chebyshev::MomentosDelta(H, psi, nMom_DOS);

    // ── 5. Write moments to file ─────────────────────────────────────────────
    // Filename: mu_01.txt  
    std::ofstream f("mu_01.txt");
    if (!f) throw std::runtime_error(std::string("Cannot open mu_01.txt") );

    // Moments (1-based index in file to match Fortran output)
    for (int i = 0; i < nMom_DOS; ++i) {
        f << (i+1)
          << " " << std::real(mu_0[i]) << " " << std::imag(mu_0[i])
	  << "\n";
    }

    f.close();
    std::cout << "  Written mu_01.txt   \n";


    
      
    double inorm2 = 0.0;
    for (int i = 0; i < N; ++i) inorm2 += std::norm(psi[i]);
    std::cout << "\n  initial psi |U(t)|  = " << std::sqrt(inorm2) << "\n";

		
    // ── 6. Time-evolution loop ────────────────────────────────────────────────
    for (int it = 1; it <= nT; ++it) {
        std::cout << "\n\nTime step = " << it << "\n";

        // Evolve psi, xupsi, yupsi by one step T:
        //   psi   →  U(T)|psi>
        //   xupsi →  U(T)|xupsi> + [X, U(T)]|psi>
        //   yupsi →  U(T)|yupsi> + [Y, U(T)]|psi>
        evolution(H, c, psi, xupsi, yupsi);

        // Compute Chebyshev moments of ||xupsi||² and ||yupsi||² and write
        // them to  muxy_NNNNN.txt  (KPM input for post-processing dX²/dY²)
        get_dXY2(H, xupsi, yupsi, nMom_DOS, it);
    }
  return 0;
}





std::vector<cdouble> chebyshev::MomentosDelta(chebyshev::MomentsTD &  H,
                                          const std::vector<cdouble>&  psi,
                                          int                   nMom)
{
  const int N = H.SystemSize();

    std::vector<cdouble> mu(nMom + 1, cdouble{0.0, 0.0});
    std::vector<cdouble> r_n(N), r_nm1(N);

    // ── mu[0] = <psi|psi> = 1 ────────────────────────────────────────────────
    mu[0] = cdouble{1.0, 0.0};

    // ── Build r_nm1 = |psi>,  r_n = H_bar|psi> ───────────────────────────────
    // and collect mu[1] = <psi|H_bar|psi>,  mu[2] = 2<r_n|r_n> - mu[0]


    for (int i = 0; i < N; ++i) r_nm1[i] = psi[i];

    H.Hamiltonian().Multiply ( 1.0, r_nm1.data(), 0.0, r_n.data());

    mu[1] = linalg::vdot( psi, r_n ) - mu[1];
    mu[2] = 2.0 * linalg::vdot( r_n, r_n ) - mu[0];

    // ── Main recursion: n = 3 .. nMom/2+1  (0-based: n-1 = 2 .. nMom/2) ─────
    for (int n = 3; n <= nMom / 2 + 1; ++n) {
        // One Chebyshev step + moment accumulation
      
        H.Hamiltonian().Multiply ( 2.0, r_nm1.data(), 0.0, r_n.data());
        mu[2*n - 3] = linalg::vdot( psi, r_n ) - mu[1];
	mu[2*n - 2] = 2.0 * linalg::vdot( r_n, r_n ) - mu[0];

        std::swap(r_n, r_nm1);
    }

    return mu;
}







// ─────────────────────────────────────────────────────────────────────────────
// get_dXY2
//
// Computes the mean-square displacement along x and y by running the
// Chebyshev recursion on the commutator vectors xupsi and yupsi.
//
// Workflow:
//   1. Normalise xupsi and yupsi (save norms for later restoration).
//   2. Run MomentosDelta on each to get mux and muy.
//   3. Restore xupsi and yupsi to their original (un-normalised) values.
//   4. Write norms + moments to file  muxy_NNNNN.txt
//
// Arguments:
//   H      : Graphene_NNN (adimensionalised)
//   xupsi  : [X, U(T)]|psi>  vector (modified in place, restored on return)
//   yupsi  : [Y, U(T)]|psi>  vector (modified in place, restored on return)
//   nMom   : number of Chebyshev moments
//   it     : time-step index (used for filename)
// ─────────────────────────────────────────────────────────────────────────────


void chebyshev::get_dXY2(chebyshev::MomentsTD & H,
                     std::vector<cdouble> &            xupsi,
                     std::vector<cdouble> &            yupsi,
                     int                 nMom,
                     int                 it)
{
  const int N = H.SystemSize();

    // ── 1. Compute norms and normalise ───────────────────────────────────────
    double normdx = 0.0, normdy = 0.0;
    for (int i = 0; i < N; ++i) {
        normdx += std::norm(xupsi[i]);
        normdy += std::norm(yupsi[i]);
    }
    normdx = std::sqrt(normdx);
    normdy = std::sqrt(normdy);

    for (int i = 0; i < N; ++i) {
        xupsi[i] /= normdx;
        yupsi[i] /= normdy;
    }

    // ── 2. Run Chebyshev recursion ───────────────────────────────────────────
    std::cout << "\nRunning Chebyshev recursion (dX2) ...\n";

    auto mux = chebyshev::MomentosDelta(H, xupsi, nMom);
    auto muy = chebyshev::MomentosDelta(H, yupsi, nMom);

    std::cout << "...done\n\n";

    // ── 3. Restore original (un-normalised) vectors ──────────────────────────
    for (int i = 0; i < N; ++i) {
        xupsi[i] *= normdx;
        yupsi[i] *= normdy;
    }



    
    // ── 4. Write moments to file ─────────────────────────────────────────────
    // Filename: muxy_NNNNN.txt  (zero-padded 5-digit time index)
    char buf[32];
    std::snprintf(buf, sizeof(buf), "muxy_%05d.txt", it);
    std::ofstream f(buf);
    if (!f) throw std::runtime_error(std::string("Cannot open ") + buf);

    // Header: norms (mirrors Fortran: WRITE(106,*) "#",normdx,normdy)
    f << "# " << normdx << " " << normdy << "\n";

    // Moments (1-based index in file to match Fortran output)
    for (int i = 0; i < nMom; ++i) {
        f << (i+1)
          << " " << std::real(mux[i]) << " " << std::imag(mux[i])
          << " " << std::real(muy[i]) << " " << std::imag(muy[i])
          << "\n";
    }

    f.close();
    std::cout << "  Written " << buf << "\n";
}







// ─────────────────────────────────────────────────────────────────────────────
// evolution
//
// Evolves the wave packet one time step and tracks position via commutators.
// All three vectors are updated in-place:
//
//   psi    →  U(T)|psi>
//   xupsi  →  U(T)|xupsi> + [X, U(T)]|psi>
//   yupsi  →  U(T)|yupsi> + [Y, U(T)]|psi>
//
// PART 1 — U(T)|xupsi>, U(T)|yupsi>:
//   Standard Chebyshev sum using H.H_ket() and H.update_cheb().
//
// PART 2 — [X/Y, U(T)]|psi> and U(T)|psi>:
//   Coupled alpha/beta recursion:
//     alpha_{n+1} = 2*H_bar*alpha_n - alpha_{n-1}
//     beta_{n+1}  = 2*H_bar*beta_n + 2*[X,H_bar]*alpha_n - beta_{n-1}
//   where [X,H_bar]*alpha_n is computed cleanly via H.XYH_ket().
// ─────────────────────────────────────────────────────────────────────────────


void chebyshev::evolution(chebyshev::MomentsTD &         H,
               std::vector<cdouble>& c,
               std::vector<cdouble>&                    psi,
               std::vector<cdouble>&                    xupsi,
               std::vector<cdouble>&                    yupsi)
{
  const int N    = H.SystemSize();
    const int npol = static_cast<int>(c.size());

    

    std::vector<cdouble> pnpsi(N),pn1psi(N), vel_tmp(N),
                         xpnpsi(N),xpn1psi = xupsi,
                         ypnpsi(N),ypn1psi = yupsi;

    
    cdouble xtemp,ytemp,ctemp, cnum, xcnum, ycnum; 

    // ═══════════════════════════════════════════════════════════════════════════
    // PART 1:  U(T)|xupsi>  and  U(T)|yupsi>
    // ═══════════════════════════════════════════════════════════════════════════
    {
      
      H.Hamiltonian().Multiply( 1.0, xpn1psi.data(), 0.0, xpnpsi.data() );
      H.Hamiltonian().Multiply(1.0,  ypn1psi.data(), 0.0, ypnpsi.data() );

  // n=0:  T_0 = I
        for (int i = 0; i < N; ++i) {
            xupsi[i] = c[0] * xupsi[i];
            yupsi[i] = c[0] * yupsi[i];
        }
	

	
        for (int n = 1; n < npol; ++n) {
            for (int i = 0; i < N; ++i) {
                xupsi[i] += c[n] * xpnpsi[i];
                yupsi[i] += c[n] * ypnpsi[i];
            }

	    // T_{n+1} = 2*H_bar*T_n - T_{n-1}
	    H.Hamiltonian().Multiply(2.0, xpnpsi.data(), -1.0, xpn1psi.data());
	    H.Hamiltonian().Multiply(2.0, ypnpsi.data(), -1.0, ypn1psi.data());
	    
            std::swap(xpnpsi, xpn1psi);
            std::swap(ypnpsi, ypn1psi);

        }
    }



    

    // ═══════════════════════════════════════════════════════════════════════════
    // PART 2:  U(T)|psi>  and  [X/Y, U(T)]|psi>
    // ═══════════════════════════════════════════════════════════════════════════

    // Initialise:
    //   app = alpha_{n-1} = T_0|psi> = |psi>
    //   ap  = alpha_n     = T_1|psi> = H_bar|psi>
    //   xbpp = beta_x_{n-1} = 0
    //   xbp  = beta_x_n     = [X,H_bar]|psi>   (first non-trivial beta)
    pn1psi = psi;
    
    std::fill(xpn1psi.begin(), xpn1psi.end(), cdouble{0,0});
    std::fill(ypn1psi.begin(), ypn1psi.end(), cdouble{0,0});
    

    H.VX().Multiply( pn1psi.data(), xpnpsi.data());
    H.VY().Multiply( pn1psi.data(), ypnpsi.data());   

    H.Hamiltonian().Multiply( pn1psi.data(), pnpsi.data());


    for (int i = 0; i < N; ++i) 
        psi[i]    = c[0] * psi[i];


    //2(H*xpnpsi+vx*pnpsi)-xpn1psi
    for (int n = 1; n < npol; ++n) {
        // ap  = T_n|psi>   xbp = beta_x_n
        for (int i = 0; i < N; ++i) {
            psi[i]   += c[n] * pnpsi[i];
            xupsi[i] += c[n] * xpnpsi[i];
            yupsi[i] += c[n] * ypnpsi[i];
        }

	H.Hamiltonian().Multiply(2.0, pnpsi.data(), -1.0, pn1psi.data());
	H.Hamiltonian().Multiply(2.0, xpnpsi.data(), -1.0, xpn1psi.data());	
	H.Hamiltonian().Multiply(2.0, ypnpsi.data(), -1.0, ypn1psi.data());

	
	H.VX().Multiply (pnpsi.data(), vel_tmp.data());
	for (int i = 0; i < N; ++i) 
           xpn1psi[i]   +=  2.0 * vel_tmp[i];


	H.VY().Multiply(pnpsi.data(), vel_tmp.data());
	for (int i = 0; i < N; ++i) 
           ypn1psi[i]   +=  2.0 * vel_tmp[i];

	
	
 
	std::swap(pnpsi, pn1psi);
        std::swap(xpnpsi, xpn1psi);
        std::swap(ypnpsi, ypn1psi);

    }

    // Unitarity check

    double norm2=0;
    for (int i = 0; i < N; ++i) norm2 += std::norm(psi[i]);
    std::cout << "\n  |U(t)|  = " << std::sqrt(norm2) << "\n";

}













std::array<double,2> utility::SpectralBounds( SparseMatrixType& HAM, const std::string& descriptor_path)
{
	// DESIGN -- spectral-bounds resolution, in priority order:
	//   1. BOUNDS file ("lo hi" in the run dir) -- explicit per-run override; wins if present.
	//   2. operator .desc sidecar (descriptor_path) -- the PRODUCER's authoritative (a,b) from
	//      Lanczos (wannier2sparse --bounds). The Phase-2 path: bounds come from the descriptor,
	//      not from a hack. Used when no BOUNDS file is supplied.
	//   3. Gershgorin enclosure + bounds_pad -- safe (loose) fallback when neither is available.
	// A BOUNDS file gives tight, exact bounds and is REQUIRED for analytic work (loose bounds
	// distort the moment matrix; thesis Fig. 4.2). We must NOT fall back to an arbitrary window
	// (the old -100..100 silently rescaled every spectrum into a tiny sub-interval).
	double highest=0, lowest=0;
	double desc_a=0, desc_b=0;
	std::ifstream bounds_file( "BOUNDS" );
	if( bounds_file.is_open() )
	{
		bounds_file>>lowest>>highest;
		bounds_file.close();
	}
	else if( !descriptor_path.empty() && lsquant::descriptor_bounds(descriptor_path, desc_a, desc_b) )
	{
		lowest  = desc_a;
		highest = desc_b;
		std::cerr << "[SpectralBounds] bounds from descriptor " << descriptor_path
		          << " -> [" << lowest << ", " << highest << "] (producer Lanczos estimate)."
		          << std::endl;
	}
	else
	{
		const double PAD = chebyshev::safety_factors().bounds_pad; // Gershgorin enclosure margin (recon_grid.hpp)
		auto& mat = HAM.eigen_matrix();
		double rho = 0.0;        // max absolute row sum = ||H||_inf >= spectral radius
		for (int r = 0; r < mat.outerSize(); ++r)
		{
			double rowsum = 0.0;
			for (Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType>::InnerIterator it(mat, r); it; ++it)
				rowsum += std::abs(it.value());
			if (rowsum > rho) rho = rowsum;
		}
		if (rho <= 0.0) rho = 1.0; // degenerate guard
		highest =  rho*(1.0+PAD);
		lowest  = -rho*(1.0+PAD);
		std::cerr << "[SpectralBounds] no BOUNDS file: using Gershgorin estimate |E|<="
		          << rho << " with " << (int)(PAD*100) << "% pad -> [" << lowest << ", "
		          << highest << "]. Supply a BOUNDS file ('lo hi') for tight/exact bounds."
		          << std::endl;
	}
	return { lowest, highest};
};


int chebyshev::ComputeMomTable( chebyshev::Vectors &chevVecL, chebyshev::Vectors& chevVecR ,  chebyshev::Moments2D& sub)
{
	const size_t maxMR = sub.HighestMomentNumber(0);
	const size_t maxML = sub.HighestMomentNumber(1);
		
//	mkl_set_num_threads_local(1); 
	for( auto m0 = 0; m0 < maxML; m0++)
	{
	  //#pragma omp parallel for default(none) shared(chevVecL,chevVecR,m0,sub)
	  for( auto m1 = 0; m1 < maxMR; m1++){

 		sub(m0,m1) = linalg::vdot( chevVecL.Vector(m0) , chevVecR.Vector(m1) );
	  }
	}
//	mkl_set_num_threads_local(nthreads); 

return 0;
};


int chebyshev::CorrelationExpansionMoments(	const vector_t& PhiR, const vector_t& PhiL,
											SparseMatrixType &OPL,
											SparseMatrixType &OPR,  
											chebyshev::Vectors &chevVecL,
											chebyshev::Vectors &chevVecR,
											chebyshev::Moments2D &chebMoms
											)
{
	const size_t numRVecs = chevVecR.NumberOfVectors();
	const size_t numLVecs = chevVecL.NumberOfVectors();
	const size_t NumMomsR = chevVecR.HighestMomentNumber();
	const size_t NumMomsL = chevVecL.HighestMomentNumber();
	const size_t momvecSize = (size_t)( numRVecs*numLVecs );

	auto start = std::chrono::high_resolution_clock::now();
	chebyshev::Moments2D sub(numLVecs,numRVecs);

	
	chevVecL.SetInitVectors( OPL, PhiL );
	for(int  mL = 0 ; mL <  NumMomsL ; mL+=numLVecs )
	{
		chevVecL.IterateAll();

		chevVecR.SetInitVectors( PhiR );
		for(int  mR = 0 ; mR <  NumMomsR ; mR+=numRVecs)
		{
			chevVecR.IterateAll();
			chevVecR.Multiply( OPR );
			chebyshev::ComputeMomTable(chevVecL,chevVecR, sub );		
			chebMoms.AddSubMatrix(sub,mL,mR);
		}
	}	   
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Calculation of the moments in sequential solver took: " << elapsed.count() << " seconds\n\n";

	return 0;
};



int chebyshev::CorrelationExpansionMoments( SparseMatrixType &OPL, SparseMatrixType &OPR,  chebyshev::Moments2D &chebMoms, qstates::generator& gen )
{
	//Set Batch Behavior
	int batchSize;
    if( !chebyshev::GetBatchSize(batchSize)  )
    {
		batchSize = chebMoms.HighestMomentNumber(0);
		std::cout<<"Using default BATCH_SIZE = "<<batchSize<<std::endl;
	}
	else
		std::cout<<"Using user's BATCH_SIZE = "<<batchSize<<std::endl;


	int 
	or_mom0=chebMoms.HighestMomentNumber(0),
	or_mom1=chebMoms.HighestMomentNumber(1);
	if( or_mom0%batchSize != 0 ||
		or_mom1%batchSize != 0 )
	{
		std::cout<<"\nWARNING: This implementation need to the batchsize to be a multiple of ";
		std::cout<<or_mom0<<" and "<<or_mom1<<std::endl;
		std::cout<<"The moments requiremenets will be increased to ";
		
		int mom0= ( (or_mom0+batchSize-1)/batchSize )*batchSize;
		int mom1= ( (or_mom1+batchSize-1)/batchSize )*batchSize;
		std::cout<<mom0<<"x"<<mom1<<std::endl;
		std::cout<<" which is suboptimal "<<std::endl<<std::endl;
		chebMoms = chebyshev::Moments2D( chebMoms, mom0,mom1);
	}

	
	const int DIM  = chebMoms.SystemSize(); 
	
	chebyshev::Vectors 
	chevVecL( chebMoms ), chevVecR( chebMoms );

	//Allocate the memory
	
	chevVecL.SetNumberOfVectors( batchSize );
	chevVecR.SetNumberOfVectors( batchSize );
	printf("Chebyshev::CorrelationExpansionMoments will used %f GB\n", chevVecL.MemoryConsumptionInGB() + chevVecR.MemoryConsumptionInGB() );

	//This operation is memory intensive
	std::cout<<"Initializing chevVecL"<<std::endl;
	chevVecL.CreateVectorSet( );
	std::cout<<"Initialize chevVecR"<<std::endl;
	chevVecR.CreateVectorSet( );

	gen.SystemSize(DIM);
	while( gen.getQuantumState() )
	{
		std::cout<<"Computing with ID: "<<gen.count<<" states" <<std::endl;

		//SELECT RUNNING TYPE
		chebyshev::CorrelationExpansionMoments(gen.State(),gen.State(), OPL, OPR, chevVecL,chevVecR, chebMoms);
	}

	//Fix the scaling of the moments
    const int NumMomsL = chebMoms.HighestMomentNumber(0);
    const int NumMomsR = chebMoms.HighestMomentNumber(1);
	for (int mL = 0 ; mL < NumMomsL; mL++)				  
	for (int mR = mL; mR < NumMomsR; mR++)
	{
		double scal=4.0/gen.NumberOfStates();
		if( mL==0) scal*=0.5;
		if( mR==0) scal*=0.5;

		const value_t tmp = scal*( chebMoms(mL,mR) + std::conj(chebMoms(mR,mL)) )	/2.0;
		chebMoms(mL,mR)= tmp;
		chebMoms(mR,mL)= std::conj(tmp);
	}
	
	
	if( NumMomsL != or_mom0 or NumMomsR != or_mom1 )
	{
		std::cout<<"Changing from "<<NumMomsL<<"x"<<NumMomsR<<" to "<<or_mom0<<"x"<<or_mom1<<std::endl;
		chebMoms.MomentNumber(or_mom0,or_mom1);
	}
	 
return 0;
};



int chebyshev::SpectralMoments( SparseMatrixType &OP,  chebyshev::Moments1D &chebMoms, qstates::generator& gen )
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();

	gen.SystemSize(Dim);
	while( gen.getQuantumState() )
	{		
		auto Phi = gen.State();
		//Set the evolved vector as initial vector of the chebyshev iterations
		if (OP.isIdentity() )
			chebMoms.SetInitVectors( Phi );
		else
			chebMoms.SetInitVectors( OP,Phi );
			
		for(int m = 0 ; m < NumMoms ; m++ )
		{
			double scal=2.0/gen.NumberOfStates();
			if( m==0) scal*=0.5;
			chebMoms(m) += scal*linalg::vdot( Phi, chebMoms.Chebyshev0() ) ;
			chebMoms.Iterate();
		}
	}
	return 0;
};





int chebyshev::SpectralMoments_S_inv( SparseMatrixType &S,  chebyshev::Moments1D &chebMoms, qstates::generator& gen, std::array<double,2> spectral_bounds )
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();


	
	Eigen::VectorXcd acc_result = Eigen::VectorXcd::Zero(Dim); 
	Eigen::MatrixXd convergence = Eigen::MatrixXd::Zero(Dim,2);
	
	double x0 = 2.0 * ( 0.0 - spectral_bounds[0] ) / (spectral_bounds[1] - spectral_bounds[0] )  - 1.0 ;


	
	  
	gen.SystemSize(Dim);
	while( gen.getQuantumState() )
	{		
		auto Phi = gen.State();
		//Set the evolved vector as initial vector of the chebyshev iterations
		chebMoms.SetInitVectors( S,Phi );
			
		for(int m = 0 ; m < NumMoms ; m++ )
		{
			double scal=2.0/gen.NumberOfStates();
			if( m==0) scal*=0.5;
			chebMoms(m) += scal*linalg::vdot( Phi, chebMoms.Chebyshev0() ) ;
			
			chebMoms.Iterate();

			if( m == 0 )
			  acc_result(0) = greenR_chebF(x0, 0) * chebMoms.JacksonKernel(0,  NumMoms ) * chebMoms(0); 
			else
			  acc_result(m) = acc_result(m - 1) + greenR_chebF(x0, m) * chebMoms.JacksonKernel(m,  NumMoms ) * chebMoms(m); 

			std::cout<<m<<":   "<<std::abs(acc_result(m))<< "    " <<linalg::vdot( Phi, Phi )<<"   "<<x0<<std::endl;

			convergence(m,0)=m;
			convergence(m,1)=std::abs(std::abs(acc_result(m))-1.0);
			
		}
	}

       std::ofstream dataP;
       dataP.open("S_convergence.dat");


       dataP<< convergence<<std::endl;
       dataP.close();

       return 0;
};


int chebyshev::ComputeDeltaPhi2( SparseMatrixType &OP,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen ,double Energy)
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();

	gen.SystemSize(Dim);
	while( gen.getQuantumState() )
	{		
		auto Phi = gen.State();
	  	chebMoms.resetCurrentTm();

		//Set the evolved vector as initial vector of the chebyshev iterations
		if (OP.isIdentity() )
			chebMoms.SetInitVectors( Phi );
		else
			chebMoms.SetInitVectors( OP,Phi );
		
			
		for(int m = 0 ; m < NumMoms ; m++ )
		{
  			double g_D_m=chebMoms.JacksonKernel(m,  NumMoms );
			double scal=2.0/gen.NumberOfStates();
			if( m==0) scal*=0.5;
			linalg::axpy( scal*delta_chebF(Energy,m)*g_D_m /chebMoms.HalfWidth(), chebMoms.Chebyshev0(), chebMoms.DeltaPhi());
			//chebMoms(m) += scal*linalg::vdot( Phi, chebMoms.Chebyshev0() ) ;
			chebMoms.Iterate();
		}
	}
	return 0;
};

int chebyshev::SpectralMoments_nonOrth( SparseMatrixType &OP,  chebyshev::Moments1D_nonOrth &chebMoms, qstates::generator& gen )
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();


	gen.SystemSize(Dim);
	while( gen.getQuantumState() )
	{		
		auto Phi = gen.State();
		//Set the evolved vector as initial vector of the chebyshev iterations
		if (OP.isIdentity() )
			chebMoms.SetInitVectors_nonOrthogonal( Phi );
		else
			chebMoms.SetInitVectors_nonOrthogonal( OP,Phi );

		double scal=2.0/gen.NumberOfStates();
		chebMoms(0) += 0.5*scal*linalg::vdot( Phi, chebMoms.Chebyshev0() ) ;
		
		
		chebMoms(1) +=scal*linalg::vdot( Phi, chebMoms.Chebyshev1() ) ;
		
		
		for(int m = 2 ; m < NumMoms ; m++ )
		{
		        chebMoms.Iterate_nonOrthogonal();
			chebMoms(m) += scal*linalg::vdot( Phi, chebMoms.Chebyshev1() ) ;
			

		}
	}
	return 0;
};

int chebyshev::SpectralMoments_nonOrth_test( SparseMatrixType &OP, SparseMatrixType &orth_Ham, chebyshev::Moments1D_nonOrth &chebMoms, qstates::generator& gen )
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();

	std::vector<double> errors(NumMoms,0);

	gen.SystemSize(Dim);
	while( gen.getQuantumState() )
	{		
		auto Phi = gen.State();
		//Set the evolved vector as initial vector of the chebyshev iterations
		if (OP.isIdentity() )
			chebMoms.SetInitVectors_nonOrthogonal( Phi );
		else
			chebMoms.SetInitVectors_nonOrthogonal( OP,Phi );

		double scal=2.0/gen.NumberOfStates();
		chebMoms(0) += 0.5*scal*linalg::vdot( Phi, chebMoms.Chebyshev0() ) ;
		
		
		chebMoms(1) +=scal*linalg::vdot( Phi, chebMoms.Chebyshev1() ) ;
		
		
		for(int m = 2 ; m < NumMoms ; m++ )
		{
		     std::cout<<m+1<<"/"<<NumMoms<<std::endl;
		     errors[m] = chebMoms.Iterate_nonOrthogonal_test(orth_Ham);
		     chebMoms(m) += scal*linalg::vdot( Phi, chebMoms.Chebyshev1() ) ;
		     std::cout<<errors[m]<<std::endl<<std::endl<<std::endl;
		}
	}

	std::ofstream data;
        data.open("orthogonalization_errors.dat");

	for(int m = 0; m < NumMoms; m++)  
          data<< m <<"  "<< errors[m] <<std::endl;
	  
        data.close();
  

	return 0;
};

int chebyshev::ComputeDeltaPhi( SparseMatrixType &OP,  chebyshev::MomentsLocal &chebMoms, qstates::generator& gen)
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();
	const auto Norb= chebMoms.NumberOfOrbitals();

	gen.SystemSize(Dim);
	while( gen.getQuantumState() )
	{		
		auto Phi = gen.State();

		//Set the evolved vector as initial vector of the chebyshev iterations
		if (OP.isIdentity() )
			chebMoms.SetInitVectors( Phi );
		else
			chebMoms.SetInitVectors( OP,Phi );
		
			
		for(size_t m = 0 ; m < NumMoms ; m++ )
		{
			//std::cout<<m<<" "<<std::endl;//chebMoms(1073,1480+m)<<std::endl;
			double scal=2.0/gen.NumberOfStates();
			if( m==0) scal*=0.5;
			for (size_t n = 0; n<Norb;n++)
			{
				chebMoms(m,n) += scal*conj(Phi[n])*( chebMoms.Chebyshev0()[n]);
			}
			//chebMoms(m) += scal*linalg::vdot( Phi, chebMoms.Chebyshev0() ) ;
			chebMoms.Iterate();
		}
	}
	return 0;
};

int chebyshev::TimeDependentCorrelations(SparseMatrixType &OPL, SparseMatrixType &OPR,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen  )
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();
	const auto NumTimes= chebMoms.MaxTimeStep();
	
	//Initialize the Random Phase vector used for the Trace approximation
	gen.SystemSize(Dim);	
	while( gen.getQuantumState() )
	{
		chebMoms.ResetTime();

		auto PhiR =gen.State();
		auto PhiL =gen.State();
		 
		//Multiply right operator its operator
		OPL.Multiply(PhiR,PhiL); //Defines <Phi| OPL 
		
		//Evolve state vector from t=0 to Tmax
		while ( chebMoms.CurrentTimeStep() !=  chebMoms.MaxTimeStep()  )
		{
			const auto n = chebMoms.CurrentTimeStep();

			//Set the evolved vector as initial vector of the chebyshev iterations
			chebMoms.SetInitVectors( OPR , PhiR );

			for(int m = 0 ; m < NumMoms ; m++ )
			{
				double scal=2.0/gen.NumberOfStates();
				if( m==0) scal*=0.5;
				chebMoms(m,n) += scal*linalg::vdot( PhiL, chebMoms.Chebyshev0() ) ;
				chebMoms.Iterate();
			}
			
			chebMoms.IncreaseTimeStep();
			//evolve PhiL ---> PhiLt , PhiR ---> PhiRt 
			chebMoms.Evolve(PhiL) ;
			chebMoms.Evolve(PhiR) ;
		}
	
	}
	
	return 0;
};

int chebyshev::MeanSquareDisplacement(chebyshev::MomentsTD &chebMoms, qstates::generator& gen  )
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();
	const auto NumTimes= chebMoms.MaxTimeStep();

	//Initialize the Random Phase vector used for the Trace approximation
	gen.SystemSize(Dim);	
	while( gen.getQuantumState() )
	{
		chebMoms.ResetTime();
		chebMoms.ResetGaussian();

		auto PhiR = gen.State();
		auto PhiL = PhiR;
		 
		//Multiply right operator its operator
		//{
		//auto PhiT = PhiR;
		// OP.Multiply(PhiL,tempPhiL); //Defines <Phi| OP
		//OPPRJ.Multiply(PhiR, PhiL); //Defines <Phi| OPPRJ
		//linalg::copy(PhiL, PhiR);
		//}
		
		linalg::copy(PhiR, PhiL);
		//Evolve state vector from t=0 to Tmax
		while ( chebMoms.CurrentTimeStep() !=  chebMoms.MaxTimeStep()  )
		{
			const auto n = chebMoms.CurrentTimeStep();

			//Set the evolved vector as initial vector of the chebyshev iterations
			//chebMoms.SetInitVectors( PhiR );
			chebMoms.SetInitVectors( chebMoms.MSDWF() );

			for(int m = 0 ; m < NumMoms ; m++ )
			{
				double scal=2.0/gen.NumberOfStates();
				if( m==0) scal*=0.5;
				//OP.Multiply( chebMoms.Chebyshev0(), PhiT );
				//linalg::copy(chebMoms.Chebyshev0(),PhiT);
				//PhiT=chebMoms.Chebyshev0();
				chebMoms(m,n) += scal*linalg::vdot( chebMoms.MSDWF(), chebMoms.Chebyshev0() ) ;
			//	chebMoms(m,n) += scal*linalg::vdot( PhiL, chebMoms.Chebyshev0() ) ;
				//std::cout<<chebMoms(m,n)<<std::endl;
				chebMoms.Iterate();
			}
			
			//SaveWFatEachTimeSpace
			std::string WFfilename= "WavefunctionatT"+std::to_string(chebMoms.CurrentTimeStep())+".dat";
  			//ofstream outputfile(WFfilename.c_str());
  			//for ( auto wfcoef : PhiL )
    			//	outputfile << wfcoef.real() << " " << wfcoef.imag() << std::endl;
  			//outputfile.close();
			
			
			chebMoms.IncreaseTimeStep();
			//evolve PhiL ---> PhiLt , PhiR ---> PhiRt 
			//chebMoms.MSD_Evolve(PhiL) ;
			chebMoms.MSD_Evolve(PhiR) ;
			linalg::copy(PhiR,PhiL)   ;
			//chebMoms.Evolve(PhiL) ;
			//chebMoms.Evolve(PhiR) ;
		}
	
	}
	
	return 0;
};



int chebyshev::TimeEvolvedProjectedOperatorWF(SparseMatrixType &OP, SparseMatrixType &OPRJ,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen  )
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();
	const auto NumTimes= chebMoms.MaxTimeStep();
	
	//Initialize the Random Phase vector used for the Trace approximation
	gen.SystemSize(Dim);	
	while( gen.getQuantumState() )
	{
		chebMoms.ResetTime();

		auto PhiR =gen.State();
		auto PhiL =gen.State();
		 
		//Multiply right operator its operator
		OPRJ.Multiply(PhiR,PhiL); //Defines <Phi| OPL 
		OP.Multiply(PhiL,PhiL); //Defines <Phi| OPL OPR
		
		//Evolve state vector from t=0 to Tmax
		while ( chebMoms.CurrentTimeStep() !=  chebMoms.MaxTimeStep()  )
		{
			const auto n = chebMoms.CurrentTimeStep();

			//Set the evolved vector as initial vector of the chebyshev iterations
			chebMoms.SetInitVectors( OPRJ , PhiR );

			for(int m = 0 ; m < NumMoms ; m++ )
			{
				double scal=2.0/gen.NumberOfStates();
				if( m==0) scal*=0.5;
				chebMoms(m,n) += scal*linalg::vdot( PhiL, chebMoms.Chebyshev0() ) ;
				chebMoms.Iterate();
			}

			chebMoms.IncreaseTimeStep();
			//evolve PhiL ---> PhiLt , PhiR ---> PhiRt 
			chebMoms.Evolve(PhiL) ;
			chebMoms.Evolve(PhiR) ;
		}
	
	}
	
	return 0;
};


int chebyshev::TimeEvolvedProjectedOperator(SparseMatrixType &OP, SparseMatrixType &OPPRJ,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen  )
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();
	const auto NumTimes= chebMoms.MaxTimeStep();
	
	//Initialize the Random Phase vector used for the Trace approximation
	gen.SystemSize(Dim);	
	while( gen.getQuantumState() )
	{
		chebMoms.ResetTime();

		auto PhiR = gen.State();
		auto PhiL = PhiR;
		 
		//Multiply right operator its operator
		//{
		auto PhiT = PhiR;
		// OP.Multiply(PhiL,tempPhiL); //Defines <Phi| OP
		OPPRJ.Multiply(PhiR, PhiL); //Defines <Phi| OPPRJ
		linalg::copy(PhiL, PhiR);
		//}
		
		//Evolve state vector from t=0 to Tmax
		while ( chebMoms.CurrentTimeStep() !=  chebMoms.MaxTimeStep()  )
		{
			const auto n = chebMoms.CurrentTimeStep();

			//Set the evolved vector as initial vector of the chebyshev iterations
			chebMoms.SetInitVectors( PhiR );

			for(int m = 0 ; m < NumMoms ; m++ )
			{
				double scal=2.0/gen.NumberOfStates();
				if( m==0) scal*=0.5;
				//scal*=2;//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////REMOVE
				OP.Multiply( chebMoms.Chebyshev0(), PhiT );
				chebMoms(m,n) += scal*linalg::vdot( PhiL, PhiT ) ;
				chebMoms.Iterate();
			}
			
			chebMoms.IncreaseTimeStep();
			//evolve PhiL ---> PhiLt , PhiR ---> PhiRt 
			chebMoms.Evolve(PhiL) ;
			chebMoms.Evolve(PhiR) ;
		}
	
	}
	
	return 0;
};







int chebyshev::TimeEvolvedOperator(SparseMatrixType &OP,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen  )
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();
	const auto NumTimes= chebMoms.MaxTimeStep();
	
	//Initialize the Random Phase vector used for the Trace approximation
	gen.SystemSize(Dim);	
	while( gen.getQuantumState() )
	{
		chebMoms.ResetTime();

		auto PhiR = gen.State();
		auto PhiL = PhiR;
		 
		//Multiply right operator its operator
		//{
		auto PhiT = PhiR;
		// OP.Multiply(PhiL,tempPhiL); //Defines <Phi| OP
		linalg::copy(PhiR, PhiL);
		//}

		int t=0;
		//Evolve state vector from t=0 to Tmax
		while ( chebMoms.CurrentTimeStep() !=  chebMoms.MaxTimeStep()  )
		{
			const auto n = chebMoms.CurrentTimeStep();

			//Set the evolved vector as initial vector of the chebyshev iterations
			chebMoms.SetInitVectors( PhiR );

			for(int m = 0 ; m < NumMoms ; m++ )
			{
				double scal=2.0/gen.NumberOfStates();
				if( m==0) scal*=0.5;
				//scal*=2;//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////REMOVE
				OP.Multiply( chebMoms.Chebyshev0(), PhiT );
				chebMoms(m,n) += scal*linalg::vdot( PhiL, PhiT ) ;
				chebMoms.Iterate();
			}

			
			t++;
			cout<< t<<"/"<< chebMoms.MaxTimeStep()<<"  timestep:  " <<std::endl;
			chebMoms.IncreaseTimeStep();
			//evolve PhiL ---> PhiLt , PhiR ---> PhiRt 
			cout<<"  First time evolution:"<<std::endl;
			chebMoms.Evolve(PhiL) ;
			cout<<"  Second time evolution:"<<std::endl;
			chebMoms.Evolve(PhiR) ;

			
		        std::string outputfilename="TimeEvol_currentResult.chebmomTD";	
			chebMoms.saveIn(outputfilename);
			std::cout<<std::endl<<std::endl;

		}
	
	}
	
	return 0;
};




int chebyshev::TimeEvolvedOperatorWithWF(chebyshev::MomentsTD &chebMoms, qstates::generator& gen  )
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();
	const auto NumTimes= chebMoms.MaxTimeStep();
	
	//Initialize the Random Phase vector used for the Trace approximation
	gen.SystemSize(Dim);	
	while( gen.getQuantumState() )
	{
		chebMoms.ResetTime();

		auto PhiR = gen.State();
		auto PhiL = PhiR;
		 
		//Multiply right operator its operator
		//{
		auto PhiT = PhiR;
		// OP.Multiply(PhiL,tempPhiL); //Defines <Phi| OP
		linalg::copy(PhiR, PhiL);
		//}
		
		//Evolve state vector from t=0 to Tmax
		while ( chebMoms.CurrentTimeStep() !=  chebMoms.MaxTimeStep()  )
		{
			const auto n = chebMoms.CurrentTimeStep();

			//Set the evolved vector as initial vector of the chebyshev iterations
			chebMoms.SetInitVectors( PhiR );

			chebMoms.IncreaseTimeStep();
			//SaveWFatEachTimeSpace
			std::string WFfilename= "WavefunctionatT"+std::to_string(chebMoms.CurrentTimeStep())+".dat";
  			typedef std::numeric_limits<double> dbl;
  			ofstream outputfile(WFfilename.c_str());
  			outputfile.precision(dbl::digits10);
  			for ( auto wfcoef : PhiL )
    				outputfile << wfcoef.real() << " " << wfcoef.imag() << std::endl;
  			outputfile.close();
			//evolve PhiL ---> PhiLt , PhiR ---> PhiRt 
			chebMoms.Evolve(PhiL) ;
			chebMoms.Evolve(PhiR) ;
		}
	
	}
	
	return 0;
};
