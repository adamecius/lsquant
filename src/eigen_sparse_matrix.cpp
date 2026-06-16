#include "sparse_matrix.hpp"
#include "util/timer.hpp"   // LSQ_SCOPED -- time the SpMV hot loop (CLAUDE.md reporting)


const int MKL_SPMVMUL = 1000;

// Below this many nonzeros the SpMV/SpMM run SERIALLY: the OpenMP fork/join cost of a
// parallel region dwarfs the work on a tiny matrix, and small matrices are called in
// tight inner loops (e.g. the exact-trace spin tutorials -> tens of thousands of SpMVs).
// Measured pathology before this guard: a tiny N=128 exact-trace run took 339 s on 256
// threads vs 2.0 s on 1 thread (a 170x regression). Serial and parallel give bit-identical
// results (each row is an independent reduction), so the threshold never changes output.
// Tunable; ~3e4 nnz is well above the per-region overhead and below any DRAM-bound case.
static const long LSQ_SPMV_MIN_PARALLEL_NNZ = 30000;

void SparseMatrixType::ConvertFromCOO(vector<indexType> &rows, vector<indexType> &cols, vector<complex<double> > &vals)
{
	rows_ = vector<indexType>(rows);
	cols_ = vector<indexType>(cols);
	vals_ = vector<complex<double> >(vals);

	std::size_t NNZ = vals_.size();
	

        std::vector<Eigen::Triplet<std::complex<double>,indexType>> triplets;
        triplets.reserve(NNZ);


        for (std::size_t i = 0; i < NNZ; ++i) {
	  triplets.push_back(Eigen::Triplet<std::complex<double>, indexType>(rows_[i], cols_[i], vals_[i]));
        }


        int num_rows = *std::max_element(rows_.begin(), rows_.end()) + 1;
        int num_cols = *std::max_element(cols_.begin(), cols_.end()) + 1;



        Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor, indexType> matrix(num_rows, num_cols);
        matrix.setFromTriplets(triplets.begin(), triplets.end());

	
	setDimensions(num_rows, num_cols);
	std::cout<<"ERROR: Incomplete function for constructing Hamiltonian from COO input matrix. I dont want to be generating extra copies of a huge input."<<std::endl;
        std::terminate();
	//new (&matrix_) Eigen::Map<Eigen::SparseMatrix<complex<double>, Eigen::RowMajor> >(rows_.size(), cols_.size(), NNZ, rows_.data(), cols_.data(), vals_.data());
};


#include <cmath>
#include <numeric> // For std::accumulate
#include <stdexcept> // std::runtime_error (Hermiticity guard)
#include <algorithm> // std::max


double compute_norm(const std::vector<complex<double>>& vec) {
    // Compute the L2 norm (Euclidean norm)
    double sum_of_squares = std::accumulate(vec.begin(), vec.end(), 0.0,
        [](double sum, complex<double> val) {
	  return sum + abs(val) * abs(val); // Add square of the element
        });
    return std::sqrt(sum_of_squares); // Return square root of the sum of squares
}


double compute_norm(const std::vector<indexType>& vec) {
    // Compute the L2 norm (Euclidean norm)
    double sum_of_squares = std::accumulate(vec.begin(), vec.end(), 0.0,
        [](double sum, double val) {
	  return sum + val * val; // Add square of the element
        });
    return std::sqrt(sum_of_squares); // Return square root of the sum of squares
}

void SparseMatrixType::ConvertFromCSR(vector<indexType> &cols, vector<indexType> &rowIndex, vector<complex<double> > &vals)
{
  //============================================================================================================================================//
  //
  // 
  //WARNING: proper implementation involves generating the rows_ cols_ and vals_ internal variables DIRECTLY, not passing them from elsewhere.  //
  //
  //
  //============================================================================================================================================//
  
        rows_ = vector<indexType>(rowIndex);
        cols_ = vector<indexType>(cols);
        vals_ = vector<complex<double> >(vals);

        indexType NNZ = vals.size();


	matrix_ = Eigen::Map<Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType> >(rows_.size()-1, rows_.size()-1, NNZ, rows_.data(), cols_.data(), vals_.data());

	//if(matrix_ =)
	/*
	//Buntcha testing
	std::cout<<"NNZ: "<<NNZ<<"  rows size"<<rows_.size()<<"  cols size: "<<cols_.size()<<" vals_ size"<<vals_.size()  <<std::endl;
        Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType> matrix_adj = matrix_.adjoint();
	Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType> result = matrix_adj - matrix_;
	
	std::cout<<"matrix self adjointedness: "<<result.norm()  <<std::endl;
	std::cout<<"rows: "<<compute_norm(rows_)  <<std::endl;
	std::cout<<"cols: "<<compute_norm(cols_)  <<std::endl;
	std::cout<<"vals: "<<compute_norm(vals_)  <<std::endl;




	std::cout<<matrix_.block(1999989,1999989,10,10)<<std::endl;
	
        const indexType* outer_indices = matrix_.outerIndexPtr(); // Column pointers
        const indexType* inner_indices = matrix_.innerIndexPtr(); // Row indices
        const complex<double>* values = matrix_.valuePtr();       // Non-zero values


	vector<indexType> rows_2 = vector<indexType>(outer_indices, outer_indices + rows_.size()-1);
        vector<indexType> cols_2 = vector<indexType>(inner_indices, inner_indices + NNZ );
        vector<complex<double> > vals_2 = vector<complex<double> >(values, values + NNZ);


	std::cout<<"NNZ: "<<NNZ<<"  rows size:"<<rows_2.size()<<"  cols size: "<<cols_2.size()<<" vals_ size"<<vals_2.size()  <<std::endl;	
	std::cout<<"Second: "<<std::endl;
	std::cout<<"rows_2: "<<compute_norm(rows_2)-compute_norm(rows_)   <<std::endl;
	std::cout<<"cols_2: "<<compute_norm(cols_2)-compute_norm(cols_)   <<std::endl;
	std::cout<<"vals_2: "<<compute_norm(vals_2)-compute_norm(vals_)  <<std::endl;
	*/	




	
	
	
	//Eigen::SelfAdjointView<Eigen::SparseMatrix<double>, Eigen::Upper> symmetric_view(sparse_matrix);// ->Maybe saves half the memory? Couldnt see the documentation

	// --- Hermiticity guard: validate the RAW input before folding ---
	// Computed on matrix_ as read, before the H = 0.5*(M + M^dagger) fold below.
	const double num = (matrix_ - Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType>(matrix_.adjoint())).norm();
	const double den = std::max(matrix_.norm(), 1e-300);
	const double herm_res = num / den;
	const double herm_tol = 1e-8;
	if (herm_res > herm_tol) {
		std::cerr << "[CSR] non-Hermitian input: ||M-M^H||/||M|| = " << herm_res
		          << " (tol " << herm_tol << "). Likely an upper-triangle file read "
		             "under the full-matrix convention; off-diagonals would be halved."
		          << std::endl;
		throw std::runtime_error("ConvertFromCSR: non-Hermitian CSR input under full-matrix convention");
	}

	// Input CSR is the FULL Hermitian operator. Fold to enforce/clean Hermiticity:
	//   H = 0.5*(M + M^dagger)   (identity for an exactly-Hermitian M).
	// A large self-adjointness residual ||M - M^dagger|| / ||M|| means the producer
	// wrote a half/triangle matrix under the full convention (off-diagonals would be
	// silently halved) -> fail rather than corrupt; see the assertion above.
	// Do NOT reintroduce a plain "U + U^dagger" path: on a diagonal-bearing triangle it
	// doubles the on-site terms. A real triangle mode must use U + U^dagger - diag(U)
	// and be selected by an explicit storage tag, never assumed.
	Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType> transp = Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType>( matrix_.adjoint() );
	matrix_ =  0.5 * (transp + matrix_);

	

	/*
	std::cout<<Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>(matrix_).block(0,0,10,10)<<std::endl;

	std::cout<<"HEERE 3: "<<matrix_.innerIndexPtr()[0]<<"  "<<matrix_.innerIndexPtr()[1]<<"  "<<matrix_.innerIndexPtr()[2]<<std::endl;
	std::cout<<"HEERE 3: "<<matrix_.outerIndexPtr()[0]<<"  "<<matrix_.outerIndexPtr()[1]<<"  "<<matrix_.outerIndexPtr()[2]<<std::endl;
	std::cout<<"HEERE 3: "<<matrix_.valuePtr()[0]<<"  "<<matrix_.valuePtr()[1]<<"  "<<matrix_.valuePtr()[2]<<std::endl;

	std::complex<double> errVal=0.0;
	double errCol=0.0, errRow=0.0;
	
	for(indexType i=0;i<(long int)NNZ;i++){
	  errVal += matrix_.valuePtr()[i] - vals_[i];
	  errRow += cols_[i] - matrix_.innerIndexPtr()[i];  
	}

	
	for(indexType i = 0; i < rows_.size(); i++ ){
	  errCol += rows_[i] - matrix_.outerIndexPtr()[i];  
	}

	
	std::cout<<"Errors: val-"<<errVal<<" col-"<<errCol<<" row-"<< errRow<<std::endl;
	*/

	
	setDimensions(rows_.size()-1, rows_.size()-1);

}


// SparseMatrixType::Multiply -- the leaf SpMV  y = a*M*x + b*y  (raw-pointer form).
//
// PROBLEM. This is the single O(nnz) hot loop of the whole KPM machinery: the
// Chebyshev moment recursion reaches it once per moment through op().multiply()
// (HamiltonianOp -> here), so on large systems essentially all wall time is spent
// here. The stock body was `eig_y = a*matrix_*eig_x + b*eig_y`, an Eigen sparse x
// dense product that runs on ONE thread (Eigen does not parallelise sparse*vector),
// leaving 255 of 256 EPYC cores idle on the recursion. Measured: the production
// recursion was flat past ~4 threads; a row-parallel kernel reaches ~92% of the
// per-socket read-bandwidth roofline and ~30x the production driver (see
// docs/perf/PERF_ANALYSIS.md).
//
// APPROACH. Apply the matrix one CSR row at a time with an OpenMP parallel-for over
// rows. For RowMajor compressed storage, row i owns the contiguous nonzero range
// [outerIndexPtr[i], outerIndexPtr[i+1]); each y[i] is the sum over that range,
// accumulated by a single thread in stored column order, then scaled: a*acc + b*y[i].
//
// WHY THIS IS BIT-IDENTICAL TO THE STOCK EIGEN PRODUCT (golden safety). Each output
// y[i] is an independent reduction over its own row, summed in the SAME stored order
// Eigen uses, with the SAME scalar complex<double> arithmetic. Distinct rows are
// independent writes, so partitioning rows across threads cannot change any single
// sum -> the result is independent of the thread count and equals the old expression
// term-for-term. The exact-match moment goldens stay byte-identical; the gate
// `spmv_threaded_bitexact` asserts this at 1/2/4/8 threads.
//
// COMPLEXITY. O(nnz) work, O(rows) parallel tasks; memory-bandwidth bound (the
// matrix stream dominates; x is cache-resident for banded/local operators).
//
// Pitfall: this is bit-exact ONLY because each row is reduced by a single thread.
//   Do NOT parallelise a cross-vector reduction (the moment dot-product, vdot) on
//   the golden path -- that changes summation order and moves the byte-exact goldens.
// Warning: assumes COMPRESSED RowMajor storage (outer/inner/value pointers). The
//   isCompressed() guard below makes it robust if an upstream op left gaps; it is a
//   no-op on the already-compressed matrices the producers hand us.
void SparseMatrixType::Multiply(const complex<double> a, const complex<double> *x, const complex<double> b, complex<double> *y)
{
	LSQ_SCOPED("spmv");

	if (!matrix_.isCompressed())
		matrix_.makeCompressed();

	const indexType        n   = static_cast<indexType>( numRows() );
	const indexType* const rp  = matrix_.outerIndexPtr();   // row start offsets [n+1]
	const indexType* const cl  = matrix_.innerIndexPtr();   // column indices    [nnz]
	const complex<double>* const vl = matrix_.valuePtr();   // nonzero values    [nnz]

	// if(): serial on small matrices (no fork/join), parallel once the work is worth it.
	#pragma omp parallel for schedule(static) if(rp[n] > LSQ_SPMV_MIN_PARALLEL_NNZ)
	for (indexType i = 0; i < n; ++i)
	{
		complex<double> acc(0.0, 0.0);
		for (indexType p = rp[i]; p < rp[i + 1]; ++p)
			acc += vl[p] * x[cl[p]];          // stored-order sum == Eigen's order
		y[i] = a * acc + b * y[i];                // reads old y[i]; rows independent
	}

	return ;
};

void SparseMatrixType::Multiply(const complex<double> a, const vector< complex<double> >& x, const complex<double> b, vector< complex<double> >& y)
{
        Multiply(a, x.data(), b, y.data() );
        return ;
};



void SparseMatrixType::Rescale(const complex<double> a, const complex<double> b)
{
	//Create Identity Matrix
        Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType> bID(numRows(), numCols()); 


	bID.setIdentity();
        bID *= b;
	
	matrix_ = a * matrix_ + bID;


	return ;
}




// SparseMatrixType::BatchMultiply -- block SpMM  Y = a*M*X + b*Y  for R = batchSize
// right-hand sides at once. Previously a std::terminate() stub (the Eigen backend never
// had an SpMM); implemented here so multi-vector callers (block/stochastic-trace KPM)
// can amortise the matrix read across R vectors.
//
// PROBLEM. Running R independent SpMVs re-streams the whole matrix R times. When the
// recursion is bandwidth-bound (it is -- see Phase A), reading M once and applying it
// to all R vectors trades that redundant matrix traffic for arithmetic, lifting the
// arithmetic intensity (measured ~1.6x moment throughput, saturating by R~8).
//
// LAYOUT. Column-major dense blocks, leading dimension numCols() (the MKL convention the
// old stub referenced): right-hand-side r is the contiguous N-vector x + r*N, and output
// r is y + r*N. So Y[:,r] = a*M*X[:,r] + b*Y[:,r] for r in [0,R).
//
// APPROACH. One OpenMP parallel-for over rows; each thread owns whole rows. For row i it
// reads the row's nonzeros ONCE and accumulates into R running sums (one per RHS), then
// writes a*acc[r] + b*y[r*N+i]. Bit-identical to R separate Multiply() calls: the per
// (row,RHS) arithmetic and stored summation order are exactly those of the single-vector
// kernel, so column r of BatchMultiply == Multiply on column r, term-for-term. The gate
// `spmm_batchmultiply_equiv` asserts this for R in {1,8}.
//
// Pitfall: this is NOT used on any byte-exact golden path today (no caller); R=1 MUST stay
//   bit-identical to Multiply so that if a future caller adopts it the goldens cannot move.
// Complexity: O(R*nnz) work, O(rows) parallel tasks; same memory-bound character as the
//   SpMV but with R-fold reuse of each streamed matrix element.
void SparseMatrixType::BatchMultiply(const int batchSize, const complex<double> a, const complex<double> *x, const complex<double> b, complex<double> *y)
{
	LSQ_SCOPED("spmm");

	if (batchSize <= 0) return;
	if (!matrix_.isCompressed())
		matrix_.makeCompressed();

	const indexType        n   = static_cast<indexType>( numRows() );
	const indexType        R   = static_cast<indexType>( batchSize );
	const indexType        ld  = static_cast<indexType>( numCols() );   // column leading dim
	const indexType* const rp  = matrix_.outerIndexPtr();
	const indexType* const cl  = matrix_.innerIndexPtr();
	const complex<double>* const vl = matrix_.valuePtr();

	// if(): serial on small blocks (no fork/join); R-fold work so scale the nnz threshold by R.
	#pragma omp parallel if(rp[n] * R > LSQ_SPMV_MIN_PARALLEL_NNZ)
	{
		std::vector<complex<double> > acc(R);   // per-thread R accumulators
		#pragma omp for schedule(static)
		for (indexType i = 0; i < n; ++i)
		{
			for (indexType r = 0; r < R; ++r) acc[r] = complex<double>(0.0, 0.0);
			for (indexType p = rp[i]; p < rp[i + 1]; ++p)
			{
				const complex<double> v = vl[p];
				const indexType       j = cl[p];
				for (indexType r = 0; r < R; ++r)
					acc[r] += v * x[r * ld + j];     // same stored-order sum as Multiply, per RHS
			}
			for (indexType r = 0; r < R; ++r)
				y[r * ld + i] = a * acc[r] + b * y[r * ld + i];
		}
	}
}

