#include "sparse_matrix.hpp"


const int MKL_SPMVMUL = 1000;

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


void SparseMatrixType::Multiply(const complex<double> a, const complex<double> *x, const complex<double> b, complex<double> *y)
{
        Eigen::Map<const Eigen::Vector<std::complex<double>, -1>>
	  eig_x(x, numRows());

	Eigen::Map<Eigen::Vector<std::complex<double>, -1>>
	  eig_y(y, numRows());




	eig_y = a * matrix_ * eig_x + b * eig_y; 


	  
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




void SparseMatrixType::BatchMultiply(const int batchSize, const complex<double> a, const complex<double> *x, const complex<double> b, complex<double> *y)
{
  	std::cout<<"ERROR: SparseMatrixType:BatchMultiply function not implemented in eigen yet."<<std::endl;
	std::terminate();
	//auto status = mkl_sparse_z_mm(SPARSE_OPERATION_NON_TRANSPOSE,a,Matrix, descr,SPARSE_LAYOUT_COLUMN_MAJOR,x, batchSize, numCols(), b, y, numCols() );
  //assert( status == SPARSE_STATUS_SUCCESS);
}

