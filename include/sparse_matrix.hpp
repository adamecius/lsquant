/**
 * @file sparse_matrix.hpp
 * @brief Sparse operator type used throughout LinQT.
 *
 * Defines @ref SparseMatrixType, the matrix-free complex sparse operator (CSR, Eigen
 * RowMajor backend) on which all KPM recursions act, and @ref SparseMatrixBuilder, which
 * loads an operator from a CSR text file. Operators are the Hamiltonian and the
 * observables (velocity, spin, ...) read from the input tool's `.CSR` files.
 */
#ifndef SPARSE_MATRIX
#define SPARSE_MATRIX

#include <assert.h> /* assert */
#include <iostream>
#include <vector>
#include <complex>
#include <fstream>

#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Sparse>

using std::complex;
using std::vector;
using std::string;

typedef long int indexType;

namespace Sparse
{
bool OPERATOR_FromCSRFile(const std::string input, int &dim, vector<indexType> &columns, vector<indexType> &rowIndex, vector<complex<double> > &values);
};

class SparseMatrixType_BASE
{
	
	
public:
	typedef complex<double> value_t;
	typedef vector< value_t > vector_t;
	
	  int numRows() { return numRows_; };
	  int numCols() { return numCols_; };
	  int rank() { return ((this->numRows() > this->numCols()) ? this->numCols() : this->numRows()); };
	  void setDimensions(const int numRows, const int numCols)
	  {
		numRows_ = numRows;
		numCols_ = numCols;
	  };
	  void SetID(string id) { id_ = id; }
	  string ID() const { return id_; }
		
	  bool isIdentity(){ return (bool)( ID()=="1"); };

private:
  int numRows_, numCols_;
  string id_;
};


/**
 * @brief Complex sparse operator (Eigen RowMajor / CSR) with the matrix-free product
 *        used by the KPM recursion.
 *
 * Holds an operator (Hamiltonian or observable) and provides the sparse matrix-vector
 * product `y = a*M*x + b*y` that drives every Chebyshev iteration. Built from a CSR text
 * file via @ref SparseMatrixBuilder; `Rescale` maps the spectrum into `[-1,1]` for the
 * Chebyshev expansion.
 */
class SparseMatrixType  : public SparseMatrixType_BASE
{
public:
  SparseMatrixType()
  {}
  ~SparseMatrixType() {}


  string matrixType() const { return "CSR Matrix from Eigen."; };
  /// @brief Sparse matrix-vector product `y = a*M*x + b*y` (raw pointers).
  void Multiply(const value_t a, const value_t *x, const value_t b, value_t *y);
  /// @brief Sparse matrix-vector product `y = a*M*x + b*y` (std::vector).
  void Multiply(const value_t a, const vector_t& x, const value_t b, vector_t& y);
  /// @brief Affine rescale `M <- a*M + b*I` (used to map the spectrum into [-1,1]).
  void Rescale(const value_t a,const value_t b);
  inline 
  void Multiply(const value_t *x, value_t *y){ Multiply(value_t(1.0,0),x,value_t(0,0),y);};
  inline 
  void Multiply(const vector_t& x, vector_t& y){ Multiply(value_t(1.0,0),x,value_t(0,0),y);};


  void BatchMultiply(const int batchSize, const value_t a, const value_t *x, const value_t b, value_t *y);

  /// @brief Assemble the operator from coordinate (COO) triplets.
  void ConvertFromCOO(vector<indexType> &rows, vector<indexType> &cols, vector<complex<double> > &vals);
  /// @brief Assemble the operator from CSR arrays; the FULL Hermitian matrix is expected
  ///        and folded as H = 0.5*(M + M^dagger) (a non-Hermitian input trips a guard).
  void ConvertFromCSR(vector<indexType> &rowIndex, vector<indexType> &cols, vector<complex<double> > &vals);

  vector<indexType>* rows() {return &rows_;};
  vector<indexType>* cols() {return &cols_;};
  vector<complex<double> >* vals(){return &vals_;};
  Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType>& eigen_matrix(){return matrix_;};

  void set_eigen_matrix( Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType>& new_matrix ){
    matrix_ = new_matrix;
    setDimensions(new_matrix.rows(),new_matrix.cols());
  };
  
private:
  Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType> matrix_;
  vector<indexType> rows_;
  vector<indexType> cols_;
  vector<complex<double> > vals_;
};



class SparseMatrixBuilder
{
public:
  void setSparseMatrix(SparseMatrixType *b)
  {
    _matrix_type = b;
  };

public:
  void BuildOPFromCSRFile(const std::string input)
  {
    vector<indexType> columns, rowIndex;
    vector<complex<double> > values;
    int dim;

    if(!Sparse::OPERATOR_FromCSRFile(input, dim, columns, rowIndex, values))
      return;
    
    _matrix_type->setDimensions(dim, dim);
    _matrix_type->ConvertFromCSR( columns,rowIndex, values);
    std::cout << "OPERATOR SUCCESSFULLY BUILT" << std::endl;
  }
  SparseMatrixType *_matrix_type;
};

#endif
