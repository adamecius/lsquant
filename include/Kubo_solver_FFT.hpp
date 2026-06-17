/**
 * @file Kubo_solver_FFT.hpp
 * @brief FFT-accelerated Kubo–Bastin conductivity solver.
 *
 * Reconstructs the Kubo–Bastin conductivity from the 2D Chebyshev moments using an
 * FFT over the Chebyshev/energy grid (see `Bastin_FFTs` / `Greenwood_FFTs`), which turns
 * the double Chebyshev sum into FFTs for large moment counts.
 */
#ifndef KUBO_BASTIN_SOLVER_HPP
#define KUBO_BASTIN_SOLVER_HPP

#include<string>
#include <cstring>

#include<iostream>
#include<complex>
#include"sparse_matrix.hpp"
#include"chebyshev_moments.hpp"
#include"chebyshev_solver.hpp"


namespace chebyshev 
{
  
enum formula{
  NIL = -1,
  KUBO_GREENWOOD = 0,
  KUBO_BASTIN = 1
};

  


  
class Kubo_solver_FFT{

        public:
        typedef double  r_value_t;


        Kubo_solver_FFT(int M, int num_sections, int nump, formula sym_formula, chebyshev::Vectors_sliced& chebVecL, chebyshev::Vectors_sliced& chebVecR, std::string& outputfilename)
	  : M_(M), num_sections_(num_sections), nump_(nump), section_size_(0),
	    sym_formula_(sym_formula), chebVecL_(&chebVecL), chebVecR_(&chebVecR), outputfilename_(outputfilename) {};

        ~Kubo_solver_FFT(){};
  
        void compute(SparseMatrixType &, SparseMatrixType &,  qstates::generator& );


  
        //GETTERS
        inline
	int nump(){return nump_; };

        inline
	int num_sections(){return num_sections_; };

        inline
	size_t section_size(){return section_size_; };

        inline
	int num_vectors(){return M_; };


        inline
	chebyshev::Vectors_sliced Hamiltonian(){return *chebVecL_; };

        inline
	std::string& OutputFilename(){return outputfilename_;};

        inline
	formula SolverFormula(){return sym_formula_;};


        double relativeL2Norm(const value_t*, const value_t* );
  
        protected:  
        void allocate(Vectors_sliced& , Vectors_sliced& );

        void polynomial_recursion(const vector_t& , const vector_t& ,
				  SparseMatrixType &,
				  SparseMatrixType &,  
				  Vectors_sliced *,
				  Vectors_sliced *,
				  value_t*);


  
        void Greenwood_FFTs ( Vectors_sliced&, Vectors_sliced&,  value_t*, int);
        void Bastin_FFTs ( Vectors_sliced&, Vectors_sliced&, value_t*, int);


  
	//Initializers
        inline
        void reset_data(value_t data[]){
          for(int i = 0; i < 2 * nump_; i++)
	    data[i] = 0.0;
	};

        inline
        void update_data(value_t final_data[], const value_t new_r_data[], int  r){
          for(int i = 0; i < 2 * nump_; i++)
	    final_data[i] = ( final_data[i] * value_t( r - 1 ) + new_r_data[i] ) / value_t(r);
	};

       
        r_value_t Jackson_kernel(const int m, const int M){
          return 1/(M+1.0) *
           (
   	     (M-m+1.0) *
             cos( M_PI * m / (M + 1.0)) +
             sin( M_PI * m / (M + 1.0)) /
             tan( M_PI     / (M + 1.0))
	   );
        };


  

        private:
        int M_, num_sections_, nump_;
	size_t section_size_;
        formula sym_formula_;
        chebyshev::Moments Hamiltonian_dummyMoms_;
        chebyshev::Vectors_sliced *chebVecL_, *chebVecR_;

        std::string outputfilename_;

        friend class Kubo_solver_FFT_postProcess;
};



  
  

class Kubo_solver_FFT_postProcess{//will interpret data_set of points k=0,...,nump-1 as associated to  the energies e_k=a_*cos(MPI*(2*k+0.5))-b_; 

        public:
        typedef double  r_value_t;
        typedef std::complex< r_value_t > value_t;

        Kubo_solver_FFT_postProcess(Kubo_solver_FFT&);
        ~Kubo_solver_FFT_postProcess(){ delete [] E_points_;/*, delete [] conv_R_;*/};

  
        void operator()(const value_t*, const value_t*, int);
  
        void Greenwood_postProcess (const value_t*, const value_t*, int  );
        void Bastin_postProcess (const value_t*, const value_t*, int  );

        void integration ( const r_value_t*, const r_value_t*, r_value_t* );
        void integration_lsquant ( const r_value_t*, const r_value_t*, r_value_t* );
        void partial_integration ( const r_value_t*, const r_value_t*, r_value_t* );

        void rearrange_crescent_order( r_value_t* );
        void eta_CAP_correct(r_value_t*, r_value_t* );  
        void plot_data   ( const std::string&, const std::string& );

        private:
        Kubo_solver_FFT parent_solver_;
        r_value_t *E_points_;
        //           *conv_R_;

    
  };


};//Ends namespace chebyshev


#endif //KUBO_BASTIN_SOLVER_HPP
