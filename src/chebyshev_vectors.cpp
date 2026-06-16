#include "chebyshev_moments.hpp"



void chebyshev::Vectors_sliced_nonOrthogonal::SetInitVectors_2( const  Moments::vector_t& T0 )
{
  


	assert( T0.size() == this->SystemSize() );
	const auto dim = this->SystemSize();

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()

	
        vector_t tmp_(dim,0.0);

	
	
	linalg::copy ( T0, this->Chebyshev0() );	
	this->op().multiply( 1.0, this->Chebyshev0(), 0.0, tmp_ );
        linalg::orthogonalize(*S_,  tmp_, this->Chebyshev1());



	/*
	double b =  this->ShiftFactor();
	
 #pragma omp parallel for
	for(std::size_t i=0; i< dim;i++)	  
	  Chebyshev1()[i] += b * Chebyshev0()[i];
	*/
};

void chebyshev::Vectors_sliced_nonOrthogonal::SetInitVectors_2( SparseMatrixType &OP , const Moments::vector_t& T0 )
{
	const auto dim = this->SystemSize();
	assert( OP.rank() == this->SystemSize() && T0.size() == this->SystemSize() );

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()


	assert( OP.rank() == this->SystemSize() );
	if( OPV().size()!= OP.rank() )
	       OPV() = Moments::vector_t ( OP.rank() );

	
	vector_t tmp_(dim,0.0);
        Moments::vector_t tmp2(this->SystemSize());


	double b =  this->ShiftFactor(),
	  a =  this->ScaleFactor();

	linalg::copy ( T0, this->Chebyshev1() );
        

	  dHdK_->Multiply(1.0,this->Chebyshev1(),0.0, tmp_); //Multiply(  OPV(), tmp2 );
	  linalg::orthogonalize(*S_, tmp_, this->Chebyshev0());

	
	  
	  this->op().multiply(1.0,this->Chebyshev1(),0.0, tmp_);
	  
#pragma omp parallel for
	  for(std::size_t i=0; i< dim;i++){
	    tmp_[i] /= a;
	    //tmp_[i] += b * Chebyshev0()[i];
	  }
	  
	  linalg::orthogonalize(*S_, tmp_, tmp2);
	  dSdK_->Multiply(1.0,tmp2,0.0, tmp_); //Multiply(  OPV(), tmp2 );
	  linalg::orthogonalize(*S_, tmp_, tmp2);
	  
	  
	  linalg::axpy( -1.0, tmp2, this->Chebyshev0());
	  

	  
	
	this->op().multiply( 1.0, this->Chebyshev0(), 0.0, tmp_ );
        linalg::orthogonalize(*S_, tmp_ , this->Chebyshev1());


	/*
	double b =  this->ShiftFactor();
	
 #pragma omp parallel for
	for(std::size_t i=0; i< dim;i++)	  
	  Chebyshev1()[i] += b * Chebyshev0()[i];
	*/




	return ;
};



int chebyshev::Vectors_sliced_nonOrthogonal::MultiplySliced( SparseMatrixType &OP , int s )
{
  size_t segment_size = ( s == this->NumberOfSections()-1 ? this->LastSectionSize() : this->SectionSize() ),
    segment_start = s * this->SectionSize(),
    DIM = this->SystemSize();

  Moments::vector_t tmp2(this->SystemSize()),
    tmp3(this->SystemSize()),
    tmp4(this->SystemSize());

  double b =  this->ShiftFactor(),
	 a =  this->ScaleFactor();

  



	assert( OP.rank() == this->SystemSize() );
	if( OPV().size()!= OP.rank() )
	       OPV() = Moments::vector_t ( OP.rank() );

#pragma omp parallel for
	for(int i=0; i<OP.rank(); i++){//not parallelized; with omp/ eigen this is straightforward;
	  OPV()[i] = 0.0;
	  tmp2[i] = 0.0;
	  tmp3[i] = 0.0;
	  tmp4[i] = 0.0;
	}


	for(int m=0; m < this->NumberOfVectors(); m++ )
	{

	  	  
	  linalg::introduce_segment(this->Chebmu().ListElem(m), tmp4, segment_start);

	  
	  dHdK_->Multiply(1.0,tmp4,0.0, tmp2); //Multiply(  OPV(), tmp2 );
	  linalg::orthogonalize(*S_, tmp2, tmp3);

	  
	  this->op().multiply(1.0,tmp4,0.0, tmp2);


#pragma omp parallel for
	  for(std::size_t i=0; i< DIM;i++){
	    tmp2[i] /= a;
	    //Chebyshev1()[i] += b * Chebyshev0()[i];
	  }
	  	  

	  linalg::orthogonalize(*S_, tmp2, tmp4);
	  dSdK_->Multiply(1.0,tmp4,0.0, tmp2); //Multiply(  OPV(), tmp2 );	  
	  linalg::orthogonalize(*S_, tmp2, tmp4);
	  
	  
	  linalg::axpy( -1.0, tmp4,tmp3);
	  



	  linalg::extract_segment(tmp3, segment_start,  this->Chebmu().ListElem(m));
	  
	}

	return 0;

};



int chebyshev::Vectors_sliced_nonOrthogonal::IterateAllSliced( int s )
{

        std::size_t dim = Chebyshev0().size();
        vector_t tmp_( dim, 0.0 ),
	         tmp_2_( dim, 0.0);

	size_t segment_size = ( s == this->NumberOfSections()-1 ? this->LastSectionSize() : this->SectionSize() ),
	  segment_start = s * this->SectionSize(),
         DIM = this->SystemSize();

  
	//The vectorss Chebyshev0() and Chebyshev1() are assumed to have
	// been initialized

        linalg::extract_segment( Chebyshev0(),  segment_start, Vector(0));
	for(int m=1; m < this->NumberOfVectors(); m++ )
	{

	  linalg::extract_segment( Chebyshev1(),  segment_start, Vector(m));
	  this->op().multiply( this->Chebyshev1(),tmp_);
	
       	  linalg::orthogonalize(*S_, tmp_, tmp_2_);
	
	  #pragma omp parallel for
	  for(std::size_t i=0; i< dim;i++){
	    Chebyshev0()[i] = 2.0 * tmp_2_[i] - Chebyshev0()[i];
	  }
	
	  this->Chebyshev0().swap(this->Chebyshev1());


	}
	

   return 0;

};













double chebyshev::Vectors_sliced_nonOrthogonal::IterateAllSliced_test(int s,  SparseMatrixType &orth_Ham )
{

        double error=0;  
        std::size_t dim = Chebyshev0().size();
        vector_t tmp_( dim, 0.0 ),
	  tmp_2_( dim, 0.0),
	  compare( dim, 0.0 );

	double b =  this->ShiftFactor();


	size_t segment_size = ( s == this->NumberOfSections()-1 ? this->LastSectionSize() : this->SectionSize() ),
	  segment_start = s * this->SectionSize(),
         DIM = this->SystemSize();

	
	this->op().multiply( this->Chebyshev1(),tmp_);
       	linalg::orthogonalize(*S_, tmp_, tmp_2_);

	
	orth_Ham.Multiply(this->Chebyshev1(), compare);
	

        linalg::extract_segment( Chebyshev0(),  segment_start, Vector(0));	
	#pragma omp parallel for
	for(std::size_t i=0; i< dim;i++){
	  std::complex<double> tmp_2 = tmp_2_[i] + b * Chebyshev1()[i];
	  
	  Chebyshev0()[i] = 2.0 * tmp_2 - Chebyshev0()[i];

	  error += std::abs( ( compare[i] - tmp_2 ) / compare[i] );

	  linalg::extract_segment( Chebyshev1(),  segment_start, Vector(0));
	  this->op().multiply(2.0,Chebyshev1(),-1.0,Chebyshev0());
	  Chebyshev0().swap(Chebyshev1());

	}
	
	this->Chebyshev0().swap(this->Chebyshev1());

	

	return error/dim;

};
































void chebyshev::Vectors_sliced::SetInitVectors_2( SparseMatrixType &OP , const Moments::vector_t& T0 )
{

        const auto dim = this->SystemSize();
	assert( OP.rank() == this->SystemSize() && T0.size() == this->SystemSize() );

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()


	linalg::copy ( T0, this->Chebyshev1() );
	OP.Multiply( 1.0, this->Chebyshev1(), 0.0, this->Chebyshev0() );
	this->op().multiply( 1.0, this->Chebyshev0(), 0.0, this->Chebyshev1() );

	return ;
};


void chebyshev::Vectors_sliced::SetInitVectors_2( const Moments::vector_t& T0 )
{

	assert( T0.size() == this->SystemSize() );
	const auto dim = this->SystemSize();

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()

	linalg::copy ( T0, this->Chebyshev0() );
	this->op().multiply( 1.0, this->Chebyshev0(), 0.0, this->Chebyshev1() );
};




int chebyshev::Vectors::IterateAll( )
{	
	//The vectorss Chebyshev0() and Chebyshev1() are assumed to have
	// been initialized
	linalg::copy( this->Chebyshev0() ,this->Vector(0) );
	for(int m=1; m < this->NumberOfVectors(); m++ )
	{
		linalg::copy( Chebyshev1() , this->Vector(m) );
		this->op().multiply(2.0,Chebyshev1(),-1.0,Chebyshev0());
		Chebyshev0().swap(Chebyshev1());
	}
	return 0;
};


int chebyshev::Vectors::Multiply( SparseMatrixType &OP )
{
	assert( OP.rank() == this->SystemSize() );
	if( this->OPV.size()!= OP.rank() )
		this->OPV = Moments::vector_t ( OP.rank() );
	
	for(size_t m=0; m < this->NumberOfVectors(); m++ )
	{
		linalg::copy( this->Chebmu.ListElem(m), this->OPV ); 
		OP.Multiply(  this->OPV, this->Chebmu.ListElem(m) );
	}

	return 0;
};






int chebyshev::Vectors_sliced::IterateAllSliced(int s )
{

  size_t segment_size = ( s == num_sections_-1 ? last_section_size_ : section_size_ ),
         segment_start = s * section_size_,
         DIM = this->SystemSize();

  
	//The vectorss Chebyshev0() and Chebyshev1() are assumed to have
	// been initialized

        linalg::extract_segment( Chebyshev0(),  segment_start, Vector(0));
	//linalg::copy( this->Chebyshev0() ,this->Vector(0) );

	for(int m=1; m < this->NumberOfVectors(); m++ )
	{
	  
	  linalg::extract_segment( Chebyshev1(),  segment_start, Vector(m));
	  this->op().multiply(2.0,Chebyshev1(),-1.0,Chebyshev0());
	  Chebyshev0().swap(Chebyshev1());
	  	
}
	return 0;
};


int chebyshev::Vectors_sliced::MultiplySliced( SparseMatrixType &OP, int s)
{
  size_t segment_size = ( s == num_sections_-1 ? last_section_size_ : section_size_ ),
    segment_start = s * section_size_,
    DIM = this->SystemSize();

  Moments::vector_t tmp2(this->SystemSize());


	assert( OP.rank() == this->SystemSize() );
	if( OPV().size()!= OP.rank() )
	       OPV() = Moments::vector_t ( OP.rank() );

#pragma omp parallel for
	for(int i=0; i<OP.rank(); i++){//not parallelized; with omp/ eigen this is straightforward;
	  OPV()[i] = 0.0;
	  tmp2[i] = 0.0;
	}

	
	for(int m=0; m < this->NumberOfVectors(); m++ )
	{
	  linalg::introduce_segment(Chebmu_.ListElem(m), OPV(), segment_start);
	  OP.Multiply(1.0,OPV(),0.0, tmp2); //Multiply(  OPV(), tmp2 );
	  linalg::extract_segment(tmp2, segment_start,  Chebmu_.ListElem(m));	  
	}

	return 0;
};
/*
int chebyshev::Vectors_sliced::MultiplySliced( SparseMatrixType &OP, int s)
{
  size_t segment_size = ( s == num_sections_-1 ? last_section_size_ : section_size_ ),
    segment_start = s * section_size_,
    DIM = this->SystemSize();

  
  Moments::vector_t tmp2(this->SystemSize());

	assert( OP.rank() == this->SystemSize() );
	if( OPV().size()!= OP.rank() )
	       OPV() = Moments::vector_t ( OP.rank() );

#pragma omp parallel for
	for(int i=0; i<OP.rank(); i++){//not parallelized; with omp/ eigen this is straightforward;
	  OPV()[i] = 0.0;
	  tmp2[i] = 0.0;
	}

	
	for(int m=0; m < this->NumberOfVectors(); m++ )
	{
	  	  
          linalg::introduce_segment(Chebmu_.ListElem(m), OPV(), segment_start);
	  OP.Multiply(1.0,OPV(),0.0, tmp2); //Multiply(  OPV(), tmp2 );
	  linalg::extract_segment(tmp2, segment_start,  Chebmu_.ListElem(m));	  
	  	
        }

	return 0;
};
*/

int chebyshev::Vectors::EvolveAll(const double DeltaT, const double Omega0)
{
	const auto dim = this->SystemSize();
	const auto numVecs = this->NumberOfVectors();

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()

	const auto I = Moments::value_t(0, 1);
	const double x = Omega0*DeltaT;
	for(size_t m=0; m < this->NumberOfVectors(); m++ )
	{
		auto& myVec = this->Vector(m);
		
		int n = 0;
		double Jn = besselJ(n,x);
		linalg::copy(myVec , Chebyshev0());
		linalg::scal(0, myVec); //Set to zero
		linalg::axpy( Jn , Chebyshev0(), myVec);

		double Jn1 = besselJ(n+1,x);	
		this->op().multiply(Chebyshev0(), Chebyshev1());
		linalg::axpy(-value_t(2) * I * Jn1, Chebyshev1(), myVec);
		
		auto nIp =-I;
		while( 0.5*(std::abs(Jn)+std::abs(Jn1) ) > 1e-15)
		{
			nIp*=-I ;
			Jn  = Jn1;
			Jn1 = besselJ(n, x);
			this->op().multiply(2.0, Chebyshev1(), -1.0, Chebyshev0());
			linalg::axpy(value_t(2) * nIp * value_t(Jn1), Chebyshev0(), myVec);
			Chebyshev0().swap(Chebyshev1());
			n++;
		}
	}
  return 0;
};


double chebyshev::Vectors::MemoryConsumptionInGB()
{
	return SizeInGB()+2.0*( (double)this->SystemSize() )*pow(2.0,-30.0) ;
}
	
	
