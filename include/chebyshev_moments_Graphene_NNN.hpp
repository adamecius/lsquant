// With contributions made by Angel D. Prieto S.
#ifndef CHEBYSHEV_MOMENTS_GRAPHENE_NNN
#define CHEBYSHEV_MOMENTS_GRAPHENE_NNN


#include <complex>
#include <vector>
#include <string>
#include "units.hpp"
#include <array>

#include "sparse_matrix.hpp" //contain SparseMatrixType
#include <cassert>			 //needed for assert
#include <fstream>   		 //For ifstream and ofstream
#include <limits>    		 //Needed for dbl::digits10
#include "linear_algebra.hpp"
#include "vector_list.hpp"
#include "special_functions.hpp"
#include "chebyshev_coefficients.hpp"
#include "../src/Device/Device.hpp"
#include "../src/Graphene_NNN.hpp"

namespace chebyshev 
{
	const double CUTOFF_2 = 1.00;

class Moments_Graphene_NNN
{
	public:
	typedef std::complex<double>  value_t;
	typedef std::vector< value_t > vector_t;

	//default constructor
	Moments_Graphene_NNN(Graphene_NNN& device):
	  device_(device),_pNHAM(0),system_label(""),system_size(device.Natoms()),
	band_width(0),band_center(0){};

  	~Moments_Graphene_NNN(){};

	//GETTERS

        Graphene_NNN& device(){return device_;};
  
  void getMomentsParams( Moments_Graphene_NNN& mom)
	{
		this->SetHamiltonian( mom.Hamiltonian() ) ; 
		this->SystemLabel( mom.SystemLabel());
		this->BandWidth( mom.BandWidth() );
		this->BandCenter( mom.BandCenter() );
	};	
	
	inline
	size_t SystemSize() const { return system_size; };

	inline
	string SystemLabel() const { return system_label; };

	inline
	double BandWidth() const { return band_width; };

	inline
	double HalfWidth() const { return BandWidth()/2.0; };

	inline
	double BandCenter() const { return band_center; };

	inline
	double ScaleFactor() const { return chebyshev::CUTOFF_2/HalfWidth(); };

	inline
	double ShiftFactor() const { return -BandCenter()/HalfWidth()*chebyshev::CUTOFF_2; };

	inline 
	vector_t& MomentVector() { return mu ;}

	inline
	value_t& MomentVector(const size_t i){return  mu[i]; };

	inline
	Moments_Graphene_NNN::vector_t& Chebyshev0(){ return ChebV0; } 

	inline
	Moments_Graphene_NNN::vector_t& Chebyshev1(){ return ChebV1; } 

	inline
	Moments_Graphene_NNN::vector_t& Chebyshev2(){ return ChebV2; } 

  
	inline
	SparseMatrixType& Hamiltonian()
	{ 
		return *_pNHAM; 
	};



	//SETTERS
	inline
	void SetHamiltonian( SparseMatrixType& NHAM )
	{ 
		if ( this->SystemSize() == 0 ) //Use the rank of the hamiltonian as system size
			this->SystemSize( NHAM.rank() );	
		assert( NHAM.rank() == this->SystemSize()  );
		_pNHAM = &NHAM; 
	};



	//Heavy functions
	int  Rescale2ChebyshevDomain()
	{
		this->Hamiltonian().Rescale(this->ScaleFactor(),this->ShiftFactor());
		return 0;
	};


	inline
	void SetAndRescaleHamiltonian(SparseMatrixType& NHAM)
	{ 
		this->SetHamiltonian(NHAM );
		this->Rescale2ChebyshevDomain();
	};


	inline
	void SystemSize(const int dim)  { system_size = dim; };

	inline
	void SystemLabel(string label)  { system_label = label; };

	inline
	void BandWidth( const double x)  { band_width = x; };

	inline
	void BandCenter(const double x) { band_center = x; };

	inline 
	void MomentVector(const vector_t _mu ) { mu= _mu;}

        virtual void SetInitVectors( const vector_t& );

        virtual void SetInitVectors( SparseMatrixType & ,const vector_t&  );


	inline
	double Rescale2ChebyshevDomain(const double energ)
	{ 
		return (energ - this->BandCenter() )/this->HalfWidth(); 
	};

	int Iterate( );

	//light functions
    int JacksonKernelMomCutOff( const double broad );
	
	//light functions
    double JacksonKernel(const double m,  const double Mom );


	private:
        Graphene_NNN& device_;  
        SparseMatrixType* _pNHAM;

   
        Moments_Graphene_NNN::vector_t ChebV0,ChebV1, ChebV2,OPV;
	std::string system_label;
	size_t system_size;
	double band_width,band_center;
	vector_t mu;	
};


class Moments1D_Graphene_NNN: public Moments_Graphene_NNN
{
	public: 
  
        Moments1D_Graphene_NNN( Graphene_NNN& device):Moments_Graphene_NNN(device),_numMoms(0){};

        Moments1D_Graphene_NNN( Graphene_NNN& device, const size_t m0):Moments_Graphene_NNN(device),_numMoms(m0){ this->MomentVector( Moments_Graphene_NNN::vector_t(_numMoms, 0.0) ); };

        Moments1D_Graphene_NNN( Graphene_NNN& device, std::string momfilename );

	//GETTERS
	inline
	size_t MomentNumber() const { return _numMoms; };

	inline
	size_t HighestMomentNumber() const { return  _numMoms; };


	//OPERATORS
	inline
	Moments_Graphene_NNN::value_t& operator()(const size_t m0) { return this->MomentVector(m0); };

	void MomentNumber(const size_t _numMoms );

	void saveIn(std::string filename);

	//Transformation
	void ApplyJacksonKernel( const double broad );

	//Transformation
        void ApplyLorentzKernel( const double broad, const double lambda );
  
	// Input/Output
	void Print();

	private:
	size_t _numMoms;
};

  

class Moments2D_Graphene_NNN: public Moments_Graphene_NNN
{
	public:
        Moments2D_Graphene_NNN(Graphene_NNN& device):Moments_Graphene_NNN(device), numMoms({0,0}){};
        ~Moments2D_Graphene_NNN(){};


  Moments2D_Graphene_NNN(Graphene_NNN& device, const size_t m0,const size_t m1):Moments_Graphene_NNN(device), numMoms({int(m0),int(m1)}){ this->MomentVector( Moments_Graphene_NNN::vector_t(numMoms[1]*numMoms[0], 0.0) );    };

  Moments2D_Graphene_NNN( Graphene_NNN& device, std::string momfilename ):Moments_Graphene_NNN(device){};


  Moments2D_Graphene_NNN( Graphene_NNN& device, Moments2D_Graphene_NNN mom, const size_t m0,const size_t m1 ):Moments_Graphene_NNN(device)
	{
		this->getMomentsParams(mom);
		this->numMoms={m0,m1};
		this->MomentVector( Moments_Graphene_NNN::vector_t(numMoms[1]*numMoms[0], 0.0) );    
	};

	//GETTERS

	array<int, 2> MomentNumber() const { return numMoms; };

	int HighestMomentNumber(const int i) const { return numMoms[i]; };

	inline
	int HighestMomentNumber() const { return  (numMoms[1] > numMoms[0]) ? numMoms[1] : numMoms[0]; };

	//SETTERS
	void MomentNumber(const int mom0, const int mom1 );

	//OPERATORS
	inline
	Moments_Graphene_NNN::value_t& operator()(const int m0,const int m1) { return this->MomentVector(m0*numMoms[1] + m1 ); };

	//Transformation
	void ApplyJacksonKernel( const double b0, const double b1 );

	// Input/Output 
	//COSTFUL FUNCTIONS
	void saveIn(std::string filename);

	
	void AddSubMatrix( Moments2D_Graphene_NNN& sub , const int mL, const int mR)
	{
		for(int m0=0; m0<sub.HighestMomentNumber(0); m0++)
		for(int m1=0; m1<sub.HighestMomentNumber(1); m1++)
			this->operator()(mL+m0,mR+m1) += sub(m0,m1);
	} 

	void InsertSubMatrix( Moments2D_Graphene_NNN& sub , const int mL, const int mR)
	{
		for(int m0=0; m0<sub.HighestMomentNumber(0); m0++)
		for(int m1=0; m1<sub.HighestMomentNumber(1); m1++)
			this->operator()(mL+m0,mR+m1) = sub(m0,m1);
	} 
	
	

	void Print();

	protected:
	array<int, 2> numMoms;
};

  class MomentsTD_Graphene_NNN : public Moments_Graphene_NNN
  {
	  private:
		SparseMatrixType* _pNCON;
	  public:
	 


	  MomentsTD_Graphene_NNN(Graphene_NNN& device):
	  Moments_Graphene_NNN(device),_numMoms(1), _maxTimeStep(1),
	  _timeStep(0),
	  _dt(0)
	  {};
	  ~MomentsTD_Graphene_NNN(){};
	  
          MomentsTD_Graphene_NNN( Graphene_NNN& device, const size_t m, const size_t n ): 
	  Moments_Graphene_NNN(device),_numMoms(m), _maxTimeStep(n),
	  _timeStep(0),
	  _dt(0)
	  { this->MomentVector( Moments_Graphene_NNN::vector_t(m*n, 0.0) );    };
	  
	  MomentsTD_Graphene_NNN( Graphene_NNN& , std::string momfilename );
	  
	  //GETTERS
	  inline
	  size_t MomentNumber() const { return _numMoms;};
	  
	  inline
	  size_t HighestMomentNumber() const { return _numMoms;};
	  
	  inline
	  size_t CurrentTimeStep() const { return _timeStep;};

	  inline 
	  size_t MaxTimeStep() const { return _maxTimeStep; };
	  
	  inline
	  double TimeDiff() const   { return _dt; };
	  	  
	  inline
	  double ChebyshevFreq() const   { return HalfWidth()/chebyshev::CUTOFF_2/HBAR; };

	  inline
	  double ChebyshevFreq_0() const   { return BandCenter()/HBAR; };
	  
	  inline
	  Moments_Graphene_NNN::vector_t& Chebyshevx0(){ return ChebVx0; } 

	  inline
	  Moments_Graphene_NNN::vector_t& Chebyshevx1(){ return ChebVx1; } 

	  inline
	  Moments_Graphene_NNN::vector_t& DeltaPhi(){ return DelPhi; } 

	  inline
	  Moments_Graphene_NNN::vector_t& MSDWF(){ return WFGaussian; } 

	  void SetConmuter( SparseMatrixType& NCON )
	  { 
		if ( this->SystemSize() == 0 ) //Use the rank of the hamiltonian as system size
			this->SystemSize( NCON.rank() );	
		assert( NCON.rank() == this->SystemSize()  );
		_pNCON = &NCON; 
	  };

	  int  RescaleCon2ChebyshevDomain()
	  {
		const auto I = Moments_Graphene_NNN::value_t(0, 1);
		this->Conmuter().Rescale(this->ScaleFactor(),this->ShiftFactor());
		std::cout<<"Warning only works for symmetric boundaries in energy"<<std::endl;
		return 0;
          };

	  int ResetGaussian()
	  {
		const auto dim = this->SystemSize();
		if( this->MSDWF().size()!= dim )
			this->MSDWF() = Moments_Graphene_NNN::vector_t(dim,Moments_Graphene_NNN::value_t(0)); 
		linalg::scal(0.0,this->MSDWF());
		return 0;
	  };
 
	  void SetAndRescaleConmuter(SparseMatrixType& NCON)
	  { 
		this->SetConmuter(NCON );
		this->RescaleCon2ChebyshevDomain();
	  };



	  SparseMatrixType& Conmuter()
	  { 
		return *_pNCON; 
  	  };



	  //SETTERS
	
	  void resetCurrentTm(){
		  const auto dim = this->SystemSize();
		  if( this->DeltaPhi().size()!= dim )
			  this->DeltaPhi() = Moments_Graphene_NNN::vector_t(dim,Moments_Graphene_NNN::value_t(0)); 
		  currentTm=0;
	  }

          inline
	  void reset_MomentNumber(int new_numMoms)  {  _numMoms = new_numMoms; };

	  void MomentNumber(const size_t mom);
	  
	  void MaxTimeStep(const  size_t maxTimeStep )  {  _maxTimeStep = maxTimeStep; };

	  inline
	  void IncreaseTimeStep(){ _timeStep++; };

	  inline
	  void ResetTime(){ _timeStep=0; };
	  
	  inline
	  void TimeDiff(const double dt ) { _dt = dt; };
	  
	  
	  int Evolve(  vector_t& Phi);

	  int Conv_Evolve(  vector_t& Phi);

	  int Deltaiter(double E);

	  int MSD_Evolve(  vector_t& Phi);

	  //OPERATORS
	  inline
	  Moments_Graphene_NNN::value_t& operator()(const size_t m, const size_t n)
	  {
		  return this->MomentVector( m*MaxTimeStep() + n );
	  };
	  
	  //Transformation
	  void ApplyJacksonKernel( const double broad );
	  
	  //COSTFUL FUNCTIONS
	  void saveIn(std::string filename);

          int dT(){return _dt;};

    
	  void Print();

  private:
    int currentTm;
    Moments_Graphene_NNN::vector_t ChebVx0,ChebVx1,WFGaussian,DelPhi;
    size_t _numMoms, _maxTimeStep, _timeStep;
    double _dt;
  };




class Vectors_sliced_Graphene_NNN : public Moments_Graphene_NNN
{
	public: 
	typedef VectorList< Moments_Graphene_NNN::value_t > vectorList_t;



        ~Vectors_sliced_Graphene_NNN( ){};

  //virtual void SetInitVectors( const vector_t& );

  //virtual void SetInitVectors( SparseMatrixType & ,const vector_t&  );
  virtual void SetInitVectors_2( const vector_t& );

  	virtual void SetInitVectors_2( SparseMatrixType & ,const vector_t&  );

	  
	int NumberOfVectors() const
	{ return numVecs_;}

	int NumberOfSections() const
	{ return num_sections_;}
  
  
        Moments_Graphene_NNN::vector_t& OPV()
        {return OPV_;};
  
  	int SectionSize() const
	{ return section_size_;}

        int LastSectionSize() const
	{ return last_section_size_;}

        vectorList_t& Chebmu(){return Chebmu_;};
  
  	void SetNumberOfVectors( const int x)
	{
	  numVecs_ = x;
	}
  
	int SetNumSections( const int x)
	{

	    num_sections_ = x;
	    section_size_ = SystemSize()/num_sections_; // All of this assuming */* is larger than *%*. Otherwise, i'd change num_sections until it is.
	    last_section_size_ = SystemSize()/num_sections_ + SystemSize() % num_sections_;
	  
	    return 0;
	}
  
        void SetChebmu(const int num_vecs, const size_t last_section_size)
        {
	     Chebmu_(num_vecs, last_section_size);
        }



  Vectors_sliced_Graphene_NNN(Graphene_NNN& device):Moments_Graphene_NNN(device)
        {
                SetChebmu(0,0);
                SetNumSections(1);
	};
  

  Vectors_sliced_Graphene_NNN(Graphene_NNN& device, chebyshev::Vectors_sliced_Graphene_NNN& vec):Moments_Graphene_NNN(device){ *this = vec; };
  
  Vectors_sliced_Graphene_NNN( Graphene_NNN& device,const int numVecs, const size_t num_sections): Moments_Graphene_NNN(device), numVecs_(numVecs), num_sections_(num_sections), Chebmu_(0,0) {};
  
  Vectors_sliced_Graphene_NNN(Graphene_NNN& device, Moments_Graphene_NNN& mom, const int numVecs, const size_t num_sections): Moments_Graphene_NNN(device), numVecs_(numVecs), num_sections_(num_sections), Chebmu_(numVecs, mom.SystemSize()/num_sections + mom.SystemSize() % num_sections ) { 

	  section_size_ = mom.SystemSize()/num_sections; // All of this assuming */* is larger than *%*. Otherwise, i'd change num_sections until it is.
	  last_section_size_ = mom.SystemSize()/num_sections + mom.SystemSize() % num_sections;
	  
	  this->getMomentsParams(mom);	  
        };

  Vectors_sliced_Graphene_NNN(Graphene_NNN& device, Moments2D_Graphene_NNN& mom, const int numVecs, const size_t num_sections): Moments_Graphene_NNN(device), numVecs_(numVecs), num_sections_(num_sections), Chebmu_(numVecs, mom.SystemSize()/num_sections + mom.SystemSize() % num_sections ) { 

	  section_size_ = mom.SystemSize()/num_sections; // All of this assuming */* is larger than *%*. Otherwise, i'd change num_sections until it is.
	  last_section_size_ = mom.SystemSize()/num_sections + mom.SystemSize() % num_sections;
	  
	  this->getMomentsParams(mom);	  
        };

  
  
        int CreateVectorSet()
        {
		try
		{

		        this->SystemSize( this->Hamiltonian().rank() );
			section_size_ = SystemSize()/num_sections_; // All of this assuming */* is larger than *%*. Otherwise, i'd change num_sections until it is.
	                last_section_size_ = SystemSize()/num_sections_ + SystemSize() % num_sections_;

		        const int vec_size  = last_section_size_;
			const int list_size = this->NumberOfVectors();
			Chebmu_.Resize(list_size, vec_size );
		}
		catch (...)
		{ std::cerr<<"Failed to initilize the vector list."<<std::endl;}

		
		return 0;
	}
    
    
	void getMomentsParams( Moments_Graphene_NNN& mom)
	{
		this->SetHamiltonian( mom.Hamiltonian() ) ; 
		this->SystemLabel( mom.SystemLabel());
		this->BandWidth( mom.BandWidth() );
		this->BandCenter( mom.BandCenter() );
		if ( mom.Chebyshev0().size() == this->SystemSize() )
			this->SetInitVectors( mom.Chebyshev0() );
	}

	inline
	size_t Size() const
	{
		return  (long unsigned int)section_size_*
				(long unsigned int)this->NumberOfVectors();
	}

	inline 
	double SizeInGB() const
	{
		const double list_size = this->NumberOfVectors();
		return  sizeof(value_t)*section_size_*list_size*pow(2.0,-30.0);
	}

	inline
	size_t HighestMomentNumber() const { return  this->Chebmu_.ListSize(); };

	inline
	vectorList_t& List() { return this->Chebmu_; };

	inline
	Moments_Graphene_NNN::vector_t& Vector(const size_t m0) { return this->Chebmu_.ListElem(m0); };

	inline
	Moments_Graphene_NNN::value_t& operator()(const size_t m0) { return this->Chebmu_(m0,0); };

        void operator=( Vectors_sliced_Graphene_NNN& vec)
	{
		this->SetHamiltonian( vec.Hamiltonian() ) ; 
		this->SystemLabel( vec.SystemLabel());
		this->BandWidth( vec.BandWidth() );
		this->BandCenter( vec.BandCenter() );
		if ( vec.Chebyshev0().size() == this->SystemSize() )
			this->SetInitVectors( vec.Chebyshev0() );

		num_sections_ = vec.NumberOfSections();
		numVecs_ = vec.NumberOfVectors();
		section_size_ = vec.SectionSize();
	        last_section_size_ = vec.LastSectionSize();

	}

  
	virtual int IterateAllSliced( int );

        virtual int MultiplySliced( SparseMatrixType & , int );

	double MemoryConsumptionInGB();


  

	private:
	Moments_Graphene_NNN::vector_t OPV_;
        int numVecs_, num_sections_;
        size_t section_size_, last_section_size_;
        vectorList_t Chebmu_;	
  };




class Vectors_Graphene_NNN : public Moments_Graphene_NNN
{
	public: 
	typedef VectorList< Moments_Graphene_NNN::value_t > vectorList_t;


	int NumberOfVectors() const
	{ return numVecs;}


	int SetNumberOfVectors( const int x)
	{ numVecs = x; return 0;}


  Vectors_Graphene_NNN(Graphene_NNN& device):Moments_Graphene_NNN(device),Chebmu(0,0){};

	~Vectors_Graphene_NNN(){};

  
	Vectors_Graphene_NNN(Graphene_NNN& device,const size_t nMoms,const size_t dim ):Moments_Graphene_NNN(device),Chebmu(nMoms,dim) {  };

	Vectors_Graphene_NNN( Graphene_NNN& device, Moments1D_Graphene_NNN& mom ): Moments_Graphene_NNN(device),Chebmu(mom.HighestMomentNumber(), mom.SystemSize() )
	{ 
		this->getMomentsParams(mom);
	};
	
	Vectors_Graphene_NNN( Graphene_NNN& device, Moments2D_Graphene_NNN& mom ): Moments_Graphene_NNN(device), Chebmu(mom.HighestMomentNumber(), mom.SystemSize() )
	{ 
		this->getMomentsParams(mom);
	};

	
	Vectors_Graphene_NNN( Graphene_NNN& device, Moments2D_Graphene_NNN& mom, size_t i  ): Moments_Graphene_NNN(device), Chebmu(mom.HighestMomentNumber(), mom.SystemSize() )
	{ 
		this->getMomentsParams(mom);
	};
	
	
   
    int CreateVectorSet()
    {
		try
		{
			const int vec_size  = this->SystemSize();
			const int list_size = this->NumberOfVectors();
			Chebmu(list_size, vec_size );
		}
		catch (...)
		{ std::cerr<<"Failed to initilize the vector list."<<std::endl;}

		
		return 0;
	}
    
    
    
  Vectors_Graphene_NNN( Graphene_NNN& device, MomentsTD_Graphene_NNN& mom ):Moments_Graphene_NNN(device), Chebmu(mom.HighestMomentNumber(), mom.SystemSize() )
	{ 
		this->getMomentsParams(mom);
	};

	void getMomentsParams( Moments_Graphene_NNN& mom)
	{
		this->SetHamiltonian( mom.Hamiltonian() ) ; 
		this->SystemLabel( mom.SystemLabel());
		this->BandWidth( mom.BandWidth() );
		this->BandCenter( mom.BandCenter() );
		if ( mom.Chebyshev0().size() == this->SystemSize() )
			this->SetInitVectors( mom.Chebyshev0() );
	}

	inline
	size_t Size() const
	{
		return  (long unsigned int)this->SystemSize()*
				(long unsigned int)this->NumberOfVectors();
	}

	inline 
	double SizeInGB() const
	{
		const double vec_size  = this->SystemSize();
		const double list_size = this->NumberOfVectors();
		return  sizeof(value_t)*vec_size*list_size*pow(2.0,-30.0);
	}

	inline
	size_t HighestMomentNumber() const { return  this->Chebmu.ListSize(); };


	inline
	vectorList_t& List() { return this->Chebmu; };

	inline
	Moments_Graphene_NNN::vector_t& Vector(const size_t m0) { return this->Chebmu.ListElem(m0); };


	inline
	Moments_Graphene_NNN::value_t& operator()(const size_t m0) { return this->Chebmu(m0,0); };

	int IterateAll( );

	int EvolveAll( const double DeltaT, const double Omega0);

	int Multiply( SparseMatrixType &OP );


	double MemoryConsumptionInGB();


	private:
	Moments_Graphene_NNN::vector_t OPV;
	vectorList_t Chebmu;	
	int numVecs;
};



};

#endif 


