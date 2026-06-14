/**
 * @file quantum_states.hpp
 * @brief Source states for the stochastic-trace KPM evaluation.
 *
 * Defines `qstates::generator`, which yields the vectors whose expectation values
 * approximate operator traces: random-phase states (the default, for the stochastic
 * trace), local states, or user-supplied states. Random phases are drawn from a
 * per-instance `std::mt19937` seeded deterministically from `KPM_SEED` (default 12345),
 * so moment outputs are bitwise reproducible.
 */
#ifndef QUANTUM_STATES
#define QUANTUM_STATES

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <fstream>
#include <vector>
#include <complex>
#include <numeric>
#include <iostream>
#include <random>
using namespace std; // for std::complex<double> , and std::vector





enum StateType
{
	LOCAL_STATE,
	USER_STATE,
	RANDOM_STATE
};

typedef std::complex<double>  Complex;
typedef std::vector<Complex>  Vector;

namespace qstates
{
	typedef std::complex<double> scalar;

	/**
	 * @brief Iterator-style source-state generator for the stochastic trace.
	 *
	 * Call `getQuantumState()` repeatedly to fill the internal state vector (accessed via
	 * `State()`) with successive source vectors until `NumberOfStates()` is exhausted.
	 * `kind` selects random-phase, local, or user states; the random-phase engine is
	 * seeded from `KPM_SEED` (else 12345) for reproducibility.
	 */
	struct generator
	{
		public:

	  generator():count(0),dim(0),label("default"), num_states(1),kind(RANDOM_STATE){
		// Deterministic random-phase RNG: per-instance std::mt19937 (NOT the global
		// C rand()), so multiple generators never clobber each other's stream.
		// Seed = KPM_SEED env when set, else the fixed compile-time default 12345.
		// Never time(0) -> goldens are bitwise reproducible from a single seed.
		unsigned long kpm_seed = 12345UL;
		if( getenv("KPM_SEED") )
			kpm_seed = std::stoul(string(getenv("KPM_SEED")));
		rng_.seed( static_cast<std::mt19937::result_type>(kpm_seed) );
	  }

		int SystemSize() const
		{
			return dim;
		}

		int NumberOfStates() const
		{
			return num_states;
		}


		void NumberOfStates(const int n)
		{
			num_states = n;
			spos.reserve(n);
		}

		void SystemSize( const int n)
		{
			if( SystemSize() == 0 )
			{
				dim = n;
				out = std::vector< scalar >(dim);
				return ;
			}
			assert( SystemSize() == n);
			return ;
			
		}

		bool getQuantumState()
		{
			if( count < this->NumberOfStates() )
			{
				switch (kind)
				{
					case LOCAL_STATE:
						std::fill(out.begin(), out.end(), 0.0);
						out[ spos[count] ] = 1.0;
						break;
					case USER_STATE:
						std::cout<<"Using vector at"<<spos[count]<<std::endl;
						std::cout<<"With dimension: "<<count<<std::endl;
						out.assign(data.begin()+spos[count] , data.begin()+spos[count]+ out.size() ); 

						break;
						
						
					case RANDOM_STATE:
					{
						// Random phases drawn from the per-instance mt19937 engine.
						// Phase convention unchanged: x = (cos f + i sin f)/sqrt(N), f = 2*pi*u, u in [0,1).
						const double norm = sqrt(out.size());
						std::uniform_real_distribution<double> uni(0.0, 1.0);
						for( auto&x : out )
						{
							auto phi = 2.0*M_PI*uni(rng_);
							x = Complex(cos(phi)/norm,sin(phi)/norm );
						}
						break;
					}
				}
				count++;
				return true;
			}
			return false;
		}
		
		std::string StateLabel()
		{return label;
		}
	
		
		std::vector< scalar >& State()
		{
			return out;
		}
		

		std::string label;
		int count;
		std::vector< scalar > out;
		std::vector< scalar > data;
		std::vector<int> spos;
		int dim;
		int num_states;
		StateType kind;
		std::mt19937 rng_; // per-instance random-phase engine (deterministic; see ctor)

};
	
	inline 
	StateType String2State( std::string kind )
	{
		if( kind == "local" )
			return LOCAL_STATE;
		if( kind == "vector" )
			return USER_STATE;
		
		return RANDOM_STATE;
	}

	
	inline 
	void readLocalState(std::ifstream& file, generator& data)
	{
		int ibuff;
		file>>ibuff;
		data.NumberOfStates(ibuff);
		
		for ( int i = 0 ; i < data.num_states; i++ )
		{
			int ibuff;
			file>>ibuff;
			data.spos.push_back(ibuff);
		}
		assert( data.spos.size() == data.num_states  );	
		
		return ;
	};


	inline 
	void readCustomState(std::ifstream& file, generator& data)
	{
		int ibuff;
		file>>ibuff;
		std::cout<<"NumberOfStates"<<ibuff<<std::endl;
		data.NumberOfStates(ibuff);

		file>>ibuff;
		data.SystemSize(ibuff);
		
		const int size= data.SystemSize();
		const int num = data.NumberOfStates();
		data.data = std::vector< scalar >( size*num )  ;


		
		for ( int i = 0 ; i < num; i++ )
		{
			for ( int j = 0 ; j < size; j++ )
			{
				double re,im;
				file>>re>>im;
				data.data[i*size + j] = std::complex<double>(re,im);
				
			};
			data.spos.push_back(i*size);
			
		}
	
		return ;
	};
	
	
	inline 
	generator CreateStateSet(std::ifstream& file, std::string kind)
	{
		generator data;
		data.kind = String2State( kind );

		switch( data.kind)
		{
			case(LOCAL_STATE):
			readLocalState(file,data);
			break;
			case(USER_STATE):
			readCustomState(file,data);
			break;

			default:
			break;
			
		}
		

			
	return data;
	};



	inline 
	generator CreateLocalSet(std::vector<int> spos)
	{
		generator data;
		data.kind = LOCAL_STATE;

		data.NumberOfStates(spos.size());
		data.spos=spos;
		assert( data.spos.size() == data.num_states  );	
	return data;
	};


	
	inline
	generator LoadStateFile( string filename)
	{
		std::ifstream infile(filename);
    
		std::string kind("random_phase");
		if( infile.good() )
			infile>>kind;
		else
			std::cout<<"State file not found, using default"<<std::endl;

		std::cout<<"Building "<<kind<<" states"<<std::endl;
		auto states = CreateStateSet(infile,kind);
		std::cout<<"Success"<<std::endl;
		states.label=filename;
		return states;
	}

	inline
	int FillWithRandomPhase(Vector& X)
	{
		int kpm_seed = time(0); 	
		if( getenv("KPM_SEED") ) 
			kpm_seed = std::stoi(string(getenv("KPM_SEED")));
		srand( kpm_seed );
		std::cout<<"Current seed is "<<kpm_seed<<std::endl;
		const double norm = sqrt(X.size());
		int i = 0;
		for(auto& elem : X )
		{
			auto phi = 2.0*M_PI*(double)rand() / (double)RAND_MAX ;
			elem = Complex(cos(phi)/norm,sin(phi)/norm );
			if( i < 10 )
				std::cout<<elem.real()<<" "<<elem.imag()<<std::endl;
			i++;
			
		}


		return 0;	
	}
}
#endif
