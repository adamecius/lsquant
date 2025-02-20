#ifndef  CHEBYSHEV_COEFFICIENTS
#define CHEBYSHEV_COEFFICIENTS

#include <complex>		/* for std::vector mostly class*/

using namespace std;

inline
double delta_chebF(double x, const double m)
{
	const double f0 =  sqrt(1.0 - x*x ) ;
	const double fm =  cos(m*acos(x));
	return fm/f0/ M_PI;
};

inline
complex<double> greenR_chebF(double x, const double m)
{
	const complex<double> I(0.0,1.0);
	const double f0 =  sqrt(1.0 - x*x ) ;
	const complex<double> fm =  pow( x - I*f0, m );
	return -I*fm/f0;
};

inline
complex<double> DgreenR_chebF(double x, const double m)
{
	const complex<double> I(0.0,1.0);
	const double f0 =  sqrt(1.0 - x*x ) ;
	const complex<double> fm =  pow( x - I*f0 , m )*(x + I*m*f0);
	return -I*fm/f0/f0/f0;
};

#endif
