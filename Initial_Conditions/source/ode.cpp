// ****
// *
// * Routines to calculate the source and gradient terms 
// * in the coupled energy and momentum equations
// * and to calculate n, T, P and Fc
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 07/16/2019
// *
// ****


#include <math.h>

#include "config.h"
#include "ode.h"
#include  "../../Resources/source/fitpoly.h"
#include "../../Resources/source/constants.h"
#include "../../Resources/Utils/regPoly/regpoly.h"


double CalcSolarGravity( double s, double *pfGravityCoefficients )
{
double sum = 0.0;
int i;

sum += ( 1.0 * pfGravityCoefficients[0] );
for( i=1; i<(POLY_ORDER+1); i++ ) {
	sum += ( pow(s,i) * pfGravityCoefficients[i] );
}

return sum;
}

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
double CalcCrossSection( double s, double *pfMagneticFieldCoefficients )
{
double B = 0.0;
int i;

B += ( 1.0 * pfMagneticFieldCoefficients[0] );
for( i=1; i<(POLY_ORDER+1); i++ ) {
	B += ( pfMagneticFieldCoefficients[i] * pow(s,i) );
}

// The cross-section area varies inversely with the magnetic field strength
return ( 1.0 / B );
}

double CalcdAbyds( double s, double Lfull, double *pfMagneticFieldCoefficients )
{
double dBbydx = 0.0, A;
int i;

dBbydx += ( 1.0 * pfMagneticFieldCoefficients[1] );
for( i=2; i<(POLY_ORDER+1); i++ ) {
	dBbydx += ( ((double)i) * pfMagneticFieldCoefficients[i] * pow(s,(i-1)) );
}

A = CalcCrossSection( s, pfMagneticFieldCoefficients );

return ( -((A*A)/Lfull) * dBbydx );
}

double CalcdPbyds( double s, double n, double P, double Lfull, double *pfGravityCoefficients, double *pfMagneticFieldCoefficients )
{
double term1, term2;

term1 = AVERAGE_PARTICLE_MASS * n * CalcSolarGravity( s/Lfull, pfGravityCoefficients );
term2 = - ( P / CalcCrossSection( s/Lfull, pfMagneticFieldCoefficients ) ) * CalcdAbyds( s/Lfull, Lfull, pfMagneticFieldCoefficients );

return ( term1 + term2 );
}

double CalcdFcbyds( double s, double EH, double ER, double Fc, double Lfull, double *pfMagneticFieldCoefficients )
{
double term1;

term1 = - ( Fc / CalcCrossSection( s/Lfull, pfMagneticFieldCoefficients ) ) * CalcdAbyds( s/Lfull, Lfull, pfMagneticFieldCoefficients );

return ( EH + ER + term1 );
}
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
double CalcdPbyds( double s, double n, double Lfull, double *pfGravityCoefficients )
{
return AVERAGE_PARTICLE_MASS * n * CalcSolarGravity( s/Lfull, pfGravityCoefficients );
}

double CalcdFcbyds( double EH, double ER )
{
return ( EH + ER );
}
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

double CalcdTbyds( double Fc, double T )
{
double term1;

term1 = SPITZER_SINGLE_FLUID_CONDUCTIVITY * pow( T, (2.5) );

return ( - Fc / term1 );
}

double Eheat( double s, double H0, double sH0, double sH )
{
double term1;

term1 = s - sH0;
term1 *= term1;

return ( H0 * exp( - term1 / (2.0*sH*sH) ) );
}