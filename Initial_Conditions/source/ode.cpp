// ****
// *
// * Routines to calculate the source and gradient terms 
// * in the coupled energy and momentum equations
// * and to calculate n, T, P and Fc
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 01/12/2016
// *
// ****


#include <math.h>

#include "config.h"
#include  "../../Resources/source/fitpoly.h"
#include "../../Resources/source/constants.h"
#include "../../Resources/Utils/regPoly/regpoly.h"


double CalcdPbyds( double s, double n, double Lfull, double *pfGravityCoefficients )
{
double CalcSolarGravity( double s, double *a );

return AVERAGE_PARTICLE_MASS * n * CalcSolarGravity( s/Lfull, pfGravityCoefficients );
}

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

double CalcdFcbyds( double EH, double ER )
{
return ( EH + ER );
}

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