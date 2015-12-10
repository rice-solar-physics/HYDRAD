// ****
// *
// * Routines to calculate the source and gradient terms 
// * in the coupled energy and momentum equations
// * and to calculate n, T, P and Fc
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 19/12/2012
// *
// ****


#include <math.h>

#include "config.h"
#include  "../../Resources/source/fitpoly.h"
#include "../../Resources/source/constants.h"


double CalcdPbyds( double s, double n, int igdp, double **ppGravity )
{
double CalcSolarGravity( double s, int igdp, double **ppGravity );

return AVERAGE_PARTICLE_MASS * n * CalcSolarGravity( s, igdp, ppGravity );
}

double CalcSolarGravity( double s, int igdp, double **ppGravity )
{
double x[3], y[3], g_parallel;
int i;

for( i=0; i<igdp; i++ )
	if( s < ppGravity[i][0] ) break;

if( i == 0 ) i = 1;
if( i == igdp ) i = igdp - 1;

x[1] = ppGravity[i-1][0];
x[2] = ppGravity[i][0];

y[1] = ppGravity[i-1][1];
y[2] = ppGravity[i][1];

LinearFit( x, y, s, &g_parallel );

return g_parallel;
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