// ****
// *
// * Routines to calculate the source and gradient terms 
// * in the coupled energy and momentum equations
// * and to calculate n, T, P and Fc
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 07/20/2020
// *
// ****


#include <math.h>

#include "config.h"
#include "../../Resources/Utils/generatePieceWiseFit/source/piecewisefit.h"
#include "ode.h"
#include  "../../Resources/source/fitpoly.h"
#include "../../Resources/source/constants.h"


double CalcSolarGravity( double x, PPIECEWISEFIT pGravity )
{
return pGravity->GetPieceWiseFit( x );
}

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
double CalcCrossSection( double x, PPIECEWISEFIT pMagneticField )
{
// The cross-section area varies inversely with the magnetic field strength
return ( 1.0 / pMagneticField->GetPieceWiseFit( x ) );
}

double CalcdAbyds( double x, double Lfull, PPIECEWISEFIT pMagneticField )
{
double dBbydx, A;

dBbydx = pMagneticField->GetDerivative( x, 1 );
A = CalcCrossSection( x, pMagneticField );

return ( -((A*A)/Lfull) * dBbydx );
}

double CalcdFcbyds( double s, double EH, double ER, double Fc, double Lfull, PPIECEWISEFIT pMagneticField )
{
double term1;

term1 = - ( Fc / CalcCrossSection( s/Lfull, pMagneticField ) ) * CalcdAbyds( s/Lfull, Lfull, pMagneticField );

return ( EH + ER + term1 );
}
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
double CalcdFcbyds( double EH, double ER )
{
return ( EH + ER );
}
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

double CalcdPbyds( double s, double n, double Lfull, PPIECEWISEFIT pGravity )
{
return AVERAGE_PARTICLE_MASS * n * CalcSolarGravity( s/Lfull, pGravity );
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