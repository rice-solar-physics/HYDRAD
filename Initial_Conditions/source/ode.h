// ****
// *
// * Include file for routines to calculate the source and
// * gradient terms in the coupled energy and momentum equations 
// * and to calculate n, T, P and Fc
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 07/20/2020
// *
// ****

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
double CalcCrossSection( double x, PPIECEWISEFIT pMagneticField );
double CalcdAbyds( double x, double Lfull, PPIECEWISEFIT pMagneticField );
double CalcdFcbyds( double s, double EH, double ER, double Fc, double Lfull, PPIECEWISEFIT pMagneticField );
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
double CalcdFcbyds( double EH, double ER );
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

double CalcdPbyds( double s, double n, double Lfull, PPIECEWISEFIT pGravity );
double CalcdTbyds( double Fc, double T );

double Eheat( double s, double EH0, double sH0, double sH );