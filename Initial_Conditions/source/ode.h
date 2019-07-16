// ****
// *
// * Include file for routines to calculate the source and
// * gradient terms in the coupled energy and momentum equations 
// * and to calculate n, T, P and Fc
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 07/16/2019
// *
// ****

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
double CalcCrossSection( double s, double *pfMagneticFieldCoefficients );
double CalcdAbyds( double s, double Lfull, double *pfMagneticFieldCoefficients );
double CalcdPbyds( double s, double n, double P, double Lfull, double *pfGravityCoefficients, double *pfMagneticFieldCoefficients );
double CalcdFcbyds( double s, double EH, double ER, double Fc, double Lfull, double *pfMagneticFieldCoefficients );
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
double CalcdPbyds( double s, double n, double Lfull, double *pfGravityCoefficients );
double CalcdFcbyds( double EH, double ER );
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

double CalcdTbyds( double Fc, double T );

double Eheat( double s, double EH0, double sH0, double sH );