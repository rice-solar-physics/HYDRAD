// ****
// *
// * Include file for routines to calculate the source and
// * gradient terms in the coupled energy and momentum equations 
// * and to calculate n, T, P and Fc
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 01/12/2016
// *
// ****


double CalcdPbyds( double s, double n, double Lfull, double *pfGravityCoefficients );
double CalcdFcbyds( double EH, double ER );
double CalcdTbyds( double Fc, double T );

double Eheat( double s, double EH0, double sH0, double sH );