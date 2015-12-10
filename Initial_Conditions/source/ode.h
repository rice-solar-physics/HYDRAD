// ****
// *
// * Include file for routines to calculate the source and
// * gradient terms in the coupled energy and momentum equations 
// * and to calculate n, T, P and Fc
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 19/12/2012
// *
// ****


double CalcdPbyds( double s, double n, int igdp, double **ppGravity );
double CalcdFcbyds( double EH, double ER );
double CalcdTbyds( double Fc, double T );

double Eheat( double s, double EH0, double sH0, double sH );