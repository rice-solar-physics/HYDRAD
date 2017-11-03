// ****
// *
// * Header file for the routines from Numerical Recipes in C to 
// * perform a polynomial fit to a set of data points using general 
// * linear least squares fitting (Ch. 15.4) and solution by use of 
// * the normal equations
// *
// * Dr. Stephen J. Bradshaw
// *
// * Date last modified: 01/12/2016
// *
// ****


#define POLY_ORDER	6


// ** Routine to return the values of the polynomial coefficients a[1..ma] **
void lfit(double x[], double y[], double sig[], int ndat, double a[], int ia[], int ma, double **covar, double *chisq, void (*funcs)(double, double [], int));


// ** Routine to return the values of the polynomial coefficients a[1..ma] (using Single Value Decomposition) **
void svdfit(double x[], double y[], double sig[], int ndata, double a[], int ma, double **u, double **v, double w[], double *chisq, void (*funcs)(double, double [], int));

// ** To evaluate the covariance matrix cvm[1..ma][1..ma] of the fit for ma parameters obtained by svdfit, call this routine with matrices v[1..ma][1..ma], w[1..ma] as returned from svdfit **
void svdvar(double **v, int ma, double w[], double **cvm);
