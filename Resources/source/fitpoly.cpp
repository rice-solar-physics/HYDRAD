// ****
// *
// * Routines from Numerical Recipes in C to perform polynomial
// * interpolation in one and two dimensions. Also contains some
// * of my own functions for linear interpolation and finding
// * maxima and minima
// *
// * Dr. Stephen J. Bradshaw
// *
// * Date last modified: 20/12/2012
// *
// ****


#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fitpoly.h"


// ** Numerical Recipes standard error handler **
void nrerror(char error_text[])
{
fprintf(stderr,"\nNumerical Recipes run-time error...\n");
fprintf(stderr,"%s\n",error_text);
fprintf(stderr,"...now exiting to system...\n");
exit(1);
}

// ** Routine to perform a linear interpolation / extrapolation **
void LinearFit( double xa[], double ya[], double x, double *y )
{
double fGrad;

fGrad = ( ya[2] - ya[1] ) / ( xa[2] - xa[1] );

*y = ( fGrad * ( x - xa[2] ) ) + ya[2];
}

// ** Routine to calculate the degree N-1 polynomial fit to N data points **
void FitPolynomial(double xa[], double ya[], int n, double x, double *y, double *dy)
{
int i,m,ns=1;
double den,dif,dift,ho,hp,w;
double *c,*d;
    
dif=fabs(x-xa[1]);
c=(double*)alloca((size_t) ((n+1)*sizeof(double)));
d=(double*)alloca((size_t) ((n+1)*sizeof(double)));
for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
        ns=i;
    	dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
}
*y=ya[ns--];
for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
        ho=xa[i]-x;
    	hp=xa[i+m]-x;
    	w=c[i+1]-d[i];
    	if ( (den=ho-hp) == 0.0) nrerror((char*)"Error in routine FitPolynomial");
    	den=w/den;
    	d[i]=hp*den;
    	c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
}
}

// ** Routine to calculate the degree 3 polynomial fit to 4 data points **
void FitPolynomial4(double xa[], double ya[], double x, double *y, double *dy)
{
int i,m,ns=1;
double den,dif,dift,ho,hp,w;
double *c,*d;
    
dif=fabs(x-xa[1]);
c=(double*)alloca((size_t) (5*sizeof(double)));
d=(double*)alloca((size_t) (5*sizeof(double)));
for (i=1;i<=4;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
        ns=i;
        dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
}
*y=ya[ns--];
for (m=1;m<4;m++) {
    for (i=1;i<=4-m;i++) {
        ho=xa[i]-x;
    	hp=xa[i+m]-x;
    	w=c[i+1]-d[i];
    	if ( (den=ho-hp) == 0.0) nrerror((char*)"Error in routine FitPolynomial");
    	den=w/den;
    	d[i]=hp*den;
    	c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (4-m) ? c[ns+1] : d[ns--]));
}

// Ensure a stable interpolation for n=4 (introduce no new maxima or minima)
if( (*y) > min(ya[2],ya[3]) && (*y) < max(ya[2],ya[3]) )
    return;

LinearFit( &(xa[1]), &(ya[1]), x, y );
}

// ** Routine to perform a two dimensional interpolation **
void FitPolynomial2D( double x1a[], double x2a[], double **ya, int m, int n, double x1, double x2, double *y, double *dy )
{
int j;
double *ymtmp;

ymtmp=(double*)alloca((size_t) ((m+1)*sizeof(double)));

for( j=1; j<=m; j++ )
    FitPolynomial( x2a, ya[j], n, x2, &ymtmp[j], dy );
	
FitPolynomial( x1a, ymtmp, m, x1, y, dy );
}

// ** Routine to find the larger of two values **
double max( double a, double b )
{
if( a > b ) return a;

return b;
}

// ** Routine to find the smaller of two values **
double min( double a, double b )
{
if( a < b ) return a;

return b;
}
