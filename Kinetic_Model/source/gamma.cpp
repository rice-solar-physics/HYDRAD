// ****
// *
// * Routines from Numerical Recipes in C to
// * calculate the lower incomplete gamma function
// *
// * Dr. Stephen J. Bradshaw
// *
// * Date last modified: 12/08/2022
// *
// ****

// #include <iostream>
#include <cmath>

// #include "message.hpp"

// Lower incomplete gamma function
double gammp( double a, double x )
{
    void gcf( double *gammcf, double a, double x, double *gln );
    void gser( double *gamser, double a, double x, double *gln );
    double gamser, gammcf, gln;

    if( x < 0.0 || a <= 0.0 ) {
//        message( "\nInvalid arguments in routine gammp.\n" );
        exit( 1 );
    }

    if( x < ( a + 1.0 ) ) {
        gser( &gamser, a, x, &gln );
        return gamser * exp( gln );
    } else {
        gcf( &gammcf, a, x, &gln );
        return ( 1.0 - gammcf ) * exp( gln );
    }
}

// Derivative (w.r.t x) of the lower incomplete gamma function
double gammpPrime( double a, double x )
{
    return pow( x, ( a - 1.0 ) ) * exp( -x );
}

double gammln( double xx )
{
    double x, y, tmp, ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
                          24.01409824083091,-1.231739572450155,
                          0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y = x = xx;
    tmp = x + 5.5;
    tmp -= ( x + 0.5 ) * log( tmp );
    ser = 1.000000000190015;
    for ( j=0; j<=5; j++ )
        ser += cof[j] / ++y;

    return -tmp + log( 2.5066282746310005 * ser / x );
}

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
void gcf( double *gammcf, double a, double x, double *gln )
{
    double an, b, c, d, del, h;
    int i;

    *gln = gammln( a );
    b = x + 1.0 - a;
    c = 1.0 / FPMIN;
    d = 1.0 / b;
    h = d;
    for ( i=1; i<=ITMAX; i++ ) {
        an = -i * ( i - a );
    	b += 2.0;
	    d = an * d + b;
	    if( fabs( d ) < FPMIN ) d = FPMIN;
	    c = b + an / c;
	    if( fabs( c ) < FPMIN ) c = FPMIN;
	    d = 1.0 / d;
	    del = d * c;
	    h *= del;
	    if( fabs( del - 1.0 ) < EPS ) break;
    }

    if( i > ITMAX ) {
//        message( "\na too large, ITMAX too small in gcf.\n" );
        exit( 1 );
    }

    *gammcf = exp( -x + a * log( x ) - (*gln) ) * h;
}

void gser(double *gamser, double a, double x, double *gln)
{
    double sum, del, ap;
    int n;

    *gln = gammln( a );

    if( x <= 0.0 ) {
        if( x < 0.0 ) {
//            message( "\nx less than 0 in routine gser.\n" );
            exit( 1 );
        }
        *gamser = 0.0;
        return;
    } else {
        ap = a;
        del = sum = 1.0 / a;
        for( n=1; n<=ITMAX; n++) {
            ++ap;
	        del *= x / ap;
	        sum += del;
	        if( fabs( del ) < fabs( sum ) * EPS ) {
                *gamser = sum * exp( -x +a * log( x ) - (*gln) );
                return;
	        }
        }
//        message( "\na too large, ITMAX too small in routine gser.\n" );
        exit( 1 );
    }
}