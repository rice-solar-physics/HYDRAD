/* Taken from "Numerical Recipes in C" by Press et al.
	Evaluates Euler's Gamma function and the complete and incomplete Beta functions */


#include <math.h>
#include <stdio.h>
#include "gammabeta.h"


// Return the natural log of Euler's Gamma function: ln Gamma(xx)
double lnGamma(double xx)
{
	double x,y,tmp,ser;
	double cof[6];
	cof[0]=76.18009172947146; cof[1]=-86.50532032941677; cof[2]=24.01409824083091; cof[3]=-1.231739572450155; cof[4]=0.1208650973866179e-2; cof[5]=-0.5395239384953e-5;
	int j;

	y = xx;
	x = xx;
	tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser = 1.000000000190015;
	for(j=0; j<=5; j++) ser += cof[j]/++y;
	
	return -tmp+log(2.5066282746310005*ser/x);
}


// Return the beta function: B(z,w) = Gamma(z)Gamma(w) / (Gamma(z+w)
double Beta(double z, double w)
{
	if (exp( lnGamma(z) + lnGamma(w) - lnGamma(z+w) ) < 0) printf("Routine Beta returning negative value.\n");
	return exp( lnGamma(z) + lnGamma(w) - lnGamma(z+w) );
}


/* 
    Two incomplete beta functions:
	(1) a continued fraction evaluation routine
	(2) evaluation of the incomplete beta function
*/

#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

// Evaluates continued fraction for incomplete beta function by modified Lentz's method
double betacf(double a, double b, double x)
{
	int m, m2;
	double aa,c,d,del,h,qab,qam,qap;

	qab = a + b;	
	qap = a + 1.0;
	qam = a - 1.0;
	c = 1.0;
	d = 1.0 - qab * x / qap;
	if(fabs(d) < FPMIN) d = FPMIN;
	d = 1.0 / d;
	h = d;
	
	for(m=1; m<=MAXIT; m++) 
	{
		m2 = 2 * m;
		aa = m*(b-m)*x/((qam + m2)*(a + m2));
		d = 1.0 + aa * d;
		if(fabs(d) < FPMIN) d = FPMIN;
		c = 1.0 + aa / c;
		if(fabs(c) < FPMIN) c = FPMIN;
		d = 1.0 / d;
		h *= d*c;
		aa = -(a+m) * (qab+m) * x / ((a+m2)*(qap+m2));
		
		d = 1.0 + aa * d;
		if(fabs(d) < FPMIN) d = FPMIN;
		c = 1.0 + aa / c;
		if(fabs(c) < FPMIN) c = FPMIN;
		d = 1.0 / d;
		del = d*c;
		h *= del;
		if( fabs(del-1.0) < EPS) break;
	}
	
	if (m > MAXIT) printf("a or b too big, or MAXIT too small in betacf\na %e\tb %e\tx %e\n", a, b, x);
	return h;	
}

// Returns the (normalized) incomplete beta function
//	To obtain the un-normalized function, multiply by the complete beta function: B(a,b)
double incompleteBeta( double a, double b, double x)
{
	double bt;
	
	if( x < 0.0 || x > 1.0) printf("Bad x value in routine incompleteBeta\n");	

	if( x == 0.0 || x == 1.0) bt = 0.0;

	else
		bt = exp(lnGamma(a+b) - lnGamma(a) - lnGamma(b) + a * log(x) + b * log(1.0 - x));

	if( x <  (a+1.0)/(a+b+2.0) )
	{
		// if (bt*betacf(a,b,x)/a < 0) printf("Routine incompleteBeta returning negative value. (1)\n");
		return bt*betacf(a,b,x)/a;
	}
	else 
	{
		// if (1.0 - bt*betacf(b,a,1.0-x)/b < 0) printf("Routine incompleteBeta returning negative value. (2)\n");
		return 1.0 - bt*betacf(b,a,1.0-x)/b;
	}
}

/*
// Returns the incomplete gamma function: P(a,x)
double gamma_p(double a, double x)
{
    //void gcf(double *gammcf, double a, double x, double *gln);
    //void gser(double *gammser, double a, double x, double *gln);
    //void nrerror(char error_text[]);
    double gamser, gammcf, gln;
    
    if ( x < 0.0 || a <= 0.0 ) printf("Invalid arguments in routine gammp\n");
    if ( x < (a + 1.0)) {
        gamma_ser(&gamser,a,x,&gln);
        return gamser;
    } else {
        gamma_cf(&gammcf,a,x,&gln);
        return 1.0-gammcf;
    }
}

// Returns the incomplete gamma function Q(a,x) = 1 - P(a,x)
double gamma_q(double a, double x)
{
    //void gcf(double *gammcf, double a, double x, double *gln);
    //void gser(double *gammser, double a, double x, double *gln);
    //void nrerror(char error_text[]);
    double gamser, gammcf, gln;
    
    if( x< 0.0 || a <= 0.0) printf("Invalid arguments in routine gammq\n");
    if( x < (a+1.0)) {
        gamma_ser(&gamser,a,x,&gln);
        return 1.0-gamser;
    } else {
        gamma_cf(&gammcf,a,x,&gln);
        return gammcf;
    }
}

// Returns the incomplete gamma function P(a,x) evaluated by its series representation as gamser
// Also returns ln Gamma(a) as gln
void gamma_ser(double *gamser, double a, double x, double *gln)
{
    double lnGamma(double xx);
    void nrerror(char error_text[]);
    
    int n;
    double sum, del, ap;
    
    *gln = lnGamma(a);
    if( x <= 0.0) {
        if( x < 0.0) printf("x less than 0 in routine gser\n");
        *gamser=0.0;
        return;
    } else {
        ap = a;
        del=sum=1.0/a;
        for( n=1; n<=MAXIT; ++n) {
            ++ap;
            del *= x/ap;
            sum += del;
            if( fabs(del) < fabs(sum)*EPS) {
                *gamser=sum*exp(-x+a*log(x)-(*gln));
                return;
            }
        }
        printf("a too large, MAXIT too small in routine gser\n");
        return;
    }
}

// Returns the incomplete gamma function Q(a,x) evaluated by its continued fraction representation as gammcf
// Also returns ln Gamma(a) as gln
void gamma_cf(double *gammcf, double a, double x, double *gln)
{
    double lnGamma(double xx);
    void nrerror(char error_text[]);
    
    int i;
    double an,b,c,d,del,h;
    
    *gln=lnGamma(a);
    b=x+1.0-a;
    c=1.0/FPMIN;
    d=1.0/b;
    h=d;
    for( i=1;i<=MAXIT;++i) {
        an= -i*(i-a);
        b += 2.0;
        d=an*d+b;
        if (fabs(d) < FPMIN) d=FPMIN;
        c=b+an/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        del=d*c;
        h *= del;
        if (fabs(del-1.0) < EPS) break;
    }
    if( i > MAXIT) printf("a too large, MAXIT too small in gcf\n");
    *gammcf=exp(-x+a*log(x)-(*gln))*h;
}
*/