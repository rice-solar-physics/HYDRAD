// ****
// *
// * Routines from Numerical Recipes in C to
// * calculate the lower incomplete gamma function
// *
// * Dr. Stephen J. Bradshaw
// *
// * Date last modified: 21/12/2012
// *
// ****


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gamma.h"


double gammp(double a, double x)
{
void gcf(double *gammcf, double a, double x, double *gln);
void gser( double *gamser, double a, double x, double *gln);
void nrerror2(char error_text[]);	// Defined in fitpoly.cpp
double gamser,gammcf,gln;

if (x < 0.0 || a <= 0.0) nrerror2((char *)"Invalid arguments in routine gammp");
if(x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
#ifdef REGULARISED
    return gamser;
#endif // REGULARISED
#ifdef LOWER
    return gamser * exp( gln );
#endif // LOWER
} else {
    gcf(&gammcf,a,x,&gln);
#ifdef REGULARISED
    return 1.0-gammcf;
#endif // REGULARISED
#ifdef LOWER
    return (1.0-gammcf) * exp( gln );
#endif // LOWER
}
}


double gammln(double xx)
{
double x,y,tmp,ser;
static double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5};
int j;

y=x=xx;
tmp=x+5.5;
tmp -= (x+0.5)*log(tmp);
ser=1.000000000190015;
for (j=0;j<=5;j++)
    ser += cof[j]/++y;

return -tmp+log(2.5066282746310005*ser/x);
}

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(double *gammcf, double a, double x, double *gln)
{
    double gammln(double xx);
    void nrerror2(char error_text[]);	// Defined in fitpoly.cpp
    int i;
    double an,b,c,d,del,h;

    *gln=gammln(a);
    b=x+1.0-a;
    c=1.0/FPMIN;
    d=1.0/b;
    h=d;
    for (i=1;i<=ITMAX;i++) {
        an = -i*(i-a);
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
    if (i > ITMAX) nrerror2((char *)"a too large, ITMAX too small in gcf");
    *gammcf=exp(-x+a*log(x)-(*gln))*h;
}

void gser(double *gamser, double a, double x, double *gln)
{
double gammln(double xx);
void nrerror2(char error_text[]);	// Defined in fitpoly.cpp
int n;
double sum,del,ap;

*gln=gammln(a);
if (x <= 0.0) {
    if (x < 0.0) nrerror2((char *)"x less than 0 in routine gser");
    *gamser=0.0;
    return;
} else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
        ++ap;
	del *= x/ap;
	sum += del;
	if (fabs(del) < fabs(sum)*EPS) {
            *gamser=sum*exp(-x+a*log(x)-(*gln));
            return;
	}
    }
    nrerror2((char *)"a too large, ITMAX too small in routine gser");
    return;
}
}
