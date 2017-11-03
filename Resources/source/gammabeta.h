/* Taken from "Numerical Recipes in C" by Press et al.
	Evaluates Euler's Gamma function and the complete and incomplete Beta functions */

double lnGamma(double xx);
double Beta(double z, double w);
double betacf(double a, double b, double x);
double incompleteBeta(double a, double b, double x);

/*
double gamma_p(double a, double x);
double gamma_q(double a, double x);
void gamma_ser(double *gamser, double a, double x, double *gln);
void gamma_cf(double *gammcf, double a, double x, double *gln);
*/