// ****
// *
// * Function bodies for the piece-wise polynomial fitting method
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 07/20/2020
// *
// ****


#include <stdio.h>
#include <stdlib.h>

#include "piecewisefit.h"
#include "../../regPoly/regpoly.h"
#include "../../regPoly/nrutil.h"
#include "../../../source/file.h"


// **** PIECEWISE FIT CLASS ****

// Constructor
CPieceWiseFit::CPieceWiseFit( char *pszInputFilename )
{
	OpenPieceWiseFit( pszInputFilename );
}

CPieceWiseFit::CPieceWiseFit( char *pszInputFilename, char *pszOutputFilename )
{
	GeneratePieceWiseFit( pszInputFilename, pszOutputFilename );
	OpenPieceWiseFit( pszOutputFilename );
}

// Destructor
CPieceWiseFit::~CPieceWiseFit( void )
{
	FreeAll();
}

void CPieceWiseFit::OpenPieceWiseFit( char *pszInputFilename )
{
	FILE *pInputFile;
	int i, j;
	
	// Open the input file
	pInputFile = fopen( pszInputFilename, "r" );
		// Get the order of the polynomial fit
		fscanf( pInputFile, "%i", &iPolyOrder );
		// Get the number of sub-domains
		fscanf( pInputFile, "%i", &inumSD );
		// There are inumSD + 1 sub-domain boundaries
		pfSDBoundary = (double*)malloc( sizeof(double) * ( inumSD + 1 ) );
		// Get the sub-domain boundaries
		for( i=0; i<(inumSD+1); i++ )
			ReadDouble( pInputFile, &(pfSDBoundary[i]) );
		// Allocate memory to store the coefficients and min/max 
		// values for the fitted quantity in each sub-domain
		ppfSDCoefficient = (double**)malloc( sizeof(double*) * inumSD );
		ppfSDMinMax = (double**)malloc( sizeof(double*) * inumSD );
		for( i=0; i<inumSD; i++ ) {
			ppfSDCoefficient[i] = (double*)malloc( sizeof(double) * ( iPolyOrder + 1 ) );
			ppfSDMinMax[i] = (double*)malloc( sizeof(double) * 2 );
		}
		// Get the coefficients and the min/max values 
		// for the fitted quantity in each sub-domain
		for( i=0; i<inumSD; i++ ) {
			for( j=0; j<(iPolyOrder+1); j++ )
				ReadDouble( pInputFile, &(ppfSDCoefficient[i][j]) );
			ReadDouble( pInputFile, &(ppfSDMinMax[i][0]) );
			ReadDouble( pInputFile, &(ppfSDMinMax[i][1]) );
		}
	fclose( pInputFile );
}

void CPieceWiseFit::GeneratePieceWiseFit( char *pszInputFilename, char *pszOutputFilename )
{
	FILE *pInputFile, *pOutputFile;
	double fMin, fMax;
	int inumPoints, iStart;
	int i, j, k;
	// **** FUNCTIONS AND VARIABLES FOR THE CURVE-FITTING CODE ****
	void BasisFuncs( double x, double *bfunc, int ma );
	double *x, *y, *sig, *a, chisq;
	double **u, **v, *w;
	int ndat, ma;
	// **** FUNCTIONS AND VARIABLES FOR THE CURVE-FITTING CODE ****

	// Open the input file
	pInputFile = fopen( pszInputFilename, "r" );
		// Get the order of the polynomial fit
		fscanf( pInputFile, "%i", &iPolyOrder );
		// Allocate memory for the coefficients
		ma = iPolyOrder + 1;
		a = vector( 1, ma );
		v = matrix( 1, ma, 1, ma );
		w = vector( 1, ma );
		// Get the number of sub-domains
		fscanf( pInputFile, "%i", &inumSD );
		// There are inumSD + 1 sub-domain boundaries
		pfSDBoundary = (double*)malloc( sizeof(double) * ( inumSD + 1 ) );
		// Get the sub-domain boundaries
		for( i=0; i<(inumSD+1); i++ )
			ReadDouble( pInputFile, &(pfSDBoundary[i]) );
		// Get the total number of data points
		fscanf( pInputFile, "%i", &inumPoints );
		// Allocate memory to store the data points
		x = vector( 1, (inumPoints+1) );
		y = vector( 1, (inumPoints+1) );
		sig = vector( 1, (inumPoints+1) );
		// Open the output file
		pOutputFile = fopen( pszOutputFilename, "w" );
			fprintf( pOutputFile, "%i\n", iPolyOrder );
			fprintf( pOutputFile, "%i\n", inumSD );
			fprintf( pOutputFile, "%.16e", pfSDBoundary[0] );
			for( i=1; i<(inumSD+1); i++ )
				fprintf( pOutputFile, "\t%.16e", pfSDBoundary[i] );
			// Generate a polynomial fit within each sub-domain
			j = 1;
			for( i=0; i<inumSD; i++ ) {
				iStart = j;
				if( !i ) {
					// Get the first pair of data points in the first sub-domain
					ReadDouble( pInputFile, &(x[j]) );
					ReadDouble( pInputFile, &(y[j]) );
					sig[j] = 1.0;
				}
				// Keep track of the minimum and maximum values in the sub-domain
				fMin = fMax = y[j];
				for( ;; ) {
					j++;
					// Get the next pair of data points
					ReadDouble( pInputFile, &(x[j]) );
					ReadDouble( pInputFile, &(y[j]) );
					sig[j] = 1.0;
					// Keep track of the minimum and maximum values in the sub-domain
					if( fMin > y[j] ) fMin = y[j];
					if( fMax < y[j] ) fMax = y[j];
					// If the upper boundary of the sub-domain is reached then break
					if( x[j] >= pfSDBoundary[i+1] ) break;
				}
				// Set the number of data points in the sub-domain
				ndat = j - iStart + 1;
#ifdef VERBOSE
				printf( "\niStart = %i, ndat = %i\n", iStart, ndat );
				for( k=iStart;k<=j; k++ )
					printf( "x[%i] = %g, y[%i] = %g\n", k, x[k], k, y[k] );
				printf( "fMin = %g, fMax = %g\n", fMin, fMax );
#endif // VERBOSE
				// Call the function that returns the coefficients of the best-fit polynomial (Single Value Decomposition)
				u = matrix( 1, ndat, 1, ma );
				svdfit( &(x[iStart-1]), &(y[iStart-1]), &(sig[iStart-1]), ndat, a, ma, u, v, w, &chisq, BasisFuncs );
				// Write the fitting coefficients and the min/max values for the sub-domain to the output file
				fprintf( pOutputFile, "\n%.16e", a[1] );
				for( k=2; k<=ma; k++ )
					fprintf( pOutputFile, "\t%.16e", a[k] );
				fprintf( pOutputFile, "\n%.16e\t%.16e", fMin, fMax );	
				free_matrix( u, 1, ndat, 1, ma );
			}
		fclose( pOutputFile );
		// **** FUNCTIONS AND VARIABLES FOR THE CURVE-FITTING CODE ****
		free_vector( sig, 1, (inumPoints+1) );
		free_vector( y, 1, (inumPoints+1) );
		free_vector( x, 1, (inumPoints+1) );
		// **** FUNCTIONS AND VARIABLES FOR THE CURVE-FITTING CODE ****
		free( pfSDBoundary );
		// **** FUNCTIONS AND VARIABLES FOR THE CURVE-FITTING CODE ****
		free_vector( w, 1, ma );
		free_matrix( v, 1, ma, 1, ma );
		free_vector( a, 1, ma );
		// **** FUNCTIONS AND VARIABLES FOR THE CURVE-FITTING CODE ****
	fclose( pInputFile );
}
// Define the basis function to use for the curve fit
void BasisFuncs( double x, double *bfunc, int ma )
{
	int i;
	// [ 1.0, x, x^2, x^3, ..., x^ma ]
	bfunc[1] = 1.0;
	for( i=2; i<=ma; i++ )
		bfunc[i] = bfunc[i-1] * x;
}

void CPieceWiseFit::FreeAll( void )
{
	int i;
	
	for( i=0; i<inumSD; i++ ) {
		free( ppfSDMinMax[i] );
		free( ppfSDCoefficient[i] );
	}
	free( ppfSDMinMax );
	free( ppfSDCoefficient );
	free( pfSDBoundary );
}

void CPieceWiseFit::ShowPieceWiseFit( void )
{
	int i, j;

	printf( "\n%i sub-domains\n", inumSD );
	printf( "Sub-domain boundaries: %g", pfSDBoundary[0] );
	for( i=1; i<(inumSD+1); i++ )
		printf( ", %g", pfSDBoundary[i] );
	printf( "\n" );
	for( i=0; i<inumSD; i++ ) {
		printf( "Sub-domain %i: %g", i, ppfSDCoefficient[i][0] );
		for( j=1; j<(iPolyOrder+1); j++ )
			printf( ", %g", ppfSDCoefficient[i][j] );
		printf( "\t[%g,%g]\n", ppfSDMinMax[i][0], ppfSDMinMax[i][1] );
	}
}

double CPieceWiseFit::GetPieceWiseFit( double fx )
{
	double fu, fQ;
	int iSubDomain;
	int i;
	
	// Find the sub-domain in which to perform this polynomial fit
	for( i=1; i<(inumSD+0); i++ )
		if( fx < pfSDBoundary[i] ) break;
		
	iSubDomain = i - 1;
#ifdef VERBOSE
	printf( "\nFitting within sub-domain %i\n", iSubDomain );
	printf( "Using coefficients: %g", ppfSDCoefficient[iSubDomain][0] );
	for( i=1; i<(iPolyOrder+1); i++ )
		printf( ", %g", ppfSDCoefficient[iSubDomain][i] );
	printf( "\t[%g,%g]\n", ppfSDMinMax[iSubDomain][0], ppfSDMinMax[iSubDomain][1] );
#endif // VERBOSE
	// Calculate the fitted quantity using the coefficients for the sub-domain
	fu = fx;
	fQ = ppfSDCoefficient[iSubDomain][0];
#ifdef VERBOSE
	printf( "fQ[0] = %g x 1.0 = %g\n", ppfSDCoefficient[iSubDomain][0], fQ );
#endif // VERBOSE
	for( i=1; i<(iPolyOrder+1); i++ ) {
		fQ += ppfSDCoefficient[iSubDomain][i] * fu;
#ifdef VERBOSE
		printf( "fQ[%i] = %g x %g = %g\n", i, ppfSDCoefficient[iSubDomain][i], fu, ppfSDCoefficient[iSubDomain][i] * fu );
#endif // VERBOSE
		fu *= fx;
	}
#ifdef VERBOSE
	printf( "fQ = %g\n", fQ );
#endif // VERBOSE
	// Ensure the fitted quantity is constrained within the sub-domain limits
	if( fQ < ppfSDMinMax[iSubDomain][0] ) fQ = ppfSDMinMax[iSubDomain][0];
	else if( fQ > ppfSDMinMax[iSubDomain][1] ) fQ = ppfSDMinMax[iSubDomain][1];
#ifdef VERBOSE
	printf( "f(%g) = %g\n", fx, fQ );
#endif // VERBOSE
	return fQ;
}

double CPieceWiseFit::GetDerivative( double fx, int iDerivative )
{
	if( iDerivative > iPolyOrder ) return 0.0;

	// Function to find factorial of given number 
	unsigned int factorial(unsigned int n);
	double fu, fQ;
	int iSubDomain;
	int i;
	
	// Find the sub-domain in which to perform this polynomial fit
	for( i=1; i<(inumSD+0); i++ )
		if( fx < pfSDBoundary[i] ) break;
		
	iSubDomain = i - 1;
#ifdef VERBOSE
	printf( "\nFitting within sub-domain %i\n", iSubDomain );
	printf( "Using coefficients: %g", ppfSDCoefficient[iSubDomain][0] );
	for( i=1; i<(iPolyOrder+1); i++ )
		printf( ", %g", ppfSDCoefficient[iSubDomain][i] );
	printf( "\t[%g,%g]\n", ppfSDMinMax[iSubDomain][0], ppfSDMinMax[iSubDomain][1] );
#endif // VERBOSE
	// Calculate the fitted quantity using the coefficients for the sub-domain
	fu = fx;
	fQ = (double)factorial(iDerivative) * ppfSDCoefficient[iSubDomain][iDerivative];
#ifdef VERBOSE
	printf( "fQ[%i] = %g x %g = %g\n", iDerivative, (double)factorial(iDerivative), ppfSDCoefficient[iSubDomain][iDerivative], fQ );
#endif // VERBOSE
	for( i=(iDerivative+1); i<(iPolyOrder+1); i++ ) {
		fQ += ( (double)(factorial(i)/factorial(i-iDerivative)) * ppfSDCoefficient[iSubDomain][i] * fu );
#ifdef VERBOSE
		printf( "fQ[%i] = %g x %g x %g = %g\n", i, (double)(factorial(i)/factorial(i-iDerivative)), ppfSDCoefficient[iSubDomain][i], fu, (double)(factorial(i)/factorial(i-iDerivative)) * ppfSDCoefficient[iSubDomain][i] * fu );
#endif // VERBOSE
		fu *= fx;
	}
#ifdef VERBOSE
	printf( "f(%g) = %g\n", fx, fQ );
#endif // VERBOSE
	return fQ;
}
// Function to find factorial of given number 
unsigned int factorial(unsigned int n) 
{ 
	// Single line of code to find factorial 
    return (n == 1 || n == 0) 
               ? 1 
               : n * factorial(n - 1); 
}