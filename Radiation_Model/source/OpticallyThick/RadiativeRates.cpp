// ****
// *
// * Function bodies for the class definition of the radiative rates
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/06/2020
// *
// ****

#include "RadiativeRates.h"

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE

// **** RADIATIVE TRANSITION RATES CLASS ****

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <math.h>

#include "../../../Resources/source/fitpoly.h"
#include "../../../Resources/source/file.h"
// #include "../../../Resources/source/svd.h"


#define NR_END 1
#define FREE_ARG char*
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


double pow10( double index )
{
    return pow(10., index);
}


// Constructor
CRadiativeRates::CRadiativeRates( char *pszRatesFiles )
{
	Initialise( pszRatesFiles );
}

// Destructor
CRadiativeRates::~CRadiativeRates( void )
{
	FreeAll();
}

// Initialisation routine
void CRadiativeRates::Initialise( char *pszRatesFiles )
{
FILE *pRatesFiles;
char szRatesFile[256];
FILE *pTRANSFile, *pTERMSFile;
char szBuffer[256];
int iBuffer, i;
    
// Open the file containing the list of radiative transition rate files to use
pRatesFiles = fopen( pszRatesFiles, "r" );
	// Get the bound-bound rates
	fscanf( pRatesFiles, "%s", szRatesFile );
	GetBBRates( szRatesFile );
	// Get the bound-free rates
	fscanf( pRatesFiles, "%s", szRatesFile );
	GetBFRates( szRatesFile );
	// Get the free-bound rates
	fscanf( pRatesFiles, "%s", szRatesFile );
	GetFBRates( szRatesFile );
	// Get the collisional rates
	fscanf( pRatesFiles, "%s", szRatesFile );
	GetCollRates( szRatesFile );
fclose( pRatesFiles );
    
   
    // Get the information about the transitions
    pTRANSFile = fopen( "Radiation_Model/atomic_data/OpticallyThick/radiative_rates/transitions.txt", "r" );
    
    fscanf( pTRANSFile, "%i", &iNBBT );		// Number of bound-bound transitions
    fscanf( pTRANSFile, "%i", &iNBFT );		// Number of bound-free transitions
    
    // Allocate memory for the transition data
    pTrt =          (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
    ppScaledTrt = (double**)malloc( sizeof(double) * 2 );
		for( i=0; i<2; i++ )
			ppScaledTrt[i] = (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
    pnu0 =          (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
    pTeZ_c  =       (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
    pZ_c_LEFT =     (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
    pZ_c_RIGHT =    (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
    pMcZ_c_LEFT =   (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
    pMcZ_c_RIGHT =  (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
   
    // Column labels
    fscanf( pTRANSFile, "%s", szBuffer );
    fscanf( pTRANSFile, "%s", szBuffer );
    fscanf( pTRANSFile, "%s", szBuffer );
    fscanf( pTRANSFile, "%s", szBuffer );
    fscanf( pTRANSFile, "%s", szBuffer );
    
    for( i=0; i<iNBBT+iNBFT; i++ )
    {
        fscanf( pTRANSFile, "%i", &iBuffer );	// i
        fscanf( pTRANSFile, "%i", &iBuffer );	// j
        ReadDouble( pTRANSFile, &(pTrt[i]) );	// T_rad^top
			ppScaledTrt[0][i] = ppScaledTrt[1][i] = pTrt[i];
        ReadDouble( pTRANSFile, &(pnu0[i]) );	// Transition rest frequency (wavelength)
        ReadDouble( pTRANSFile, &(pTeZ_c[i]) );	// Temperature at Z_c
    }
    
    fclose( pTRANSFile );
    
    // Get the equation terms needed to calculate T_rad for each transition as a function of 's'
    pTERMSFile = fopen( "Radiation_Model/atomic_data/OpticallyThick/radiative_rates/terms.txt", "r" );
    
    fscanf( pTERMSFile, "%i", &iNBBT );		// Number of bound-bound transitions
    fscanf( pTERMSFile, "%i", &iNBFT );		// Number of bound-free transitions
    
    // Allocate memory for the transition data
    pterm1 = (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
    pterm2 = (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
    pH = (double*)malloc( sizeof(double) * (iNBBT+iNBFT) );
    
    // Column labels
    fscanf( pTERMSFile, "%s", szBuffer );
    fscanf( pTERMSFile, "%s", szBuffer );
    fscanf( pTERMSFile, "%s", szBuffer );
    
    for( i=0; i<iNBBT+iNBFT; i++ )
    {
        ReadDouble( pTERMSFile, &(pterm1[i]) );
        ReadDouble( pTERMSFile, &(pterm2[i]) );
        ReadDouble( pTERMSFile, &(pH[i]) );
    }
    
    fclose( pTERMSFile );
}

void CRadiativeRates::FreeAll( void )
{
int i, j;

// Free the memory allocated to the bound-bound rates
free( pfBB_logT );
for( j=0; j<iBB_TRvals; j++ )
{
	free( ppfBB_lu[j] );
	free( ppfBB_ul[j] );
}
free( ppfBB_lu );
free( ppfBB_ul );

// Free the memory allocated to the bound-free rates
free( pfBF_logT );
for( j=0; j<iBF_TRvals; j++ )
	free( ppfBF[j] );
free( ppfBF );

// Free the memory allocated to the free-bound rates
free( pfFB_logradT );
free( pfFB_logeT );
for( j=0; j<iFB_TRvals; j++ )
{
	for( i=0; i<iFB_radTvals; i++ )
		free( pppfFB[j][i] );
	free( pppfFB[j] );
}
free( pppfFB );

// Free the memory allocated to the collisional rates
free( pfColl_logT );
for( j=0; j<iColl_exTRvals; j++ )
{
	free( ppfColl_ex_lu[j] );
	free( ppfColl_ex_ul[j] );
}
free( ppfColl_ex_lu );
free( ppfColl_ex_ul );
for( j=0; j<iColl_ionTRvals; j++ )
{
	free( ppfColl_ion[j] );
	free( ppfColl_rec[j] );
}
free( ppfColl_ion );
free( ppfColl_rec );
    
free( pH );
free( pterm2 );
free( pterm1 );
    
free( pMcZ_c_RIGHT );
free( pMcZ_c_LEFT );
free( pZ_c_RIGHT );
free( pZ_c_LEFT );
free( pTeZ_c );
free( pnu0 );
for( i=0; i<2; i++ )
	free( ppScaledTrt[i] );
free( ppScaledTrt );
free( pTrt );
}

#ifdef USE_LOWER_UPPER_DECOMPOSITION
// Solves a linear set of equations for a matrix equation A x = b using LU decomposition
// Based on Numerical Recipes in C, but adapted from the version in the RH code written by Han Uitenbroek
void CRadiativeRates::SolveLinearEq(double **A, int n, double *b, bool improve)
{
    int i, j;
    int *index;
    double d, **A_copy, *b_copy, *residual;
/*  
	If improve == TRUE improve the solution of the set of linear equations by evaluating 
	the residual and correcting the initial solution
	See: Press, Flannery, Teukolsky and Vetterling, Numerical Recipes: The Art of Scientific Computing 1986, p. 41
*/
	index = (int*)alloca(n * sizeof(int));
    
    // Copy matrix and source vector    
    if (improve)
    {
        residual = (double*)alloca(n * sizeof(double));
        b_copy = (double*)alloca(n * sizeof(double));
        A_copy = (double**)alloca(n * sizeof(double*));
        for( i = 0; i < n; i++ )
        {
            A_copy[i] = (double*)alloca(n * sizeof(double));
            b_copy[i] = b[i];
            for( j = 0; j < n; j++ )
                A_copy[i][j] = A[i][j];
        }
    }
    
    // Initial solution
    LUdcmp(A, n, index, &d);
    LUbksb(A, n, index, b);
    
    if (improve)
    {
        for( i = 0; i < n; i++ )
        {
            residual[i] = b_copy[i];
            for( j = 0; j < n; j++ )
                residual[i] -= A_copy[i][j] * b[j];
        }
        
        LUbksb(A, n, index, residual);
        
        // Correct the initial solution
        for( i = 0; i < n; i++ )
        {
            b[i] += residual[i];
        }
    }
}

// LU decomposition
// Decomposes an n x n matrix into lower(L) and upper(U) diagonal matrices (LU)
void CRadiativeRates::LUdcmp(double **A, int n, int *index, double *d)
{
    int i, j, k;
    int imax = 0;
    double big, dum, sum, temp, *vv;
    const double TINY = 1.0e-20;
    
    vv = (double*)alloca(n * sizeof(double));
    *d = 1.0;
    
    for( i = 0; i < n; i++ )
    {
        big = 0.0;
        for( j = 0; j < n; j++ )
            if((temp = fabs(A[i][j])) > big)
                big = temp;
        if( big == 0.0 )
        {
            // Throw an error here!
            fprintf(stderr, "Singular matrix in LUdcmp!\n");
            return;
        }
        vv[i] = 1.0 / big;
    }
    
    for( j = 0; j < n; j++ )
    {
        for( i = 0; i < j; i++ )
        {
            sum = A[i][j];
            for( k = 0; k < i; k++ )
                sum -= A[i][k] * A[k][j];
            A[i][j] = sum;
        }
        
        big = 0.0;
        for( i = j; i < n; i++ )
        {
            sum = A[i][j];
            for( k = 0; k < j; k++ )
                sum -= A[i][k] * A[k][j];
            A[i][j] = sum;
            if((dum = vv[i]*fabs(sum)) >= big)
            {
                big = dum;
                imax = i;
            }
        }
        
        if( j != imax )
        {
            for( k = 0; k < n; k++ )
            {
                dum = A[imax][k];
                A[imax][k] = A[j][k];
                A[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        
        index[j] = imax;
        if( A[j][j] == 0.0 )
            A[j][j] = TINY;
        if( j != n )
        {
            dum = 1.0 / A[j][j];
            for( i = j+1; i < n; i++ )
                A[i][j] *= dum;
        }
    }
}

// LU back substitution
// Uses the LU decomposition to then solve the matrix equation directly
void CRadiativeRates::LUbksb(double **A, int n, int *index, double *b)
{
    int i, j, ip;
    int ii = -1;
    double sum;
    
    for( i = 0; i < n; i++ )
    {
        ip = index[i];
        sum = b[ip];
        b[ip] = b[i];
        if( ii >= 0 )
        {
            for( j = ii; j < i; j++ )
                sum -= A[i][j] * b[j];
        }
        else if( sum )
        {
            ii = i;
        }
        b[i] = sum;
    }
    
    for( i = n-1; i >= 0; i-- )
    {
        sum = b[i];
        for( j = i+1; j < n; j++ )
            sum -= A[i][j] * b[j];
        b[i] = sum / A[i][i];
    }
}
#endif // USE_LOWER_UPPER_DECOMPOSITION

#ifdef USE_SINGLE_VALUE_DECOMPOSITION
// Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A = U W V T . The matrix U replaces a on output
// The diagonal matrix of singular values W is output as a vector w[1..n]. The matrix V (not the transpose V T ) is output as v[1..n][1..n]
int CRadiativeRates::svdcmp(double **a, int m, int n, double *w, double **v)
{
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
    
    if (m < n)
    {
        fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }
    
    rv1 = (double*)alloca((unsigned int) n*sizeof(double));
    
    // Householder reduction to bidiagonal form
    for (i = 0; i < n; i++)
    {
        // Left-hand reduction
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs((double)a[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    a[k][i] = (double)((double)a[k][i]/scale);
                    s += ((double)a[k][i] * (double)a[k][i]);
                }
                f = (double)a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (double)(f - g);
                if (i != n - 1)
                {
                    for (j = l; j < n; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += ((double)a[k][i] * (double)a[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            a[k][j] += (double)(f * (double)a[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    a[k][i] = (double)((double)a[k][i]*scale);
            }
        }
        w[i] = (float)(scale * g);
        
        // Right-hand reduction
        g = s = scale = 0.0;
        if (i < m && i != n - 1)
        {
            for (k = l; k < n; k++)
                scale += fabs((double)a[i][k]);
            if (scale)
            {
                for (k = l; k < n; k++)
                {
                    a[i][k] = (double)((double)a[i][k]/scale);
                    s += ((double)a[i][k] * (double)a[i][k]);
                }
                f = (double)a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = (double)(f - g);
                for (k = l; k < n; k++)
                    rv1[k] = (double)a[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < n; k++)
                            s += ((double)a[j][k] * (double)a[i][k]);
                        for (k = l; k < n; k++)
                            a[j][k] += (double)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++)
                    a[i][k] = (double)((double)a[i][k]*scale);
            }
        }
        anorm = fmax(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
    }
    
    // Accumulate the right-hand transformation
    for (i = n - 1; i >= 0; i--)
    {
        if (i < n - 1)
        {
            if (g)
            {
                for (j = l; j < n; j++)
                    v[j][i] = (double)(((double)a[i][j] / (double)a[i][l]) / g);
                // Double division to avoid underflow
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < n; k++)
                        s += ((double)a[i][k] * (double)v[k][j]);
                    for (k = l; k < n; k++)
                        v[k][j] += (double)(s * (double)v[k][i]);
                }
            }
            for (j = l; j < n; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    
    // Accumulate the left-hand transformation
    for (i = n - 1; i >= 0; i--)
    {
        l = i + 1;
        g = (double)w[i];
        if (i < n - 1)
            for (j = l; j < n; j++)
                a[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != n - 1)
            {
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += ((double)a[k][i] * (double)a[k][j]);
                    f = (s / (double)a[i][i]) * g;
                    for (k = i; k < m; k++)
                        a[k][j] += (double)(f * (double)a[k][i]);
                }
            }
            for (j = i; j < m; j++)
                a[j][i] = (double)((double)a[j][i]*g);
        }
        else
        {
            for (j = i; j < m; j++)
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }
    
    // Diagonalize the bidiagonal form
    for (k = n - 1; k >= 0; k--)
    {                             		// Loop over singular values
        for (its = 0; its < 30; its++)
        {                         		// Loop over allowed iterations
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     		// Test for splitting
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs((double)w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = (double)w[i];
                        h = pythag(f, g);
                        w[i] = (double)h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = (double)a[j][nm];
                            z = (double)a[j][i];
                            a[j][nm] = (double)(y * c + z * s);
                            a[j][i] = (double)(z * c - y * s);
                        }
                    }
                }
            }
            z = (double)w[k];
            if (l == k)
            {                  // Convergence
                if (z < 0.0)
                {              // Make singular value nonnegative
                    w[k] = (double)(-z);
                    for (j = 0; j < n; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                free((void*) rv1);
                fprintf(stderr, "No convergence after 30,000! iterations \n");
                return(0);
            }
            
            // Shift from bottom 2 x 2 minor
            x = (double)w[l];
            nm = k - 1;
            y = (double)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
            
            // Next QR transformation
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = (double)w[i];
                h = s * g;
                g = c * g;
                z = pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++)
                {
                    x = (double)v[jj][j];
                    z = (double)v[jj][i];
                    v[jj][j] = (double)(x * c + z * s);
                    v[jj][i] = (double)(z * c - x * s);
                }
                z = pythag(f, h);
                w[j] = (double)z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = (double)a[jj][j];
                    z = (double)a[jj][i];
                    a[jj][j] = (double)(y * c + z * s);
                    a[jj][i] = (double)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (double)x;
        }
    }
    return 0;
}

// Solves A x = B for a vector x, where A is specified by the arrays u[1..m][1..n], w[1..n], v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and will be equal for square matrices
// b[1..m] is the input right-hand side. x[1..n] is the output solution vector. No input quantities are destroyed, so the routine may be called sequentially with different b's
void CRadiativeRates::svbksb(double **u, double *w, double **v, int m, int n, double *b, double *x)
{
    int jj,j,i;
    double s,*tmp;

    tmp=(double*)alloca((unsigned int) n*sizeof(double));
    
    for( j=0; j<n; j++ )
    {
        s=0.0;
        if( w[j] )
        {
            for( i=0; i<m; i++)
            {
                s+= u[i][j]*b[i];
            }
                
            s /= w[j];
        }
        tmp[j] = s;
    }
    
    for( j=0; j<n; j++ )
    {
        s = 0.0;
        for( jj=0; jj<n; jj++ )
            s += v[j][jj]*tmp[jj];
        
        x[j] = s;
    }  
}

// Computes (a2 + b2)1/2 without destructive underflow or overflow
double CRadiativeRates::pythag(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;
    
    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}
#endif // USE_SINGLE_VALUE_DECOMPOSITION

void CRadiativeRates::GetBBRates( char *pszBBRatesFile )
{
FILE *pBBRatesFile;
int i, j;

// Open the bound-bound radiative transition rates file
pBBRatesFile = fopen( pszBBRatesFile, "r" );
	// Get the number of transitions stored in the file
	fscanf( pBBRatesFile, "%i", &iBB_TRvals );
	// Get the number of temperature values in the file
	fscanf( pBBRatesFile, "%i", &iBB_Tvals );

	// Allocate sufficient memory
	// Temperature range
	pfBB_logT = (double*)malloc( sizeof(double) * iBB_Tvals );
	// Photoexcitation
	ppfBB_lu = (double**)malloc( sizeof(double*) * iBB_TRvals );
	// Radiative decay
	ppfBB_ul = (double**)malloc( sizeof(double*) * iBB_TRvals );
	for( j=0; j<iBB_TRvals; j++ )
	{
		ppfBB_lu[j] = (double*)malloc( sizeof(double) * iBB_Tvals );
		ppfBB_ul[j] = (double*)malloc( sizeof(double) * iBB_Tvals );
	}

	// Read the rates
	for( i=0; i<iBB_Tvals; i++ )
	{
		ReadDouble( pBBRatesFile, &(pfBB_logT[i]) );
		for( j=0; j<iBB_TRvals; j++ )
			ReadDouble( pBBRatesFile, &(ppfBB_lu[j][i]) );
		for( j=0; j<iBB_TRvals; j++ )
			ReadDouble( pBBRatesFile, &(ppfBB_ul[j][i]) );	
	}
fclose( pBBRatesFile );
}

void CRadiativeRates::GetBFRates( char *pszBFRatesFile )
{
    FILE *pBFRatesFile;
    int i, j;

    // Open the bound-free radiative transition rates file
    pBFRatesFile = fopen( pszBFRatesFile, "r" );
	// Get the number of transitions stored in the file
	fscanf( pBFRatesFile, "%i", &iBF_TRvals );
	// Get the number of temperature values in the file
	fscanf( pBFRatesFile, "%i", &iBF_Tvals );

	// Allocate sufficient memory
	// Temperature range
	pfBF_logT = (double*)malloc( sizeof(double) * iBF_Tvals );
	// Photoionization
	ppfBF = (double**)malloc( sizeof(double*) * iBF_TRvals );
	for( j=0; j<iBF_TRvals; j++ )
		ppfBF[j] = (double*)malloc( sizeof(double) * iBF_Tvals );

	// Read the rates
	for( i=0; i<iBF_Tvals; i++ )
	{
		ReadDouble( pBFRatesFile, &(pfBF_logT[i]) );
		for( j=0; j<iBF_TRvals; j++ )
			ReadDouble( pBFRatesFile, &(ppfBF[j][i]) );
	}
fclose( pBFRatesFile );
/*
for( i=0; i<iBF_Tvals; i++ )
{
	printf( "logT=%.3g", pfBF_logT[i] );
	for( j=0; j<iBF_TRvals; j++ )
		printf( " BF[%i]=%.2e", j, ppfBF[j][i] );
	printf( "\n" );
}
*/
}

void CRadiativeRates::GetFBRates( char *pszFBRatesFile )
{
FILE *pFBRatesFile;
int h, i, j;

// Open the free-bound radiative transition rates file
pFBRatesFile = fopen( pszFBRatesFile, "r" );
	// Get the number of transitions stored in the file
	fscanf( pFBRatesFile, "%i", &iFB_TRvals );
	// Get the number of radiation temperature and electron temperature values in the file
	fscanf( pFBRatesFile, "%i", &iFB_radTvals );
	fscanf( pFBRatesFile, "%i", &iFB_eTvals );

	// Allocate sufficient memory
	// Temperature range
	pfFB_logradT = (double*)malloc( sizeof(double) * iFB_radTvals );
	pfFB_logeT = (double*)malloc( sizeof(double) * iFB_eTvals );
	// Radiative de-excitation
	pppfFB = (double***)malloc( sizeof(double**) * iFB_TRvals );
	for( j=0; j<iFB_TRvals; j++ )
	{
		pppfFB[j] = (double**)malloc( sizeof(double*) * iFB_radTvals );
		for( i=0; i<iFB_radTvals; i++ )
			pppfFB[j][i] = (double*)malloc( sizeof(double) * iFB_eTvals );
	}
	// j = transition index; i = radiation temperature index; h = electron temperature index
	// The correct way to index this array is: pppfFB[j][i][h]

	// Read the temperature ranges
	for( i=0; i<iFB_radTvals; i++ )
		ReadDouble( pFBRatesFile, &(pfFB_logradT[i]) );
	for( h=0; h<iFB_eTvals; h++ )
		ReadDouble( pFBRatesFile, &(pfFB_logeT[h]) );

	// Read the 2D rate arrays
	for( j=0; j<iFB_TRvals; j++ )
		for( h=0; h<iFB_eTvals; h++ )
			for( i=0; i<iFB_radTvals; i++ )
				ReadDouble( pFBRatesFile, &(pppfFB[j][i][h]) );	
fclose( pFBRatesFile );
}

void CRadiativeRates::GetCollRates( char *pszCollRatesFile )
{
FILE *pCollRatesFile;
int i, j;

    // Open the collisional transition rates file
    pCollRatesFile = fopen( pszCollRatesFile, "r" );
	// Get the number of collisional excitation and collisional ionization transitions stored in the file
	fscanf( pCollRatesFile, "%i", &iColl_exTRvals );
	fscanf( pCollRatesFile, "%i", &iColl_ionTRvals );
	// Get the number of temperature values in the file
	fscanf( pCollRatesFile, "%i", &iColl_Tvals );

	// Allocate sufficient memory
	// Temperature range
	pfColl_logT = (double*)malloc( sizeof(double) * iColl_Tvals );
	// Collisional excitation
	ppfColl_ex_lu = (double**)malloc( sizeof(double*) * iColl_exTRvals );
	ppfColl_ex_ul = (double**)malloc( sizeof(double*) * iColl_exTRvals );
	for( j=0; j<iColl_exTRvals; j++ )
	{
		ppfColl_ex_lu[j] = (double*)malloc( sizeof(double) * iColl_Tvals );
		ppfColl_ex_ul[j] = (double*)malloc( sizeof(double) * iColl_Tvals );
	}
	// Collisional ionization
	ppfColl_ion = (double**)malloc( sizeof(double*) * iColl_ionTRvals );
	ppfColl_rec = (double**)malloc( sizeof(double*) * iColl_ionTRvals );
	for( j=0; j<iColl_ionTRvals; j++ )
	{
		ppfColl_ion[j] = (double*)malloc( sizeof(double) * iColl_Tvals );
		ppfColl_rec[j] = (double*)malloc( sizeof(double) * iColl_Tvals );
	}

	// Read the rates
	// Collisional excitation
	for( i=0; i<iColl_Tvals; i++ )
	{
		ReadDouble( pCollRatesFile, &(pfColl_logT[i]) );
		for( j=0; j<iColl_exTRvals; j++ )
			ReadDouble( pCollRatesFile, &(ppfColl_ex_lu[j][i]) );
	}
	for( i=0; i<iColl_Tvals; i++ )
	{
		ReadDouble( pCollRatesFile, &(pfColl_logT[i]) );
		for( j=0; j<iColl_exTRvals; j++ )
			ReadDouble( pCollRatesFile, &(ppfColl_ex_ul[j][i]) );
	}
	// Collisional ionization	
	for( i=0; i<iColl_Tvals; i++ )
	{
		ReadDouble( pCollRatesFile, &(pfColl_logT[i]) );
		for( j=0; j<iColl_ionTRvals; j++ )
			ReadDouble( pCollRatesFile, &(ppfColl_ion[j][i]) );
	}
	for( i=0; i<iColl_Tvals; i++ )
	{
		ReadDouble( pCollRatesFile, &(pfColl_logT[i]) );
		for( j=0; j<iColl_ionTRvals; j++ )
			ReadDouble( pCollRatesFile, &(ppfColl_rec[j][i]) );
	}
fclose( pCollRatesFile );
}

void CRadiativeRates::GetBoundBoundRates( double *pfBB_lu, double *pfBB_ul, double *pflog10T )
{
double x[3], y[3];
int i, j;

for( j=0; j<iBB_TRvals; j++ )
{
	// The first task is the locate the temperature range within which to perform the interpolation of the rate for each transition
	// Catch out-of-range values
	if( pflog10T[j] < pfBB_logT[0] ) pflog10T[j] = pfBB_logT[0];
	else if( pflog10T[j] > pfBB_logT[iBB_Tvals-1] ) pflog10T[j] = pfBB_logT[iBB_Tvals-1];
	// Find the temperature range
	for( i=0; i<iBB_Tvals-1; i++ )
		// if( pflog10T[j] >= pfBB_logT[i] && pflog10T[j] <= pfBB_logT[i+1] )
		if( pflog10T[j] <= pfBB_logT[i+1] )
			break;

	x[1] = pfBB_logT[i];
	x[2] = pfBB_logT[i+1];

	// Photoexcitation
	y[1] = log10( ppfBB_lu[j][i] );
	y[2] = log10( ppfBB_lu[j][i+1] );
	LinearFit( x, y, pflog10T[j], &(pfBB_lu[j]) );
	pfBB_lu[j] = pow10( pfBB_lu[j] );

	// Radiative decay
	y[1] = log10( ppfBB_ul[j][i] );
	y[2] = log10( ppfBB_ul[j][i+1] );
	LinearFit( x, y, pflog10T[j], &(pfBB_ul[j]) );
	pfBB_ul[j] = pow10( pfBB_ul[j] );
}
}

void CRadiativeRates::GetBoundFreeRates( double *pfBF, double *pflog10T )
{
double x[3], y[3];
int i, j;

for( j=0; j<iBF_TRvals; j++ )
{
	// The first task is the locate the temperature range within which to perform the interpolation of the rate for each transition
	// Catch out-of-range values
	if( pflog10T[j] < pfBF_logT[0] ) pflog10T[j] = pfBF_logT[0];
	else if( pflog10T[j] > pfBF_logT[iBF_Tvals-1] ) pflog10T[j] = pfBF_logT[iBF_Tvals-1];
	// Find the temperature range
	for( i=0; i<iBF_Tvals-1; i++ )
		// if( pflog10T[j] >= pfBF_logT[i] && pflog10T[j] <= pfBF_logT[i+1] )
		if( pflog10T[j] <= pfBF_logT[i+1] )
			break;

	x[1] = pfBF_logT[i];
	x[2] = pfBF_logT[i+1];
	y[1] = log10( ppfBF[j][i] );
	y[2] = log10( ppfBF[j][i+1] );
	LinearFit( x, y, pflog10T[j], &(pfBF[j]) );
	pfBF[j] = pow10( pfBF[j] );
}
}

void CRadiativeRates::GetFreeBoundRates( double *pfFB, double *pflog10radT, double flog10eT, double fne )
{
double x1[3], x2[3], y1[3], y2[3];
int h, i, j;

// Set the minimum density limit to ensure well-behaved solutions at low densities (e.g. in the solar corona)
if( fne < MIN_DENSITY_LIMIT )
	fne = MIN_DENSITY_LIMIT;

// The first task is to locate the temperature range within which to perform the interpolation of the rate for each transition
// Catch out-of-range values
if( flog10eT < pfFB_logeT[0] ) flog10eT = pfFB_logeT[0];
else if( flog10eT > pfFB_logeT[iFB_eTvals-1] ) flog10eT = pfFB_logeT[iFB_eTvals-1];
// Find the temperature range
for( h=0; h<iFB_eTvals-1; h++ )
	// if( flog10eT >= pfFB_logeT[h] && flog10eT <= pfFB_logeT[h+1] )
	if( flog10eT <= pfFB_logeT[h+1] )
		break;

x2[1] = pfFB_logeT[h];
x2[2] = pfFB_logeT[h+1];

for( j=0; j<iFB_TRvals; j++ )
{
	// The first task is the locate the temperature range within which to perform the interpolation of the rate for each transition
	// Catch out-of-range values
	if( pflog10radT[j] < pfFB_logradT[0] ) pflog10radT[j] = pfFB_logradT[0];
	else if( pflog10radT[j] > pfFB_logradT[iFB_radTvals-1] ) pflog10radT[j] = pfFB_logradT[iFB_radTvals-1];
	// Find the temperature range
	for( i=0; i<iFB_radTvals-1; i++ )
		// if( pflog10radT[j] >= pfFB_logradT[i] && pflog10radT[j] <= pfFB_logradT[i+1] )
		if( pflog10radT[j] <= pfFB_logradT[i+1] )
			break;

	x1[1] = pfFB_logradT[i];
	x1[2] = pfFB_logradT[i+1];

	y1[1] = log10( pppfFB[j][i][h] );
	y1[2] = log10( pppfFB[j][i+1][h] );
	LinearFit( x1, y1, pflog10radT[j], &(y2[1]) );

	y1[1] = log10( pppfFB[j][i][h+1] );
	y1[2] = log10( pppfFB[j][i+1][h+1] );
	LinearFit( x1, y1, pflog10radT[j], &(y2[2]) );

	LinearFit( x2, y2, flog10eT, &(pfFB[j]) );
	pfFB[j] = fne * pow10( pfFB[j] );
}
}

void CRadiativeRates::GetFreeBoundRatesRH( double *pfFB, double *pflog10radT, double flog10eT, double fne, double *pfLevel_Ratio )
{
double x1[3], x2[3], y1[3], y2[3];
int h, i, j;
double fgij[] = {8.0,18.0,32.0,50.0};    // Ratio of statistical weights of levels i to j (j>i) for the transitions
double fdE[] = {39464.221,17536.454,9867.2552,6311.0119};   // = h * nu / k_B for the each transition (K)

// Set the minimum density limit to ensure well-behaved solutions at low densities (e.g. in the solar corona)
if( fne < MIN_DENSITY_LIMIT )
	fne = MIN_DENSITY_LIMIT;

// The first task is to locate the temperature range within which to perform the interpolation of the rate for each transition
// Catch out-of-range values
if( flog10eT < pfFB_logeT[0] ) flog10eT = pfFB_logeT[0];
else if( flog10eT > pfFB_logeT[iFB_eTvals-1] ) flog10eT = pfFB_logeT[iFB_eTvals-1];
// Find the temperature range
for( h=0; h<iFB_eTvals-1; h++ )
	// if( flog10eT >= pfFB_logeT[h] && flog10eT <= pfFB_logeT[h+1] )
	if( flog10eT <= pfFB_logeT[h+1] )
        	break;
    
x2[1] = pfFB_logeT[h];
x2[2] = pfFB_logeT[h+1];
    
for( j=0; j<iFB_TRvals; j++ )
{
	// The first task is the locate the temperature range within which to perform the interpolation of the rate for each transition
        // Catch out-of-range values
        if( pflog10radT[j] < pfFB_logradT[0] ) pflog10radT[j] = pfFB_logradT[0];
        else if( pflog10radT[j] > pfFB_logradT[iFB_radTvals-1] ) pflog10radT[j] = pfFB_logradT[iFB_radTvals-1];
        // Find the temperature range
        for( i=0; i<iFB_radTvals-1; i++ )
            // if( pflog10radT[j] >= pfFB_logradT[i] && pflog10radT[j] <= pfFB_logradT[i+1] )
	    if( pflog10radT[j] <= pfFB_logradT[i+1] )
                break;
        
        x1[1] = pfFB_logradT[i];
        x1[2] = pfFB_logradT[i+1];
        
        y1[1] = log10( pppfFB[j][i][h] );
        y1[2] = log10( pppfFB[j][i+1][h] );
        LinearFit( x1, y1, pflog10radT[j], &(y2[1]) );
        
        y1[1] = log10( pppfFB[j][i][h+1] );
        y1[2] = log10( pppfFB[j][i+1][h+1] );
        LinearFit( x1, y1, pflog10radT[j], &(y2[2]) );
        
        LinearFit( x2, y2, flog10eT, &(pfFB[j]) );
        pfFB[j] = fne * pow10( pfFB[j] );
        
        if( !(pfLevel_Ratio[j+11] == 0.0) )
        {
        	// 4.829e15 = 1 / 2.071e-16
        	// 2.071e-16 = 1/2 * (2 pi m_e k_b / h^2 )^(-1.5)
        	// Divide by the LTE population ratio
        	pfFB[j] *= (4.829e15 * exp( -fdE[j]/pow10(flog10eT) ) * pow( pow10(flog10eT), 1.5 ) ) / ( fgij[j] * fne ) ;
        
        	// Multiply by the NLTE ratio:
        	pfFB[j] /= pfLevel_Ratio[j+11];
        }
}
}

void CRadiativeRates::GetCollisionalRates( double *pfColl_ex_lu, double *pfColl_ex_ul, double *pfColl_ion, double *pfColl_rec, double flog10T, double fne )
{
double x[3], y[3];
int i, j;

// Set the minimum density limit to ensure well-behaved solutions at low densities (e.g. in the solar corona)
if( fne < MIN_DENSITY_LIMIT )
	fne = MIN_DENSITY_LIMIT;

// The first task is to locate the temperature range within which to perform the interpolation of the rate for each transition
// Catch out-of-range values
if( flog10T < pfColl_logT[0] ) flog10T = pfColl_logT[0];
else if( flog10T > pfColl_logT[iColl_Tvals-1] ) flog10T = pfColl_logT[iColl_Tvals-1];
// Find the temperature range
for( i=0; i<iColl_Tvals-1; i++ )
	// if( flog10T >= pfColl_logT[i] && flog10T <= pfColl_logT[i+1] )
	if( flog10T <= pfColl_logT[i+1] )
		break;

x[1] = pfColl_logT[i];
x[2] = pfColl_logT[i+1];
// Collisional excitation
for( j=0; j<iColl_exTRvals; j++ )
{
	y[1] = log10( ppfColl_ex_lu[j][i] );
	y[2] = log10( ppfColl_ex_lu[j][i+1] );
	LinearFit( x, y, flog10T, &(pfColl_ex_lu[j]) );
	pfColl_ex_lu[j] = fne * pow10( pfColl_ex_lu[j] );

	y[1] = log10( ppfColl_ex_ul[j][i] );
	y[2] = log10( ppfColl_ex_ul[j][i+1] );
	LinearFit( x, y, flog10T, &(pfColl_ex_ul[j]) );
	pfColl_ex_ul[j] = fne * pow10( pfColl_ex_ul[j] );   
}
// Collisional ionization
for( j=0; j<iColl_ionTRvals; j++ )
{
	y[1] = log10( ppfColl_ion[j][i] );
	y[2] = log10( ppfColl_ion[j][i+1] );
	LinearFit( x, y, flog10T, &(pfColl_ion[j]) );
	pfColl_ion[j] = fne * pow10( pfColl_ion[j] );

	y[1] = log10( ppfColl_rec[j][i] );
	y[2] = log10( ppfColl_rec[j][i+1] );
	LinearFit( x, y, flog10T, &(pfColl_rec[j]) );
	pfColl_rec[j] = fne * fne * pow10( pfColl_rec[j] );
}
}

void CRadiativeRates::GetCollisionalRatesRH( double *pfColl_ex_lu, double *pfColl_ex_ul, double *pfColl_ion, double *pfColl_rec, double flog10T, double fne )
{
double x[3], y[3];
double fgij[] = {0.25,0.11111111,0.0625,0.04,0.44444444,0.25,0.16,0.5625,0.36,0.64,2.0,8.0,18.0,32.0,50.0};    // Ratio of statistical weights of levels i to j (j>i) for the transitions
double fdE[] = {118301.48,140234.04,147864.85,151464.29,21918.168,29587.368,33148.41,
                7673.9985,11225.443,3551.444,157895.28,39464.221,17536.454,9867.2552,6311.0119};   // = h * nu / k_B for the each transition (K)
int i, j;
    
// Set the minimum density limit to ensure well-behaved solutions at low densities (e.g. in the solar corona)
if( fne < MIN_DENSITY_LIMIT )
	fne = MIN_DENSITY_LIMIT;

// The first task is the locate the temperature range within which to perform the interpolation of the rate for each transition
// Catch out-of-range values
if( flog10T < pfColl_logT[0] ) flog10T = pfColl_logT[0];
else if( flog10T > pfColl_logT[iColl_Tvals-1] ) flog10T = pfColl_logT[iColl_Tvals-1];
// Find the temperature range
for( i=0; i<iColl_Tvals-1; i++ )
	// if( flog10T >= pfColl_logT[i] && flog10T <= pfColl_logT[i+1] )
	if( flog10T <= pfColl_logT[i+1] )
        	break;
    
x[1] = pfColl_logT[i];
x[2] = pfColl_logT[i+1];

// Collisional excitation
for( j=0; j<iColl_exTRvals; j++ )
{
        y[1] = log10( ppfColl_ex_ul[j][i] );
        y[2] = log10( ppfColl_ex_ul[j][i+1] );
        LinearFit( x, y, flog10T, &(pfColl_ex_ul[j]) );
        
        pfColl_ex_ul[j] = fne * fgij[j] * sqrt( pow10( flog10T ) ) * pow10( pfColl_ex_ul[j] ) * 1.e6;
        pfColl_ex_lu[j] = (pfColl_ex_ul[j] * exp( -fdE[j] / pow10(flog10T) ) ) / fgij[j] ;
}
// Collisional ionization
for( j=0; j<iColl_ionTRvals; j++ )
{
	y[1] = log10( ppfColl_ion[j][i] );
        y[2] = log10( ppfColl_ion[j][i+1] );
        LinearFit( x, y, flog10T, &(pfColl_ion[j]) );
        
        pfColl_ion[j] = fne * sqrt( pow10( flog10T ) ) * pow10( pfColl_ion[j] ) * exp( -fdE[j+iColl_exTRvals] / pow10(flog10T) ) * 1.e6;
	// 2.071e-16 = 1/2 * (2 pi m_e k_b / h^2 )^(-1.5)
        pfColl_rec[j] = 2.071e-16 * pfColl_ion[j] * fne * fgij[j+iColl_exTRvals] * exp( fdE[j+iColl_exTRvals] / pow10(flog10T) ) * pow( pow10(flog10T), -1.5 );        
}
}

#ifdef USE_LOWER_UPPER_DECOMPOSITION
// Solves the 6x6 matrix equation A x = b, for the level populations of hydrogen x (where the vector x contains the level populations)
void CRadiativeRates::SolveHIIFraction( double *pfHstate, double *pfColl_ex_lu, double *pfColl_ex_ul, double *pfColl_ion, double *pfColl_rec, double *pfBB_lu, double *pfBB_ul, double *pfBF, double *pfFB )
{
    int i, j;
    int imax = 0;
    double Hmax = 0.0;
    double fSum;
    double **A, *b;

    A = (double**)alloca( sizeof(double*) * 6 );
    for( i=0; i<6; i++ )
        A[i] = (double*)alloca( sizeof(double) * 6 );

    b = (double*)alloca( sizeof(double) * 6 );
    
    // Row 1:
    fSum = 0.0;
    fSum += pfColl_ex_lu[0] + pfColl_ex_lu[1] + pfColl_ex_lu[2] + pfColl_ex_lu[3];
    fSum += pfColl_ion[0];
    A[0][0] = -fSum;
    
    A[0][1] = pfColl_ex_ul[0];
    A[0][2] = pfColl_ex_ul[1];
    A[0][3] = pfColl_ex_ul[2];
    A[0][4] = pfColl_ex_ul[3];
    A[0][5] = pfColl_rec[0];
    
    // Row 2:
    fSum = 0.0;
    fSum += pfColl_ex_ul[0];
    fSum += pfBB_lu[0] + pfBB_lu[1] + pfBB_lu[2];
    fSum += pfColl_ex_lu[4] + pfColl_ex_lu[5] + pfColl_ex_lu[6];
    fSum += pfBF[0];
    fSum += pfColl_ion[1];
    A[1][1] = -fSum;
    
    A[1][0] = pfColl_ex_lu[0];
    A[1][2] = pfBB_ul[0] + pfColl_ex_ul[4];
    A[1][3] = pfBB_ul[1] + pfColl_ex_ul[5];
    A[1][4] = pfBB_ul[2] + pfColl_ex_ul[6];
    A[1][5] = pfFB[0] + pfColl_rec[1];
    
    // Row 3:
    fSum = 0.0;
    fSum += pfColl_ex_ul[1] + pfColl_ex_ul[4];
    fSum += pfBB_ul[0];
    fSum += pfBB_lu[3] + pfBB_lu[4];
    fSum += pfColl_ex_lu[7] + pfColl_ex_lu[8];
    fSum += pfBF[1];
    fSum += pfColl_ion[2];
    A[2][2] = -fSum;
    
    A[2][0] = pfColl_ex_lu[1];
    A[2][1] = pfColl_ex_lu[4] + pfBB_lu[0];
    A[2][3] = pfColl_ex_ul[7] + pfBB_ul[3];
    A[2][4] = pfColl_ex_ul[8] + pfBB_ul[4];
    A[2][5] = pfFB[1] + pfColl_rec[2];
    
    // Row 4:
    fSum = 0.0;
    fSum += pfColl_ex_ul[2] + pfColl_ex_ul[5] + pfColl_ex_ul[7];
    fSum += pfBB_ul[1] + pfBB_ul[3];
    fSum += pfBB_lu[5];
    fSum += pfColl_ex_lu[9];
    fSum += pfBF[2];
    fSum += pfColl_ion[3];
    A[3][3] = -fSum;
    
    A[3][0] = pfColl_ex_lu[2];
    A[3][1] = pfColl_ex_lu[5] + pfBB_lu[1];
    A[3][2] = pfColl_ex_lu[7] + pfBB_lu[3];
    A[3][4] = pfColl_ex_ul[9] + pfBB_ul[5];
    A[3][5] = pfFB[2] + pfColl_rec[3];
    
    // Row 5:
    fSum = 0.0;
    fSum += pfColl_ex_ul[3] + pfColl_ex_ul[6] + pfColl_ex_ul[8] + pfColl_ex_ul[9];
    fSum += pfBB_ul[2] + pfBB_ul[4] + pfBB_ul[5];
    fSum += pfBF[3];
    fSum += pfColl_ion[4];
    A[4][4] = -fSum;
    
    A[4][0] = pfColl_ex_lu[3];
    A[4][1] = pfColl_ex_lu[6] + pfBB_lu[2];
    A[4][2] = pfColl_ex_lu[8] + pfBB_lu[4];
    A[4][3] = pfColl_ex_lu[9] + pfBB_lu[5];
    A[4][5] = pfFB[3] + pfColl_rec[4];
    
    // Row 6:
    fSum = 0.0;
    fSum += pfColl_rec[0] + pfColl_rec[1] + pfColl_rec[2] + pfColl_rec[3] + pfColl_rec[4];
    fSum += pfFB[0] + pfFB[1] + pfFB[2] + pfFB[3];
    A[5][5] = -fSum;
    
    A[5][0] = pfColl_ion[0];
    A[5][1] = pfColl_ion[1] + pfBF[0];
    A[5][2] = pfColl_ion[2] + pfBF[1];
    A[5][3] = pfColl_ion[3] + pfBF[2];
    A[5][4] = pfColl_ion[4] + pfBF[3];
    
    // Initialize the vector b
    for( i = 0; i < 6; i++ )
    {
        // Find the row with the maximum level population which we will eliminate from the matrix
        if( pfHstate[i] > Hmax )
        {
            Hmax = pfHstate[i];
            imax = i;
        }
        b[i] = 0.0;
    }
    b[imax] = 1.0;

    for( j = 0; j < 6; j++ )
        A[imax][j] = 1.0;

    SolveLinearEq(A, 6, b, true);
    
    for( i = 0; i < 6; i++ )
        pfHstate[i] = b[i];
    
    // Normalize the sum of all levels to 1
    fSum = 0.0;
    for( i=0; i<6; i++)
    {
		if( pfHstate[i] <= 0.0 ) pfHstate[i] = 1E-300;
        fSum += pfHstate[i];
    }
    for( i=0; i<6; i++)
        pfHstate[i] /= fSum;
}
#endif // USE_LOWER_UPPER_DECOMPOSITION

#ifdef USE_SINGLE_VALUE_DECOMPOSITION
// Solves the 7x6 matrix equation M x = b, for the level populations of hydrogen x (where the vector x contains the level populations)
void CRadiativeRates::SolveHIIFraction( double *pfHstate, double *pfColl_ex_lu, double *pfColl_ex_ul, double *pfColl_ion, double *pfColl_rec, double *pfBB_lu, double *pfBB_ul, double *pfBF, double *pfFB )
{
    int i;
    double fSum;
    double **M, *w, **V, *b;

    M = (double**)alloca( sizeof(double*) * 7 );
    for( i=0; i<7; i++ )
        M[i] = (double*)alloca( sizeof(double) * 6 );
    
    V = (double**)alloca( sizeof(double*) * 6 );
    for( i=0; i<6; i++ )
        V[i] = (double*)alloca( sizeof(double) * 6 );
    
    w = (double*)alloca( sizeof(double) * 6 );
    b = (double*)alloca( sizeof(double) * 7 );
    
    // Row 1:
    fSum = 0.0;
    fSum += pfColl_ex_lu[0] + pfColl_ex_lu[1] + pfColl_ex_lu[2] + pfColl_ex_lu[3];
    fSum += pfColl_ion[0];
    M[0][0] = -fSum;
    
    M[0][1] = pfColl_ex_ul[0];
    M[0][2] = pfColl_ex_ul[1];
    M[0][3] = pfColl_ex_ul[2];
    M[0][4] = pfColl_ex_ul[3];
    M[0][5] = pfColl_rec[0];
    
    // Row 2:
    fSum = 0.0;
    fSum += pfColl_ex_ul[0];
    fSum += pfBB_lu[0] + pfBB_lu[1] + pfBB_lu[2];
    fSum += pfColl_ex_lu[4] + pfColl_ex_lu[5] + pfColl_ex_lu[6];
    fSum += pfBF[0];
    fSum += pfColl_ion[1];
    M[1][1] = -fSum;
    
    M[1][0] = pfColl_ex_lu[0];
    M[1][2] = pfBB_ul[0] + pfColl_ex_ul[4];
    M[1][3] = pfBB_ul[1] + pfColl_ex_ul[5];
    M[1][4] = pfBB_ul[2] + pfColl_ex_ul[6];
    M[1][5] = pfFB[0] + pfColl_rec[1];
    
    // Row 3:
    fSum = 0.0;
    fSum += pfColl_ex_ul[1] + pfColl_ex_ul[4];
    fSum += pfBB_ul[0];
    fSum += pfBB_lu[3] + pfBB_lu[4];
    fSum += pfColl_ex_lu[7] + pfColl_ex_lu[8];
    fSum += pfBF[1];
    fSum += pfColl_ion[2];
    M[2][2] = -fSum;
    
    M[2][0] = pfColl_ex_lu[1];
    M[2][1] = pfColl_ex_lu[4] + pfBB_lu[0];
    M[2][3] = pfColl_ex_ul[7] + pfBB_ul[3];
    M[2][4] = pfColl_ex_ul[8] + pfBB_ul[4];
    M[2][5] = pfFB[1] + pfColl_rec[2];
    
    // Row 4:
    fSum = 0.0;
    fSum += pfColl_ex_ul[2] + pfColl_ex_ul[5] + pfColl_ex_ul[7];
    fSum += pfBB_ul[1] + pfBB_ul[3];
    fSum += pfBB_lu[5];
    fSum += pfColl_ex_lu[9];
    fSum += pfBF[2];
    fSum += pfColl_ion[3];
    M[3][3] = -fSum;
    
    M[3][0] = pfColl_ex_lu[2];
    M[3][1] = pfColl_ex_lu[5] + pfBB_lu[1];
    M[3][2] = pfColl_ex_lu[7] + pfBB_lu[3];
    M[3][4] = pfColl_ex_ul[9] + pfBB_ul[5];
    M[3][5] = pfFB[2] + pfColl_rec[3];
    
    // Row 5:
    fSum = 0.0;
    fSum += pfColl_ex_ul[3] + pfColl_ex_ul[6] + pfColl_ex_ul[8] + pfColl_ex_ul[9];
    fSum += pfBB_ul[2] + pfBB_ul[4] + pfBB_ul[5];
    fSum += pfBF[3];
    fSum += pfColl_ion[4];
    M[4][4] = -fSum;
    
    M[4][0] = pfColl_ex_lu[3];
    M[4][1] = pfColl_ex_lu[6] + pfBB_lu[2];
    M[4][2] = pfColl_ex_lu[8] + pfBB_lu[4];
    M[4][3] = pfColl_ex_lu[9] + pfBB_lu[5];
    M[4][5] = pfFB[3] + pfColl_rec[4];
    
    // Row 6:
    fSum = 0.0;
    fSum += pfColl_rec[0] + pfColl_rec[1] + pfColl_rec[2] + pfColl_rec[3] + pfColl_rec[4];
    fSum += pfFB[0] + pfFB[1] + pfFB[2] + pfFB[3];
    M[5][5] = -fSum;
    
    M[5][0] = pfColl_ion[0];
    M[5][1] = pfColl_ion[1] + pfBF[0];
    M[5][2] = pfColl_ion[2] + pfBF[1];
    M[5][3] = pfColl_ion[3] + pfBF[2];
    M[5][4] = pfColl_ion[4] + pfBF[3];
    
    // Row 7:
    for( i=0; i<6; i++ )
        M[6][i] = 1.0;

    // Singular value decomposition of the matrix
    i = svdcmp( M, 7, 6, w, V );
    
    // Initialize the vector b
    for( i = 0; i<6; i++)
        b[i] = 0.0;
    b[6] = 1.0;

    // Back substitution
    svbksb( M, w, V, 7, 6, b, pfHstate );
    
    // Normalize the sum of all levels to 1
    fSum = 0.0;
    for( i=0; i<6; i++)
    {
		if( pfHstate[i] <= 0.0 ) pfHstate[i] = 1E-300;
    		fSum += pfHstate[i];
    }
    for( i=0; i<6; i++)
        pfHstate[i] /= fSum;
}
#endif // USE_SINGLE_VALUE_DECOMPOSITION

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
void CRadiativeRates::GetAllDel_Hstate_dot_v( double *pHstate0, double *pHstate1, double *pHstate2, double *pHstate3, double *pHstate4, double *s, double *s_pos, double *pv, double *cross_section, double cell_volume, double *pDel_Hstate_dot_v )
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
void CRadiativeRates::GetAllDel_Hstate_dot_v( double *pHstate0, double *pHstate1, double *pHstate2, double *pHstate3, double *pHstate4, double *s, double *s_pos, double *pv, double delta_s, double *pDel_Hstate_dot_v )
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
{
// Variables used for interpolation
double x[3], y[3];

// Variables used for transport (flux) calculation
double Q1, Q2, Q3, QT, Hstate0, Hstate2;

int iLevel;

// for( iLevel=0; iLevel<6; iLevel++ )
iLevel = 5;
{
    // Calculate the fluxes to be used with Barton's Method

    if( pv[0] > 0.0 )
    {
        // Calculate the level population fraction
        x[1] = s[0];
		x[2] = s[1];
		y[1] = pHstate0[iLevel];
		y[2] = pHstate1[iLevel];
		LinearFit( x, y, s_pos[0], &Q1 );

		x[1] = s[1];
		x[2] = s[2];
		y[1] = pHstate1[iLevel];
		y[2] = pHstate2[iLevel];
		LinearFit( x, y, s_pos[0], &Q2 );

		Q3 = pHstate1[iLevel];

		if( pHstate2[iLevel] <= pHstate1[iLevel] )
		{
            QT = max( Q1, Q2 );
            if( Q3 < QT )
                Hstate0 = Q3;
            else
                Hstate0 = QT;
		}
		else
		{
            QT = min( Q1, Q2 );
            if( Q3 > QT )
                Hstate0 = Q3;
            else
                Hstate0 = QT;
		}
    }
    else
    {
        // Calculate the level population fraction
        x[1] = s[2];
        x[2] = s[3];
        y[1] = pHstate2[iLevel];
        y[2] = pHstate3[iLevel];
        LinearFit( x, y, s_pos[0], &Q1 );

        x[1] = s[1];
        x[2] = s[2];
        y[1] = pHstate1[iLevel];
        y[2] = pHstate2[iLevel];
        LinearFit( x, y, s_pos[0], &Q2 );

        Q3 = pHstate2[iLevel];

        if( pHstate2[iLevel] <= pHstate1[iLevel] )
        {
            QT = min( Q1, Q2 );
            if( Q3 > QT )
                Hstate0 = Q3;
            else
                Hstate0 = QT;
        }
        else
        {
            QT = max( Q1, Q2 );
            if( Q3 < QT )
                Hstate0 = Q3;
            else
                Hstate0 = QT;
        }
    }

    if( pv[2] > 0.0 )
    {
        // Calculate the level population fraction
		x[1] = s[1];
		x[2] = s[2];
		y[1] = pHstate1[iLevel];
		y[2] = pHstate2[iLevel];
		LinearFit( x, y, s_pos[2], &Q1 );

		x[1] = s[2];
		x[2] = s[3];
		y[1] = pHstate2[iLevel];
		y[2] = pHstate3[iLevel];
		LinearFit( x, y, s_pos[2], &Q2 );

		Q3 = pHstate2[iLevel];

		if( pHstate3[iLevel] <= pHstate2[iLevel] )
		{
            QT = max( Q1, Q2 );
            if( Q3 < QT )
                Hstate2 = Q3;
            else
                Hstate2 = QT;
		}
		else
		{
            QT = min( Q1, Q2 );
            if( Q3 > QT )
                Hstate2 = Q3;
            else
                Hstate2 = QT;
		}
    }
    else
    {
        // Calculate the level population fraction
		x[1] = s[3];
        x[2] = s[4];
        y[1] = pHstate3[iLevel];
        y[2] = pHstate4[iLevel];
        LinearFit( x, y, s_pos[2], &Q1 );

        x[1] = s[2];
        x[2] = s[3];
        y[1] = pHstate2[iLevel];
        y[2] = pHstate3[iLevel];
        LinearFit( x, y, s_pos[2], &Q2 );

        Q3 = pHstate3[iLevel];

        if( pHstate3[iLevel] <= pHstate2[iLevel] )
        {
            QT = min( Q1, Q2 );
            if( Q3 > QT )
                Hstate2 = Q3;
            else
                Hstate2 = QT;
        }
        else
        {
            QT = max( Q1, Q2 );
            if( Q3 < QT )
                Hstate2 = Q3;
            else
                Hstate2 = QT;
        }
    }
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
    pDel_Hstate_dot_v[iLevel] = - ( ( Hstate2 * pv[2] * cross_section[2] ) - ( Hstate0 * pv[0] * cross_section[0] ) ) / cell_volume;
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
    pDel_Hstate_dot_v[iLevel] = - ( ( Hstate2 * pv[2] ) - ( Hstate0 * pv[0] ) ) / delta_s;
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
}
}

void CRadiativeRates::Normalize( double *pHstate )
{
double fSum = 0.0;
int i;

for( i=0; i<6; i++ )
	fSum += pHstate[i];

for( i=0; i<6; i++ )
	pHstate[i] /= fSum;
}

int CRadiativeRates::GetNumberOfTransitions( void )
{
    return (iNBBT+iNBFT);
}

double CRadiativeRates::GetZ_c_LEFT( int i )
{
    return (pZ_c_LEFT[i]);
}

double CRadiativeRates::GetZ_c_RIGHT( int i )
{
    return (pZ_c_RIGHT[i]);
}

double CRadiativeRates::GetMcZ_c_LEFT( int i )
{
    return (pMcZ_c_LEFT[i]);
}

double CRadiativeRates::GetMcZ_c_RIGHT( int i )
{
    return (pMcZ_c_RIGHT[i]);
}


double CRadiativeRates::GetTeZ_c( int i )
{
    return (pTeZ_c[i]);
}

double CRadiativeRates::GetNu0( int i )
{
    return (pnu0[i]);
}

double CRadiativeRates::GetH( int i )
{
    return (pH[i]);
}

double CRadiativeRates::GetTrt( int i )
{
	return (pTrt[i]);
}

double CRadiativeRates::GetScaledTrt( int i, bool iFlag )
{
	return(ppScaledTrt[iFlag][i]);
}

void CRadiativeRates::SetZ_c_LEFT( int i, double Z_c )
{
    pZ_c_LEFT[i] = Z_c;
}

void CRadiativeRates::SetZ_c_RIGHT( int i, double Z_c )
{
    pZ_c_RIGHT[i] = Z_c;
}

void CRadiativeRates::SetMcZ_c_LEFT( int i, double McZ_c )
{
    pMcZ_c_LEFT[i] = McZ_c;
}

void CRadiativeRates::SetMcZ_c_RIGHT( int i, double McZ_c )
{
    pMcZ_c_RIGHT[i] = McZ_c;
}

void CRadiativeRates::SetScaledTrt( int i, double Trt, bool iFlag )
{
	ppScaledTrt[iFlag][i] = Trt;
}

#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
