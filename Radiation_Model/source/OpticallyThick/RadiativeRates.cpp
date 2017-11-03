// ****
// *
// * Function bodies for the class definition of the radiative rates
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 11/02/2017
// *
// ****

// **** RADIATIVE TRANSITION RATES CLASS ****

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "RadiativeRates.h"
#include "../../../Resources/source/fitpoly.h"
#include "../../../Resources/source/file.h"
#include "../../../Resources/source/svd.h"

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
free( pTrt );

}

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
		if( pflog10T[j] >= pfBB_logT[i] && pflog10T[j] <= pfBB_logT[i+1] )
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
		if( pflog10T[j] >= pfBF_logT[i] && pflog10T[j] <= pfBF_logT[i+1] )
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

// The first task is the locate the temperature range within which to perform the interpolation of the rate for each transition
// Catch out-of-range values
if( flog10eT < pfFB_logeT[0] ) flog10eT = pfFB_logeT[0];
else if( flog10eT > pfFB_logeT[iFB_eTvals-1] ) flog10eT = pfFB_logeT[iFB_eTvals-1];
// Find the temperature range
for( h=0; h<iFB_eTvals-1; h++ )
	if( flog10eT >= pfFB_logeT[h] && flog10eT <= pfFB_logeT[h+1] )
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
		if( pflog10radT[j] >= pfFB_logradT[i] && pflog10radT[j] <= pfFB_logradT[i+1] )
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

    
    // The first task is the locate the temperature range within which to perform the interpolation of the rate for each transition
    // Catch out-of-range values
    if( flog10eT < pfFB_logeT[0] ) flog10eT = pfFB_logeT[0];
    else if( flog10eT > pfFB_logeT[iFB_eTvals-1] ) flog10eT = pfFB_logeT[iFB_eTvals-1];
    // Find the temperature range
    for( h=0; h<iFB_eTvals-1; h++ )
        if( flog10eT >= pfFB_logeT[h] && flog10eT <= pfFB_logeT[h+1] )
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
            if( pflog10radT[j] >= pfFB_logradT[i] && pflog10radT[j] <= pfFB_logradT[i+1] )
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

// The first task is the locate the temperature range within which to perform the interpolation of the rate for each transition
// Catch out-of-range values
if( flog10T < pfColl_logT[0] ) flog10T = pfColl_logT[0];
else if( flog10T > pfColl_logT[iColl_Tvals-1] ) flog10T = pfColl_logT[iColl_Tvals-1];
// Find the temperature range
for( i=0; i<iColl_Tvals-1; i++ )
	if( flog10T >= pfColl_logT[i] && flog10T <= pfColl_logT[i+1] )
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
    
//    printf("Original j %i lu %.2e ul %.2e\n",j,pfColl_ex_ul[j],pfColl_ex_lu[j]);
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
//    printf("Original j %i ion %.2e rec %.2e\n",j,pfColl_ion[j],pfColl_rec[j]);

}
}


void CRadiativeRates::GetCollisionalRatesRH( double *pfColl_ex_lu, double *pfColl_ex_ul, double *pfColl_ion, double *pfColl_rec, double flog10T, double fne )
{
    double x[3], y[3];
    double fgij[] = {0.25,0.11111111,0.0625,0.04,0.44444444,0.25,0.16,0.5625,0.36,0.64,2.0,8.0,18.0,32.0,50.0};    // Ratio of statistical weights of levels i to j (j>i) for the transitions
    double fdE[] = {118301.48,140234.04,147864.85,151464.29,21918.168,29587.368,33148.41,
                    7673.9985,11225.443,3551.444,157895.28,39464.221,17536.454,9867.2552,6311.0119};   // = h * nu / k_B for the each transition (K)
    int i, j;
    
    // The first task is the locate the temperature range within which to perform the interpolation of the rate for each transition
    // Catch out-of-range values
    if( flog10T < pfColl_logT[0] ) flog10T = pfColl_logT[0];
    else if( flog10T > pfColl_logT[iColl_Tvals-1] ) flog10T = pfColl_logT[iColl_Tvals-1];
    // Find the temperature range
    for( i=0; i<iColl_Tvals-1; i++ )
        if( flog10T >= pfColl_logT[i] && flog10T <= pfColl_logT[i+1] )
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
        
        pfColl_rec[j] = 2.071e-16 * pfColl_ion[j] * fne * fgij[j+iColl_exTRvals] * exp( fdE[j+iColl_exTRvals] / pow10(flog10T) ) * pow( pow10(flog10T), -1.5 );
        // 2.071e-16 = 1/2 * (2 pi m_e k_b / h^2 )^(-1.5)
        
    }
    
}


void CRadiativeRates::SolveHIIFraction( double *pfHstate, double *pfColl_ex_lu, double *pfColl_ex_ul, double *pfColl_ion, double *pfColl_rec, double *pfBB_lu, double *pfBB_ul, double *pfBF, double *pfFB )
{
    // Solves the 7x6 matrix equation M x = b, for the level populations of hydrogen x
    // Returns x[5], the ionized fraction of hydrogen
    
    int i;
    double fSum;
    double **M, *w, **V, *b, *x;
    
    M = (double**)malloc( sizeof(double*) * (7) );
    for( i=0; i<7; i++ )
        M[i] = (double*)malloc( sizeof(double) * 6 );
    
    V = (double**)malloc( sizeof(double*) * (6) );
    for( i=0; i<6; i++ )
        V[i] = (double*)malloc( sizeof(double) * 6 );
    
    w = (double*)malloc( sizeof(double) * 6 );
    b = (double*)malloc( sizeof(double) * 7 );
    x = (double*)malloc( sizeof(double) * 6 );

    
    // row 1:
    fSum = 0.0;
    fSum += pfColl_ex_lu[0] + pfColl_ex_lu[1] + pfColl_ex_lu[2] + pfColl_ex_lu[3];
    fSum += pfColl_ion[0];
    M[0][0] = -fSum;
    
    M[0][1] = pfColl_ex_ul[0];
    M[0][2] = pfColl_ex_ul[1];
    M[0][3] = pfColl_ex_ul[2];
    M[0][4] = pfColl_ex_ul[3];
    M[0][5] = pfColl_rec[0];
    
    // row 2:
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
    
    // row 3:
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
    
    // row 4:
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
    
    // row 5:
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
    
    // row 6:
    fSum = 0.0;
    fSum += pfColl_rec[0] + pfColl_rec[1] + pfColl_rec[2] + pfColl_rec[3] + pfColl_rec[4];
    fSum += pfFB[0] + pfFB[1] + pfFB[2] + pfFB[3];
    M[5][5] = -fSum;
    
    M[5][0] = pfColl_ion[0];
    M[5][1] = pfColl_ion[1] + pfBF[0];
    M[5][2] = pfColl_ion[2] + pfBF[1];
    M[5][3] = pfColl_ion[3] + pfBF[2];
    M[5][4] = pfColl_ion[4] + pfBF[3];
    
    // row 7:
    for( i=0; i<6; i++ )
        M[6][i] = 1.0;
    
    /*    printf("M\n");
     printf("%.3e %.3e %.3e %.3e %.3e %.3e\n", M[0][0], M[0][1], M[0][2], M[0][3], M[0][4], M[0][5]);
     printf("%.3e %.3e %.3e %.3e %.3e %.3e\n", M[1][0], M[1][1], M[1][2], M[1][3], M[1][4], M[1][5]);
     printf("%.3e %.3e %.3e %.3e %.3e %.3e\n", M[2][0], M[2][1], M[2][2], M[2][3], M[2][4], M[2][5]);
     printf("%.3e %.3e %.3e %.3e %.3e %.3e\n", M[3][0], M[3][1], M[3][2], M[3][3], M[3][4], M[3][5]);
     printf("%.3e %.3e %.3e %.3e %.3e %.3e\n", M[4][0], M[4][1], M[4][2], M[4][3], M[4][4], M[4][5]);
     printf("%.3e %.3e %.3e %.3e %.3e %.3e\n", M[5][0], M[5][1], M[5][2], M[5][3], M[5][4], M[5][5]);
     printf("%.3e %.3e %.3e %.3e %.3e %.3e\n", M[6][0], M[6][1], M[6][2], M[6][3], M[6][4], M[6][5]);
     */

    // Singular value decomposition of the matrix:
    i = svdcmp( M, 7, 6, w, V );
    
    // Initialize the vector b
    for( i = 0; i<6; i++)
        b[i] = 0.0;
    b[6] = 1.0;
    
    
    //   NOT CLEAR IF THE FOLLOWING IS NEEDED:
    //     (from the usage in Numerical Recipes book)
    //     (gave incorrect answer when wmin was set to wmax/1e6)
    /*
     double wmax = 0.0;
     double wmin;
     for( i=0; i<7; i++ )
     if( w[i] > wmax ) wmax = w[i];
     wmin = wmax * 1.0e-20;
     for( i=0; i<7; i++ )
     if( w[i] < wmin ) w[i] = 0.0;
     */
    
    // Back substitution:
    svbksb( M, w, V, 7, 6, b, x );
    
    // Normalize the sum of all levels to 1
    fSum = 0.0;
    for( i=0; i<6; i++)
    {
	if( x[i] <= 0.0 ) x[i] = 1E-300;
        fSum += x[i];
    }
    for( i=0; i<6; i++)
        x[i] /= fSum;
    
//    printf("x[0] %.2e x[1] %.2e x[2] %.2e x[3] %.2e x[4] %.2e x[5] %.2e\n", x[0],x[1],x[2],x[3],x[4],x[5]);
    
    memcpy( pfHstate, x, sizeof(double)*6 );
    
    // Free memory
    for( i=0; i<7; i++ )
        free( M[i] );
    free( M );
    for( i=0; i<6; i++ )
        free( V[i] );
    free( V );
    free( w );
    free( x );
    free( b );
}

void CRadiativeRates::SolveHIIFraction( double *pfHstate, double *pfGradients, double *pfColl_ex_lu, double *pfColl_ex_ul, double *pfColl_ion, double *pfColl_rec, double *pfBB_lu, double *pfBB_ul, double *pfBF, double *pfFB )
{
    // Solves the 7x6 matrix equation M x = b, for the level populations of hydrogen x
    // Returns x[5], the ionized fraction of hydrogen
    
    int i;
    double fSum;
    double **M, *w, **V, *b, *x;
    
    M = (double**)malloc( sizeof(double*) * (7) );
    for( i=0; i<7; i++ )
        M[i] = (double*)malloc( sizeof(double) * 6 );
    
    V = (double**)malloc( sizeof(double*) * (6) );
    for( i=0; i<6; i++ )
        V[i] = (double*)malloc( sizeof(double) * 6 );
    
    w = (double*)malloc( sizeof(double) * 6 );
    b = (double*)malloc( sizeof(double) * 7 );
    x = (double*)malloc( sizeof(double) * 6 );
    
    
    // row 1:
    fSum = 0.0;
    fSum += pfColl_ex_lu[0] + pfColl_ex_lu[1] + pfColl_ex_lu[2] + pfColl_ex_lu[3];
    fSum += pfColl_ion[0];
    M[0][0] = -fSum;
    
    M[0][1] = pfColl_ex_ul[0];
    M[0][2] = pfColl_ex_ul[1];
    M[0][3] = pfColl_ex_ul[2];
    M[0][4] = pfColl_ex_ul[3];
    M[0][5] = pfColl_rec[0];
    
    // row 2:
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
    
    // row 3:
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
    
    // row 4:
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
    
    // row 5:
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
    
    // row 6:
    fSum = 0.0;
    fSum += pfColl_rec[0] + pfColl_rec[1] + pfColl_rec[2] + pfColl_rec[3] + pfColl_rec[4];
    fSum += pfFB[0] + pfFB[1] + pfFB[2] + pfFB[3];
    M[5][5] = -fSum;
    
    M[5][0] = pfColl_ion[0];
    M[5][1] = pfColl_ion[1] + pfBF[0];
    M[5][2] = pfColl_ion[2] + pfBF[1];
    M[5][3] = pfColl_ion[3] + pfBF[2];
    M[5][4] = pfColl_ion[4] + pfBF[3];
    
    // row 7:
    for( i=0; i<6; i++ )
        M[6][i] = 1.0;
    
    // Singular value decomposition of the matrix:
    i = svdcmp( M, 7, 6, w, V );
    
    // Initialize the vector b
    for( i = 0; i<6; i++)
        b[i] = pfGradients[i];
    b[6] = 1.;
    
    
    //   NOT CLEAR IF THE FOLLOWING IS NEEDED:
    //     (from the usage in Numerical Recipes book)
    //     (gave incorrect answer when wmin was set to wmax/1e6)
    /*
     double wmax = 0.0;
     double wmin;
     for( i=0; i<7; i++ )
     if( w[i] > wmax ) wmax = w[i];
     wmin = wmax * 1.0e-20;
     for( i=0; i<7; i++ )
     if( w[i] < wmin ) w[i] = 0.0;
     */
    
    // Back substitution:
    svbksb( M, w, V, 7, 6, b, x );
    
    // Normalize the sum of all levels to 1
    fSum = 0.0;
    for( i=0; i<6; i++)
        fSum += x[i];
    for( i=0; i<6; i++)
        x[i] /= fSum;
    
    memcpy( pfHstate, x, sizeof(double)*6 );
    
    // Free memory
    for( i=0; i<7; i++ )
        free( M[i] );
    free( M );
    for( i=0; i<6; i++ )
        free( V[i] );
    free( V );
    free( w );
    free( x );
    free( b );
}

void CRadiativeRates::GetAllDel_Hstate_dot_v( double *pHstate0, double *pHstate1, double *pHstate2, double *pHstate3, double *pHstate4, double *s, double *s_pos, double *pv, double delta_s, double *pDel_Hstate_dot_v )
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

    // pDel_Hstate_dot_v[iLevel] = - ( ( Hstate2 * pv[2] ) - ( Hstate0 * pv[0] ) ) / delta_s;
    pDel_Hstate_dot_v[iLevel] = ( ( Hstate0 * pv[0] ) - ( Hstate2 * pv[2] ) ) / delta_s;
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

double CRadiativeRates::GetTerm1( int i )
{
    return (pterm1[i]);
}

double CRadiativeRates::GetTerm2( int i )
{
    return (pterm2[i]);
}

double CRadiativeRates::GetNu0( int i )
{
    return (pnu0[i]);
}

double CRadiativeRates::GetH( int i )
{
    return (pH[i]);
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
