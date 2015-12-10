// ****
// *
// * Forward Model class function bodies
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 11/03/2015
// *
// ****


#include <stdio.h>
#include <malloc.h>

#include "forward.h"
#include "../../Resources/source/constants.h"
#include "../../Resources/source/file.h"


CForward::CForward( char *szFilename )
{
Initialise( szFilename );
}

CForward::~CForward( void )
{
FreeAll();
}

void CForward::Initialise( char *szFilename )
{
FILE *pFile;
char szBuffer[256];
int i;

// Open the forward model configuration file
pFile = fopen( szFilename, "r" );

// Get the number of virtual instruments to use
fscanf( pFile, "%i", &iNumInstruments );

// Allocate sufficient memory for the virtual instruments
ppInstrument = (PPINSTRUMENT)malloc( sizeof( CInstrument ) * iNumInstruments );

for( i=0; i<iNumInstruments; i++ )
{
    // Get the instrument name and create the virtual instrument object
    fscanf( pFile, "%s", szBuffer );
    ppInstrument[i] = new CInstrument( szBuffer );
}

fclose( pFile );
}

void CForward::FreeAll( void )
{
int i;

for( i=0; i<iNumInstruments; i++ )
    delete ppInstrument[i];
	
free( ppInstrument );
}

int CForward::GetNumInstruments( void )
{
return iNumInstruments;
}

void CForward::WriteIonLineListToFile( char *szFilename )
{
FILE *pFile;
int i;

pFile = fopen( szFilename, "w" );

for( i=0; i<iNumInstruments; i++ )
    ppInstrument[i]->WriteIonLineListToFile( pFile );

fclose( pFile );
}

void CForward::CreateVirtualDetector( PLOOP pLoop )
{
int i;

for( i=0; i<iNumInstruments; i++ )
    ppInstrument[i]->CreateVirtualDetector( pLoop );
}

void CForward::ResetVirtualDetector( void )
{
int i;

for( i=0; i<iNumInstruments; i++ )
    ppInstrument[i]->ResetVirtualDetector();
}

void CForward::WriteDetectorToFile( char *pszWorkingDirectory )
{
int i;

for( i=0; i<iNumInstruments; i++ )
    ppInstrument[i]->WriteDetectorToFile( pszWorkingDirectory );
}

void CForward::WriteDetectorToFile( char *pszWorkingDirectory, int iFileNumber )
{
int i;

for( i=0; i<iNumInstruments; i++ )
    ppInstrument[i]->WriteDetectorToFile( pszWorkingDirectory, iFileNumber );
}

void CForward::ForwardModel( PLOOP pLoop )
{
int i;

for( i=0; i<iNumInstruments; i++ )
	ppInstrument[i]->ForwardModel( pLoop );
}

void CForward::EM_Loci( double flog10_Tmin, double flog10_Tmax, double flog10_step, double flog10_n, char *pszWorkingDirectory )
{
FILE *pOutputFile, *pPottaschOutputFile = NULL, **ppInputFile;
char szOutputFilename[256], szInstrumentName[256], szInputFilename[256];
double flog10_T, fBuffer, fCoord;
double fEQEM_Loci, fNEQEM_Loci;
double *pflog10_TGmin, *pfEQEM_Loci_min, *pfNEQEM_Loci_min = NULL;
int iNumVals, iBuffer;
int i, j;
bool bNEQ, bSort;

for( i=0; i<iNumInstruments; i++ )
    ppInstrument[i]->EM_Loci( flog10_Tmin, flog10_Tmax, flog10_step, flog10_n, pszWorkingDirectory );

// Consolidate the emission measure loci plots for each individual instrument into a single file
// The file format will be:
// Number of instruments (NI)
// Number of temperature values
// Temperature values [ROW]
// Spatial coordinate
// Equilibrium (true) emission measure for instrument 1 [ROW]
// Non-equilibrium emission measure for instrument 1 [ROW]
// .
// .
// Equilibrium (true) emission measure for instrument NI [ROW]
// Non-equilibrium emission measure for instrument NI [ROW]
// Spatial coordinate
// .
// .

// Check whether non-equilibrium emission measure data exists
bNEQ = false;
for( i=0; i<iNumInstruments; i++ )
{
    bNEQ = ppInstrument[i]->IsNonEquilibriumDetector();
    if( bNEQ ) break;
}

pflog10_TGmin = (double*)malloc( sizeof(double) * iNumInstruments );
pfEQEM_Loci_min = (double*)malloc( sizeof(double) * iNumInstruments );
if( bNEQ )
{
    pfNEQEM_Loci_min = (double*)malloc( sizeof(double) * iNumInstruments );
}

// Open the data file for output
sprintf( szOutputFilename, "%s/EM_Loci.dat", pszWorkingDirectory );
pOutputFile = fopen( szOutputFilename, "w" );
sprintf( szOutputFilename, "%s/EM_Loci_Pottasch.dat", pszWorkingDirectory );
pPottaschOutputFile = fopen( szOutputFilename, "w" );

// Write the number of instruments into the output files
fprintf( pOutputFile, "%i\n", iNumInstruments );
fprintf( pPottaschOutputFile, "%i\n", iNumInstruments );

// Write the temperature range
iNumVals = (int)( ( ( flog10_Tmax - flog10_Tmin ) / flog10_step ) + 1.0 );
fprintf( pOutputFile, "%i\n", iNumVals );
for( i=0; i<iNumVals; i++ )
{
    flog10_T = flog10_Tmin+(((double)i)*flog10_step);
    fprintf( pOutputFile, "\t%g", flog10_T );
}
fprintf( pOutputFile, "\n" );

// Allocate sufficient space for the pointers to the emission measure data files
ppInputFile = (FILE**)malloc( sizeof(FILE) * iNumInstruments );

// Open the emission measure data files for each instrument
for( i=0; i<iNumInstruments; i++ )
{
    // Get the name of the instrument
    ppInstrument[i]->GetInstrumentName( szInstrumentName );
    // Construct the filename
    sprintf( szInputFilename, "%s/%s.EM_Loci", pszWorkingDirectory, szInstrumentName );
    ppInputFile[i] = fopen( szInputFilename, "r" );

    // Skip past the temperature data
    fscanf( ppInputFile[i], "%i", &iBuffer );
    for( j=0; j<iNumVals; j++ )
        ReadDouble( ppInputFile[i], &fBuffer );
}

// For each spatial coordinate, write the emission measures
// into the output file following the above format
for ( ;; )
{
    // Get the spatial coordinate
    ReadDouble( ppInputFile[0], &fCoord );
    if( feof( ppInputFile[0]) ) break;
    for( i=1; i<iNumInstruments; i++ )
        ReadDouble( ppInputFile[i], &fCoord );

    // Write the spatial coordinate into the output files
    fprintf( pOutputFile, "%.8e\n", fCoord );
    fprintf( pPottaschOutputFile, "%.8e\n", fCoord );

    for( i=0; i<iNumInstruments; i++ )
    {
        // Get the equilibrium emission measures and write them into the output file
        pfEQEM_Loci_min[i] = LARGEST_DOUBLE;
        pflog10_TGmin[i] = 0.0;
        for( j=0; j<iNumVals; j++ )
        {
            flog10_T = flog10_Tmin+(((double)j)*flog10_step);
            ReadDouble( ppInputFile[i], &fEQEM_Loci );
            fprintf( pOutputFile, "\t%.8e", fEQEM_Loci );
            if( fEQEM_Loci && fEQEM_Loci < pfEQEM_Loci_min[i] )
            {
                pfEQEM_Loci_min[i] = fEQEM_Loci;
                pflog10_TGmin[i] = flog10_T;
            }
        }
        fprintf( pOutputFile, "\n" );
        if( pfEQEM_Loci_min[i] < LARGEST_DOUBLE )
        {
            pfEQEM_Loci_min[i] /= 0.7;
        }
        else
        {
            pfEQEM_Loci_min[i] = 0.0;
        }

        if( bNEQ )
        {
            // Get the non-equilibrium emission measures and write them into the output file
            pfNEQEM_Loci_min[i] = LARGEST_DOUBLE;
            for( j=0; j<iNumVals; j++ )
            {
                ReadDouble( ppInputFile[i], &fNEQEM_Loci );
                fprintf( pOutputFile, "\t%.8e", fNEQEM_Loci );
                if( fNEQEM_Loci && fNEQEM_Loci < pfNEQEM_Loci_min[i] )
                {
                    pfNEQEM_Loci_min[i] = fNEQEM_Loci;
                }
            }
            fprintf( pOutputFile, "\n" );
            if( pfNEQEM_Loci_min[i] < LARGEST_DOUBLE )
            {
                pfNEQEM_Loci_min[i] /= 0.7;
            }
            else
            {
                pfNEQEM_Loci_min[i] = 0.0;
            }
        }
    }

    // **** SORT ROUTINE ****
    // Implement a quick and dirty switch-based sort
    do
    {
      bSort = false;
      for( i=0; i<iNumInstruments-1; i++ )
      {
          if( pflog10_TGmin[i] > pflog10_TGmin[i+1] )
          {
              fBuffer = pflog10_TGmin[i];
              pflog10_TGmin[i] = pflog10_TGmin[i+1];
              pflog10_TGmin[i+1] = fBuffer;

              fBuffer = pfEQEM_Loci_min[i];
              pfEQEM_Loci_min[i] = pfEQEM_Loci_min[i+1];
              pfEQEM_Loci_min[i+1] = fBuffer;

              if( bNEQ )
              {
                  fBuffer = pfNEQEM_Loci_min[i];
                  pfNEQEM_Loci_min[i] = pfNEQEM_Loci_min[i+1];
                  pfNEQEM_Loci_min[i+1] = fBuffer;
              }

              bSort = true;
          }
      }
    } while( bSort );
    // **** SORT ROUTINE ****

    // Write the Pottasch emission measures to the data file
    for( i=0; i<iNumInstruments; i++ )
    {
        fprintf( pPottaschOutputFile, "\t%g", pflog10_TGmin[i] );
    }
    fprintf( pPottaschOutputFile, "\n" );

    for( i=0; i<iNumInstruments; i++ )
    {
        fprintf( pPottaschOutputFile, "\t%.8e", pfEQEM_Loci_min[i] );
    }
    fprintf( pPottaschOutputFile, "\n" );

    if( bNEQ )
    {
        for( i=0; i<iNumInstruments; i++ )
        {
            fprintf( pPottaschOutputFile, "\t%.8e", pfNEQEM_Loci_min[i] );
        }
        fprintf( pPottaschOutputFile, "\n" );
    }
}

for( i=0; i<iNumInstruments; i++ )
{
    fclose( ppInputFile[i] );
}
free( ppInputFile );

fclose( pPottaschOutputFile );
fclose( pOutputFile );

free( pflog10_TGmin );
free( pfEQEM_Loci_min );
if( bNEQ )
{
    free( pfNEQEM_Loci_min );
}
}

void CForward::EM_Loci( double flog10_Tmin, double flog10_Tmax, double flog10_step, double flog10_n, char *pszWorkingDirectory, int iFileNumber )
{
FILE *pOutputFile, *pPottaschOutputFile = NULL, **ppInputFile;
char szOutputFilename[256], szInstrumentName[256], szInputFilename[256];
double flog10_T, fBuffer, fCoord;
double fEQEM_Loci, fNEQEM_Loci;
double *pflog10_TGmin, *pfEQEM_Loci_min, *pfNEQEM_Loci_min = NULL;
int iNumVals, iBuffer;
int i, j;
bool bNEQ, bSort;

for( i=0; i<iNumInstruments; i++ )
    ppInstrument[i]->EM_Loci( flog10_Tmin, flog10_Tmax, flog10_step, flog10_n, pszWorkingDirectory, iFileNumber );

// Consolidate the emission measure loci plots for each individual instrument into a single file
// The file format will be:
// Number of instruments (NI)
// Number of temperature values
// Temperature values [ROW]
// Spatial coordinate
// Equilibrium (true) emission measure for instrument 1 [ROW]
// Non-equilibrium emission measure for instrument 1 [ROW]
// .
// .
// Equilibrium (true) emission measure for instrument NI [ROW]
// Non-equilibrium emission measure for instrument NI [ROW]
// Spatial coordinate
// .
// .

// Check whether non-equilibrium emission measure data exists
bNEQ = false;
for( i=0; i<iNumInstruments; i++ )
{
    bNEQ = ppInstrument[i]->IsNonEquilibriumDetector();
    if( bNEQ ) break;
}

pflog10_TGmin = (double*)malloc( sizeof(double) * iNumInstruments );
pfEQEM_Loci_min = (double*)malloc( sizeof(double) * iNumInstruments );
if( bNEQ )
{
    pfNEQEM_Loci_min = (double*)malloc( sizeof(double) * iNumInstruments );
}

// Open the data file for output
sprintf( szOutputFilename, "%s/EM_Loci.%i.dat", pszWorkingDirectory, iFileNumber );
pOutputFile = fopen( szOutputFilename, "w" );
sprintf( szOutputFilename, "%s/EM_Loci_Pottasch.%i.dat", pszWorkingDirectory, iFileNumber );
pPottaschOutputFile = fopen( szOutputFilename, "w" );

// Write the number of instruments into the output files
fprintf( pOutputFile, "%i\n", iNumInstruments );
fprintf( pPottaschOutputFile, "%i\n", iNumInstruments );

// Write the temperature range
iNumVals = (int)( ( ( flog10_Tmax - flog10_Tmin ) / flog10_step ) + 1.0 );
fprintf( pOutputFile, "%i\n", iNumVals );
for( i=0; i<iNumVals; i++ )
{
    flog10_T = flog10_Tmin+(((double)i)*flog10_step);
    fprintf( pOutputFile, "\t%g", flog10_T );
}
fprintf( pOutputFile, "\n" );

// Allocate sufficient space for the pointers to the emission measure data files
ppInputFile = (FILE**)malloc( sizeof(FILE) * iNumInstruments );

// Open the emission measure data files for each instrument
for( i=0; i<iNumInstruments; i++ )
{
    // Get the name of the instrument
    ppInstrument[i]->GetInstrumentName( szInstrumentName );
    // Construct the filename
    sprintf( szInputFilename, "%s/%s.%i.EM_Loci", pszWorkingDirectory, szInstrumentName, iFileNumber );
    ppInputFile[i] = fopen( szInputFilename, "r" );

    // Skip past the temperature data
    fscanf( ppInputFile[i], "%i", &iBuffer );
    for( j=0; j<iNumVals; j++ )
        ReadDouble( ppInputFile[i], &fBuffer );
}

// For each spatial coordinate, write the emission measures
// into the output file following the above format
for ( ;; )
{
    // Get the spatial coordinate
    ReadDouble( ppInputFile[0], &fCoord );
    if( feof( ppInputFile[0]) ) break;
    for( i=1; i<iNumInstruments; i++ )
        ReadDouble( ppInputFile[i], &fCoord );

    // Write the spatial coordinate into the output files
    fprintf( pOutputFile, "%.8e\n", fCoord );
    fprintf( pPottaschOutputFile, "%.8e\n", fCoord );

    for( i=0; i<iNumInstruments; i++ )
    {
        // Get the equilibrium emission measures and write them into the output file
        pfEQEM_Loci_min[i] = LARGEST_DOUBLE;
        pflog10_TGmin[i] = 0.0;
        for( j=0; j<iNumVals; j++ )
        {
            flog10_T = flog10_Tmin+(((double)j)*flog10_step);
            ReadDouble( ppInputFile[i], &fEQEM_Loci );
            fprintf( pOutputFile, "\t%.8e", fEQEM_Loci );
            if( fEQEM_Loci && fEQEM_Loci < pfEQEM_Loci_min[i] )
            {
                pfEQEM_Loci_min[i] = fEQEM_Loci;
                pflog10_TGmin[i] = flog10_T;
            }
        }
        fprintf( pOutputFile, "\n" );
        if( pfEQEM_Loci_min[i] < LARGEST_DOUBLE )
        {
            pfEQEM_Loci_min[i] /= 0.7;
        }
        else
        {
            pfEQEM_Loci_min[i] = 0.0;
        }

        if( bNEQ )
        {
            // Get the non-equilibrium emission measures and write them into the output file
            pfNEQEM_Loci_min[i] = LARGEST_DOUBLE;
            for( j=0; j<iNumVals; j++ )
            {
                ReadDouble( ppInputFile[i], &fNEQEM_Loci );
                fprintf( pOutputFile, "\t%.8e", fNEQEM_Loci );
                if( fNEQEM_Loci && fNEQEM_Loci < pfNEQEM_Loci_min[i] )
                {
                    pfNEQEM_Loci_min[i] = fNEQEM_Loci;
                }
            }
            fprintf( pOutputFile, "\n" );
            if( pfNEQEM_Loci_min[i] < LARGEST_DOUBLE )
            {
                pfNEQEM_Loci_min[i] /= 0.7;
            }
            else
            {
                pfNEQEM_Loci_min[i] = 0.0;
            }
        }
    }

    // **** SORT ROUTINE ****
    // Implement a quick and dirty switch-based sort
    do
    {
      bSort = false;
      for( i=0; i<iNumInstruments-1; i++ )
      {
          if( pflog10_TGmin[i] > pflog10_TGmin[i+1] )
          {
              fBuffer = pflog10_TGmin[i];
              pflog10_TGmin[i] = pflog10_TGmin[i+1];
              pflog10_TGmin[i+1] = fBuffer;

              fBuffer = pfEQEM_Loci_min[i];
              pfEQEM_Loci_min[i] = pfEQEM_Loci_min[i+1];
              pfEQEM_Loci_min[i+1] = fBuffer;

              if( bNEQ )
              {
                  fBuffer = pfNEQEM_Loci_min[i];
                  pfNEQEM_Loci_min[i] = pfNEQEM_Loci_min[i+1];
                  pfNEQEM_Loci_min[i+1] = fBuffer;
              }

              bSort = true;
          }
      }
    } while( bSort );
    // **** SORT ROUTINE ****

    // Write the Pottasch emission measures to the data file
    for( i=0; i<iNumInstruments; i++ )
    {
        fprintf( pPottaschOutputFile, "\t%g", pflog10_TGmin[i] );
    }
    fprintf( pPottaschOutputFile, "\n" );

    for( i=0; i<iNumInstruments; i++ )
    {
        fprintf( pPottaschOutputFile, "\t%.8e", pfEQEM_Loci_min[i] );
    }
    fprintf( pPottaschOutputFile, "\n" );

    if( bNEQ )
    {
        for( i=0; i<iNumInstruments; i++ )
        {
            fprintf( pPottaschOutputFile, "\t%.8e", pfNEQEM_Loci_min[i] );
        }
        fprintf( pPottaschOutputFile, "\n" );
    }
}

for( i=0; i<iNumInstruments; i++ )
{
    fclose( ppInputFile[i] );
}
free( ppInputFile );

fclose( pPottaschOutputFile );
fclose( pOutputFile );

free( pflog10_TGmin );
free( pfEQEM_Loci_min );
if( bNEQ )
{
    free( pfNEQEM_Loci_min );
}
}