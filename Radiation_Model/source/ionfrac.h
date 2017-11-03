// ****
// *
// * Ionisation Fraction Class Definition for Radiative Emission Model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/14/2017
// *
// ****


#include "radiation.h"


class CIonFrac {

    private:

    // Pointer to a radiation object
    PRADIATION pRadiation;

    // The number of elements for which ion fractional populations are available
    int NumElements;

    // Pointer to an array containing each element's atomic number
    // The offset of the atomic number corresponds to the offset of the ionisation fraction object in the ppIonFrac array
    int *pZ;

    // Pointer to an array of pointers containing the fractional populations of the ions for each element at the current temperature
    double **ppIonFrac;

    // Pointer to an array of pointers containing the rate of change with respect to time of the fractional population of the ions for each element at the current temperature
    double **ppdnibydt;

    // Function to initialise the object
    void Initialise( CIonFrac *pIonFrac, char *szFilename, PRADIATION pRadiationObj );

    // Function to free all allocated memory
    void FreeAll( void );
	
    public:

    // Constructor
    CIonFrac( CIonFrac *pIonFrac, char *szFilename, PRADIATION pRadiationObj );
	
    // Destructor
    ~CIonFrac( void );

    // Functions to return a pointer to the fractional populations of the ions at the current temperature
    double** ppGetIonFrac( void );
    double* pGetIonFrac( int iZ );

    // Functions to write either all of the ion fractional populations or a particular set of fractional populations to a data file
    void WriteAllIonFracToFile( void *pFile );
    void WriteIonFracToFile( void *pFile, int iZ );

    // Functions to read either all of the ion fractional populations or a particular set of fractional populations from a data file
    void ReadAllIonFracFromFile( void *pFile );
    void ReadIonFracFromFile( void *pFile, int iZ );

    // Functions to return a pointer to the rate of change with respect to time of the fractional populations of the ions at the current temperature
    double** ppGetdnibydt( void );
    double* pGetdnibydt( int iZ );

    // Functions to integrate the fractional populations of the ions
    void IntegrateAllIonFrac( double delta_t );
    void IntegrateIonFrac( int iZ, double delta_t );

    // Function to return the number of elements for which ion fractional population data is stored and the atomic numbers of the elements
    int* pGetElementInfo( int *pNumElements );

    // Functions to overwrite fractional populations in the current object
    void CopyAllIonFrac( CIonFrac *pIonFrac ); 
    void CopyIonFrac( int iZ, CIonFrac *pIonFrac );

    // Functions to overwrite dnibydt's in the current object
    void CopyAlldnibydt( CIonFrac *pIonFrac );
    void Copydnibydt( int iZ, CIonFrac *pIonFrac );

    // Functions to interpolate the fractional populations in the current object from the adjacent fractional populations
    void InterpolateAllIonFrac( double *x, double ***pppIonFrac, int iPoints, double s );
    void InterpolateIonFrac( int iZ, double *x, double ***pppIonFrac, int iPoints, double s );

    // Functions to reset the fractional populations to their equilibrium values at the specified temperature and density
    void ResetIonFrac( int iZ, double flog_10T );
    void ResetIonFrac( int iZ, double flog_10T, double flog_10n );
    void ResetAllIonFrac( double flog_10T );
    void ResetAllIonFrac( double flog_10T, double flog_10n );

};

typedef CIonFrac* PIONFRAC;
