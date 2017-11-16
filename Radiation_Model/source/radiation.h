// ****
// *
// * Radiation Class Definition for Radiative Emission Model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 11/15/2017
// *
// ****


#include "element.h"


class CRadiation {

    private:
	
    // The number of elements for which ion fractional populations are available
    int NumElements;
	
    // Pointer to an array of pointers of type PELEMENT
    PPELEMENT ppElements;
	
    // Pointer to an array containing each element's atomic number
    // The offset of the atomic number corresponds to the offset of the element object in the ppElements array
    int *pZ;

    // The temperature and density values in log_10 form
    int NumTemp, NumDen;
    double *pTemp, *pDen;
		
    // Pointer to the factor total phi( n, T ) for all of the elements
    double *pTotalPhi;

    // Function to initialise the radiation object with a set of elements
    void Initialise( char *szFilename );

    // Function to open and read the ranges data file
    void OpenRangesFile( char *szRangesFilename );

    // Function to calculate the factor total phi( n, T ), which is multiplied by n^2 to calculate the radiated energy
    void CalculateTotalPhi( void );

    // Function to free all allocated memory
    void FreeAll( void );
	
    public:

    // Constructor
    CRadiation( char *szFilename );
	
    // Destructor
    ~CRadiation( void );

    // Function to return a list of the atomic numbers of the elements comprising the radiation model
    int *pGetAtomicNumbers( int *iNumElements );

    // Function to return the abundance of a specified element
    double GetAbundance( int iZ );

    // Function to return the ion fractional populations of a particular element at a specified temperature and density in equilibrium
    void GetEquilIonFrac( int iZ, double *pni, double flog_10T );
    void GetEquilIonFrac( int iZ, double *pni, double flog_10T, double flog_10n );

    // Functions to write either a given set of ion fractional populations or all fractional populations to a data file
    void WriteEquilIonFracToFile( void *pFile, int iZ, double flog_10T );
    void WriteEquilIonFracToFile( void *pFile, int iZ, double flog_10T, double flog_10n );
    void WriteEquilIonFracToFile( void *pFile, double flog_10T );
    void WriteEquilIonFracToFile( void *pFile, double flog_10T, double flog_10n );

    // Function to normalise the sum total of the ion fractional populations to 1
    void Normalise( int iZ, double *pni, double fTotal );

    // Functions to calculate the amount of energy radiated in equilibrium
    double GetRadiation( int iZ, int iIon, double flog_10T, double flog_10n );
    double GetRadiation( int iZ, double flog_10T, double flog_10n );
    double GetRadiation( double flog_10T, double flog_10n );

    // Functions to calculate the rate of change with respect to time of the fractional populations of the ions and the characteristic time-scale
    void Getdnibydt( int iZ, double flog_10T, double flog_10n, double *pni0, double *pni1, double *pni2, double *pni3, double *pni4, double *s, double *s_pos, double *pv, double delta_s, double *pdnibydt, double *pTimeScale );
    void GetAlldnibydt( double flog_10T, double flog_10n, double **ppni0, double **ppni1, double **ppni2, double **ppni3, double **ppni4, double *s, double *s_pos, double *pv, double delta_s, double **ppdnibydt, double *pTimeScale );

    // Functions to calculate the amount of energy radiated in nonequilibrium
    double GetRadiation( int iZ, int iIon, double flog_10T, double flog_10n, double ni );
    double GetRadiation( int iZ, double flog_10T, double flog_10n, double *pni );
    double GetRadiation( double flog_10T, double flog_10n, double **ppni );
	
    // Functions to calculate energy radiated based upon power-laws
    double GetPowerLawRad( double flog_10T, double fne, double fnH );
    double GetFreeFreeRad( double flog_10T, double flog_10n );

};

typedef CRadiation* PRADIATION;
