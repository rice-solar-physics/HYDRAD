// ****
// *
// * Radiation Class Definition for Radiative Emission Model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/04/2020
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
		
    // Pointer to the factor total phi( ne, T ) for all of the elements
    double *pTotalPhi;

    // Function to initialise the radiation object with a set of elements
    void Initialise( char *szFilename );

    // Function to open and read the ranges data file
    void OpenRangesFile( char *szRangesFilename );

    // Function to calculate the factor total phi( ne, T ), which is multiplied by ne x nH to calculate the radiated energy
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

	// Functions to return the ionization and recombination rates of a particular element and ion, at a specified temperature and density
	void GetRates( int iZ, int iIon, double flog_10T, double *pfIonRate, double *pfRecRate );
	void GetRates( int iZ, int iIon, double flog_10T, double flog_10n, double *pfIonRate, double *pfRecRate );

    // Functions to return the ion fractional populations of a particular element at a specified temperature and density in equilibrium
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
    double GetRadiation( int iZ, int iIon, double flog_10T, double fne, double fnH );
    double GetRadiation( int iZ, double flog_10T, double fne, double fnH );
    double GetRadiation( double flog_10T, double fne, double fnH );

    // Functions to calculate the rate of change with respect to time of the fractional populations of the ions and the characteristic time-scale
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
    void Getdnibydt( int iZ, double flog_10T, double flog_10n, double *pni0, double *pni1, double *pni2, double *pni3, double *pni4, double *s, double *s_pos, double *pv, double *cross_section, double cell_volume, double *pdnibydt, double *pTimeScale );
    void GetAlldnibydt( double flog_10T, double flog_10n, double **ppni0, double **ppni1, double **ppni2, double **ppni3, double **ppni4, double *s, double *s_pos, double *pv, double *cross_section, double cell_volume, double **ppdnibydt, double *pTimeScale );
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
    void Getdnibydt( int iZ, double flog_10T, double flog_10n, double *pni0, double *pni1, double *pni2, double *pni3, double *pni4, double *s, double *s_pos, double *pv, double delta_s, double *pdnibydt, double *pTimeScale );
    void GetAlldnibydt( double flog_10T, double flog_10n, double **ppni0, double **ppni1, double **ppni2, double **ppni3, double **ppni4, double *s, double *s_pos, double *pv, double delta_s, double **ppdnibydt, double *pTimeScale );
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

    // Functions to calculate the amount of energy radiated in nonequilibrium
    double GetRadiation( int iZ, int iIon, double flog_10T, double fne, double fnH, double ni );
    double GetRadiation( int iZ, double flog_10T, double fne, double fnH, double *pni );
    double GetRadiation( double flog_10T, double fne, double fnH, double **ppni );
	
    // Functions to calculate energy radiated based upon power-laws
    double GetPowerLawRad( double flog_10T, double fne, double fnH );
    double GetFreeFreeRad( double flog_10T, double fne, double fnH );

};

typedef CRadiation* PRADIATION;
