// ****
// *
// * Optically-thick Ion Class Definition for Radiative Emission Model
// *
// * Based on the formulation of Carlsson, M., & Leenaarts, J., 2012, A&A, 539, A39
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/06/2017
// *
// ****

#include "../config.h"

#ifdef OPTICALLY_THICK_RADIATION

#define AMU 1.660538782E-24                     // [g]
#define AMU_FILENAME "Radiation_Model/atomic_data/masses/masses.amu"

#define OPTICALLY_THICK_TEMPERATURE 2.4E4       // [K]

class COpticallyThickIon {

    private:

    // The atomic number of the element
    int Z;

    // The abundance of the element relative to hydrogen and
    // the number of hydrogen particles per gram
    double Abund, N_H;

    // Ion population fraction as a function of temperature
    int iIonFracDP;
    double **ppIonFrac;

    // Emissivity as a function of temperature
    int iEmissDP;
    double **ppOriginalEmiss, **ppEmiss;

    // Escape probability as a function of optical depth or column density
    int iEscProbDP;
    double **ppEscProb;

    // Coefficient of thermal conduction as a function of temperature
    int ikappa_0DP;
    double **ppkappa_0;

    // Functions to initialise the class and free allocated memory
    void Initialise( int iZ, char *szIon, char *szAbundFilename );
    void FreeAll( void );

    // Functions to open and read the ion data files
    void GetAbundData( char *szAbundFilename );
    void GetIonFracData( char *szIonFracFilename );
    void GetEmissData( char *szEmissFilename );
    void GetEscProbData( char *szEscProbFilename );
    void Getkappa_0Data( char *szkappa_0Filename );

    // Functions to return particular elements of ion data
    double GetEmiss( double flog_10T );
    double GetEmiss( double flog_10T, double fIonFrac );
    double GetEscProb( double fX );
   
    public:

    // Constructor
    COpticallyThickIon( int iZ, char *szIon, char *szAbundFilename );
	
    // Destructor
    ~COpticallyThickIon( void );

    // Functions to return particular elements of ion data
    double GetIonFrac( double flog_10T );
    double Getkappa_0( double flog_10T );

    // Function to calculate the optically-thick emission
    // fX can be optical depth (HI) or column density (MgII, CaII)
    double GetVolumetricLossRate( double fLog_10T, double fX, double n_e_rho );
    double GetVolumetricLossRate( double fLog_10T, double fIonFrac, double fX, double n_e_rho );
    
};

typedef COpticallyThickIon* POPTICALLYTHICKION;
typedef COpticallyThickIon** PPOPTICALLYTHICKION;

#endif // OPTICALLY_THICK_RADIATION