// ****
// *
// * Ion Class definition
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/14/2017
// *
// ****


// The cut-off temperature below which there is no emission
#define LOG_LINE_EMISSION_CUTOFF_TEMP 4.0

// Choose the nearest wavelength if the specified wavelength does
// not match any stored wavelengths
#define CHOOSE_NEAREST_WAVELENGTH

#define UNITS_DN                                // Use intensity units of DN pixel^-1 s^-1
// #define UNITS_PHOTONS                        // Use intensity units of photons cm^-2 s^-1 sr^-1
// #define UNITS_ERG                            // Use intensity units of erg cm^-2 s^-1 sr^-1


class CIon {

    private:

    // The ion label, atomic and spectroscopic numbers, and its mass (g)
    char szLabel[16];
    int Z, SpecNum;
    double mi;
	
    // The number of lines emitted by the ion in the wavelength sensitivity range of the instrument
    int NumLines;

    // Pointer to a list of the individual wavelengths
    double *pLambda;
	
    // The number of temperature and density values
    int NumTemp, NumDen, NumTempxNumDen;
	
    // The temperature and density values in log_10 form
    double *pTemp, *pDen;
		
    // Pointer to an array of pointers, each pointing to the emissivity data for a particular line in a NumTemp * NumDen size array
    double **ppEmiss;
    double *pIonEmiss;

    // Function to initialise the ion object
    void Initialise( char *pszLabel, int iZ, int iSpecNum, char *szEmissFilename, double *pRespFunc, int iNumDataPoints );

    // Function to free all allocated memory
    void FreeAll( void );
	
    // Function to open and read the emissivity data file
    void OpenEmissivityFile( char *szEmissFilename, double *pRespFunc, int iNumDataPoints );

    // Function to calculate the total line emission as a function of temperature and density (summed over the wavelength range)
    void CalculateTotalIonEmission( void );

    public:
	
    // Constructor
    CIon( char *pszLabel, int iZ, int iSpecNum, char *szEmissFilename, double *pRespFunc, int iNumDataPoints );
	
    // Destructor
    ~CIon( void );

    // Function to return the atomic number of the ion
    int GetAtomicNumber( void );

    // Function to return the spectroscopic number of the ion
    int GetSpecNumber( void );

    // Function to return the mass (g) of the ion
    double GetMass( void );

    // Function to return the number of emission lines for the ion
    int GetNumLines( void );

    // Function to write the emission line wavelength list into a pre-allocated array
    void GetLineList( double *pfLineList );

    // Function to write a list of the emission lines within the wavelength response region of the instrument to a file
    void WriteIonLineListToFile( FILE *pFile );

    // Function to return the line emission as a function of wavelength, temperature and density
    double GetLineEmission( double fLambda, double flog10T, double flog10n );
	
    // Function to return the total line emission as a function of temperature and density (summed over the wavelength range)
    double GetIonEmission( double flog10T, double flog10n );

};

typedef CIon* PION;
typedef CIon** PPION;
