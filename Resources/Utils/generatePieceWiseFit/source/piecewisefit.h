// ****
// *
// * Class definition of the piece-wise polynomial fitting method
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 07/20/2020
// *
// ****


// #define VERBOSE


// **** PIECE-WISE FIT CLASS ****

// Define the piece-wise fit class
class CPieceWiseFit {

    private:
	
	// Order of the polynomial fit (e.g. O(6): f(x) = a_0x^0 + ... + a_6x^6; 7 coefficients) 
	int iPolyOrder;
	// Number of sub-domains, each with their own piece-wise fit
	int inumSD;
	// Array of sub-domain boundaries in the range [0,1]
	double *pfSDBoundary;
	// Array of coefficients for the fitted polynomial in each sub-domain
	double **ppfSDCoefficient;
	// Minimum and maximum values of the fitted quantity in each sub-domain
	double **ppfSDMinMax;

	// Open and read the contents of the data file which defines the piece-wise fit
	void OpenPieceWiseFit( char *pszInputFilename );

	// Generate a new piece-wise fit
	void GeneratePieceWiseFit( char *pszInputFilename, char *pszOutputFilename );

    // Function to free all allocated memory
    void FreeAll( void );
	
	public:
	
    // Constructor
    CPieceWiseFit( char *pszInputFilename );
    CPieceWiseFit( char *pszInputFilename, char *pszOutputFilename );
	
    // Destructor
    ~CPieceWiseFit( void );

	// Show the parameter values for the piece-wise fit
	void ShowPieceWiseFit( void );
	// Return the piece-wise fit at a specified location
	double GetPieceWiseFit( double fx );
	// Return the nth derivative of the piece-wise fit at a specified location
	double GetDerivative( double fx, int iDerivative );
	
};

// Define a type for the piece-wise fit class
typedef CPieceWiseFit* PPIECEWISEFIT;