// ****
// *
// * Class definition of the radiative transition rates
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 03/24/2020
// *
// ****

#include "../config.h"

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE

// Use this pre-processor directive to select the method to invert the matrix equation for the hydrogen level populations
// LU decomposition seems to be faster in most cases and is the most stable
// #define USE_SINGLE_VALUE_DECOMPOSITION
#ifndef USE_SINGLE_VALUE_DECOMPOSITION
	#define USE_LOWER_UPPER_DECOMPOSITION
#endif // USE_SINGLE_VALUE_DECOMPOSITION

// Definitions for the temperature at which the foot-point density, used to scale T_b^top as a function of time, is defined in the NLTE chromosphere
// and the foot-point density itself (scaled by the quantities defined below for each transition)
#define NLTE_T_FP	9E3			// K
#define NLTE_n0     1.9e11     // cm^-3
// Indices for each transition that are used to scale the ratio of the density to the foot-point density. When multiplied by the temperature defined
// above this gives the brightness temperature for each transition according to the formula: T_b = T_0 * (n/n_0)^m
    #define NLTE_m_Ha   0.1188
    #define NLTE_m_Hb   0.1116
    #define NLTE_m_Hg   0.1061
    #define NLTE_m_Pa   0.1460
    #define NLTE_m_Pb   0.1402
    #define NLTE_m_Bra  0.1979

#define MAX_ITERATIONS			300
#define CONVERGENCE_EPSILON		0.1
#define CONVERGENCE_CONDITION	1E-6


// **** RADIATIVE TRANSITION RATES CLASS ****

// Define the radiative transition rates class
class CRadiativeRates {
	
    private:

    // Radiative bound-bound rates:
	int iBB_TRvals, iBB_Tvals;
	double *pfBB_logT;
	// Photoexcitation:
	double **ppfBB_lu;
	// Radiative decay:
	double **ppfBB_ul;

    // Radiative bound-free rates:
	int iBF_TRvals, iBF_Tvals;
	double *pfBF_logT;
	// Photoionization
	double **ppfBF;

    // Radiative free-bound rates:
	// The radiation temperatures are in columns and the electron temperatures are in rows
	int iFB_TRvals, iFB_radTvals, iFB_eTvals;
	double *pfFB_logradT, *pfFB_logeT;
	// Radiative de-excitation
	double ***pppfFB;
    
    // Collisional bound-bound and bound-free rates
	int iColl_Tvals;
	double *pfColl_logT;
	// Collisional excitation
	int iColl_exTRvals;
	double **ppfColl_ex_lu, **ppfColl_ex_ul;
	// Collisional ionization
	int iColl_ionTRvals;
	double **ppfColl_ion, **ppfColl_rec;
    
    int iNBBT, iNBFT;
    double *pTrt, **ppScaledTrt, *pnu0, *pTeZ_c, *pZ_c_LEFT, *pZ_c_RIGHT, *pMcZ_c_LEFT, *pMcZ_c_RIGHT;
    double *pterm1, *pterm2, *pH;

#ifdef USE_LOWER_UPPER_DECOMPOSITION
    // Functions to calculate the LU decomposition of a square matrix A and back-substitute to solve a matrix equation A x = b
    void SolveLinearEq(double **a, int n, double *b, bool improve);
    // LU decomposition of matrix A of size n x n
	void LUdcmp(double **a, int n, int *index, double *d);
	// Back substitution of the LU decomposition to solve for vector x
    void LUbksb(double **a, int n, int *index, double *b);
#endif // USE_LOWER_UPPER_DECOMPOSITION

#ifdef USE_SINGLE_VALUE_DECOMPOSITION
    // Functions to calculate the singular value decomposition of a matrix A and back-substitute to solve a matrix equation A x = b
    // (Used by SolveHIIFraction)
    // Singular value decomposition of matrix A of size m x n
	int svdcmp(double **a, int m, int n, double *w, double **v);
	// Back substitution of the singular value decomposition to solve for vector x
	void svbksb(double **u, double *w, double **v, int m, int n, double *b, double *x);
	// Pythagorean distance
	double pythag(double a, double b);
#endif // USE_SINGLE_VALUE_DECOMPOSITION
	
    void GetBBRates( char *pszBBRatesFile );
    void GetBFRates( char *pszBFRatesFile );
    void GetFBRates( char *pszFBRatesFile );
    void GetCollRates( char *pszCollRatesFile );

    void Initialise( char *pszRatesFiles );
    void FreeAll( void );

    public:

    // Constructor
    CRadiativeRates( char *pszRatesFiles );
	
    // Destructor
    ~CRadiativeRates( void );

    // Radiative transition rates class functions

    void GetBoundBoundRates( double *pfBB_lu, double *pfBB_ul, double *pflog10T );
    void GetBoundFreeRates( double *pfBF, double *pflog10T );
    void GetFreeBoundRates( double *pfFB, double *pflog10radT, double flog10eT, double fne );
    void GetFreeBoundRatesRH( double *pfFB, double *pflog10radT, double flog10eT, double fne, double *pfLevel_Ratio );
    void GetCollisionalRates( double *pfColl_ex_lu, double *pfColl_ex_ul, double *pfColl_ion, double *pfColl_rec, double flog10T, double fne );
    void GetCollisionalRatesRH( double *pfColl_ex_lu, double *pfColl_ex_ul, double *pfColl_ion, double *pfColl_rec, double flog10T, double fne );
    
    void SolveHIIFraction( double *pfHstate, double *pfColl_ex_lu, double *pfColl_ex_ul, double *pfColl_ion, double *pfColl_rec, double *pfBB_lu, double *pfBB_ul, double *pfBF, double *pfFB );

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
    void GetAllDel_Hstate_dot_v( double *pHstate0, double *pHstate1, double *pHstate2, double *pHstate3, double *pHstate4, double *s, double *s_pos, double *pv, double *cross_section, double cell_volume, double *pDel_Hstate_dot_v );
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
    void GetAllDel_Hstate_dot_v( double *pHstate0, double *pHstate1, double *pHstate2, double *pHstate3, double *pHstate4, double *s, double *s_pos, double *pv, double delta_s, double *pDel_Hstate_dot_v );
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
    void Normalize( double *pHstate );

    int GetNumberOfTransitions( void );
    double GetZ_c_LEFT( int i );
    double GetZ_c_RIGHT( int i );
    double GetMcZ_c_LEFT( int i );
    double GetMcZ_c_RIGHT( int i );
    double GetTeZ_c( int i );
    double GetNu0( int i );
    double GetH( int i );
	double GetTrt( int i );
	double GetScaledTrt( int i, bool iFlag );
    void SetZ_c_LEFT( int i, double Z_c );
    void SetZ_c_RIGHT( int i, double Z_c );
    void SetMcZ_c_LEFT( int i, double McZ_c );
    void SetMcZ_c_RIGHT( int i, double McZ_c );
    void SetScaledTrt( int i, double Trt, bool iFlag );
	
};

// Define a type for the radiative transition rates class
typedef CRadiativeRates* PRADIATIVERATES;

#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
