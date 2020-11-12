// ****
// *
// * Function bodies for the class definition of the 
// * time-dependent hydrodynamic equations, inherited by the 
// * adaptive mesh class
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 11/11/2020
// *
// ****


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "eqns.h"
#include "../../Resources/source/constants.h"
#include "../../Resources/source/file.h"
#include "../../Resources/source/fitpoly.h"

#ifdef OPENMP
	#include <omp.h>
#endif // OPENMP


double fLogLambda_ei( double Te, double Ti, double n )
{
double limit;

	Te *= 1.29199251618e-4;
	Ti *= 1.29199251618e-4;

	limit = Ti * ( ELECTRON_MASS / AVERAGE_PARTICLE_MASS );

	if( limit < Te && Te < 10.0 )
	{
    	return 23.0 - log( sqrt( n ) * pow( Te, (-1.5) ) );
	}
	else if( limit < 10.0 && 10.0 < Te )
	{
    	return 24.0 - log( sqrt( n ) * ( 1.0 / Te ) );
	}
	else // if( Te < limit )
	{
    	return 30.0 - log( sqrt( n ) * pow( Ti, (-1.5) ) );
	}
}

double fLogLambda_ii( double Ti, double n )
{
double term1;

	// Assume that proton-proton collisions dominate ion-ion interactions
	// sqrt(2.0) * Z^3 = 1.4142135623730950488016887242097 * 1.0 * 1.0 * 1.0
	term1 = 1.4142135623730950488016887242097 * ( sqrt( n ) / pow( Ti, 1.5 ) );
	return 23.0 - log( term1 );
}


// Constructor
CEquations::CEquations( void )
{
	Initialise();
}

// Destructor
CEquations::~CEquations( void )
{
	FreeAll();
}

void CEquations::Initialise( void )
{
FILE *pFile;
double fTemp;
int iTemp;

#if defined (OPENMP) || defined (USE_KINETIC_MODEL)
	ppCellList = NULL;
#endif // OPENMP || USE_KINETIC_MODEL

	// Get the HYDRAD initial conditions
	pFile = fopen( "HYDRAD/config/hydrad.cfg", "r" );
		// Get the initial profiles
		fscanf( pFile, "%s", Params.Profiles );
		// Get the gravity polynomial-fit coefficients filename
		fscanf( pFile, "%s", Params.GravityFilename );
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
		// Get the magnetic field polynomial-fit coefficients filename
		fscanf( pFile, "%s", Params.MagneticFieldFilename );
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
		// Get the duration
		ReadDouble( pFile, &Params.Duration );
		// Get the output period
		ReadDouble( pFile, &Params.OutputPeriod );
	fclose( pFile );

	// Get the loop length from the profiles file
	pFile = fopen( Params.Profiles, "r" );
		ReadDouble( pFile, &fTemp );
		fscanf( pFile, "%i", &iTemp );
		ReadDouble( pFile, &Params.L );
	fclose( pFile );

#ifdef USE_KINETIC_MODEL
	// Get the tabulated values from tables I and II (for Z = 1) in Spitzer & Harm, 1953, Phys. Rev., 89, 977
	Get_SH_Table();
#endif // USE_KINETIC_MODEL

	// Initialise the gravitational acceleration profile in the field-aligned direction
	pGravity = new CPieceWiseFit( Params.GravityFilename );

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
	// Initialise the magnetic field profile in the field-aligned direction
	pMagneticField = new CPieceWiseFit( Params.MagneticFieldFilename );
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

	// Create the heating object
	pHeat = new CHeat();

	// Create the radiation objects
	pRadiation = new CRadiation( (char *)"Radiation_Model/config/elements_neq.cfg" );
	pRadiation2 = new CRadiation( (char *)"Radiation_Model/config/elements_eq.cfg" );
	lower_radiation_temperature_boundary = MINIMUM_RADIATION_TEMPERATURE + ZERO_OVER_TEMPERATURE_INTERVAL;
#ifdef BEAM_HEATING
	// Identify which radiation objects contain the data for hydrogen and helium
	// If the element does not exist in an object then it will return zero for the abundance
	if( !(pRadiation->GetAbundance( 1 ) ) ) pRadiation_H = pRadiation2;
	else pRadiation_H = pRadiation;

	if( !(pRadiation->GetAbundance( 2 ) ) ) pRadiation_He = pRadiation2;
	else pRadiation_He = pRadiation;
#endif // BEAM_HEATING

#ifdef OPTICALLY_THICK_RADIATION
	// Create the optically-thick ion objects
	pHI = new COpticallyThickIon( 1, (char *)"h_1", (char *)"Radiation_Model/atomic_data/abundances/asplund.ab" );
	pMgII = new COpticallyThickIon( 12, (char *)"mg_2", (char *)"Radiation_Model/atomic_data/abundances/asplund.ab" );
	pCaII = new COpticallyThickIon( 20, (char *)"ca_2", (char *)"Radiation_Model/atomic_data/abundances/asplund.ab" );
#ifdef NLTE_CHROMOSPHERE
	pRadiativeRates = new CRadiativeRates( (char *)"Radiation_Model/atomic_data/OpticallyThick/radiative_rates/rates_files.txt" );
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
}

void CEquations::FreeAll( void )
{
#if defined (OPENMP) || defined (USE_KINETIC_MODEL)
	free( ppCellList );
#endif // OPENMP || USE_KINETIC_MODEL

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
	// Free the memory allocated to the magnetic field profile in the field-aligned direction
	delete pMagneticField;
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

	// Free the memory allocated to the gravitational acceleration profile in the field-aligned direction
	delete pGravity;

	// Delete the heating object
	delete pHeat;

	// Delete the radiation objects
	delete pRadiation;
	delete pRadiation2;

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
	delete pRadiativeRates;
#endif // NLTE_CHROMOSPHERE
	// Delete the optically-thick ion objects
	delete pHI;
	delete pMgII;
	delete pCaII;
#endif // OPTICALLY_THICK_RADIATION
}

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
void CEquations::CalculateInitialPhysicalQuantities( void )
{
PCELL pNextActiveCell;
CELLPROPERTIES CellProperties, NextCellProperties;

#ifdef USE_JB
	Tc = 0.0;
	Te_max = 0.0;
	double old_s = 0.0, old_T = 0.0, tau_JB = 1.0, L_T;
#endif // USE_JB

// General variables
double term1, term2, x[3], y[3];
int i, j;

double frho_c, fLeftFPn_H;
#ifdef OPEN_FIELD
#else // OPEN_FIELD
	double fRightFPn_H;
#endif // OPEN_FIELD
double fH, fnu, fTrt[2], fMcZ_c, fA;
double fTeZ_c, fZ_c;
double flog10_Trad[10];
double fBB_lu[6], fBB_ul[6], fBF[4], fFB[4], fColl_ex_lu[10], fColl_ex_ul[10], fColl_ion[5], fColl_rec[5];

#ifdef BEAM_HEATING
	double fQbeam, fLambda1, fLambda2, log_fAvgEE;
	log_fAvgEE = log( pHeat->GetAvgEE() );
	fLambda2 = 25.1 + log_fAvgEE;
#endif // BEAM_HEATING

	memset( &CellProperties, 0, sizeof(CELLPROPERTIES) );	// Avoids a warning that variables might be used undefined

	// We need to make T_b^top a time-dependent function, which means finding the foot-point densities at each leg of the loop prior to executing the
	// main loop. The foot-point density is defined as the density at a particular temperature (approx. the H_alpha formation temperature of 9,000 K)
	frho_c = 0.0;
	fLeftFPn_H = 0.0;
	pNextActiveCell = pCentreOfCurrentRow;
	while( pNextActiveCell )
	{
		pActiveCell = pNextActiveCell;
		pActiveCell->GetCellProperties( &CellProperties );
		// Calculate the electron density and temperature
    	CellProperties.n[ELECTRON] = CellProperties.rho_e / ELECTRON_MASS;
	    // GAMMA_MINUS_ONE / BOLTZMANN_CONSTANT = 4.830917874e15
    	CellProperties.T[ELECTRON] = (4.830917874e15) * ( CellProperties.TE_KE[1][ELECTRON] / CellProperties.n[ELECTRON] );
		// Calculate the column mass density
		frho_c += CellProperties.rho[1] * CellProperties.cell_width;
		CellProperties.rho_c = frho_c;
		// ( GAMMA_MINUS_ONE * ELECTRON_MASS ) / BOLTZMANN_CONSTANT = 4.40066838247E-12
		if( ( (4.40066838247e-12) * ( CellProperties.TE_KE[1][ELECTRON] / CellProperties.rho_e ) ) <= NLTE_T_FP && !fLeftFPn_H  )
			fLeftFPn_H = CellProperties.rho[1] / AVERAGE_PARTICLE_MASS;
	
		pActiveCell->UpdateCellProperties( &CellProperties );
	
		pNextActiveCell = pActiveCell->pGetPointer( LEFT );
	}

#ifdef OPEN_FIELD
#else // OPEN_FIELD
	frho_c = 0.0;
	fRightFPn_H = 0.0;
	pNextActiveCell = pCentreOfCurrentRow;
	while( pNextActiveCell )
	{
		pActiveCell = pNextActiveCell;
		pActiveCell->GetCellProperties( &CellProperties );	
		// Calculate the electron density and temperature
    	CellProperties.n[ELECTRON] = CellProperties.rho_e / ELECTRON_MASS;
	    // GAMMA_MINUS_ONE / BOLTZMANN_CONSTANT = 4.830917874e15
    	CellProperties.T[ELECTRON] = (4.830917874e15) * ( CellProperties.TE_KE[1][ELECTRON] / CellProperties.n[ELECTRON] );
		// Calculate the column mass density
		frho_c += CellProperties.rho[1] * CellProperties.cell_width;
		CellProperties.rho_c = frho_c;
		// ( GAMMA_MINUS_ONE * ELECTRON_MASS ) / BOLTZMANN_CONSTANT = 4.40066838247E-12
		if( ( (4.40066838247e-12) * ( CellProperties.TE_KE[1][ELECTRON] / CellProperties.rho_e ) ) <= NLTE_T_FP && !fRightFPn_H )
			fRightFPn_H = CellProperties.rho[1] / AVERAGE_PARTICLE_MASS;

		pActiveCell->UpdateCellProperties( &CellProperties );
	
		pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
	}
#endif // OPEN_FIELD

	for(i=0; i<pRadiativeRates->GetNumberOfTransitions(); i++)
    {
		fTrt[0] = fTrt[1] = pRadiativeRates->GetTrt(i);
		// Identify the current transition
		switch( i ) {
						
			case 0 :	// H-alpha
                if( fLeftFPn_H > NLTE_n0 )
                    fTrt[0] *= pow((fLeftFPn_H/NLTE_n0), NLTE_m_Ha);
#ifdef OPEN_FIELD
#else // OPEN_FIELD
                if( fRightFPn_H > NLTE_n0 )
                    fTrt[1] *= pow((fRightFPn_H/NLTE_n0), NLTE_m_Ha);
#endif // OPEN_FIELD
                break;
                    
            case 1 :	// H-beta
                if( fLeftFPn_H > NLTE_n0 )
                    fTrt[0] *= pow((fLeftFPn_H/NLTE_n0), NLTE_m_Hb);
#ifdef OPEN_FIELD
#else // OPEN_FIELD
                if( fRightFPn_H > NLTE_n0 )
                    fTrt[1] *= pow((fRightFPn_H/NLTE_n0), NLTE_m_Hb);
#endif // OPEN_FIELD
                break;
                    
            case 2 :	// H-gamma
                if( fLeftFPn_H > NLTE_n0 )
                    fTrt[0] *= pow((fLeftFPn_H/NLTE_n0), NLTE_m_Hg);
#ifdef OPEN_FIELD
#else // OPEN_FIELD
                if( fRightFPn_H > NLTE_n0 )
                    fTrt[1] *= pow((fRightFPn_H/NLTE_n0), NLTE_m_Hg);
#endif // OPEN_FIELD
                break;
                    
            case 3 :	// Paschen-alpha
                if( fLeftFPn_H > NLTE_n0 )
                    fTrt[0] *= pow((fLeftFPn_H/NLTE_n0), NLTE_m_Pa);
#ifdef OPEN_FIELD
#else // OPEN_FIELD
                if( fRightFPn_H > NLTE_n0 )
                    fTrt[1] *= pow((fRightFPn_H/NLTE_n0), NLTE_m_Pa);
#endif // OPEN_FIELD
                break;
                    
            case 4 :	// Paschen-beta
                if( fLeftFPn_H > NLTE_n0 )
                    fTrt[0] *= pow((fLeftFPn_H/NLTE_n0), NLTE_m_Pb);
#ifdef OPEN_FIELD
#else // OPEN_FIELD
                if( fRightFPn_H > NLTE_n0 )
                    fTrt[1] *= pow((fRightFPn_H/NLTE_n0), NLTE_m_Pb);
#endif // OPEN_FIELD
                break;
                    
            case 5 :	// Brackett-alpha
                if( fLeftFPn_H > NLTE_n0 )
                    fTrt[0] *= pow((fLeftFPn_H/NLTE_n0), NLTE_m_Bra);
#ifdef OPEN_FIELD
#else // OPEN_FIELD
                if( fRightFPn_H > NLTE_n0 )
                    fTrt[1] *= pow((fRightFPn_H/NLTE_n0), NLTE_m_Bra);
#endif // OPEN_FIELD
                break;
						
			default :
				break;
		}					

		// If T_rad^b has been scaled then we must update other quantities
		if( fTrt[0] != pRadiativeRates->GetScaledTrt(i,0) )
		{
			// log(0.5) + 0.5 = -0.19314718
			// h / k_B = 4.7979e-11
	        fnu = pRadiativeRates->GetNu0(i);
			fTeZ_c = ( (4.7979e-11) * fnu ) / ( -0.19314718 + ( (4.7979e-11) * (fnu/fTrt[0]) ) );
	    	if( fTeZ_c >= MINIMUM_TEMPERATURE )
			{
				pNextActiveCell = pCentreOfCurrentRow;
				pActiveCell = pNextActiveCell;
				pActiveCell->GetCellProperties( &CellProperties );
				pNextActiveCell = pActiveCell->pGetPointer( LEFT );
				while( pNextActiveCell )
				{
					pActiveCell = pNextActiveCell;
					pActiveCell->GetCellProperties( &NextCellProperties );
						
					if( ( fTeZ_c <= CellProperties.T[ELECTRON] && fTeZ_c >=NextCellProperties.T[ELECTRON] ) || ( fTeZ_c >= CellProperties.T[ELECTRON] && fTeZ_c <= NextCellProperties.T[ELECTRON] ) )
						break;
						
					CellProperties.T[ELECTRON] = NextCellProperties.T[ELECTRON];
					CellProperties.s[1] = NextCellProperties.s[1];
					CellProperties.rho_c = NextCellProperties.rho_c;
						
					pNextActiveCell = pActiveCell->pGetPointer( LEFT );
				}
				// Note that while loop must exit with break condition or x[1] = x[2] and interpolation will fail		
				x[1] = NextCellProperties.T[ELECTRON];
				x[2] = CellProperties.T[ELECTRON];
				y[1] = NextCellProperties.s[1];
				y[2] = CellProperties.s[1];
				LinearFit( x, y, fTeZ_c, &fZ_c );
				
				x[1] = NextCellProperties.rho_c;
				x[2] = CellProperties.rho_c;
				LinearFit( y, x, fZ_c, &fMcZ_c );
			
				pRadiativeRates->SetZ_c_LEFT( i, fZ_c );
				pRadiativeRates->SetMcZ_c_LEFT( i, fMcZ_c );
				pRadiativeRates->SetScaledTrt( i, fTrt[0], 0 );
			}
		}
#ifdef OPEN_FIELD
#else // OPEN_FIELD
		// If T_rad^b has been scaled then we must update other quantities
		if( fTrt[1] != pRadiativeRates->GetScaledTrt(i,1) )
		{
			// log(0.5) + 0.5 = -0.19314718
			// h / k_B = 4.7979e-11
	        fnu = pRadiativeRates->GetNu0(i);
			fTeZ_c = ( (4.7979e-11) * fnu ) / ( -0.19314718 + ( (4.7979e-11) * (fnu/fTrt[1]) ) );
	    	if( fTeZ_c >= MINIMUM_TEMPERATURE )
			{
				pNextActiveCell = pCentreOfCurrentRow;
				pActiveCell = pNextActiveCell;
				pActiveCell->GetCellProperties( &CellProperties );
				pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
				while( pNextActiveCell )
				{
					pActiveCell = pNextActiveCell;
					pActiveCell->GetCellProperties( &NextCellProperties );

					if( ( fTeZ_c <= CellProperties.T[ELECTRON] && fTeZ_c >=NextCellProperties.T[ELECTRON] ) || ( fTeZ_c >= CellProperties.T[ELECTRON] && fTeZ_c <= NextCellProperties.T[ELECTRON] ) )
						break;
							
					CellProperties.T[ELECTRON] = NextCellProperties.T[ELECTRON];
					CellProperties.s[1] = NextCellProperties.s[1];
					CellProperties.rho_c = NextCellProperties.rho_c;
						
					pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
				}
				// Note that while loop must exit with break condition or x[1] = x[2] and interpolation will fail		
				x[1] = CellProperties.T[ELECTRON];
				x[2] = NextCellProperties.T[ELECTRON];
				y[1] = CellProperties.s[1];
				y[2] = NextCellProperties.s[1];
				LinearFit( x, y, fTeZ_c, &fZ_c );
				
				x[1] = CellProperties.rho_c;
				x[2] = NextCellProperties.rho_c;
				LinearFit( y, x, fZ_c, &fMcZ_c );
			
				pRadiativeRates->SetZ_c_RIGHT( i, fZ_c );
				pRadiativeRates->SetMcZ_c_RIGHT( i, fMcZ_c );
				pRadiativeRates->SetScaledTrt( i, fTrt[1], 1 );
			}
		}
#endif // OPEN_FIELD
	}

	iMaxRL = 0;
	pNextActiveCell = pStartOfCurrentRow;
	while( pNextActiveCell )
	{
    	pActiveCell = pNextActiveCell;
#ifdef USE_JB
		if( pActiveCell->pGetPointer( LEFT ) )
		{
			old_s = CellProperties.s[1];
			old_T = CellProperties.T[ELECTRON];
		}
#endif // USE_JB
    	pActiveCell->GetCellProperties( &CellProperties );

		if( iMaxRL < CellProperties.iRefinementLevel ) iMaxRL = CellProperties.iRefinementLevel;

// *****************************************************************************
// *                                                                           *
// *    CALCULATE THE PHYSICAL QUANTITIES                                      *
// *                                                                           *
// *****************************************************************************

	    CellProperties.n[HYDROGEN] = CellProperties.rho[1] / AVERAGE_PARTICLE_MASS;
	    CellProperties.v[1] = CellProperties.rho_v[1] / CellProperties.rho[1];

	    term1 = CellProperties.TE_KE[1][HYDROGEN] / CellProperties.rho[1];
 	   	term2 = 0.5 * CellProperties.v[1] * CellProperties.v[1];
	    // GAMMA_MINUS_ONE / BOLTZMANN_CONSTANT = 4.830917874e15
	    CellProperties.T[HYDROGEN] = (4.830917874e15) * AVERAGE_PARTICLE_MASS * ( term1 - term2 );

		// Initial guess for electron density and temperature are found in the loops above this main loop when NLTE_CHROMOSPHERE is defined

		for(i=0; i<pRadiativeRates->GetNumberOfTransitions(); i++)
    	{
#ifdef OPEN_FIELD
        	if( CellProperties.s[1] <= pRadiativeRates->GetZ_c_LEFT(i) )
#else // OPEN_FIELD
        	if( CellProperties.s[1] <= pRadiativeRates->GetZ_c_LEFT(i) || CellProperties.s[1] >= pRadiativeRates->GetZ_c_RIGHT(i) )
#endif // OPEN_FIELD
	   			CellProperties.Trad[i] = CellProperties.T[ELECTRON];
			else
			{
 	        	fH = pRadiativeRates->GetH(i);
            	fnu = pRadiativeRates->GetNu0(i);
#ifdef OPEN_FIELD
				fTrt[0] = pRadiativeRates->GetScaledTrt(i,0);
            	fMcZ_c = pRadiativeRates->GetMcZ_c_LEFT(i);
#else // OPEN_FIELD
	        	if( CellProperties.s[1] <= Params.L / 2.0 )
				{
					fTrt[0] = pRadiativeRates->GetScaledTrt(i,0);
	            	fMcZ_c = pRadiativeRates->GetMcZ_c_LEFT(i);
				}
	        	else
				{
					fTrt[0] = pRadiativeRates->GetScaledTrt(i,1);
                	fMcZ_c = pRadiativeRates->GetMcZ_c_RIGHT(i);
				}
 #endif // OPEN_FIELD
 
	       		// Equation (3) in Leenaarts & Wedemeyer-Bohm, 2006, A&A, 460, 301
    	   		// Note that this equation is different to Espen Sollum's thesis (from which this treatment was taken), equation (5.3),
       			// where Te(z) is used in place of Te(Z_crit) in Sollum. Leenaarts claims the difference is small (private communication)
				term1 = exp( ( (4.7979e-11) * fnu ) / fTrt[0] ) - 1.0;
				term2 = 1.0 + pow( ( CellProperties.rho_c / fMcZ_c ), fH );
				fA = term2 / term1;
       			CellProperties.Trad[i] = ( (4.7979e-11) * fnu ) / log( (1.0/fA) + 1.0 );
			}
		}
	    
	    // Convert the radiation temperatures for each transition to log values
	    for(j=0; j<pRadiativeRates->GetNumberOfTransitions(); j++)
    	   	flog10_Trad[j] = log10( CellProperties.Trad[j]) ;

	    // CALCULATE NLTE HI FRACTION
	    pRadiativeRates->GetBoundBoundRates( fBB_lu, fBB_ul, &(flog10_Trad[0]) );
	    pRadiativeRates->GetBoundFreeRates( fBF, &(flog10_Trad[6]) );
	    pRadiativeRates->GetFreeBoundRates( fFB, &(flog10_Trad[6]), log10(CellProperties.T[ELECTRON]), CellProperties.n[ELECTRON] );
	    pRadiativeRates->GetCollisionalRatesRH( fColl_ex_lu, fColl_ex_ul, fColl_ion, fColl_rec, log10( CellProperties.T[ELECTRON] ), CellProperties.n[ELECTRON] );

#ifdef BEAM_HEATING
		fQbeam = pHeat->GetQbeam( CellProperties.s[1] );
		// Calculate the contribution to collisional excitation / ionization by non-thermal electrons
		if( fQbeam > Qbeam_CUT_OFF )
		{
			fLambda1 = 66.0 + ( 1.5 * log_fAvgEE ) - ( 0.5 * log( CellProperties.n[ELECTRON] ) );
			term1 = fQbeam * ( fLambda2 /  ( ( fLambda1 * CellProperties.n[ELECTRON] ) + ( fLambda2 * CellProperties.n[HYDROGEN] ) ) );
			fColl_ex_lu[0] += 2.94E10 * term1;
			fColl_ex_lu[1] += 5.35E9 * term1;
			fColl_ex_lu[2] += 1.91E9 * term1;
			fColl_ion[0] += 1.73E10 * term1;
		}
#endif // BEAM_HEATING

	    pRadiativeRates->SolveHIIFraction( &(CellProperties.Hstate[0]), fColl_ex_lu, fColl_ex_ul, fColl_ion, fColl_rec, fBB_lu, fBB_ul, fBF, fFB );
	    CellProperties.HI = 1.0 - CellProperties.Hstate[5];

	    // Temperature limiter
	    for( j=0; j<SPECIES; j++ )
	    {
			if( CellProperties.T[j] < MINIMUM_TEMPERATURE )
			{
            	CellProperties.T[j] = MINIMUM_TEMPERATURE;
            	// BOLTZMANN_CONSTANT / GAMMA_MINUS_ONE = 2.07e-16
            	CellProperties.TE_KE[1][j] = (2.07e-16) * CellProperties.n[j] * CellProperties.T[j];

	            if( j == HYDROGEN )
    	            CellProperties.TE_KE[1][j] += 0.5 * CellProperties.rho_v[1] * CellProperties.v[1];
        	}
    	}

	    for( j=0; j<SPECIES; j++ )
    	{
        	CellProperties.P[1][j] = BOLTZMANN_CONSTANT * CellProperties.n[j] * CellProperties.T[j];
			CellProperties.TE_KE_P[1][j] = CellProperties.TE_KE[1][j] + CellProperties.P[1][j];
    	}

	    term1 = ( GAMMA * ( CellProperties.P[1][ELECTRON] + CellProperties.P[1][HYDROGEN] ) ) / CellProperties.rho[1];
	    CellProperties.Cs = pow( term1, 0.5 );
	    CellProperties.M = fabs( CellProperties.v[1] / CellProperties.Cs );

#ifdef USE_JB
	    L_T = (0.5 * (old_T + CellProperties.T[ELECTRON])) / fabs( (CellProperties.T[ELECTRON] - old_T) / (CellProperties.s[1] - old_s) );
	    if( CellProperties.T[ELECTRON] > Te_max ) Te_max = CellProperties.T[ELECTRON];
	    if( (CellProperties.cell_width/L_T) > 0.5 && CellProperties.T[ELECTRON] > Tc ) Tc = CellProperties.T[ELECTRON];
#endif // USE_JB

// *****************************************************************************
// *                                                                           *
// *    CALCULATE THE TIMESCALES                                               *
// *                                                                           *
// *****************************************************************************

// *****************************************************************************
// *    COLLISIONS                                                             *
// *****************************************************************************

#ifdef FORCE_SINGLE_FLUID
    	CellProperties.nu_ie = 1.0 / MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE;
#else // FORCE_SINGLE_FLUID
    	// Calculate the collision frequency between the electrons and ions
    	// 4.820055089755540 * ELECTRON_MASS = 4.391070186767300e-27
    	// For derivation of 4.820055089755540 refer to research diary entry: Monday 6th February 2007 (p26)
    	CellProperties.nu_ie = ( ( (4.391070186767300e-27) / AVERAGE_PARTICLE_MASS ) * CellProperties.n[ELECTRON] * fLogLambda_ei( CellProperties.T[ELECTRON], CellProperties.T[HYDROGEN], CellProperties.n[ELECTRON] ) ) / pow( CellProperties.T[ELECTRON], 1.5 );
#endif // FORCE_SINGLE_FLUID
    	// Calculated in EvaluateTerms when the collisional term is known
    	CellProperties.collision_delta_t = Params.Duration;

// *****************************************************************************
// *    ADVECTION                                                              *
// *****************************************************************************

	    // The time-step is calculated by the CFL condition
    	CellProperties.advection_delta_t = SAFETY_ADVECTION * ( CellProperties.cell_width / ( fabs( CellProperties.v[1] ) + CellProperties.Cs ) );

// *****************************************************************************
// *    THERMAL CONDUCTION                                                     *
// *****************************************************************************

	    // Calculated in EvaluateTerms when the thermal conduction terms are known
    	for( j=0; j<SPECIES; j++ )
			CellProperties.conduction_delta_t[j] = Params.Duration;

// *****************************************************************************
// *    RADIATION                                                              *
// *****************************************************************************

	    // Calculated in EvaluateTerms when the radiative term is known
    	CellProperties.radiation_delta_t = Params.Duration;

// *****************************************************************************
// *    DYNAMIC VISCOSITY                                                      *
// *****************************************************************************

	    // Calculated in EvaluateTerms when the dynamic viscosity term is known
    	CellProperties.viscosity_delta_t = Params.Duration;

	    pActiveCell->UpdateCellProperties( &CellProperties );

    	pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
	}

#ifdef USE_JB
	// Prevent the critical temperature for the TRAC method decreasing by more than 10% over any 1 second interval
	// (See equation A.4 in Johnston et al. 2020)
	if( old_Tc && Tc < old_Tc )
	{
		term1 = TIME_STEP_LIMIT / tau_JB;
		term2 = pow( 0.1, term1 );
		Tc = old_Tc * term2;
	}
	old_Tc = Tc;
#endif // USE_JB
}
#endif // OPTICALLY_THICK_RADIATION
#endif // NLTE_CHROMOSPHERE

void CEquations::CalculatePhysicalQuantities( void )
{
#ifdef OPENMP
	#ifdef OPTICALLY_THICK_RADIATION
	#ifdef NLTE_CHROMOSPHERE
		PCELL pNextActiveCell;
	#endif // NLTE_CHROMOSPHERE
	#endif // OPTICALLY_THICK_RADIATION
	#ifdef USE_JB
		double localTe_max = 0.0, localTc = 0.0;
	#endif // USE_JB
#else // OPENMP
	PCELL pNextActiveCell;
	#ifdef USE_JB
		Te_max = 0.0;
		Tc = 0.0;
	#endif // USE_JB
#endif // OPENMP
CELLPROPERTIES CellProperties;

#ifdef USE_JB
	double old_s = 0.0, old_T = 0.0, tau_JB = 1.0, L_T;
#endif // USE_JB

// General variables
double term1, term2;
int j;

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
	CELLPROPERTIES NextCellProperties;
	double frho_c, fLeftFPn_H;
#ifdef OPEN_FIELD
#else // OPEN_FIELD
	double fRightFPn_H;
#endif // OPEN_FIELD
	double fH, fnu, fTrt[2], fMcZ_c, fA, fZ[31];
	double fTeZ_c, fZ_c;
	double fElement, fSum, flog10_Trad[10], fPreviousIteration;
	double fBB_lu[6], fBB_ul[6], fBF[4], fFB[4], fColl_ex_lu[10], fColl_ex_ul[10], fColl_ion[5], fColl_rec[5];
	double x[3], y[3];
	int iNumElements, *piA;
	int i;
#ifdef BEAM_HEATING
	double fQbeam, fLambda1, fLambda2, log_fAvgEE;
	log_fAvgEE = log( pHeat->GetAvgEE() );
	fLambda2 = 25.1 + log_fAvgEE;
#endif // BEAM_HEATING

	memset( &CellProperties, 0, sizeof(CELLPROPERTIES) );	// Avoids a warning that variables might be used undefined

	// We need to make T_b^top a time-dependent function, which means finding the foot-point densities at each leg of the loop prior to executing the
	// main loop. The foot-point density is defined as the density at a particular temperature (approx. the H_alpha formation temperature of 9,000 K)
	frho_c = 0.0;
	fLeftFPn_H = 0.0;
	pNextActiveCell = pCentreOfCurrentRow;
	while( pNextActiveCell )
	{
		pActiveCell = pNextActiveCell;
		pActiveCell->GetCellProperties( &CellProperties );
		// Calculate the electron density and temperature
		CellProperties.n[ELECTRON] = CellProperties.rho_e / ELECTRON_MASS;
	    // GAMMA_MINUS_ONE / BOLTZMANN_CONSTANT = 4.830917874e15
		CellProperties.T[ELECTRON] = (4.830917874e15) * ( CellProperties.TE_KE[1][ELECTRON] / CellProperties.n[ELECTRON] );
		// Calculate the column mass density
		frho_c += CellProperties.rho[1] * CellProperties.cell_width;
		CellProperties.rho_c = frho_c;
		// ( GAMMA_MINUS_ONE * ELECTRON_MASS ) / BOLTZMANN_CONSTANT = 4.40066838247E-12
		if( ( (4.40066838247e-12) * ( CellProperties.TE_KE[1][ELECTRON] / CellProperties.rho_e ) ) <= NLTE_T_FP && !fLeftFPn_H  )
			fLeftFPn_H = CellProperties.rho[1] / AVERAGE_PARTICLE_MASS;
	
		pActiveCell->UpdateCellProperties( &CellProperties );
		pNextActiveCell = pActiveCell->pGetPointer( LEFT );
	}
#ifdef OPEN_FIELD
#else // OPEN_FIELD
	frho_c = 0.0;
	fRightFPn_H = 0.0;
	pNextActiveCell = pCentreOfCurrentRow;
	while( pNextActiveCell )
	{
		pActiveCell = pNextActiveCell;
		pActiveCell->GetCellProperties( &CellProperties );
		// Calculate the electron density and temperature
    	CellProperties.n[ELECTRON] = CellProperties.rho_e / ELECTRON_MASS;
		// GAMMA_MINUS_ONE / BOLTZMANN_CONSTANT = 4.830917874e15
    	CellProperties.T[ELECTRON] = (4.830917874e15) * ( CellProperties.TE_KE[1][ELECTRON] / CellProperties.n[ELECTRON] );
		// Calculate the column mass density
		frho_c += CellProperties.rho[1] * CellProperties.cell_width;
		CellProperties.rho_c = frho_c;
		// ( GAMMA_MINUS_ONE * ELECTRON_MASS ) / BOLTZMANN_CONSTANT = 4.40066838247E-12
		if( ( (4.40066838247e-12) * ( CellProperties.TE_KE[1][ELECTRON] / CellProperties.rho_e ) ) <= NLTE_T_FP && !fRightFPn_H )
			fRightFPn_H = CellProperties.rho[1] / AVERAGE_PARTICLE_MASS;

		pActiveCell->UpdateCellProperties( &CellProperties );
		pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
	}
#endif // OPEN_FIELD

    for(i=0; i<pRadiativeRates->GetNumberOfTransitions(); i++)
   	{
		fTrt[0] = fTrt[1] = pRadiativeRates->GetTrt(i);
		// Identify the current transition
		switch( i ) {
						
			case 0 :	// H-alpha
				if( fLeftFPn_H > NLTE_n0 )
                    fTrt[0] *= pow((fLeftFPn_H/NLTE_n0), NLTE_m_Ha);
#ifdef OPEN_FIELD
#else // OPEN_FIELD
				if( fRightFPn_H > NLTE_n0 )
					fTrt[1] *= pow((fRightFPn_H/NLTE_n0), NLTE_m_Ha);
#endif // OPEN_FIELD
			break;
						
			case 1 :	// H-beta
                if( fLeftFPn_H > NLTE_n0 )
                    fTrt[0] *= pow((fLeftFPn_H/NLTE_n0), NLTE_m_Hb);
#ifdef OPEN_FIELD
#else // OPEN_FIELD
                if( fRightFPn_H > NLTE_n0 )
                    fTrt[1] *= pow((fRightFPn_H/NLTE_n0), NLTE_m_Hb);
#endif // OPEN_FIELD
			break;
						
			case 2 :	// H-gamma
                if( fLeftFPn_H > NLTE_n0 )
                    fTrt[0] *= pow((fLeftFPn_H/NLTE_n0), NLTE_m_Hg);
#ifdef OPEN_FIELD
#else // OPEN_FIELD
                if( fRightFPn_H > NLTE_n0 )
                    fTrt[1] *= pow((fRightFPn_H/NLTE_n0), NLTE_m_Hg);
#endif // OPEN_FIELD
			break;
						
			case 3 :	// Paschen-alpha
                if( fLeftFPn_H > NLTE_n0 )
                    fTrt[0] *= pow((fLeftFPn_H/NLTE_n0), NLTE_m_Pa);
#ifdef OPEN_FIELD
#else // OPEN_FIELD
                if( fRightFPn_H > NLTE_n0 )
                    fTrt[1] *= pow((fRightFPn_H/NLTE_n0), NLTE_m_Pa);
#endif // OPEN_FIELD
			break;
						
			case 4 :	// Paschen-beta
                if( fLeftFPn_H > NLTE_n0 )
                    fTrt[0] *= pow((fLeftFPn_H/NLTE_n0), NLTE_m_Pb);
#ifdef OPEN_FIELD
#else // OPEN_FIELD
                if( fRightFPn_H > NLTE_n0 )
                    fTrt[1] *= pow((fRightFPn_H/NLTE_n0), NLTE_m_Pb);
#endif // OPEN_FIELD
			break;
                
            case 5 :	// Brackett-alpha
                if( fLeftFPn_H > NLTE_n0 )
                    fTrt[0] *= pow((fLeftFPn_H/NLTE_n0), NLTE_m_Bra);
#ifdef OPEN_FIELD
#else // OPEN_FIELD
                if( fRightFPn_H > NLTE_n0 )
                    fTrt[1] *= pow((fRightFPn_H/NLTE_n0), NLTE_m_Bra);
#endif // OPEN_FIELD
            break;
						
			default :
			break;
		}					

		// If T_rad^b has been scaled then we must update other quantities
		if( fTrt[0] != pRadiativeRates->GetScaledTrt(i,0) )
		{
			// log(0.5) + 0.5 = -0.19314718
			// h / k_B = 4.7979e-11
            fnu = pRadiativeRates->GetNu0(i);
			fTeZ_c = ( (4.7979e-11) * fnu ) / ( -0.19314718 + ( (4.7979e-11) * (fnu/fTrt[0]) ) );
      		if( fTeZ_c >= MINIMUM_TEMPERATURE )
			{
				pNextActiveCell = pCentreOfCurrentRow;
				pActiveCell = pNextActiveCell;
				pActiveCell->GetCellProperties( &CellProperties );
				pNextActiveCell = pActiveCell->pGetPointer( LEFT );
				while( pNextActiveCell )
				{
					pActiveCell = pNextActiveCell;
					pActiveCell->GetCellProperties( &NextCellProperties );
						
					if( ( fTeZ_c <= CellProperties.T[ELECTRON] && fTeZ_c >=NextCellProperties.T[ELECTRON] ) || ( fTeZ_c >= CellProperties.T[ELECTRON] && fTeZ_c <= NextCellProperties.T[ELECTRON] ) )
						break;

					CellProperties.T[ELECTRON] = NextCellProperties.T[ELECTRON];
					CellProperties.s[1] = NextCellProperties.s[1];
					CellProperties.rho_c = NextCellProperties.rho_c;
					
					pNextActiveCell = pActiveCell->pGetPointer( LEFT );
				}
				// Note that while loop must exit with break condition or x[1] = x[2] and interpolation will fail		
				x[1] = NextCellProperties.T[ELECTRON];
				x[2] = CellProperties.T[ELECTRON];
				y[1] = NextCellProperties.s[1];
				y[2] = CellProperties.s[1];
				LinearFit( x, y, fTeZ_c, &fZ_c );
				
				x[1] = NextCellProperties.rho_c;
				x[2] = CellProperties.rho_c;
				LinearFit( y, x, fZ_c, &fMcZ_c );
			
				pRadiativeRates->SetZ_c_LEFT( i, fZ_c );
				pRadiativeRates->SetMcZ_c_LEFT( i, fMcZ_c );
				pRadiativeRates->SetScaledTrt( i, fTrt[0], 0 );
			}
		}
#ifdef OPEN_FIELD
#else // OPEN_FIELD
		// If T_rad^b has been scaled then we must update other quantities
		if( fTrt[1] != pRadiativeRates->GetScaledTrt(i,1) )
		{
			// log(0.5) + 0.5 = -0.19314718
			// h / k_B = 4.7979e-11
	        fnu = pRadiativeRates->GetNu0(i);
			fTeZ_c = ( (4.7979e-11) * fnu ) / ( -0.19314718 + ( (4.7979e-11) * (fnu/fTrt[1]) ) );
	   		if( fTeZ_c >= MINIMUM_TEMPERATURE )
			{
				pNextActiveCell = pCentreOfCurrentRow;
				pActiveCell = pNextActiveCell;
				pActiveCell->GetCellProperties( &CellProperties );
				pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
				while( pNextActiveCell )
				{
					pActiveCell = pNextActiveCell;
					pActiveCell->GetCellProperties( &NextCellProperties );

					if( ( fTeZ_c <= CellProperties.T[ELECTRON] && fTeZ_c >=NextCellProperties.T[ELECTRON] ) || ( fTeZ_c >= CellProperties.T[ELECTRON] && fTeZ_c <= NextCellProperties.T[ELECTRON] ) )
						break;

					CellProperties.T[ELECTRON] = NextCellProperties.T[ELECTRON];
					CellProperties.s[1] = NextCellProperties.s[1];
					CellProperties.rho_c = NextCellProperties.rho_c;
						
					pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
				}
				// Note that while loop must exit with break condition or x[1] = x[2] and interpolation will fail		
				x[1] = CellProperties.T[ELECTRON];
				x[2] = NextCellProperties.T[ELECTRON];
				y[1] = CellProperties.s[1];
				y[2] = NextCellProperties.s[1];
				LinearFit( x, y, fTeZ_c, &fZ_c );
				
				x[1] = CellProperties.rho_c;
				x[2] = NextCellProperties.rho_c;
				LinearFit( y, x, fZ_c, &fMcZ_c );
			
				pRadiativeRates->SetZ_c_RIGHT( i, fZ_c );
				pRadiativeRates->SetMcZ_c_RIGHT( i, fMcZ_c );
				pRadiativeRates->SetScaledTrt( i, fTrt[1], 1 );
			}
		}
#endif // OPEN_FIELD
	}
#ifdef NON_EQUILIBRIUM_RADIATION
	double *pfZ;
#endif // NON_EQUILIBRIUM_RADIATION
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

#ifdef OPENMP
	PCELL pLocalActiveCell;
	int iCounter, iLocalMaxRL = 0, iLocalNumberOfCells = Params.iNumberOfCells;
	#ifdef OPTICALLY_THICK_RADIATION
	#ifdef NLTE_CHROMOSPHERE
		#define LTE_VARS_1 iNumElements, piA, fZ, fElement, fSum,
		#ifdef NON_EQUILIBRIUM_RADIATION
			#define NLTE_VARS_1 pfZ,
		#else // NON_EQUILIBRIUM_RADIATION
			#define NLTE_VARS_1
		#endif // NON_EQUILIBRIUM_RADIATION
		#define NLTE_CHR_VARS_1 fH, fnu, fTrt, fMcZ_c, fA, flog10_Trad, fPreviousIteration, fBB_lu, fBB_ul, fBF, fFB, fColl_ex_lu, fColl_ex_ul, fColl_ion, fColl_rec, i,
		#ifdef BEAM_HEATING
			#define BEAM_HEATING_VARS_1 fQbeam, fLambda1,
		#else // BEAM_HEATING
			#define BEAM_HEATING_VARS_1
		#endif // BEAM_HEATING
	#else // NLTE_CHROMOSPHERE
		#define LTE_VARS_1
		#define NLTE_VARS_1
		#define NLTE_CHR_VARS_1
		#define BEAM_HEATING_VARS_1
	#endif //NLTE_CHROMOSPHERE
	#else // OPTICALLY_THICK_RADIATION
		#define LTE_VARS_1
		#define NLTE_VARS_1
		#define NLTE_CHR_VARS_1
		#define BEAM_HEATING_VARS_1
	#endif // OPTICALLY_THICK_RADIATION

	#ifdef USE_JB
		#define USE_JB_VARS_1 old_s, old_T, L_T,
		#define REDUCTION_VARS reduction(max: iLocalMaxRL), reduction(max: localTe_max), reduction(max: localTc)
	#else // USE_JB
		#define USE_JB_VARS_1
		#define REDUCTION_VARS reduction(max: iLocalMaxRL)
	#endif // USE_JB

#pragma omp parallel shared( ppCellList, iLocalNumberOfCells )	\
					private(	LTE_VARS_1						\
							 	NLTE_VARS_1						\
							 	NLTE_CHR_VARS_1					\
							 	BEAM_HEATING_VARS_1				\
								USE_JB_VARS_1					\
								iCounter,						\
								CellProperties,					\
								j, term1, term2 )
	{
#pragma omp for schedule(dynamic, CHUNK_SIZE)		\
    	            lastprivate(pLocalActiveCell)	\
        	        REDUCTION_VARS
	for( iCounter=0; iCounter<iLocalNumberOfCells; iCounter++ )
	{
		pLocalActiveCell = ppCellList[iCounter];
#ifdef USE_JB
		if( pLocalActiveCell->pGetPointer( LEFT ) )
		{
			old_s = CellProperties.s[1];
			old_T = CellProperties.T[ELECTRON];
		}
#endif // USE_JB
		pLocalActiveCell->GetCellProperties( &CellProperties );
		if( iLocalMaxRL < CellProperties.iRefinementLevel ) iLocalMaxRL = CellProperties.iRefinementLevel;
#else // OPENMP
	iMaxRL = 0;
	pNextActiveCell = pStartOfCurrentRow;
	while( pNextActiveCell )
	{
		pActiveCell = pNextActiveCell;
#ifdef USE_JB
		if( pActiveCell->pGetPointer( LEFT ) )
		{
			old_s = CellProperties.s[1];
			old_T = CellProperties.T[ELECTRON];
		}
#endif // USE_JB
   		pActiveCell->GetCellProperties( &CellProperties );
		if( iMaxRL < CellProperties.iRefinementLevel ) iMaxRL = CellProperties.iRefinementLevel;
#endif // OPENMP

// *****************************************************************************
// *                                                                           *
// *    CALCULATE THE PHYSICAL QUANTITIES                                      *
// *                                                                           *
// *****************************************************************************

	    CellProperties.n[HYDROGEN] = CellProperties.rho[1] / AVERAGE_PARTICLE_MASS;
    	CellProperties.v[1] = CellProperties.rho_v[1] / CellProperties.rho[1];

	    term1 = CellProperties.TE_KE[1][HYDROGEN] / CellProperties.rho[1];
    	term2 = 0.5 * CellProperties.v[1] * CellProperties.v[1];
    	// GAMMA_MINUS_ONE / BOLTZMANN_CONSTANT = 4.830917874e15
    	CellProperties.T[HYDROGEN] = (4.830917874e15) * AVERAGE_PARTICLE_MASS * ( term1 - term2 );

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE

///////////////////////////////////////////////////////////////////
// Calculate the electron density for the NLTE chromosphere   	 //
///////////////////////////////////////////////////////////////////
	
		// Initial guess for electron density and temperature are found in the loops above this main loop when NLTE_CHROMOSPHERE is defined

	    for(i=0; i<pRadiativeRates->GetNumberOfTransitions(); i++)
    	{
#ifdef OPEN_FIELD
	        if( CellProperties.s[1] <= pRadiativeRates->GetZ_c_LEFT(i) )
#else // OPEN_FIELD
	        if( CellProperties.s[1] <= pRadiativeRates->GetZ_c_LEFT(i) || CellProperties.s[1] >= pRadiativeRates->GetZ_c_RIGHT(i) )
#endif // OPEN_FIELD
		    	CellProperties.Trad[i] = CellProperties.T[ELECTRON];
			else
			{
 	        	fH = pRadiativeRates->GetH(i);
            	fnu = pRadiativeRates->GetNu0(i);
#ifdef OPEN_FIELD
				fTrt[0] = pRadiativeRates->GetScaledTrt(i,0);
    	        fMcZ_c = pRadiativeRates->GetMcZ_c_LEFT(i);
#else // OPEN_FIELD
		        if( CellProperties.s[1] <= Params.L / 2.0 )
				{
					fTrt[0] = pRadiativeRates->GetScaledTrt(i,0);
	        	    fMcZ_c = pRadiativeRates->GetMcZ_c_LEFT(i);
				}
	        	else
				{
					fTrt[0] = pRadiativeRates->GetScaledTrt(i,1);
    	            fMcZ_c = pRadiativeRates->GetMcZ_c_RIGHT(i);
				}
 #endif // OPEN_FIELD

	       		// Equation (3) in Leenaarts & Wedemeyer-Bohm, 2006, A&A, 460, 301
   	    		// Note that this equation is different to Espen Sollum's thesis (from which this treatment was taken), equation (5.3),
   	    		// where Te(z) is used in place of Te(Z_crit) in Sollum. Leenaarts claims the difference is small (private communication)
				term1 = exp( ( (4.7979e-11) * fnu ) / fTrt[0] ) - 1.0;
				term2 = 1.0 + pow( ( CellProperties.rho_c / fMcZ_c ), fH );
				fA = term2 / term1;
       			CellProperties.Trad[i] = ( (4.7979e-11) * fnu ) / log( (1.0/fA) + 1.0 );
			}
		}

	    // Calculate the contribution to the electron density from elements other than hydrogen:
 	   	fSum = 1.44e-4;		// Ensure a minimum number of electrons in the coldest part of the atmosphere

#ifdef NON_EQUILIBRIUM_RADIATION
    	piA = CellProperties.pIonFrac->pGetElementInfo( &iNumElements );
    	for( i=0; i<iNumElements; i++ )
    	{
			// Don't double count hydrogen
			if( piA[i] > 1 )
			{
    			fElement = 0.0;

	            pfZ = CellProperties.pIonFrac->pGetIonFrac( piA[i] );
 	    	    for( j=1; j<piA[i]+1; j++ )
					fElement += ((double)j) * pfZ[j];

				fElement *= pRadiation->GetAbundance( piA[i] );
				fSum += fElement;
			}		
    	}
#endif // NON_EQUILIBRIUM_RADIATION

	    piA = pRadiation2->pGetAtomicNumbers( &iNumElements );
	    for( i=0; i<iNumElements; i++ )
    	{
			// Don't double count hydrogen
			if( piA[i] > 1 )
			{
    			fElement = 0.0;
			
	    		pRadiation2->GetEquilIonFrac( piA[i], fZ, log10(CellProperties.T[ELECTRON]) );
    			for( j=1; j<piA[i]+1; j++ )
					fElement += ((double)j) * fZ[j];

            	fElement *= pRadiation2->GetAbundance( piA[i] );
            	fSum += fElement;
			}
    	}

	    // Convert the radiation temperatures for each transition to log values
    	for(j=0; j<pRadiativeRates->GetNumberOfTransitions(); j++)
       		flog10_Trad[j] = log10( CellProperties.Trad[j]) ;

#ifdef BEAM_HEATING
		fQbeam = pHeat->GetQbeam( CellProperties.s[1] );
#endif // BEAM_HEATING
    	for(i=0; i<=MAX_ITERATIONS; i++ )
    	{
			fPreviousIteration = CellProperties.n[ELECTRON];

			// CALCULATE NLTE HI FRACTION
    	   	pRadiativeRates->GetBoundBoundRates( fBB_lu, fBB_ul, &(flog10_Trad[0]) );
        	pRadiativeRates->GetBoundFreeRates( fBF, &(flog10_Trad[6]) );
       		pRadiativeRates->GetFreeBoundRates( fFB, &(flog10_Trad[6]), log10(CellProperties.T[ELECTRON]), CellProperties.n[ELECTRON] );
			pRadiativeRates->GetCollisionalRatesRH( fColl_ex_lu, fColl_ex_ul, fColl_ion, fColl_rec, log10( CellProperties.T[ELECTRON] ), CellProperties.n[ELECTRON] );

#ifdef BEAM_HEATING
			// Calculate the contribution to collisional excitation / ionization by non-thermal electrons
			if( fQbeam > Qbeam_CUT_OFF )
			{
				fLambda1 = 66.0 + ( 1.5 * log_fAvgEE ) - ( 0.5 * log( CellProperties.n[ELECTRON] ) );
				term1 = fQbeam * ( fLambda2 /  ( ( fLambda1 * CellProperties.n[ELECTRON] ) + ( fLambda2 * CellProperties.n[HYDROGEN] ) ) );
				fColl_ex_lu[0] += 2.94E10 * term1;
				fColl_ex_lu[1] += 5.35E9 * term1;
				fColl_ex_lu[2] += 1.91E9 * term1;
				fColl_ion[0] += 1.73E10 * term1;
			}
#endif // BEAM_HEATING

			pRadiativeRates->SolveHIIFraction( &(CellProperties.Hstate[0]), fColl_ex_lu, fColl_ex_ul, fColl_ion, fColl_rec, fBB_lu, fBB_ul, fBF, fFB );
			CellProperties.HI = 1.0 - CellProperties.Hstate[5];

			// Estimate the new density and temperature
			CellProperties.n[ELECTRON] = CellProperties.n[ELECTRON] - CONVERGENCE_EPSILON * ( CellProperties.n[ELECTRON] - ( CellProperties.n[HYDROGEN] * ( CellProperties.Hstate[5] + fSum ) ) );
				// If the electron density changes precipitously then prevent dramatic changes in the electron temperature by recalculating the electron energy
				CellProperties.TE_KE[1][ELECTRON] = (2.07e-16) * CellProperties.n[ELECTRON] * CellProperties.T[ELECTRON];
			CellProperties.T[ELECTRON] =  CellProperties.T[ELECTRON] - CONVERGENCE_EPSILON * ( CellProperties.T[ELECTRON] - ( ( (4.830917874e15) * CellProperties.TE_KE[1][ELECTRON] ) / CellProperties.n[ELECTRON] ) );
			// Check for convergence and exit the loop if the condition is met
    	    if( i && (fabs(fPreviousIteration-CellProperties.n[ELECTRON])/CellProperties.n[ELECTRON]) < CONVERGENCE_CONDITION ) break;
		}

	    CellProperties.rho_e = CellProperties.n[ELECTRON] * ELECTRON_MASS;
    	// BOLTZMANN_CONSTANT / GAMMA_MINUS_ONE = 2.07e-16
		CellProperties.TE_KE[1][ELECTRON] = (2.07e-16) * CellProperties.n[ELECTRON] * CellProperties.T[ELECTRON];
#else // NLTE_CHROMOSPHERE

/////////////////////////////////////////////////////////////////
// Calculate the electron density for the LTE chromosphere     //
/////////////////////////////////////////////////////////////////

	    // NOTE: We can't use the non-equilibrium population here because it's based on a completely different set of rates.
    	//       If you try to use it (or the RadiationModel object to get an equilibrium population) for CellProperties.HI
    	//	     then you end up with a completely different value for the neutral hydrogen population at T[HYDROGEN] than
    	//       pHI->GetIonFrac returns and, consequently, a different electron density. The non-equilibrium populations quickly
    	//	     evolves to a new equilibrium in the dense lower atmosphere. The different electron density leads to a different
    	//	     electron temperature (very different to the hydrogen temperature) in the lower chromosphere, particularly at the 
    	//	     boundaries which don't get updated except by refinement and quickly propagates into the domain, messing
    	//	     everything up quite spectacularly!
    	//	     The bottom line is: USE CONSISTENT RATES AND ION POPULATIONS EVERYWHERE
    
		if( CellProperties.T[HYDROGEN] < OPTICALLY_THICK_TEMPERATURE )
			CellProperties.HI = pHI->GetIonFrac( log10(CellProperties.T[HYDROGEN]) );
    	else
			CellProperties.HI = 0.0;
    
		// 1.44e-4 provides some negligible number of electrons to avoid division by zero errors
 		// 1.000144 = 1.0 + 1.44e-4
   		CellProperties.n[ELECTRON] = ( 1.000144 - CellProperties.HI ) * CellProperties.n[HYDROGEN];
    	// GAMMA_MINUS_ONE / BOLTZMANN_CONSTANT = 4.830917874e15
    	CellProperties.T[ELECTRON] = ( (4.830917874e15) * CellProperties.TE_KE[1][ELECTRON] ) / CellProperties.n[ELECTRON];

#endif // NLTE_CHROMOSPHERE
#else // OPTICALLY_THICK_RADIATION
        
//////////////////////////////////////
// Calculate the electron density   //
//////////////////////////////////////

	    CellProperties.n[ELECTRON] = CellProperties.n[HYDROGEN];
	    // GAMMA_MINUS_ONE / BOLTZMANN_CONSTANT = 4.830917874e15
    	CellProperties.T[ELECTRON] = ( (4.830917874e15) * CellProperties.TE_KE[1][ELECTRON] ) / CellProperties.n[ELECTRON];

#ifdef BEAM_HEATING
    	CellProperties.HI = 0.0;
#endif // BEAM_HEATING

#endif // OPTICALLY_THICK_RADIATION

	    // Temperature limiter
	    for( j=0; j<SPECIES; j++ )
    	{
        	if( CellProperties.T[j] < MINIMUM_TEMPERATURE )
			{
            	CellProperties.T[j] = MINIMUM_TEMPERATURE;
            	// BOLTZMANN_CONSTANT / GAMMA_MINUS_ONE = 2.07e-16
            	CellProperties.TE_KE[1][j] = (2.07e-16) * CellProperties.n[j] * CellProperties.T[j];

            	if( j == HYDROGEN )
                	CellProperties.TE_KE[1][j] += 0.5 * CellProperties.rho_v[1] * CellProperties.v[1];
        	}
    	}

	    for( j=0; j<SPECIES; j++ )
    	{
        	CellProperties.P[1][j] = BOLTZMANN_CONSTANT * CellProperties.n[j] * CellProperties.T[j];
			CellProperties.TE_KE_P[1][j] = CellProperties.TE_KE[1][j] + CellProperties.P[1][j];
    	}

	    term1 = ( GAMMA * ( CellProperties.P[1][ELECTRON] + CellProperties.P[1][HYDROGEN] ) ) / CellProperties.rho[1];
    	CellProperties.Cs = pow( term1, 0.5 );
    	CellProperties.M = fabs( CellProperties.v[1] / CellProperties.Cs );

#ifdef USE_JB
	   	L_T = (0.5 * (old_T + CellProperties.T[ELECTRON])) / fabs( (CellProperties.T[ELECTRON] - old_T) / (CellProperties.s[1] - old_s) );
		#ifdef OPENMP
	    	if( CellProperties.T[ELECTRON] > localTe_max ) localTe_max = CellProperties.T[ELECTRON];
	    	if( (CellProperties.cell_width/L_T) > 0.5 && CellProperties.T[ELECTRON] > localTc ) localTc = CellProperties.T[ELECTRON];
		#else // OPENMP
	    	if( CellProperties.T[ELECTRON] > Te_max ) Te_max = CellProperties.T[ELECTRON];
	    	if( (CellProperties.cell_width/L_T) > 0.5 && CellProperties.T[ELECTRON] > Tc ) Tc = CellProperties.T[ELECTRON];
		#endif // OPENMP
#endif // USE_JB

// *****************************************************************************
// *                                                                           *
// *    CALCULATE THE TIMESCALES                                               *
// *                                                                           *
// *****************************************************************************

// *****************************************************************************
// *    COLLISIONS                                                             *
// *****************************************************************************

#ifdef FORCE_SINGLE_FLUID
    	CellProperties.nu_ie = 1.0 / MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE;
#else // FORCE_SINGLE_FLUID
    	// Calculate the collision frequency between the electrons and ions
    	// 4.820055089755540 * ELECTRON_MASS = 4.391070186767300e-27
    	// For derivation of 4.820055089755540 refer to research diary entry: Monday 6th February 2007 (p26)
    	CellProperties.nu_ie = ( ( (4.391070186767300e-27) / AVERAGE_PARTICLE_MASS ) * CellProperties.n[ELECTRON] * fLogLambda_ei( CellProperties.T[ELECTRON], CellProperties.T[HYDROGEN], CellProperties.n[ELECTRON] ) ) / pow( CellProperties.T[ELECTRON], 1.5 );
#endif // FORCE_SINGLE_FLUID
    	// Calculated in EvaluateTerms when the collisional term is known
    	CellProperties.collision_delta_t = Params.Duration;

// *****************************************************************************
// *    ADVECTION                                                              *
// *****************************************************************************

	    // The time-step is calculated by the CFL condition
		CellProperties.advection_delta_t = SAFETY_ADVECTION * ( CellProperties.cell_width / ( fabs( CellProperties.v[1] ) + CellProperties.Cs ) );

// *****************************************************************************
// *    THERMAL CONDUCTION                                                     *
// *****************************************************************************

	    // Calculated in EvaluateTerms when the thermal conduction terms are known
	    for( j=0; j<SPECIES; j++ )
			CellProperties.conduction_delta_t[j] = Params.Duration;

// *****************************************************************************
// *    RADIATION                                                              *
// *****************************************************************************

	    // Calculated in EvaluateTerms when the radiative term is known
	    CellProperties.radiation_delta_t = Params.Duration;

// *****************************************************************************
// *    DYNAMIC VISCOSITY                                                      *
// *****************************************************************************

	    // Calculated in EvaluateTerms when the dynamic viscosity term is known
	    CellProperties.viscosity_delta_t = Params.Duration;

#ifdef OPENMP
		pLocalActiveCell->UpdateCellProperties( &CellProperties );
	}
#else // OPENMP
   	pActiveCell->UpdateCellProperties( &CellProperties );
	pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
#endif // OPENMP
	}

#ifdef OPENMP
	iMaxRL = iLocalMaxRL;
#endif // OPENMP
	
#ifdef USE_JB
	#ifdef OPENMP
		Te_max = localTe_max;
		Tc = localTc;
	#endif // OPENMP
	// Prevent the critical temperature for the TRAC method decreasing by more than 10% over any 1 second interval
	// (See equation A.4 in Johnston et al. 2020)
	if( old_Tc && Tc < old_Tc )
	{
		term1 = TIME_STEP_LIMIT / tau_JB;
		term2 = pow( 0.1, term1 );
		Tc = old_Tc * term2;
	}
	old_Tc = Tc;
#endif // USE_JB
}

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
void CEquations::InitialiseRadiativeRates( void )
{
FILE *pVALAtmosphere;
PCELL pNextActiveCell, pLeftActiveCell, pLeftCell;
CELLPROPERTIES CellProperties, LeftCellProperties;
double **ppVALAtmosphere, TeZ_c, fZ_c, fMcZ_c;
double x[3], y[3], fMcZ_c_LEFT;
int iVALDataPoints;
int i, j;
#ifdef OPEN_FIELD
#else // OPEN_FIELD
	PCELL pRightActiveCell, pRightCell;
	CELLPROPERTIES RightCellProperties;
	double fMcZ_c_RIGHT;
#endif // OPEN_FIELD

	// Open and read the tabulated VAL atmosphere
	pVALAtmosphere = fopen( "Radiation_Model/atomic_data/OpticallyThick/VAL_atmospheres/VAL.T", "r" );
		fscanf( pVALAtmosphere, "%i", &iVALDataPoints );
		// Allocate memory to the tabulated VAL atmosphere
		ppVALAtmosphere = (double**)malloc( sizeof(double*) * 2 );
		ppVALAtmosphere[0] = (double*)malloc( sizeof(double) * iVALDataPoints );
		ppVALAtmosphere[1] = (double*)malloc( sizeof(double) * iVALDataPoints );
		for( j=0; j<iVALDataPoints; j++ )
		{
			ReadDouble( pVALAtmosphere, &(ppVALAtmosphere[0][j]) );
			ReadDouble( pVALAtmosphere, &(ppVALAtmosphere[1][j]) );
		}
	fclose( pVALAtmosphere );

	for(i=0; i<pRadiativeRates->GetNumberOfTransitions(); i++)
	{
		TeZ_c = pRadiativeRates->GetTeZ_c( i );
            
///////////////////////////////////
// STEP 1:                       //
// Find Z_c for each transition  //
///////////////////////////////////

		for( j=0; j<iVALDataPoints-1; j++ )
			if( ( TeZ_c <= ppVALAtmosphere[1][j] && TeZ_c >= ppVALAtmosphere[1][j+1] ) || ( TeZ_c >= ppVALAtmosphere[1][j] && TeZ_c <= ppVALAtmosphere[1][j+1] ) )
        	    break;
		
		x[1] = ppVALAtmosphere[1][j];
		x[2] = ppVALAtmosphere[1][j+1];
		y[1] = ppVALAtmosphere[0][j];
		y[2] = ppVALAtmosphere[0][j+1];
		LinearFit( x, y, TeZ_c, &fZ_c );

		pRadiativeRates->SetZ_c_LEFT( i, fZ_c );
		pRadiativeRates->SetZ_c_RIGHT( i, Params.L - fZ_c );

////////////////////////////////////////
// STEP 2:                            //
// Find McZ_c for each transition     //
////////////////////////////////////////

		// Find the centre of the grid (corresponds approximately to the apex of the loop)
		pNextActiveCell = pCentreOfCurrentRow;
		pActiveCell = pNextActiveCell;
		pActiveCell->GetCellProperties( &CellProperties );

		CellProperties.rho_c = CellProperties.rho[1] * CellProperties.cell_width;
		pActiveCell->UpdateCellProperties( &CellProperties );

		fMcZ_c_LEFT = CellProperties.rho_c;
#ifdef OPEN_FIELD
#else // OPEN_FIELD
		fMcZ_c_RIGHT = fMcZ_c_LEFT;
#endif // OPEN_FIELD

		// Left-hand chromosphere
		pLeftActiveCell = pActiveCell;
		for( ;; )
		{
			pLeftCell = pLeftActiveCell->pGetPointer( LEFT );
			pLeftCell->GetCellProperties( &LeftCellProperties );

			if( LeftCellProperties.s[1] < fZ_c ) break;

			fMcZ_c_LEFT += LeftCellProperties.rho[1] * LeftCellProperties.cell_width;
			LeftCellProperties.rho_c = fMcZ_c_LEFT;
			pLeftCell->UpdateCellProperties( &LeftCellProperties );
			pLeftActiveCell = pLeftCell;
		}
		pActiveCell = pLeftCell->pGetPointer( RIGHT );
		pActiveCell->GetCellProperties( &CellProperties );
		x[1] = LeftCellProperties.s[1];
		x[2] = CellProperties.s[1];
		y[1] = CellProperties.rho_c + ( LeftCellProperties.rho[1] * LeftCellProperties.cell_width );
		y[2] = CellProperties.rho_c;	
		LinearFit( x, y, fZ_c, &fMcZ_c );
		pRadiativeRates->SetMcZ_c_LEFT( i, fMcZ_c );

#ifdef OPEN_FIELD
#else // OPEN_FIELD
		// Right-hand chromosphere
		pRightActiveCell = pCentreOfCurrentRow;
		for( ;; )
		{
			pRightCell = pRightActiveCell->pGetPointer( RIGHT );
			pRightCell->GetCellProperties( &RightCellProperties );

			if( RightCellProperties.s[1] > ( Params.L - fZ_c ) ) break;
			
			fMcZ_c_RIGHT += RightCellProperties.rho[1] * RightCellProperties.cell_width;
			RightCellProperties.rho_c = fMcZ_c_RIGHT;
			pRightCell->UpdateCellProperties( &RightCellProperties );
			pRightActiveCell = pRightCell;
		}
		pActiveCell = pRightCell->pGetPointer( LEFT );
		pActiveCell->GetCellProperties( &CellProperties );
		x[1] = CellProperties.s[1];
		x[2] = RightCellProperties.s[1];
		y[1] = CellProperties.rho_c;
		y[2] = CellProperties.rho_c + ( RightCellProperties.rho[1] * RightCellProperties.cell_width );
		LinearFit( x, y, ( Params.L - fZ_c ), &fMcZ_c );
		pRadiativeRates->SetMcZ_c_RIGHT( i, fMcZ_c );
#endif // OPEN_FIELD
	}

	// Free the memory allocated to the tabulated VAL atmosphere
	free( ppVALAtmosphere[1] );
	free( ppVALAtmosphere[0] );
	free( ppVALAtmosphere );
}
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

void CEquations::EvaluateTerms( double current_time, double *delta_t, int iFirstStep )
{
PCELL pNextActiveCell, pFarLeftCell, pLeftCell, pRightCell;
CELLPROPERTIES CellProperties, FarLeftCellProperties, LeftCellProperties, RightCellProperties;

// Variables used by advective flux transport algorithm
double Q1, Q2, Q3, QT;

// Variables used by thermal and viscous flux transport algorithms
double Kappa[SPECIES], max_flux_coeff[SPECIES];
	Kappa[ELECTRON] = SPITZER_ELECTRON_CONDUCTIVITY;
	Kappa[HYDROGEN] = SPITZER_ION_CONDUCTIVITY;
	max_flux_coeff[ELECTRON] = 1.5 / SQRT_ELECTRON_MASS; 
	max_flux_coeff[HYDROGEN] = 1.5 / SQRT_AVERAGE_PARTICLE_MASS;
double T[3][SPECIES], gradT, n[SPECIES], P, v[2], gradv, Kappa_B, Fc_max;

#if defined (NON_EQUILIBRIUM_RADIATION) || ( defined(OPTICALLY_THICK_RADIATION) && defined (NLTE_CHROMOSPHERE) )
	PCELL pFarRightCell;
	CELLPROPERTIES FarRightCellProperties;
	double ps[5];
#ifdef NON_EQUILIBRIUM_RADIATION
	// Variables used for time-dependent ionisation balance calculation
	double **ppni0, **ppni1, **ppni2, **ppni3, **ppni4, **ppdnibydt;
#endif // NON_EQUILIBRIUM_RADIATION
#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
	double *pHstate0, *pHstate1, *pHstate2, *pHstate3, *pHstate4;
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
#endif // NON_EQUILIBRIUM_RADIATION || ( OPTICALLY_THICK_RADIATION && NLTE_CHROMOSPHERE )

#ifdef OPTICALLY_THICK_RADIATION
	double fHI_c;
#ifndef NLTE_CHROMOSPHERE
	double frho_c;
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

#ifdef BEAM_HEATING
	double BeamParams[3];
	pHeat->CalculateBeamParameters( current_time, BeamParams );

	// Cut-off energy of the beam (erg)
	// 1.602e-9 erg / keV
	double cutoff_energy = (1.602e-9) * BeamParams[1];

	// The power-law index of the beam's energy distribution
	// delta must strictly be greater than 2
	double delta = BeamParams[2];
    
	// Variables used by the beam heating function
	double fColumnDensity;				// Column density (cm^-2)
	double fColumnDensitystar;			// Nstar (defined in Hawley & Fisher 1994, ApJ, 426, 287)
	double avg_energy, log_avg_energy;	// Average energy of the electron distribution
	double fLambda1, fLambda2;
		avg_energy = ( ( 1.0 - delta )/( 2.0 - delta ) ) * cutoff_energy;
		log_avg_energy = log(avg_energy);
		fLambda2 = 25.1 + log_avg_energy;

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
	// Evaluate terms doesn't time-step in the two left-most and two right-most grid cells
	pHeat->InitQbeam( avg_energy, Params.iNumberOfCells-4 );
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

	// Variables used for electron energy loss/gain associated with ionization and recombination
	double pHydrogen[2], pHelium[3];
	double HI_IonRate, HI_RecRate;
	double AbHe;
	double HeI_IonRate, HeI_RecRate;
	double HeII_IonRate, HeII_RecRate;
	double tau_IR;
#endif // BEAM_HEATING

#ifdef USE_KINETIC_MODEL
	CalculateKineticModel( iFirstStep );
#endif // USE_KINETIC_MODEL

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
#if defined (NON_EQUILIBRIUM_RADIATION) || ( defined(OPTICALLY_THICK_RADIATION) && defined (NLTE_CHROMOSPHERE) )
	double fCrossSection[5];
#else // NON_EQUILIBRIUM_RADIATION || ( OPTICALLY_THICK_RADIATION && NLTE_CHROMOSPHERE )
	double fCrossSection[4];
#endif //  NON_EQUILIBRIUM_RADIATION || ( OPTICALLY_THICK_RADIATION && NLTE_CHROMOSPHERE )
	double fCellVolume;
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

#ifdef USE_JB
	double Tmin;
#ifdef OPTICALLY_THICK_RADIATION
	Tmin = OPTICALLY_THICK_TEMPERATURE;
#else // OPTICALLY_THICK_RADIATION
	Tmin = MINIMUM_RADIATION_TEMPERATURE;
#endif // OPTICALLY_THICK_RADIATION
	if( Tc > (JB_TEMPERATURE_FRACTION*Te_max) ) Tc = JB_TEMPERATURE_FRACTION * Te_max;
#endif // USE_JB

#ifdef NUMERICAL_VISCOSITY
	double rho_v[2];
#endif // NUMERICAL_VISCOSITY

// Variables used for interpolation
double x[5], y[5], UpperValue, LowerValue, error;

// General variables
double term1, term2;
int j;

// Macros needed for OpenMP pragma on the loop that calculates the terms of the equations
#ifdef OPENMP
	#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
		#define USE_POLY_FIT_VARS fCrossSection, fCellVolume,
	#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
		#define USE_POLY_FIT_VARS
	#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
	#ifdef NON_EQUILIBRIUM_RADIATION
		#define NLTE_VARS_2 ppni0, ppni1, ppni2, ppni3, ppni4, ppdnibydt,
	#else // NON_EQUILIBRIUM_RADIATION
		#define NLTE_VARS_2
	#endif // NON_EQUILIBRIUM_RADIATION
	#ifdef OPTICALLY_THICK_RADIATION
	#ifdef NLTE_CHROMOSPHERE
		#define NLTE_CHR_VARS_2 pHstate0, pHstate1, pHstate2, pHstate3, pHstate4,
	#else // NLTE_CHROMOSPHERE
		#define NLTE_CHR_VARS_2
	#endif // NLTE_CHROMOSPHERE
	#else // OPTICALLY_THICK_RADIATION
		#define NLTE_CHR_VARS_2
	#endif // OPTICALLY_THICK_RADIATION
	#ifdef BEAM_HEATING
		#define BEAM_HEATING_VARS_2 pHydrogen, pHelium, HI_IonRate, HI_RecRate, AbHe, HeI_IonRate, HeI_RecRate, HeII_IonRate, HeII_RecRate, tau_IR,
	#else // BEAM_HEATING
		#define BEAM_HEATING_VARS_2
	#endif // BEAM_HEATING
	#if defined (NON_EQUILIBRIUM_RADIATION) || ( defined(OPTICALLY_THICK_RADIATION) && defined (NLTE_CHROMOSPHERE) )
		#define LEFT_CELL_VARS pLeftCell, LeftCellProperties,
		#define RIGHT_CELL_VARS pRightCell, RightCellProperties,
		#define FAR_LEFT_CELL_VARS pFarLeftCell, FarLeftCellProperties,
		#define FAR_RIGHT_CELL_VARS pFarRightCell, FarRightCellProperties, ps,
	#else // (NON_EQUILIBRIUM_RADIATION) || ( defined(OPTICALLY_THICK_RADIATION) && defined (NLTE_CHROMOSPHERE) )
		#define LEFT_CELL_VARS
		#define RIGHT_CELL_VARS
		#define FAR_LEFT_CELL_VARS
		#define FAR_RIGHT_CELL_VARS
	#endif // (NON_EQUILIBRIUM_RADIATION) || ( defined(OPTICALLY_THICK_RADIATION) && defined (NLTE_CHROMOSPHERE) )
#endif // OPENMP

// *****************************************************************************
// *                                                                           *
// *    CALCULATE THE CELL-INTERFACE TERMS                                     *
// *                                                                           *
// *****************************************************************************

	pNextActiveCell = pStartOfCurrentRow->pGetPointer( RIGHT )->pGetPointer( RIGHT );

	while( pNextActiveCell->pGetPointer( RIGHT ) )
	{
    	pActiveCell = pNextActiveCell;
    	pActiveCell->GetCellProperties( &CellProperties );

	    pLeftCell = pActiveCell->pGetPointer( LEFT );
	    pLeftCell->GetCellProperties( &LeftCellProperties );

	    pRightCell = pActiveCell->pGetPointer( RIGHT );
   		pRightCell->GetCellProperties( &RightCellProperties );

	    pFarLeftCell = pLeftCell->pGetPointer( LEFT );
	    pFarLeftCell->GetCellProperties( &FarLeftCellProperties );

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
    	fCrossSection[0] = CalculateCrossSection( FarLeftCellProperties.s[1]/Params.L );
    	fCrossSection[1] = CalculateCrossSection( LeftCellProperties.s[1]/Params.L );
    	fCrossSection[2] = CalculateCrossSection( CellProperties.s[1]/Params.L );
    	fCrossSection[3] = CalculateCrossSection( RightCellProperties.s[1]/Params.L );
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

// *****************************************************************************
// *    ADVECTIVE FLUX TRANSPORT ALGORITHM									   *
// *****************************************************************************

	    x[1] = LeftCellProperties.s[1];
	    x[2] = CellProperties.s[1];
 		y[1] = LeftCellProperties.v[1];
	    y[2] = CellProperties.v[1];
	    LinearFit( x, y, CellProperties.s[0], &(CellProperties.v[0]) );
	    LeftCellProperties.v[2] = CellProperties.v[0];

	    if( CellProperties.v[0] > 0.0 )
    	{
			// CALCULATE THE DENSITY
	 
	        x[1] = FarLeftCellProperties.s[1];
			x[2] = LeftCellProperties.s[1];
			y[1] = FarLeftCellProperties.rho[1];
			y[2] = LeftCellProperties.rho[1];
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
			y[1] *= fCrossSection[0];
			y[2] *= fCrossSection[1];
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
			LinearFit( x, y, CellProperties.s[0], &Q1 );

			x[1] = LeftCellProperties.s[1];
			x[2] = CellProperties.s[1];
			y[1] = LeftCellProperties.rho[1];
			y[2] = CellProperties.rho[1];
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
			y[1] *= fCrossSection[1];
			y[2] *= fCrossSection[2];
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
			LinearFit( x, y, CellProperties.s[0], &Q2 );

	        Q3 = y[1];

			if( y[2] <= y[1] )
			{
		    	QT = max( Q1, Q2 );
			    if( Q3 < QT )
		    	    CellProperties.rho[0] = Q3;
	    		else
	        		CellProperties.rho[0] = QT;
			}
			else
			{
	    		QT = min( Q1, Q2 );
		    	if( Q3 > QT )
		        	CellProperties.rho[0] = Q3;
	    		else
	        		CellProperties.rho[0] = QT;
			}
	        
			LeftCellProperties.rho[2] = CellProperties.rho[0];
		
			// CALCULATE THE MOMENTUM

	        x[1] = FarLeftCellProperties.s[1];
			x[2] = LeftCellProperties.s[1];
			y[1] = FarLeftCellProperties.rho_v[1];
			y[2] = LeftCellProperties.rho_v[1];
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
			y[1] *= fCrossSection[0];
			y[2] *= fCrossSection[1];
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
			LinearFit( x, y, CellProperties.s[0], &Q1 );

			x[1] = LeftCellProperties.s[1];
			x[2] = CellProperties.s[1];
			y[1] = LeftCellProperties.rho_v[1];
			y[2] = CellProperties.rho_v[1];
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
			y[1] *= fCrossSection[1];
			y[2] *= fCrossSection[2];
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
			LinearFit( x, y, CellProperties.s[0], &Q2 );

	        Q3 = y[1];

			if( y[2] <= y[1] )
			{
		    	QT = max( Q1, Q2 );
	    		if( Q3 < QT )
	        		CellProperties.rho_v[0] = Q3;
	    		else
	        		CellProperties.rho_v[0] = QT;
			}
			else
			{
	    		QT = min( Q1, Q2 );
	    		if( Q3 > QT )
	        		CellProperties.rho_v[0] = Q3;
	    		else
	        		CellProperties.rho_v[0] = QT;
			}

			LeftCellProperties.rho_v[2] = CellProperties.rho_v[0];

			// CALCULATE THE ENERGY
		
			for( j=0; j<SPECIES; j++ )
			{
        	    x[1] = FarLeftCellProperties.s[1];
	            x[2] = LeftCellProperties.s[1];
    	        y[1] = FarLeftCellProperties.TE_KE_P[1][j];
        	    y[2] = LeftCellProperties.TE_KE_P[1][j];
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
				y[1] *= fCrossSection[0];
				y[2] *= fCrossSection[1];
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
            	LinearFit( x, y, CellProperties.s[0], &Q1 );

	            x[1] = LeftCellProperties.s[1];
    	        x[2] = CellProperties.s[1];
        	    y[1] = LeftCellProperties.TE_KE_P[1][j];
            	y[2] = CellProperties.TE_KE_P[1][j];
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
				y[1] *= fCrossSection[1];
				y[2] *= fCrossSection[2];
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
            	LinearFit( x, y, CellProperties.s[0], &Q2 );

	            Q3 = y[1];

    	        if( y[2] <= y[1] )
        	    {
            	    QT = max( Q1, Q2 );
					if( Q3 < QT )
                    	CellProperties.TE_KE_P[0][j] = Q3;
					else
		    			CellProperties.TE_KE_P[0][j] = QT;
            	}
            	else
            	{
					QT = min( Q1, Q2 );
					if( Q3 > QT )
		    			CellProperties.TE_KE_P[0][j] = Q3;
					else
		    			CellProperties.TE_KE_P[0][j] = QT;
            	}

            	LeftCellProperties.TE_KE_P[2][j] = CellProperties.TE_KE_P[0][j];
			}
    	}
    	else
    	{
			// CALCULATE THE DENSITY
	        
	        x[1] = CellProperties.s[1];
			x[2] = RightCellProperties.s[1];
			y[1] = CellProperties.rho[1];
			y[2] = RightCellProperties.rho[1];
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
			y[1] *= fCrossSection[2];
			y[2] *= fCrossSection[3];
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
			LinearFit( x, y, CellProperties.s[0], &Q1 );

			x[1] = LeftCellProperties.s[1];
			x[2] = CellProperties.s[1];
			y[1] = LeftCellProperties.rho[1];
			y[2] = CellProperties.rho[1];
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
			y[1] *= fCrossSection[1];
			y[2] *= fCrossSection[2];
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
			LinearFit( x, y, CellProperties.s[0], &Q2 );

	        Q3 = y[2];

			// Note: The flow is in the opposite direction and so the conditional is switched
			if( y[1] <= y[2] )
			{
    	        QT = min( Q1, Q2 );
	    		if( Q3 > QT )
	        		CellProperties.rho[0] = Q3;
	    		else
	        		CellProperties.rho[0] = QT;
			}
			else
			{
	    		QT = max( Q1, Q2 );
		    	if( Q3 < QT )
		        	CellProperties.rho[0] = Q3;
	    		else
	        		CellProperties.rho[0] = QT;
			}

			LeftCellProperties.rho[2] = CellProperties.rho[0];
		
			// CALCULATE THE MOMENTUM
		
			x[1] = CellProperties.s[1];
			x[2] = RightCellProperties.s[1];
			y[1] = CellProperties.rho_v[1];
			y[2] = RightCellProperties.rho_v[1];
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
			y[1] *= fCrossSection[2];
			y[2] *= fCrossSection[3];
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
			LinearFit( x, y, CellProperties.s[0], &Q1 );

			x[1] = LeftCellProperties.s[1];
			x[2] = CellProperties.s[1];
			y[1] = LeftCellProperties.rho_v[1];
			y[2] = CellProperties.rho_v[1];
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
			y[1] *= fCrossSection[1];
			y[2] *= fCrossSection[2];
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
			LinearFit( x, y, CellProperties.s[0], &Q2 );

	        Q3 = y[2];

			// Note: The flow is in the opposite direction and so the conditional is switched
			if( y[1] <= y[2] )
			{
	    		QT = min( Q1, Q2 );
	    		if( Q3 > QT )
	        		CellProperties.rho_v[0] = Q3;
	    		else
	        		CellProperties.rho_v[0] = QT;
			}
			else
			{
	    		QT = max( Q1, Q2 );
	    		if( Q3 < QT )
	        		CellProperties.rho_v[0] = Q3;
	    		else
	        		CellProperties.rho_v[0] = QT;
			}

			LeftCellProperties.rho_v[2] = CellProperties.rho_v[0];

			// CALCULATE THE ENERGY
		
    	    for( j=0; j<SPECIES; j++ )
			{
	            x[1] = CellProperties.s[1];
    	        x[2] = RightCellProperties.s[1];
        	    y[1] = CellProperties.TE_KE_P[1][j];
            	y[2] = RightCellProperties.TE_KE_P[1][j];
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
				y[1] *= fCrossSection[2];
				y[2] *= fCrossSection[3];
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
            	LinearFit( x, y, CellProperties.s[0], &Q1 );

	            x[1] = LeftCellProperties.s[1];
    	        x[2] = CellProperties.s[1];
        	    y[1] = LeftCellProperties.TE_KE_P[1][j];
            	y[2] = CellProperties.TE_KE_P[1][j];
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
				y[1] *= fCrossSection[1];
				y[2] *= fCrossSection[2];
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
            	LinearFit( x, y, CellProperties.s[0], &Q2 );

	            Q3 = y[2];

				// Note: The flow is in the opposite direction and so the conditional is switched
    	        if( y[1] <= y[2] )
        	    {
					QT = min( Q1, Q2 );
					if( Q3 > QT )
		    			CellProperties.TE_KE_P[0][j] = Q3;
					else
                    	CellProperties.TE_KE_P[0][j] = QT;
            	}
            	else
            	{
					QT = max( Q1, Q2 );
					if( Q3 < QT )
		    			CellProperties.TE_KE_P[0][j] = Q3;
					else
	        			CellProperties.TE_KE_P[0][j] = QT;
            	}

	            LeftCellProperties.TE_KE_P[2][j] = CellProperties.TE_KE_P[0][j];
			}
    	}

// *****************************************************************************
// *    THERMAL FLUX TRANSPORT ALGORITHM                                       *
// *****************************************************************************

	    x[1] = FarLeftCellProperties.s[1];
	    x[2] = LeftCellProperties.s[1];
    	x[3] = CellProperties.s[1];
    	x[4] = RightCellProperties.s[1];

	    // Only need to calculate the diffusive terms on the first time-step
    	if( iFirstStep )
    	{
#ifdef USE_KINETIC_MODEL
			j = ELECTRON;

			y[1] = FarLeftCellProperties.Fc[1][j];
			y[2] = LeftCellProperties.Fc[1][j];
			y[3] = CellProperties.Fc[1][j];
			y[4] = RightCellProperties.Fc[1][j];
		
			FitPolynomial4( x, y, CellProperties.s[0], &(CellProperties.Fc[0][j]), &error );
			LeftCellProperties.Fc[2][j] = CellProperties.Fc[0][j];

			for( j=1; j<SPECIES; j++ )
#else // USE_KINETIC_MODEL
			// Calculate the conducted heat fluxes for the different species
			for( j=0; j<SPECIES; j++ )
#endif // USE_KINETIC_MODEL
			{
 	 	  		y[1] = FarLeftCellProperties.T[j];
    	    	y[2] = LeftCellProperties.T[j];
        		y[3] = CellProperties.T[j];
        		y[4] = RightCellProperties.T[j];
            	FitPolynomial4( x, y, ( CellProperties.s[1] - CellProperties.cell_width ), &(T[0][j]), &error );
            	FitPolynomial4( x, y, CellProperties.s[0], &(T[1][j]), &error );
            	T[2][j] = CellProperties.T[j];
	
			    y[1] = FarLeftCellProperties.n[j];
			    y[2] = LeftCellProperties.n[j];
	    		y[3] = CellProperties.n[j];
	    		y[4] = RightCellProperties.n[j];
	    		FitPolynomial4( x, y, CellProperties.s[0], &(n[j]), &error );

	            // Calculate the conducted heat flux at the left boundary
		    	gradT = ( T[2][j] - T[0][j] ) / CellProperties.cell_width;

#ifdef USE_JB
	    		if( j == ELECTRON && CellProperties.T[j] > Tmin && CellProperties.T[j] < Tc )
	    			T[1][j] = Tc;
#endif // USE_JB

#ifdef OPTICALLY_THICK_RADIATION
            	if( j == HYDROGEN && T[1][ELECTRON] < OPTICALLY_THICK_TEMPERATURE )
					Kappa_B = ( Kappa[j] + pHI->Getkappa_0( log10(T[1][ELECTRON]) ) ) * pow( T[1][j], 2.5 );
            	else
#endif // OPTICALLY_THICK_RADIATION
					Kappa_B = Kappa[j] * pow( T[1][j], 2.5 );

            	CellProperties.Fc[0][j] = - Kappa_B * gradT;

#ifdef USE_BKE
				if( j == ELECTRON )
				{
					// Calculate the thermal conduction suppression factor, R, due to scattering by turbulence
					// c_R / k_B = 2.28e4
					term1 = (2.28e4) * ( ( T[1][j] * T[1][j] ) / ( BKE_L_T * n[j] ) );
					term2 = CellProperties.s[0] - BKE_S_R;
					term2 *= term2;
					term1 *= exp( - term2 / BKE_2L_R2 );
					// Apply the Bian, Kontar, & Emslie thermal conduction suppression factor
					CellProperties.Fc[0][j] /= 1.0 + term1;
				}
#endif // USE_BKE

	    		// Apply the heat flux limiter
	            // BOLTZMANN_CONSTANT^1.5 = 1.621132937e-24
        	    Fc_max = HEAT_FLUX_LIMITING_COEFFICIENT * (1.621132937e-24) * max_flux_coeff[j] * n[j] * pow( T[1][j], 1.5 );
            	term1 = CellProperties.Fc[0][j] * Fc_max;
            	term2 = ( CellProperties.Fc[0][j] * CellProperties.Fc[0][j] ) + ( Fc_max * Fc_max );
            	CellProperties.Fc[0][j] =  term1 / sqrt( term2 );

// **** THERMAL CONDUCTION TIME STEP ****
				// Recalculate the coefficient of thermal conduction
				Kappa_B = fabs( CellProperties.Fc[0][j] / gradT );
		    	// The time-step is calculated by: delta_t = ( n * k_B ) * ( cell width )^2 / ( 2.0 * coefficient of thermal conduction )
		    	// 0.5 * BOLTZMANN_CONSTANT = 6.9e-17
	    		term1 = SAFETY_CONDUCTION * (6.9e-17) * CellProperties.cell_width * CellProperties.cell_width;
				CellProperties.conduction_delta_t[j] = ( term1 * n[j] ) / Kappa_B;
			    if( CellProperties.conduction_delta_t[j] < TIME_STEP_LIMIT )
			    {
	    			Kappa_B = ( term1 * n[j] ) / TIME_STEP_LIMIT;
					CellProperties.Fc[0][j] = - Kappa_B * gradT;
					CellProperties.conduction_delta_t[j] = TIME_STEP_LIMIT;
	    		}
			    LeftCellProperties.conduction_delta_t[j] += CellProperties.conduction_delta_t[j];
		    	LeftCellProperties.conduction_delta_t[j] /= 2.0;
// **** THERMAL CONDUCTION TIME STEP ****

	            LeftCellProperties.Fc[2][j] = CellProperties.Fc[0][j];
		
    	        if( pLeftCell->pGetPointer( LEFT )->pGetPointer( LEFT ) )
        	        LeftCellProperties.Fc[1][j] = 0.5 * ( LeftCellProperties.Fc[0][j] + LeftCellProperties.Fc[2][j] );
			}

// *****************************************************************************
// *    VISCOUS FLUX TRANSPORT ALGORITHM                                       *
// *****************************************************************************

			j = HYDROGEN;
			y[1] = FarLeftCellProperties.v[1];
			y[2] = LeftCellProperties.v[1];
			y[3] = CellProperties.v[1];
			y[4] = RightCellProperties.v[1];
			FitPolynomial4( x, y, ( CellProperties.s[1] - CellProperties.cell_width ), &(v[0]), &error );

			v[1] = CellProperties.v[1];

			// Calculate the viscosity flux at the left boundary
			gradv = ( v[1] - v[0] ) / CellProperties.cell_width;
			CellProperties.eta = DYNAMIC_VISCOSITY * ( pow( T[1][j], 2.5 ) / fLogLambda_ii( T[1][j], n[j] ) );

// **** VISCOUS TIME STEP ****
			// The time-step is calculated by: delta_t = ( rho ) * ( cell width )^2 / ( 2.0 * (4/3) * coefficient of dynamic viscosity )
			term1 = SAFETY_VISCOSITY * ( ( AVERAGE_PARTICLE_MASS * CellProperties.cell_width * CellProperties.cell_width ) / 2.6666666666666666666666666666667 );
    		CellProperties.viscosity_delta_t = ( term1 * n[j] ) / CellProperties.eta;

			if( CellProperties.viscosity_delta_t < TIME_STEP_LIMIT )
			{
    			CellProperties.eta = ( term1 * n[j] ) / TIME_STEP_LIMIT;
	    		CellProperties.viscosity_delta_t = TIME_STEP_LIMIT;
			}

			LeftCellProperties.viscosity_delta_t += CellProperties.viscosity_delta_t;
			LeftCellProperties.viscosity_delta_t /= 2.0;
// **** VISCOUS TIME STEP ****

			CellProperties.Feta[0] = 1.3333333333333333333333333333333 * CellProperties.eta * gradv;

			// The viscous flux forms the anisotropic part of the pressure tensor and therefore cannot exceed the total pressure
			// Hence, it must be limited to the value of the total pressure
			P = BOLTZMANN_CONSTANT * n[j] * T[1][j];
			term1 = CellProperties.Feta[0] * P;
			term2 = ( CellProperties.Feta[0] * CellProperties.Feta[0] ) + ( P * P );
			CellProperties.Feta[0] =  term1 / sqrt( term2 );

			LeftCellProperties.Feta[2] = CellProperties.Feta[0];

			if( pLeftCell->pGetPointer( LEFT )->pGetPointer( LEFT ) )
    	        LeftCellProperties.Feta[1] = 0.5 * ( LeftCellProperties.Feta[0] + LeftCellProperties.Feta[2] );

// *****************************************************************************
// *    NUMERICAL VISCOSITY                                                    *
// *****************************************************************************

#ifdef NUMERICAL_VISCOSITY
			// This calculation uses the cell boundary mass density from the advection algorithm and the velocity gradient across the boundary from the viscosity algorithm
			term1 = ( CellProperties.cell_width * CellProperties.cell_width ) / ( 2.0 * RELATIVE_VISCOUS_TIME_SCALE * ( CellProperties.advection_delta_t / SAFETY_ADVECTION ) );
			CellProperties.Fnumerical[0] = term1 * CellProperties.rho[0] * gradv;
        	LeftCellProperties.Fnumerical[2] = CellProperties.Fnumerical[0];
#endif // NUMERICAL_VISCOSITY
    	}

// *****************************************************************************
// *    CELL BOUNDARY PRESSURES                                                *
// *****************************************************************************

	    for( j=0; j<SPECIES; j++ )
    	{
			y[1] = FarLeftCellProperties.P[1][j];
			y[2] = LeftCellProperties.P[1][j];
			y[3] = CellProperties.P[1][j];
			y[4] = RightCellProperties.P[1][j];
        	FitPolynomial4( x, y, CellProperties.s[0], &(CellProperties.P[0][j]), &error );
			LeftCellProperties.P[2][j] = CellProperties.P[0][j];
    	}

	    pActiveCell->UpdateCellProperties( &CellProperties );
	    pLeftCell->UpdateCellProperties( &LeftCellProperties );

	    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
	}

// *****************************************************************************
// *                                                                           *
// *    CALCULATE THE FIELD LINE INTEGRATED TERMS                              *
// *                                                                           *
// *****************************************************************************

// *****************************************************************************
// *    COLUMN NUMBER AND MASS DENSITIES								       *
// *****************************************************************************

#if defined(OPTICALLY_THICK_RADIATION) || defined(BEAM_HEATING)
	// Left-hand leg of the loop
#ifdef OPTICALLY_THICK_RADIATION
	fHI_c = 0.0;
#ifndef NLTE_CHROMOSPHERE
	frho_c = 0.0;
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
#ifdef BEAM_HEATING
	fColumnDensity = 0.0;
	fColumnDensitystar = 0.0;
#endif // BEAM_HEATING
	pNextActiveCell = pCentreOfCurrentRow;
	while( pNextActiveCell )
	{
    	pActiveCell = pNextActiveCell;
    	pActiveCell->GetCellProperties( &CellProperties );

#ifdef OPTICALLY_THICK_RADIATION
    	if( CellProperties.T[ELECTRON] < OPTICALLY_THICK_TEMPERATURE )
    	{
        	// Calculate the neutral hydrogen column number density
        	fHI_c += ( CellProperties.HI * CellProperties.n[HYDROGEN] ) * CellProperties.cell_width;
        	CellProperties.HI_c = fHI_c;
#ifndef NLTE_CHROMOSPHERE
        	// Calculate the mass column density
        	frho_c += CellProperties.rho[1] * CellProperties.cell_width;
        	CellProperties.rho_c = frho_c;
#endif // NLTE_CHROMOSPHERE
    	}
#endif // OPTICALLY_THICK_RADIATION

#ifdef BEAM_HEATING
		fLambda1 = 66.0 + ( 1.50 * log_avg_energy ) - ( 0.5 * log(CellProperties.n[ELECTRON]) );
		fColumnDensity += CellProperties.n[HYDROGEN] * CellProperties.cell_width;
		fColumnDensitystar += ( ( (fLambda1*(1.0-CellProperties.HI)) + (fLambda2*CellProperties.HI) ) / fLambda1 ) * CellProperties.n[HYDROGEN] * CellProperties.cell_width;
		CellProperties.nH_c = fColumnDensity;
		CellProperties.nH_star_c = fColumnDensitystar;
#endif // BEAM_HEATING

	    pActiveCell->UpdateCellProperties( &CellProperties );

    	pNextActiveCell = pActiveCell->pGetPointer( LEFT );
	}

#ifdef OPEN_FIELD
#else // OPEN_FIELD
	// Right-hand leg of the loop
#ifdef OPTICALLY_THICK_RADIATION
	fHI_c = 0.0;
#ifndef NLTE_CHROMOSPHERE
	frho_c = 0.0;
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
#ifdef BEAM_HEATING
	fColumnDensity = 0.0;
	fColumnDensitystar = 0.0;
#endif // BEAM_HEATING
	pNextActiveCell = pCentreOfCurrentRow;
	while( pNextActiveCell )
	{
    	pActiveCell = pNextActiveCell;
    	pActiveCell->GetCellProperties( &CellProperties );

#ifdef OPTICALLY_THICK_RADIATION
    	if( CellProperties.T[ELECTRON] < OPTICALLY_THICK_TEMPERATURE )
    	{
			// Calculate the neutral hydrogen column number density
        	fHI_c += ( CellProperties.HI * CellProperties.n[HYDROGEN] ) * CellProperties.cell_width;;
        	CellProperties.HI_c = fHI_c;
#ifndef NLTE_CHROMOSPHERE
        	// Calculate the mass column density
        	frho_c += CellProperties.rho[1] * CellProperties.cell_width;
        	CellProperties.rho_c = frho_c;
#endif // NLTE_CHROMOSPHERE
	    }
#endif // OPTICALLY_THICK_RADIATION

#ifdef BEAM_HEATING
		fLambda1 = 66.0 + ( 1.50 * log_avg_energy ) - ( 0.5 * log(CellProperties.n[ELECTRON]) );
		fColumnDensity += CellProperties.n[HYDROGEN] * CellProperties.cell_width;
		fColumnDensitystar += ( ( (fLambda1*(1.0-CellProperties.HI)) + (fLambda2*CellProperties.HI) ) / fLambda1 ) * CellProperties.n[HYDROGEN] * CellProperties.cell_width;
		CellProperties.nH_c = fColumnDensity;
		CellProperties.nH_star_c = fColumnDensitystar;
#endif // BEAM_HEATING

    	pActiveCell->UpdateCellProperties( &CellProperties );

	    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
	}
#endif // OPEN_FIELD
#endif // OPTICALLY_THICK_RADIATION || BEAM_HEATING

// *****************************************************************************
// *                                                                           *
// *    CALCULATE THE CELL-CENTERED TERMS                                      *
// *                                                                           *
// *****************************************************************************

// *****************************************************************************
// *    BEAM HEATING								       					   *
// *****************************************************************************
#ifdef BEAM_HEATING
	pNextActiveCell = pStartOfCurrentRow->pGetPointer( RIGHT )->pGetPointer( RIGHT );
	while( pNextActiveCell->pGetPointer( RIGHT )->pGetPointer( RIGHT ) )
	{
    	pActiveCell = pNextActiveCell;
    	pActiveCell->GetCellProperties( &CellProperties );

    	CellProperties.TE_KE_term[4][ELECTRON] = pHeat->CalculateBeamHeating( current_time, BeamParams, CellProperties.nH_c, CellProperties.nH_star_c, CellProperties.n[ELECTRON], CellProperties.n[HYDROGEN], (1.0-CellProperties.HI) );
		#ifdef OPTICALLY_THICK_RADIATION
			#ifdef NLTE_CHROMOSPHERE
				pHeat->SetQbeam( CellProperties.s[1], CellProperties.TE_KE_term[4][ELECTRON] );
			#endif // NLTE_CHROMOSPHERE
		#endif // OPTICALLY_THICK_RADIATION

    	pActiveCell->UpdateCellProperties( &CellProperties );

	    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
	}
#endif // BEAM_HEATING

// *****************************************************************************
// *    TERMS OF THE CONSERVATION EQUATIONS                                    *
// *****************************************************************************

#ifdef OPENMP
	PCELL pLocalActiveCell;
	int iCounter, iLocalNumberOfCells = Params.iNumberOfCells;
#pragma omp parallel shared( ppCellList, iLocalNumberOfCells )		\
   	                 private(	USE_POLY_FIT_VARS					\
								NLTE_VARS_2							\
								NLTE_CHR_VARS_2						\
								BEAM_HEATING_VARS_2					\
                             	LEFT_CELL_VARS					    \
                             	RIGHT_CELL_VARS					    \
                             	FAR_LEFT_CELL_VARS					\
								FAR_RIGHT_CELL_VARS					\
								iCounter,							\
								CellProperties,				  		\
                             	LowerValue, UpperValue,             \
                             	j, term1, term2 ) 
	{
#pragma omp for schedule(dynamic, CHUNK_SIZE)						\
   	            lastprivate(pLocalActiveCell)
	for( iCounter=2; iCounter<iLocalNumberOfCells-2; iCounter++ )
	{
		pLocalActiveCell = ppCellList[iCounter];
		pLocalActiveCell->GetCellProperties( &CellProperties );

		#if defined (NON_EQUILIBRIUM_RADIATION) || ( defined(OPTICALLY_THICK_RADIATION) && defined (NLTE_CHROMOSPHERE) )
	    	pLeftCell = pLocalActiveCell->pGetPointer( LEFT );
    		pLeftCell->GetCellProperties( &LeftCellProperties );

	    	pRightCell = pLocalActiveCell->pGetPointer( RIGHT );
	    	pRightCell->GetCellProperties( &RightCellProperties );

	    	pFarLeftCell = pLeftCell->pGetPointer( LEFT );
    		pFarLeftCell->GetCellProperties( &FarLeftCellProperties );

	    	pFarRightCell = pRightCell->pGetPointer( RIGHT );
    		pFarRightCell->GetCellProperties( &FarRightCellProperties );
		#endif // NON_EQUILIBRIUM_RADIATION || ( OPTICALLY_THICK_RADIATION && NLTE_CHROMOSPHERE )
#else // OPENMP
	pNextActiveCell = pStartOfCurrentRow->pGetPointer( RIGHT )->pGetPointer( RIGHT );
	while( pNextActiveCell->pGetPointer( RIGHT )->pGetPointer( RIGHT ) )
	{
	    pActiveCell = pNextActiveCell;
		pActiveCell->GetCellProperties( &CellProperties );

		#if defined (NON_EQUILIBRIUM_RADIATION) || ( defined(OPTICALLY_THICK_RADIATION) && defined (NLTE_CHROMOSPHERE) )
	    	pLeftCell = pActiveCell->pGetPointer( LEFT );
    		pLeftCell->GetCellProperties( &LeftCellProperties );

	    	pRightCell = pActiveCell->pGetPointer( RIGHT );
	    	pRightCell->GetCellProperties( &RightCellProperties );

	    	pFarLeftCell = pLeftCell->pGetPointer( LEFT );
    		pFarLeftCell->GetCellProperties( &FarLeftCellProperties );

	    	pFarRightCell = pRightCell->pGetPointer( RIGHT );
    		pFarRightCell->GetCellProperties( &FarRightCellProperties );
		#endif // NON_EQUILIBRIUM_RADIATION || ( OPTICALLY_THICK_RADIATION && NLTE_CHROMOSPHERE )
#endif // OPENMP

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
    	fCrossSection[0] = CalculateCrossSection( CellProperties.s[0]/Params.L );
    	fCrossSection[1] = CalculateCrossSection( CellProperties.s[1]/Params.L );
    	fCrossSection[2] = CalculateCrossSection( CellProperties.s[2]/Params.L );
    	fCellVolume = fCrossSection[1] * CellProperties.cell_width;
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

// *****************************************************************************
// *                                                                           *
// *    MASS CONSERVATION                                                      *
// *                                                                           *
// *****************************************************************************

// *****************************************************************************
// *    ADVECTION                                                              *
// *****************************************************************************

		// The cell interface quantities calculated with Barton's algorithm account for
		// a non-uniform magnetic field / cross-section in the field-aligned direction
    	LowerValue = CellProperties.rho[0] * CellProperties.v[0];
    	UpperValue = CellProperties.rho[2] * CellProperties.v[2];
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
		CellProperties.rho_term[0] = - ( UpperValue - LowerValue ) / fCellVolume;
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
    	CellProperties.rho_term[0] = - ( UpperValue - LowerValue ) / CellProperties.cell_width;
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

		// Calculate the time derivative
	    CellProperties.drhobydt = CellProperties.rho_term[0];

// *****************************************************************************
// *                                                                           *
// *    MOMENTUM CONSERVATION                                                  *
// *                                                                           *
// *****************************************************************************

// *****************************************************************************
// *    ADVECTION                                                              *
// *****************************************************************************

		// The cell interface quantities calculated with Barton's algorithm account for
		// a non-uniform magnetic field / cross-section in the field-aligned direction
	    LowerValue = CellProperties.rho_v[0] * CellProperties.v[0];
	    UpperValue = CellProperties.rho_v[2] * CellProperties.v[2];
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
	    CellProperties.rho_v_term[0] = - ( UpperValue - LowerValue ) / fCellVolume;
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
	    CellProperties.rho_v_term[0] = - ( UpperValue - LowerValue ) / CellProperties.cell_width;
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

// *****************************************************************************
// *    PRESSURE GRADIENT                                                      *
// *****************************************************************************

	    LowerValue = CellProperties.P[0][ELECTRON] + CellProperties.P[0][HYDROGEN];
	    UpperValue = CellProperties.P[2][ELECTRON] + CellProperties.P[2][HYDROGEN];
	    CellProperties.rho_v_term[1] = - ( UpperValue - LowerValue ) / CellProperties.cell_width;

// *****************************************************************************
// *    GRAVITY                                                                *
// *****************************************************************************

	    CellProperties.rho_v_term[2] = CellProperties.rho[1] * CalculateGravity( CellProperties.s[1]/Params.L );

	    // Terms that must be integrated to first order in time only, otherwise they're unconditionally unstable
	    if( iFirstStep )
    	{
// *****************************************************************************
// *    VISCOUS STRESS                                                         *
// *****************************************************************************

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
			CellProperties.rho_v_term[3] = ( ( CellProperties.Feta[2] * fCrossSection[2] ) - ( CellProperties.Feta[0] * fCrossSection[0] ) ) / fCellVolume;
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
			CellProperties.rho_v_term[3] = ( CellProperties.Feta[2] - CellProperties.Feta[0] ) / CellProperties.cell_width;
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

// *****************************************************************************
// *    NUMERICAL VISCOSITY                                                    *
// *****************************************************************************

#ifdef NUMERICAL_VISCOSITY
			// Numerical viscosity is used to stabilise the solutions as in the Lax scheme
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
			CellProperties.rho_v_term[4] = ( ( CellProperties.Fnumerical[2] * fCrossSection[2] ) - ( CellProperties.Fnumerical[0] * fCrossSection[0] ) ) / fCellVolume;
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
			CellProperties.rho_v_term[4] = ( CellProperties.Fnumerical[2] - CellProperties.Fnumerical[0] ) / CellProperties.cell_width;
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
#endif // NUMERICAL_VISCOSITY
    	}

		// Calculate the time derivative
    	CellProperties.drho_vbydt = CellProperties.rho_v_term[0] + CellProperties.rho_v_term[1] + CellProperties.rho_v_term[2] + CellProperties.rho_v_term[3] + CellProperties.rho_v_term[4];

// *****************************************************************************
// *    																	   *
// *    ENERGY CONSERVATION                                                    *
// *                                                                           *
// *****************************************************************************

// *****************************************************************************
// *    ADVECTION                                                              *
// *****************************************************************************

		for( j=0; j<SPECIES; j++ )
		{
			// The cell interface quantities calculated with Barton's algorithm account for
			// a non-uniform magnetic field / cross-section in the field-aligned direction
			LowerValue = CellProperties.TE_KE_P[0][j] * CellProperties.v[0];
			UpperValue = CellProperties.TE_KE_P[2][j] * CellProperties.v[2];
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
			CellProperties.TE_KE_term[0][j] = - ( UpperValue - LowerValue ) / fCellVolume;
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
			CellProperties.TE_KE_term[0][j] = - ( UpperValue - LowerValue ) / CellProperties.cell_width;
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
		}

	    // Terms that must be integrated to first order in time only, otherwise they're unconditionally unstable
    	if( iFirstStep )
    	{
// *****************************************************************************
// *    THERMAL CONDUCTION                                                     *
// *****************************************************************************

			for( j=0; j<SPECIES; j++ )
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
            	CellProperties.TE_KE_term[1][j] = - ( ( CellProperties.Fc[2][j] * fCrossSection[2] ) - ( CellProperties.Fc[0][j] * fCrossSection[0] ) ) / fCellVolume;
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
            	CellProperties.TE_KE_term[1][j] = - ( CellProperties.Fc[2][j] - CellProperties.Fc[0][j] ) / CellProperties.cell_width;
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

// *****************************************************************************
// *    VISCOUS STRESS                                                         *
// *****************************************************************************

	    	// Heating due to the viscous stress and work done on (by) the flow by (on) the viscous stress
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
			CellProperties.TE_KE_term[7][HYDROGEN] = ( ( CellProperties.Feta[2] * CellProperties.v[2] * fCrossSection[2] ) - ( CellProperties.Feta[0] * CellProperties.v[0] * fCrossSection[0] ) ) / fCellVolume;
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
			CellProperties.TE_KE_term[7][HYDROGEN] = ( ( CellProperties.Feta[2] * CellProperties.v[2] ) - ( CellProperties.Feta[0] * CellProperties.v[0] ) ) / CellProperties.cell_width;
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

// *****************************************************************************
// *    NUMERICAL VISCOSITY                                                    *
// *****************************************************************************

#ifdef NUMERICAL_VISCOSITY
			// Numerical viscosity is used to stabilise the solutions as in the Lax scheme
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
			CellProperties.TE_KE_term[8][HYDROGEN] = ( ( CellProperties.Fnumerical[2] * CellProperties.v[2] * fCrossSection[2] ) - ( CellProperties.Fnumerical[0] * CellProperties.v[0] * fCrossSection[0] ) ) / fCellVolume;
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
			CellProperties.TE_KE_term[8][HYDROGEN] = ( ( CellProperties.Fnumerical[2] * CellProperties.v[2] ) - ( CellProperties.Fnumerical[0] * CellProperties.v[0] ) ) / CellProperties.cell_width;
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
#endif // NUMERICAL_VISCOSITY
    	}

// *****************************************************************************
// *    GRAVITY                                                                *
// *****************************************************************************

	    CellProperties.TE_KE_term[2][HYDROGEN] = CellProperties.rho_v_term[2] * CellProperties.v[1];

// *****************************************************************************
// *    COLLISIONS                                                             *
// *****************************************************************************

	    // The collision frequency is calculated using n_e, therefore the collisional coupling depends on (n_e)(n_H)
    	CellProperties.TE_KE_term[3][ELECTRON] = (2.07e-16) * CellProperties.n[HYDROGEN] * CellProperties.nu_ie * ( CellProperties.T[HYDROGEN] - CellProperties.T[ELECTRON] );

// **** COLLISIONAL TIME STEP ****
    	// Set the collisional timescale to 1% of the timescale for order of magnitude changes in the electron energy due to collisional energy exchange
    	CellProperties.collision_delta_t = ( 0.01 * CellProperties.TE_KE[1][ELECTRON] ) / fabs( CellProperties.TE_KE_term[3][ELECTRON] );
    	// If the collisional timescale is less than the minimum specified collisional timescale then scale the rate of energy exchange so that tiny timesteps can be avoided
    	if( CellProperties.collision_delta_t < MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE )
        	CellProperties.TE_KE_term[3][ELECTRON] *= CellProperties.collision_delta_t / MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE;
// **** COLLISIONAL TIME STEP ****

	    CellProperties.TE_KE_term[3][HYDROGEN] = - CellProperties.TE_KE_term[3][ELECTRON];

// *****************************************************************************
// *    HEATING                                                                *
// *****************************************************************************

#ifdef BEAM_HEATING
	// The beam heating function is called in the section above where the cell-centered terms are calculated
#if defined (ELECTRON_HEATING_ONLY) || defined(HYDROGEN_HEATING_ONLY)
	#ifdef ELECTRON_HEATING_ONLY
    	CellProperties.TE_KE_term[4][ELECTRON] += pHeat->CalculateHeating( CellProperties.s[1], current_time );
    	CellProperties.TE_KE_term[4][HYDROGEN] = SMALLEST_DOUBLE;
	#endif // ELECTRON_HEATING_ONLY
	#ifdef HYDROGEN_HEATING_ONLY
    	CellProperties.TE_KE_term[4][HYDROGEN] = pHeat->CalculateHeating( CellProperties.s[1], current_time );
	#endif // HYDROGTEN_HEATING_ONLY
#else // ELECTRON_HEATING_ONLY || HYDROGEN_HEATING_ONLY
	term1 = pHeat->CalculateHeating( CellProperties.s[1], current_time );
   	CellProperties.TE_KE_term[4][ELECTRON] += ELECTRON_HEATING * term1;
	CellProperties.TE_KE_term[4][HYDROGEN] = HYDROGEN_HEATING * term1;
#endif // // ELECTRON_HEATING_ONLY || HYDROGEN_HEATING_ONLY
#else // BEAM_HEATING
#if defined (ELECTRON_HEATING_ONLY) || defined(HYDROGEN_HEATING_ONLY)
	#ifdef ELECTRON_HEATING_ONLY
    	CellProperties.TE_KE_term[4][ELECTRON] = pHeat->CalculateHeating( CellProperties.s[1], current_time );
    	CellProperties.TE_KE_term[4][HYDROGEN] = SMALLEST_DOUBLE;
	#endif // ELECTRON_HEATING_ONLY
	#ifdef HYDROGEN_HEATING_ONLY
    	CellProperties.TE_KE_term[4][ELECTRON] = SMALLEST_DOUBLE;
    	CellProperties.TE_KE_term[4][HYDROGEN] = pHeat->CalculateHeating( CellProperties.s[1], current_time );
	#endif // HYDROGTEN_HEATING_ONLY
#else // ELECTRON_HEATING_ONLY || HYDROGEN_HEATING_ONLY
	term1 = pHeat->CalculateHeating( CellProperties.s[1], current_time );
   	CellProperties.TE_KE_term[4][ELECTRON] = ELECTRON_HEATING * term1;
	CellProperties.TE_KE_term[4][HYDROGEN] = HYDROGEN_HEATING * term1;
#endif // // ELECTRON_HEATING_ONLY || HYDROGEN_HEATING_ONLY
#endif // BEAM_HEATING

// *****************************************************************************
// *    RADIATION                                                              *
// *****************************************************************************

		CellProperties.TE_KE_term[5][ELECTRON] = - SMALLEST_DOUBLE;
#ifdef OPTICALLY_THICK_RADIATION
    	if( CellProperties.T[ELECTRON] < OPTICALLY_THICK_TEMPERATURE )
    	{
        	CellProperties.TE_KE_term[4][ELECTRON] += pHeat->CalculateVALHeating( log10(CellProperties.rho_c ) );
        	CellProperties.TE_KE_term[5][ELECTRON] -= pHI->GetVolumetricLossRate( log10(CellProperties.T[ELECTRON]), log10((4e-14)*CellProperties.HI_c), CellProperties.n[ELECTRON] * CellProperties.rho[1] );
        	CellProperties.TE_KE_term[5][ELECTRON] -= pMgII->GetVolumetricLossRate( log10(CellProperties.T[ELECTRON]), log10(CellProperties.rho_c), CellProperties.n[ELECTRON] * CellProperties.rho[1] );
        	CellProperties.TE_KE_term[5][ELECTRON] -= pCaII->GetVolumetricLossRate( log10(CellProperties.T[ELECTRON]), log10(CellProperties.rho_c), CellProperties.n[ELECTRON] * CellProperties.rho[1] );
// **** RADIATION TIME STEP ****
        	CellProperties.radiation_delta_t = ( SAFETY_RADIATION * CellProperties.TE_KE[1][ELECTRON] ) / fabs( CellProperties.TE_KE_term[5][ELECTRON] );
// **** RADIATION TIME STEP ****
    	}
    	else
    	{
			term1 = 1.0;
#else // OPTICALLY_THICK_RADIATION
   		if( CellProperties.T[ELECTRON] < MINIMUM_RADIATION_TEMPERATURE )
   		{
       		// Provide some additional heating to the chromosphere if the temperature drops below the specified isothermal temperature
       		CellProperties.TE_KE_term[5][ELECTRON] = ( ( ( MINIMUM_RADIATION_TEMPERATURE / CellProperties.T[ELECTRON] ) - 1.0 ) * CellProperties.TE_KE[1][ELECTRON] ) / MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE;
// **** RADIATION TIME STEP ****
			CellProperties.radiation_delta_t = MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE;
// **** RADIATION TIME STEP ****
    	}
   		else
   		{   
			term1 = 1.0;
			// Decrease the radiation smoothly to zero in the chromosphere
			if( CellProperties.T[ELECTRON] < lower_radiation_temperature_boundary )
    			term1 = ( CellProperties.T[ELECTRON] - MINIMUM_RADIATION_TEMPERATURE ) / ZERO_OVER_TEMPERATURE_INTERVAL;
#endif // OPTICALLY_THICK_RADIATION

#ifdef DECOUPLE_IONISATION_STATE_SOLVER
	#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
			CellProperties.TE_KE_term[5][ELECTRON] -= term1 * pRadiation2->GetPowerLawRad( log10( CellProperties.T[ELECTRON] ), CellProperties.n[ELECTRON], CellProperties.n[HYDROGEN] );
	#else // USE_POWER_LAW_RADIATIVE_LOSSES
			// If USE_POWER_LAW_RADIATIVE_LOSSES hasn't been defined then the default must be to use the NON_EQUILIBRIUM_RADIATION method for calculating the radiative losses, because
			// NON_EQUILIBRIUM_RADIATION is *ALWAYS* defined when DECOUPLE_IONISATION_STATE_SOLVER is defined
			CellProperties.TE_KE_term[5][ELECTRON] -= term1 * ( pRadiation->GetRadiation( log10( CellProperties.T[ELECTRON] ), CellProperties.n[ELECTRON], CellProperties.n[HYDROGEN] ) + pRadiation2->GetRadiation( log10( CellProperties.T[ELECTRON] ), CellProperties.n[ELECTRON], CellProperties.n[HYDROGEN] ) + pRadiation2->GetFreeFreeRad( log10( CellProperties.T[ELECTRON] ), CellProperties.n[ELECTRON], CellProperties.n[HYDROGEN] ) );
	#endif // USE_POWER_LAW_RADIATIVE_LOSSES
#else // DECOUPLE_IONISATION_STATE_SOLVER
	#ifdef NON_EQUILIBRIUM_RADIATION
       		ppni2 = CellProperties.pIonFrac->ppGetIonFrac();
        	CellProperties.TE_KE_term[5][ELECTRON] -= term1 * ( pRadiation->GetRadiation( log10( CellProperties.T[ELECTRON] ), CellProperties.n[ELECTRON], CellProperties.n[HYDROGEN], ppni2 ) + pRadiation2->GetRadiation( log10( CellProperties.T[ELECTRON] ), CellProperties.n[ELECTRON], CellProperties.n[HYDROGEN] ) + pRadiation2->GetFreeFreeRad( log10( CellProperties.T[ELECTRON] ), CellProperties.n[ELECTRON], CellProperties.n[HYDROGEN] ) );
	#else // NON_EQUILIBRIUM_RADIATION
		#ifdef USE_POWER_LAW_RADIATIVE_LOSSES
			CellProperties.TE_KE_term[5][ELECTRON] -= term1 * pRadiation2->GetPowerLawRad( log10( CellProperties.T[ELECTRON] ), CellProperties.n[ELECTRON], CellProperties.n[HYDROGEN] );
		#else // USE_POWER_LAW_RADIATIVE_LOSSES
       		CellProperties.TE_KE_term[5][ELECTRON] -= term1 * ( pRadiation2->GetRadiation( log10( CellProperties.T[ELECTRON] ), CellProperties.n[ELECTRON], CellProperties.n[HYDROGEN] ) + pRadiation2->GetFreeFreeRad( log10( CellProperties.T[ELECTRON] ), CellProperties.n[ELECTRON], CellProperties.n[HYDROGEN] ) );
		#endif // USE_POWER_LAW_RADIATIVE_LOSSES
	#endif // NON_EQUILIBRIUM_RADIATION
#endif // DECOUPLE_IONISATION_STATE_SOLVER

#ifdef USE_JB
			if( CellProperties.T[ELECTRON] > Tmin && CellProperties.T[ELECTRON] < Tc )
			{
				term1 = CellProperties.T[ELECTRON] / Tc;
				term2 = pow( term1, 2.5 );
				// Modify heating (only when applied to electrons) per the TRAC method
				CellProperties.TE_KE_term[4][ELECTRON] *= term2;
				// Modify radiation (applies only to electrons) per the TRAC method
				CellProperties.TE_KE_term[5][ELECTRON] *= term2;
			}
#endif // USE_JB

// **** RADIATION TIME STEP ****
	    	CellProperties.radiation_delta_t = ( SAFETY_RADIATION * CellProperties.TE_KE[1][ELECTRON] ) / fabs( CellProperties.TE_KE_term[5][ELECTRON] );
// **** RADIATION TIME STEP ****
		}

// *****************************************************************************
// *    SMALL-SCALE ELECTRIC FIELDS                                            *
// *****************************************************************************

		// Derived from qnEv = v dP/ds
   		// The term added to the electron energy equation is (e)nEv = v dPe/ds
	    // The term added to the hydrogen energy equation is (-e)nEv = -v dPe/ds
    	CellProperties.TE_KE_term[6][ELECTRON] = CellProperties.v[1] * ( ( CellProperties.P[2][ELECTRON] - CellProperties.P[0][ELECTRON] ) / CellProperties.cell_width );
    	CellProperties.TE_KE_term[6][HYDROGEN] = -CellProperties.TE_KE_term[6][ELECTRON];

// *****************************************************************************
// *    IONIZATION / RECOMBINATION ENERGY                                      *
// *****************************************************************************

		term1 = 0.0;
#ifdef BEAM_HEATING
		// Only calculate the energy loss due to hydrogen and helium ionization where hydrogen is not fully ionized
		// Assume hydrogen is fully ionized above 0.04 MK
		if( CellProperties.T[ELECTRON] < 4E4 ) {
#ifdef DENSITY_DEPENDENT_RATES
			// Hydrogen
			pRadiation_H->GetEquilIonFrac( 1, pHydrogen, log10(CellProperties.T[ELECTRON]), log10(CellProperties.n[ELECTRON]) );
			pRadiation_H->GetRates( 1, 1, log10(CellProperties.T[ELECTRON]), log10(CellProperties.n[ELECTRON]), &HI_IonRate, &HI_RecRate );
			// Helium
			pRadiation_He->GetEquilIonFrac( 2, pHelium, log10(CellProperties.T[ELECTRON]), log10(CellProperties.n[ELECTRON]) );
			pRadiation_He->GetRates( 2, 1, log10(CellProperties.T[ELECTRON]), log10(CellProperties.n[ELECTRON]), &HeI_IonRate, &HeI_RecRate );
			pRadiation_He->GetRates( 2, 2, log10(CellProperties.T[ELECTRON]), log10(CellProperties.n[ELECTRON]), &HeII_IonRate, &HeII_RecRate );
#else // DENSITY_DEPENDENT_RATES
			// Hydrogen
			pRadiation_H->GetEquilIonFrac( 1, pHydrogen, log10(CellProperties.T[ELECTRON]) );
			pRadiation_H->GetRates( 1, 1, log10(CellProperties.T[ELECTRON]), &HI_IonRate, &HI_RecRate );
			// Helium
			pRadiation_He->GetEquilIonFrac( 2, pHelium, log10(CellProperties.T[ELECTRON]) );
			pRadiation_He->GetRates( 2, 1, log10(CellProperties.T[ELECTRON]), &HeI_IonRate, &HeI_RecRate );
			pRadiation_He->GetRates( 2, 2, log10(CellProperties.T[ELECTRON]), &HeII_IonRate, &HeII_RecRate );
#endif // DENSITY_DEPENDENT_RATES
			// Hydrogen
			term1 = - (2.179E-11) * ( ( pHydrogen[0] * HI_IonRate ) - ( pHydrogen[1] * HI_RecRate ) );
			// Helium
			AbHe = pRadiation_He->GetAbundance( 2 );
			term1 += - (3.940E-11) * AbHe * ( ( pHelium[0] * HeI_IonRate ) - ( pHelium[1] * HeI_RecRate ) );
			term1 += - (8.716E-11) * AbHe * ( ( pHelium[1] * HeII_IonRate ) - ( pHelium[2] * HeII_RecRate ) );
			term1 *= ( CellProperties.n[HYDROGEN] * CellProperties.n[ELECTRON] );

			// Scale the rate of energy loss for stability such that its timescale cannot be shorter than for collisions
			tau_IR = fabs( CellProperties.TE_KE[1][ELECTRON] / term1 );
			if( tau_IR < MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE ) {
				term1 *= ( tau_IR / MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE );
			}

			// Be sure not to spuriously add energy to the electrons
			if( isnan(term1) || term1 > 0.0 )
				term1 = 0.0;
		}
#endif // BEAM_HEATING
		CellProperties.TE_KE_term[9][ELECTRON] = term1;

		// Calculate the time derivative
    	for( j=0; j<SPECIES; j++ )
        	CellProperties.dTE_KEbydt[j] = CellProperties.TE_KE_term[0][j] + CellProperties.TE_KE_term[1][j] + CellProperties.TE_KE_term[2][j] + CellProperties.TE_KE_term[3][j] + CellProperties.TE_KE_term[4][j] + CellProperties.TE_KE_term[5][j] + CellProperties.TE_KE_term[6][j] + CellProperties.TE_KE_term[7][j] + CellProperties.TE_KE_term[8][j] + CellProperties.TE_KE_term[9][j];

// *****************************************************************************
// *                                                                           *
// *    TIME-DEPENDENT IONISATION                                              *
// *                                                                           *
// *****************************************************************************

#if defined (NON_EQUILIBRIUM_RADIATION) || ( defined(OPTICALLY_THICK_RADIATION) && defined (NLTE_CHROMOSPHERE) )
	    ps[0] = FarLeftCellProperties.s[1];
    	ps[1] = LeftCellProperties.s[1];
    	ps[2] = CellProperties.s[1];
    	ps[3] = RightCellProperties.s[1];
    	ps[4] = FarRightCellProperties.s[1];

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
		// Calculate the grid cell cross-sections corresponding to the cells containing the ion fractions
		// Note: It's safe to do this here because these variables are not used again in the current cell 
		//       and are reset to the cell center and interface cross-sections when the time-derivatives 
		//       are calculated in the next grid cell
    	fCrossSection[0] = CalculateCrossSection( FarLeftCellProperties.s[1]/Params.L );
    	fCrossSection[1] = CalculateCrossSection( LeftCellProperties.s[1]/Params.L );
    	fCrossSection[2] = CalculateCrossSection( CellProperties.s[1]/Params.L );
    	fCrossSection[3] = CalculateCrossSection( RightCellProperties.s[1]/Params.L );
    	fCrossSection[4] = CalculateCrossSection( FarRightCellProperties.s[1]/Params.L );
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

#ifdef NON_EQUILIBRIUM_RADIATION
    	ppni0 = FarLeftCellProperties.pIonFrac->ppGetIonFrac();
    	ppni1 = LeftCellProperties.pIonFrac->ppGetIonFrac();
    	ppni2 = CellProperties.pIonFrac->ppGetIonFrac();
    	ppni3 = RightCellProperties.pIonFrac->ppGetIonFrac();
    	ppni4 = FarRightCellProperties.pIonFrac->ppGetIonFrac();
   
    	ppdnibydt = CellProperties.pIonFrac->ppGetdnibydt();
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
    	pRadiation->GetAlldnibydt( log10( CellProperties.T[ELECTRON] ), log10( CellProperties.n[ELECTRON] ), ppni0, ppni1, ppni2, ppni3, ppni4, ps, CellProperties.s, CellProperties.v, fCrossSection, fCellVolume, ppdnibydt, &(CellProperties.atomic_delta_t) );
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
		pRadiation->GetAlldnibydt( log10( CellProperties.T[ELECTRON] ), log10( CellProperties.n[ELECTRON] ), ppni0, ppni1, ppni2, ppni3, ppni4, ps, CellProperties.s, CellProperties.v, CellProperties.cell_width, ppdnibydt, &(CellProperties.atomic_delta_t) );
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
#endif // NON_EQUILIBRIUM_RADIATION

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
   		pHstate0 = FarLeftCellProperties.Hstate;
   		pHstate1 = LeftCellProperties.Hstate;
   		pHstate2 = CellProperties.Hstate;
   		pHstate3 = RightCellProperties.Hstate;
   		pHstate4 = FarRightCellProperties.Hstate;
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
   		pRadiativeRates->GetAllDel_Hstate_dot_v( pHstate0, pHstate1, pHstate2, pHstate3, pHstate4, ps, CellProperties.s, CellProperties.v, fCrossSection, fCellVolume, CellProperties.Del_Hstate_dot_v );
#else // USE_POLY_FIT_TO_MAGNETIC_FIELD
   		pRadiativeRates->GetAllDel_Hstate_dot_v( pHstate0, pHstate1, pHstate2, pHstate3, pHstate4, ps, CellProperties.s, CellProperties.v, CellProperties.cell_width, CellProperties.Del_Hstate_dot_v );
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

#endif // NON_EQUILIBRIUM_RADIATION || ( OPTICALLY_THICK_RADIATION && NLTE_CHROMOSPHERE )

#ifdef OPENMP
		pLocalActiveCell->UpdateCellProperties( &CellProperties );
	}
	}
#else // OPENMP
	    pActiveCell->UpdateCellProperties( &CellProperties );
		pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
	}
#endif // OPENMP

	// Find the smallest characteristic time-scale
	GetSmallestTimeScale( delta_t, iFirstStep );
}

double CEquations::CalculateGravity( double x )
{
	return pGravity->GetPieceWiseFit( x );
}

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
double CEquations::CalculateMagneticField( double x )
{
	return pMagneticField->GetPieceWiseFit( x );
}

double CEquations::CalculateCrossSection( double x )
{
	// The cross-section area varies inversely with the magnetic field strength
	return ( 1.0 / pMagneticField->GetPieceWiseFit( x ) );
}
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

void CEquations::GetSmallestTimeScale( double *delta_t, int iFirstStep )
{
if( !iFirstStep )
    return;
	
PCELL pNextActiveCell;
CELLPROPERTIES CellProperties;
int j;

	*delta_t = Params.Duration;

	pNextActiveCell = pStartOfCurrentRow->pGetPointer( RIGHT )->pGetPointer( RIGHT );

	while( pNextActiveCell->pGetPointer( RIGHT )->pGetPointer( RIGHT ) )
	{
    	pActiveCell = pNextActiveCell;
    	pActiveCell->GetCellProperties( &CellProperties );

	    // Advection timescale
    	if( CellProperties.advection_delta_t < *delta_t )
			*delta_t = CellProperties.advection_delta_t;

	    // Thermal conduction timescale
    	for( j=0; j<SPECIES; j++ )
    	{
			if( CellProperties.conduction_delta_t[j] < *delta_t )
            	*delta_t = CellProperties.conduction_delta_t[j];
    	}

	    // Viscous timescale
    	if( CellProperties.viscosity_delta_t < *delta_t )
			*delta_t = CellProperties.viscosity_delta_t;

	    // Collisional timescale
    	if( CellProperties.collision_delta_t >= MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE && CellProperties.collision_delta_t < *delta_t )
        	*delta_t = CellProperties.collision_delta_t;

	    // Radiation timescale
    	if( CellProperties.radiation_delta_t < *delta_t )
			*delta_t = CellProperties.radiation_delta_t;

#ifdef NON_EQUILIBRIUM_RADIATION
	    // Ionisation and recombination timescale
    	if( CellProperties.atomic_delta_t < *delta_t )
        	*delta_t = CellProperties.atomic_delta_t;
#endif // NON_EQUILIBRIUM_RADIATION

	    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
	}
}

void CEquations::Half_Time_Step( PCELLPROPERTIES pNewCellProperties, double delta_t )
{
CELLPROPERTIES CellProperties;
int j;

pActiveCell->GetCellProperties( &CellProperties );

#ifdef NON_EQUILIBRIUM_RADIATION
	// Integrate the ion fractions
	pNewCellProperties->pIonFrac->IntegrateAllIonFrac( delta_t );
#endif // NON_EQUILIBRIUM_RADIATION

	pNewCellProperties->rho[1] = CellProperties.rho[1] + ( delta_t * CellProperties.drhobydt );
	pNewCellProperties->rho_v[1] = CellProperties.rho_v[1] + ( delta_t * CellProperties.drho_vbydt );

	for( j=0; j<SPECIES; j++ )
    	pNewCellProperties->TE_KE[1][j] = CellProperties.TE_KE[1][j] + ( delta_t * CellProperties.dTE_KEbydt[j] );

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
	double fdne;
	// Now calculate the electron mass density in the new grid cell
	fdne = CellProperties.n[HYDROGEN] * ( delta_t * CellProperties.Del_Hstate_dot_v[5] );
	pNewCellProperties->rho_e = CellProperties.rho_e + ( fdne * ELECTRON_MASS );
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
}

void CEquations::Full_Time_Step( PCELLPROPERTIES pCellProperties, double delta_t )
{
PCELL pBottomCell;
CELLPROPERTIES BottomCellProperties;
int j;

	pBottomCell = pActiveCell->pGetPointer( BOTTOM );
	pBottomCell->GetCellProperties( &BottomCellProperties );

#ifdef NON_EQUILIBRIUM_RADIATION
	// Integrate the ion fractions using the previous set of values
	pCellProperties->pIonFrac->CopyAllIonFrac( BottomCellProperties.pIonFrac );
	pCellProperties->pIonFrac->IntegrateAllIonFrac( delta_t );
#endif // NON_EQUILIBRIUM_RADIATION

	pCellProperties->rho[1] = BottomCellProperties.rho[1] + ( delta_t * pCellProperties->drhobydt );
	pCellProperties->rho_v[1] = BottomCellProperties.rho_v[1] + ( delta_t * pCellProperties->drho_vbydt );

	for( j=0; j<SPECIES; j++ )
    	pCellProperties->TE_KE[1][j] = BottomCellProperties.TE_KE[1][j] + ( delta_t * pCellProperties->dTE_KEbydt[j] );

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
	double fdne;
	// Now calculate the electron mass density in the new grid cell
	fdne = BottomCellProperties.n[HYDROGEN] * ( delta_t * BottomCellProperties.Del_Hstate_dot_v[5] );
	pCellProperties->rho_e = BottomCellProperties.rho_e + ( fdne * ELECTRON_MASS );
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
}

#if defined (OPENMP) || defined (USE_KINETIC_MODEL)
void CEquations::CreateIndexedCellList( void )
{
PCELL pNextActiveCell;
int iIndex;

	// Free the previous indexed cell list and allocate sufficient memory for the next one
	if( ppCellList )
    	free( ppCellList );

	ppCellList = (PCELL*)malloc( Params.iNumberOfCells * sizeof(PCELL) );

	// Compile an indexed list of the cells
	iIndex = 0;
	pNextActiveCell = pStartOfCurrentRow;
	while( pNextActiveCell )
	{
    	pActiveCell = pNextActiveCell;

	    ppCellList[iIndex] = pActiveCell;
    	iIndex++;
	
	    pNextActiveCell=pActiveCell->pGetPointer( RIGHT );
	}
}
#endif // OPENMP || USE_KINETIC_MODEL

// *****************************************************************************
// *                                                                           *
// *    KINETIC COMPONENT FOR DISTRIBUTION FUNCTION CALCULATIONS               *
// *                                                                           *
// *****************************************************************************

#ifdef USE_KINETIC_MODEL
// Index labels for each column of the table in Spitzer & Harm, 1953, Phys. Rev., 89, 977
#define U							0		// Relative (to u_th) velocity values in column 0
#define X_E							1		// Coefficients when an electric field is present in column 1
#define X_T							2		// Coefficients when a temperature gradient is present in column 2

void CEquations::Get_SH_Table( void )
{
FILE *pFile;
int i;

pFile = fopen( "Kinetic_Model/config/spitzer_harm_table.txt" ,"r" );
for( i=0; i<51; i++ )
{
    ReadDouble( pFile, &(SH_Table[i][U]) );
    ReadDouble( pFile, &(SH_Table[i][X_E]) );
    ReadDouble( pFile, &(SH_Table[i][X_T]) );
}
fclose( pFile );
}

void CEquations::CalculateKineticModel( int iFirstStep )
{
if( !iFirstStep )
    return;

PCELL pNextActiveCell;
CELLPROPERTIES CellProperties;
double Tmax = 0.0;

// Calculate the physical quantities required for the kinetic model
pNextActiveCell = pStartOfCurrentRow;

while( pNextActiveCell )
{
    pActiveCell = pNextActiveCell;
    pActiveCell->GetCellProperties( &CellProperties );

    // Find the maximum temperature to determine the velocity range
    if( CellProperties.T[ELECTRON] > Tmax )
		Tmax = CellProperties.T[ELECTRON];

    pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
}

// Calculate the quantities required for the kinetic model
#ifdef OPENMP
	PCELL pLocalActiveCell;
	int iCounter, iLocalNumberOfCells = Params.iNumberOfCells;
#pragma omp parallel shared( ppCellList, iLocalNumberOfCells )		\
   	                 private(	iCounter, CellProperties )
	{
#pragma omp for schedule(dynamic, CHUNK_SIZE)						\
   	            lastprivate(pLocalActiveCell)
	for( iCounter=0; iCounter<iLocalNumberOfCells; iCounter++ )
	{
		pLocalActiveCell = ppCellList[iCounter];
		pLocalActiveCell->GetCellProperties( &CellProperties );
#else // OPENMP
	pNextActiveCell = pStartOfCurrentRow;
	while( pNextActiveCell )
	{
    	pActiveCell = pNextActiveCell;
    	pActiveCell->GetCellProperties( &CellProperties );
#endif // OPENMP

	    CellProperties.pKinetic->CalculateVelocityRange( CellProperties.T[ELECTRON], Tmax );
    	CellProperties.pKinetic->CalculateCollisionalProperties( CellProperties.T[ELECTRON], CellProperties.T[HYDROGEN], CellProperties.n[ELECTRON] );
    	CellProperties.pKinetic->CalculateMaxDFN( CellProperties.T[ELECTRON], CellProperties.T[HYDROGEN], CellProperties.n[ELECTRON] );

#ifdef OPENMP
	}
#else // OPENMP
    	pNextActiveCell = pActiveCell->pGetPointer( RIGHT );
#endif // OPENMP
	}

// Calculate the non-Maxwellian distribution function and its moments in each grid cell
CalculateNonMaxDFN();
}

void CEquations::CalculateNonMaxDFN( void )
{
// General variables
double *pupsilon, *pnu_ee, *pnu_ei, *plambda_ei, *pMaxDFN_ee, *pMaxDFN_ei, *pNonMaxDFN;
double dTbyds, K_SH, Kappa;
int iIndex, i, j, half_data_points = DISTRIBUTION_DATA_POINTS / 2;

// Variables for the BGK part of the solution
CELLPROPERTIES CellProperties;
double ds, previous_cell_width, *previous_pNonMaxDFN;
double term1, term2, term3;

// Variables for the Spitzer-Harm part of the solution
CELLPROPERTIES FarLeftCellProperties, LeftCellProperties, RightCellProperties, FarRightCellProperties;
double u_th, u_n;
double fScaleLength, E, lambda_ei, x_e, x_t, fp;
double x[6], y[6], fLowerValue, fUpperValue, fError;

// *****************************************************************************
// *    INITIALISE THE BOUNDARY CONDITIONS                                      													*
// *****************************************************************************

iIndex = 0;
pActiveCell = ppCellList[iIndex];
pActiveCell->GetCellProperties( &CellProperties );
pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();
for( i=0; i<DISTRIBUTION_DATA_POINTS; i++ )
    pNonMaxDFN[i] = pMaxDFN_ee[i];

iIndex = 1;
pActiveCell = ppCellList[iIndex];
pActiveCell->GetCellProperties( &CellProperties );
pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();
for( i=0; i<DISTRIBUTION_DATA_POINTS; i++ )
    pNonMaxDFN[i] = pMaxDFN_ee[i];

iIndex = Params.iNumberOfCells-2;
pActiveCell = ppCellList[iIndex];
pActiveCell->GetCellProperties( &CellProperties );
pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();
for( i=0; i<DISTRIBUTION_DATA_POINTS; i++ )
    pNonMaxDFN[i] = pMaxDFN_ee[i];

iIndex = Params.iNumberOfCells-1;
pActiveCell = ppCellList[iIndex];
pActiveCell->GetCellProperties( &CellProperties );
pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();
for( i=0; i<DISTRIBUTION_DATA_POINTS; i++ )
    pNonMaxDFN[i] = pMaxDFN_ee[i];

// *****************************************************************************
// *    BGK PART OF THE SOLUTION                                                																*
// *****************************************************************************

// Solve for positive velocities
iIndex = 1;
pActiveCell = ppCellList[iIndex];
pActiveCell->GetCellProperties( &CellProperties );
pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();

for( iIndex=2; iIndex<Params.iNumberOfCells; iIndex++ )
{
   	previous_cell_width = CellProperties.cell_width;
   	previous_pNonMaxDFN = pNonMaxDFN;

    pActiveCell = ppCellList[iIndex];
   	pActiveCell->GetCellProperties( &CellProperties );
    pupsilon = CellProperties.pKinetic->Get_pupsilon();
    pnu_ee = CellProperties.pKinetic->Get_pnu_ee();
    pnu_ei = CellProperties.pKinetic->Get_pnu_ei();
    pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
    pMaxDFN_ei = CellProperties.pKinetic->Get_pMaxDFN_ei();
    pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();

    ds = 0.5 * ( previous_cell_width + CellProperties.cell_width );

    for( i=half_data_points; i<DISTRIBUTION_DATA_POINTS; i++ )
    {
		term1 = ( ( pnu_ee[i] * pMaxDFN_ee[i] ) + ( 4.0 * pnu_ei[i] * pMaxDFN_ei[i] ) ) * ds;
		term2 = pupsilon[i] * previous_pNonMaxDFN[i];
		term3 = pupsilon[i] + ( ( pnu_ee[i] + ( 4.0 * pnu_ei[i] ) ) * ds );

		pNonMaxDFN[i] = ( term1 + term2 ) / term3;
   	}
}

// Solve for negative velocities

iIndex = Params.iNumberOfCells-2;
pActiveCell = ppCellList[iIndex];
pActiveCell->GetCellProperties( &CellProperties );
pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();

for( iIndex=Params.iNumberOfCells-3; iIndex>=0; iIndex-- )
{
    previous_cell_width = CellProperties.cell_width;
    previous_pNonMaxDFN = pNonMaxDFN;

    pActiveCell = ppCellList[iIndex];
    pActiveCell->GetCellProperties( &CellProperties );
    pupsilon = CellProperties.pKinetic->Get_pupsilon();
    pnu_ee = CellProperties.pKinetic->Get_pnu_ee();
    pnu_ei = CellProperties.pKinetic->Get_pnu_ei();
    pMaxDFN_ee = CellProperties.pKinetic->Get_pMaxDFN_ee();
    pMaxDFN_ei = CellProperties.pKinetic->Get_pMaxDFN_ei();
    pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();

    ds = 0.5 * ( previous_cell_width + CellProperties.cell_width );

    for( i=0; i<half_data_points; i++ )
    {
        term1 = ( ( pnu_ee[i] * pMaxDFN_ee[i] ) + ( 4.0 * pnu_ei[i] * pMaxDFN_ei[i] ) ) * ds;
		term2 = pupsilon[i] * previous_pNonMaxDFN[i];
		term3 = pupsilon[i] + ( ( pnu_ee[i] + ( 4.0 * pnu_ei[i] ) ) * ds );

		pNonMaxDFN[i] = ( term1 + term2 ) / term3;
    }
}

// *****************************************************************************
// *    SPITZER PART OF THE SOLUTION                                            															*
// *****************************************************************************

#ifdef OPENMP
	PCELL pLocalActiveCell;
	int iCounter, iLocalNumberOfCells = Params.iNumberOfCells;
#pragma omp parallel shared( ppCellList, iLocalNumberOfCells )					\
   	                 private(	iCounter, CellProperties,						\
					 			LeftCellProperties, FarLeftCellProperties, 		\
								RightCellProperties, FarRightCellProperties,	\
					 			pupsilon, plambda_ei, pNonMaxDFN, u_th,			\
								x, y, fLowerValue, fUpperValue, fError, 		\
								dTbyds, fScaleLength, E, lambda_ei,				\
								i, j, u_n, x_e, x_t, term1, term2, fp,			\
								K_SH, Kappa )
{
#pragma omp for schedule(dynamic, CHUNK_SIZE)									\
   	            lastprivate(pLocalActiveCell)
	for( iCounter=2; iCounter<iLocalNumberOfCells-2; iCounter++ )
	{
		pLocalActiveCell = ppCellList[iCounter];
		pLocalActiveCell->GetCellProperties( &CellProperties );

	    ppCellList[iCounter-2]->GetCellProperties( &FarLeftCellProperties );
    	ppCellList[iCounter-1]->GetCellProperties( &LeftCellProperties );
    	ppCellList[iCounter+1]->GetCellProperties( &RightCellProperties );
    	ppCellList[iCounter+2]->GetCellProperties( &FarRightCellProperties );
#else // OPENMP
	for( iIndex=2; iIndex<Params.iNumberOfCells-2; iIndex++ )
	{
    	pActiveCell = ppCellList[iIndex];
    	pActiveCell->GetCellProperties( &CellProperties );

	    ppCellList[iIndex-2]->GetCellProperties( &FarLeftCellProperties );
    	ppCellList[iIndex-1]->GetCellProperties( &LeftCellProperties );
    	ppCellList[iIndex+1]->GetCellProperties( &RightCellProperties );
    	ppCellList[iIndex+2]->GetCellProperties( &FarRightCellProperties );
#endif // OPENMP

	    pupsilon = CellProperties.pKinetic->Get_pupsilon();
	    plambda_ei = CellProperties.pKinetic->Get_plambda_ei();
	    pNonMaxDFN = CellProperties.pKinetic->Get_pNonMaxDFN();

    	u_th = CellProperties.pKinetic->Get_u_th();

	    x[1] = FarLeftCellProperties.s[1];
    	x[2] = LeftCellProperties.s[1];
    	x[3] = CellProperties.s[1];
    	x[4] = RightCellProperties.s[1];
    	x[5] = FarRightCellProperties.s[1];

	    y[1] = FarLeftCellProperties.T[ELECTRON];
	    y[2] = LeftCellProperties.T[ELECTRON];
	    y[3] = CellProperties.T[ELECTRON];
	    y[4] = RightCellProperties.T[ELECTRON];
	    y[5] = FarRightCellProperties.T[ELECTRON];

	    FitPolynomial4( x, y, CellProperties.s[0], &fLowerValue, &fError );
    	FitPolynomial4( &(x[1]), &(y[1]), CellProperties.s[2], &fUpperValue, &fError );

	    dTbyds = ( fUpperValue - fLowerValue ) / CellProperties.cell_width;
    	fScaleLength = KNUDSEN_NUMBER * fabs( CellProperties.T[ELECTRON] / dTbyds );
    	E = -0.703 * ( BOLTZMANN_CONSTANT / ELECTRON_CHARGE ) * dTbyds;

	    lambda_ei = fSHMeanFreePath( u_th, CellProperties.T[ELECTRON], CellProperties.T[HYDROGEN], CellProperties.n[ELECTRON] );

	    // Solve for positive velocities

    	for( i=half_data_points; i<DISTRIBUTION_DATA_POINTS; i++ )
		{
    	    // Normalise the velocity to the thermal speed
			u_n = pupsilon[i] / u_th;

			if( plambda_ei[i] > fScaleLength || u_n > 3.20 )
	            break;

			// Calculate X_E and X_T
        	if( u_n > 0.10 )
			{
            	// Find the data points to interpolate at the current velocity
            	for( j=0; j<51; j++ )
            	{
                	if( u_n <= SH_Table[j][U] )
                    	break;
            	}

	            x[1] = SH_Table[j-1][U];
    	        x[2] = SH_Table[j][U];

        	    y[1] = SH_Table[j-1][X_E];
            	y[2] = SH_Table[j][X_E];
            	LinearFit( x, y, u_n, &x_e );

	            y[1] = SH_Table[j-1][X_T];
    	        y[2] = SH_Table[j][X_T];
        	    LinearFit( x, y, u_n, &x_t );
			}
			else
			{
            	x_e = 0.0;
            	x_t = 0.0;
			}

			term1 = -x_e * ( ( ELECTRON_CHARGE * E ) / ( BOLTZMANN_CONSTANT * CellProperties.T[ELECTRON] ) );
			term2 = 2.0 * x_t * ( 1.0 / CellProperties.T[ELECTRON] ) * dTbyds; 
			fp = lambda_ei * ( term1 + term2 );

			pNonMaxDFN[i] *= ( 1.0 + fp );
    	}

	    // Solve for negative velocities

    	for( i=half_data_points-1; i>=0; i-- )
    	{
        	// Normalise the velocity to the thermal speed
			u_n = -pupsilon[i] / u_th;

			if( plambda_ei[i] > fScaleLength || u_n > 3.20 )
    	        break;

			// Calculate X_E and X_T
			if( u_n > 0.10 )
			{
            	// Find the data points to interpolate at the current velocity
            	for( j=0; j<51; j++ )
            	{
                	if( u_n <= SH_Table[j][U] )
                    	break;
            	}

	            x[1] = SH_Table[j-1][U];
    	        x[2] = SH_Table[j][U];

	            y[1] = SH_Table[j-1][X_E];
    	        y[2] = SH_Table[j][X_E];
        	    LinearFit( x, y, u_n, &x_e );

	            y[1] = SH_Table[j-1][X_T];
    	        y[2] = SH_Table[j][X_T];
        	    LinearFit( x, y, u_n, &x_t );
			}
			else
			{
    	        x_e = 0.0;
        	    x_t = 0.0;
			}

			term1 = -x_e * ( ( ELECTRON_CHARGE * E ) / ( BOLTZMANN_CONSTANT * CellProperties.T[ELECTRON] ) );
			term2 = 2.0 * x_t * ( 1.0 / CellProperties.T[ELECTRON] ) * dTbyds; 
			fp = - lambda_ei * ( term1 + term2 );

			pNonMaxDFN[i] *= ( 1.0 + fp );
    	}

	    CellProperties.Fc[1][ELECTRON] = CellProperties.pKinetic->Get_Fc();

	    // Calculate the thermal conduction time-scale
    	K_SH = SPITZER_ELECTRON_CONDUCTIVITY * pow( CellProperties.T[ELECTRON], 2.5 );
    	Kappa = fabs( CellProperties.Fc[1][ELECTRON] / dTbyds );
	
	    if( Kappa > K_SH )
    	    Kappa = K_SH;

	    // The time-step is calculated by: delta_t = ( n * k_B ) * ( cell width )^2 / ( 2.0 * coefficient of thermal conduction )
    	// 0.5 * BOLTZMANN_CONSTANT = 6.9e-17
    	term1 = SAFETY_CONDUCTION * (6.9e-17) * CellProperties.cell_width * CellProperties.cell_width;
    	CellProperties.conduction_delta_t[ELECTRON] = ( term1 * CellProperties.n[ELECTRON] ) / Kappa;

#ifdef OPENMP
		pLocalActiveCell->UpdateCellProperties( &CellProperties );
#else // OPENMP
	    pActiveCell->UpdateCellProperties( &CellProperties );
#endif // OPENMP
	}
#ifdef OPENMP
}
#endif // OPENMP
}
#endif // USE_KINETIC_MODEL