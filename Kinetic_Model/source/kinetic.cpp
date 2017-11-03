// ****
// *
// * Code to calculate the electron distribution function, using the Spitzer-Harm method for
// * the bulk of the distribution function and the BGK approximation for the tail,
// * given a spatial temperature and density distribution in one dimension
// *
// * Class function bodies
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/14/2017
// *
// ****

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "kinetic.h"
#include "gamma.h"
#include "../../Resources/source/constants.h"
#include "../../Resources/source/file.h"
#include "../../Resources/source/fitpoly.h"


// The Coulomb logarithm for electron-electron collisions
// (The Coulomb logarithms for electron-ion and ion-ion collisions are already defined in eqns.cpp)

double fLogLambda_ee( double Te, double n )
{
    Te *= 1.29199251618e-4;

    if( Te < 10.0 )
    {
        return 23.0 - log( sqrt( n ) * pow( Te, (-1.5) ) );
    }
    else
    {
    	return 24.0 - log( sqrt( n ) * ( 1.0 / Te ) );
    }
}

// This function calculates the electron collisional mean free path at near thermal velocities
double fSHMeanFreePath( double u_th, double Te, double Ti, double n )
{
double upsilon2, upsilon3;
double term1_ei;
double nu_ei;

// Numerical values used to calculate nu_ei
// ( 4.0 * _PI_ * ELECTRON_CHARGE_POWER_4 ) / ELECTRON_MASS_SQUARED = 8.037815988653620e17
// 2.0 / sqrt( _PI_ ) = 1.12837916709551
// 8.037815988653620e17 * 1.12837916709551 = 9.069704110543940e17

upsilon2 = u_th * u_th;
upsilon3 = upsilon2 * u_th;

#ifdef USE_APPROXIMATE_NU_EI
term1_ei = (8.037815988653620e17) * n * fLogLambda_ei( Te, Ti, n );

nu_ei = term1_ei / upsilon3;
#else // USE_APPROXIMATE_NU_EI
double term2_ei, x_ei;
// 2.0 * BOLTZMANN_CONSTANT = 2.76e-16
term1_ei = ( AVERAGE_PARTICLE_MASS / (2.76e-16) ) / Ti;
term2_ei = (9.069704110543940e17) * ( 1.0 + ( ELECTRON_MASS / AVERAGE_PARTICLE_MASS ) ) * n * fLogLambda_ei( Te, Ti, n );

x_ei = term1_ei * upsilon2;
nu_ei = ( term2_ei / upsilon3 ) * gammp( 1.5, x_ei );
#endif // USE_APPROXIMATE_NU_EI

return ( u_th / nu_ei );
}


CKinetic::CKinetic( void )
{
Initialise();
}

CKinetic::~CKinetic( void )
{
FreeAll();
}

void CKinetic::Initialise( void )
{
pupsilon = (double*)malloc( sizeof(double) * DISTRIBUTION_DATA_POINTS );
pnu_ee = (double*)malloc( sizeof(double) * DISTRIBUTION_DATA_POINTS );
pnu_ei = (double*)malloc( sizeof(double) * DISTRIBUTION_DATA_POINTS );
plambda_ei = (double*)malloc( sizeof(double) * DISTRIBUTION_DATA_POINTS );
pMaxDFN_ee = (double*)malloc( sizeof(double) * DISTRIBUTION_DATA_POINTS );
pMaxDFN_ei = (double*)malloc( sizeof(double) * DISTRIBUTION_DATA_POINTS );
pNonMaxDFN = (double*)malloc( sizeof(double) * DISTRIBUTION_DATA_POINTS );
}

void CKinetic::FreeAll( void )
{
free( pupsilon );
free( pnu_ee );
free( pnu_ei );
free( plambda_ei );
free( pMaxDFN_ee );
free( pMaxDFN_ei );
free( pNonMaxDFN );
}

double CKinetic::CalculateThermalSpeed( double T )
{
// ( 2.0 * BOLTZMANN_CONSTANT ) / ELECTRON_MASS = 3.029637760702530e11
return sqrt( (3.029637760702530e11) * T );
}

void CKinetic::CalculateVelocityRange( double T, double Tmax )
{
double fv, fdv, u_th_max;
int i, j, half_data_points = ( DISTRIBUTION_DATA_POINTS / 2 );

u_th = CalculateThermalSpeed( T );

fdv = V_RANGE_UPPER / half_data_points;
fv = fdv;

u_th_max = CalculateThermalSpeed( Tmax );

for( i=half_data_points-1; i>=0; i-- )
{
    pupsilon[i] = - fv * u_th_max;
    fv += fdv;
}

j = half_data_points - 1;
for( i=half_data_points; i<DISTRIBUTION_DATA_POINTS; i++ )
{
    pupsilon[i] = -pupsilon[j];
    j--;
}
}

// This function calculates the collision frequencies and mean-free-paths for the current grid cell, using the equations given on page 31 of the NRL Plasma Formulary (revised 1984)
void CKinetic::CalculateCollisionalProperties( double Te, double Ti, double n )
{
double term1_ee, term2_ee, term1_ei;
double upsilon2, upsilon3;
double x_ee;
int i, j, half_data_points = ( DISTRIBUTION_DATA_POINTS / 2 );

// Numerical values used to calculate nu_ee and nu_ei
// ( 4.0 * _PI_ * ELECTRON_CHARGE_POWER_4 ) / ELECTRON_MASS_SQUARED = 8.037815988653620e17
// 2.0 / sqrt( _PI_ ) = 1.12837916709551
// 8.037815988653620e17 * 1.12837916709551 = 9.069704110543940e17

// Values used to calculate nu_ee
// ELECTRON_MASS / ( 2.0 * BOLTZMANN_CONSTANT ) = 3.30072463768116e-12
term1_ee = (3.30072463768116e-12) / Te;
// 1.0 + ( ELECTRON_MASS / ELECTRON_MASS ) = 2.0
// 9.069704110543940e17 * 2.0 = 1.813940822108790e18
term2_ee = (1.813940822108790e18) * n * fLogLambda_ee( Te, n );

// Values used to calculate nu_ei
#ifdef USE_APPROXIMATE_NU_EI
term1_ei = (8.037815988653620e17) * n * fLogLambda_ei( Te, Ti, n );
#else // USE_APPROXIMATE_NU_EI
double term2_ei, x_ei;
// 2.0 * BOLTZMANN_CONSTANT = 2.76e-16
term1_ei = ( AVERAGE_PARTICLE_MASS / (2.76e-16) ) / Ti;
term2_ei = (9.069704110543940e17) * ( 1.0 + ( ELECTRON_MASS / AVERAGE_PARTICLE_MASS ) ) * n * fLogLambda_ei( Te, Ti, n );
#endif // USE_APPROXIMATE_NU_EI

for( i=0; i<half_data_points; i++ )
{
    upsilon2 = pupsilon[i] * pupsilon[i];
    upsilon3 = upsilon2 * pupsilon[i];

    // Calculate nu_ee
    x_ee = term1_ee * upsilon2;
    pnu_ee[i] = ( term2_ee / upsilon3 ) * gammp( 1.5, x_ee );
		
    // Calculate nu_ei
#ifdef USE_APPROXIMATE_NU_EI
    pnu_ei[i] = term1_ei / upsilon3;
#else // USE_APPROXIMATE_NU_EI
    x_ei = term1_ei * upsilon2;
    pnu_ei[i] = ( term2_ei / upsilon3 ) * gammp( 1.5, x_ei );
#endif // USE_APPROXIMATE_NU_EI

    plambda_ei[i] = pupsilon[i] / pnu_ei[i];
}

j = half_data_points - 1;
for( i=half_data_points; i<DISTRIBUTION_DATA_POINTS; i++ )
{
    pnu_ee[i] = -pnu_ee[j];
    pnu_ei[i] = -pnu_ei[j];
    plambda_ei[i] = plambda_ei[j];
    j--;
}
}

void CKinetic::CalculateMaxDFN( double Te, double Ti, double n )
{
int i;

double term1e, term1i, term2e, term2i, term3, term4e, term4i;

// Simplified Greene, 1973, equation (16)
// For details of simplification procedure refer to research diary entry: Wednesday 31st January 2007 (p18)
double alpha_E, Te_bar;
double temp0, temp1, temp2;
// Omit extra factor of n so that alpha_E is in units of [T]^-1
temp0 = n * fLogLambda_ei( Te, Ti, n );
temp1 = pow( ( ( ELECTRON_MASS * Ti ) + ( AVERAGE_PARTICLE_MASS * Te ) ), (1.5) );
alpha_E = (6.16e-41) * ( temp0 / temp1 );

temp2 = alpha_E * ( Ti - Te );

// ELECTRON_MASS / ( 2.0 * BOLTZMANN_CONSTANT ) = 3.30072463768116e-12
term1e = (3.30072463768116e-12) / Te;
term2e = pow( ( term1e / _PI_ ), (1.5) );
// End of Simplified Greene, 1973, equation (16)

for( i=0; i<DISTRIBUTION_DATA_POINTS; i++ )
{
    // Simplified Greene, 1973, equation (16)
    // For details of simplification procedure refer to research diary entry: Wednesday 31st January 2007 (p18)
    Te_bar = Te + ( temp2 / fabs( pnu_ei[i] ) );
    if( Te_bar < 0.0 ) Te_bar = Ti; // Seems only necessary for highest velocity (least influential) particles
    // ELECTRON_MASS / ( 2.0 * BOLTZMANN_CONSTANT ) = 3.30072463768116e-12
    term1i = (3.30072463768116e-12) / Te_bar;
    term2i = pow( ( term1i / _PI_ ), (1.5) );
    // End of Simplified Greene, 1973, equation (16)
	
    term3 = pupsilon[i];
    term3 *= term3;

    term4e = term1e * term3;
    term4i = term1i * term3;

    pMaxDFN_ee[i] = n * term2e * exp( -term4e );
    if( pMaxDFN_ee[i] < F_MIN )
        pMaxDFN_ee[i] = F_MIN;

    pMaxDFN_ei[i] = n * term2i * exp( -term4i );
    if( pMaxDFN_ei[i] < F_MIN )
        pMaxDFN_ei[i] = F_MIN;
}
}

double CKinetic::Get_u_th( void )
{
return u_th;
}

double* CKinetic::Get_pupsilon( void )
{
return pupsilon;
}

double* CKinetic::Get_pnu_ee( void )
{
return pnu_ee;
}

double* CKinetic::Get_pnu_ei( void )
{
return pnu_ei;
}

double* CKinetic::Get_plambda_ei( void )
{
return plambda_ei;
}

double* CKinetic::Get_pMaxDFN_ee( void )
{
return pMaxDFN_ee;
}

double* CKinetic::Get_pMaxDFN_ei( void )
{
return pMaxDFN_ei;
}

double* CKinetic::Get_pNonMaxDFN( void )
{
return pNonMaxDFN;
}

double CKinetic::Get_rho( void )
{
int j;
double fIntegral = 0.0, term1, term2;

// Use the Trapezium rule to perform the integration
term1 = pNonMaxDFN[0] * pupsilon[0] * pupsilon[0];
for( j=1; j<DISTRIBUTION_DATA_POINTS; j++ )
{
    term2 = pNonMaxDFN[j] * pupsilon[j] * pupsilon[j];
    fIntegral += ( 0.5 * ( pupsilon[j] - pupsilon[j-1] ) * ( term1 + term2 ) );
    term1 = term2;
}

// 2.0 * _PI_ = 6.283185307179586
return (6.283185307179586) * AVERAGE_PARTICLE_MASS * fIntegral;
}

double CKinetic::Get_TE_KE( void )
{
int j;
double fIntegral = 0.0, term1, term2;

// Use the Trapezium rule to perform the integration
term1 = pNonMaxDFN[0] * pupsilon[0] * pupsilon[0] * pupsilon[0] * pupsilon[0];
for( j=1; j<DISTRIBUTION_DATA_POINTS; j++ )
{
    term2 = pNonMaxDFN[j] * pupsilon[j] * pupsilon[j] * pupsilon[j] * pupsilon[j];
    fIntegral += ( 0.5 * ( pupsilon[j] - pupsilon[j-1] ) * ( term1 + term2 ) );
    term1 = term2;
}

// _PI_ * ELECTRON_MASS = 2.86199090742030e-27
return (2.86199090742030e-27) * fIntegral;
}

double CKinetic::Get_Fc( void )
{
int j;
double fIntegral = 0.0, term1, term2;

// Use the Trapezium rule to perform the integration
term1 = pNonMaxDFN[0] * pupsilon[0] * pupsilon[0] * pupsilon[0] * pupsilon[0] * pupsilon[0];
for( j=1; j<DISTRIBUTION_DATA_POINTS; j++ )
{
    term2 = pNonMaxDFN[j] * pupsilon[j] * pupsilon[j] * pupsilon[j] * pupsilon[j] * pupsilon[j];
    fIntegral += ( 0.5 * ( pupsilon[j] - pupsilon[j-1] ) * ( term1 + term2 ) );
    term1 = term2;
}

// ( _PI_ / 3.0 ) * ELECTRON_MASS = 9.53996969140101e-28
return (9.53996969140101e-28) * fIntegral;
}
