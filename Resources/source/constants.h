// ****
// *
// * Physical constants
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 31/12/2012
// *
// ****


// All constants are given in cgs units

#define _PI_					3.1415926535897932384626433832795
#define SQRT_PI                                 1.7724538509055160272981674833411
#define GAMMA					1.666666667                         // Adiabatic index
#define GAMMA_MINUS_ONE				0.666666667
#define SPEED_OF_LIGHT                          3e10                                // (cm)(s^-1)

#define SPITZER_ELECTRON_CONDUCTIVITY  		7.8e-7                              // (erg) (cm^-1) (s^-1) (K^-7/2)
#define SPITZER_ION_CONDUCTIVITY		3.2e-8                              // (erg) (cm^-1) (s^-1) (K^-7/2)
// #define SPITZER_SINGLE_FLUID_CONDUCTIVITY	9.2e-7                              // (erg) (cm^-1) (s^-1) (K^-7/2)
#define SPITZER_SINGLE_FLUID_CONDUCTIVITY	8.12e-7                             // (erg) (cm^-1) (s^-1) (K^-7/2)

#define	DYNAMIC_VISCOSITY			2.5224364780829600e-15              // (g) (cm^-1) (s^-1)
                                                                                    // Also denoted eta_0

#define BOLTZMANN_CONSTANT    			1.38e-16                            // (erg) (K^-1)

#define ELECTRON_MASS				9.11e-28                            // (g)
#define SQRT_ELECTRON_MASS			3.018e-14                           // (g^0.5)
#define ELECTRON_MASS_SQUARED			8.29921e-55                         // (g^2)
#define PROTON_MASS                             1.6726231e-24                       // (g)
#define NEUTRON_MASS                            1.6749286e-24                       // (g)
#define AVERAGE_PARTICLE_MASS  			2.171e-24                           // (g)
#define SQRT_AVERAGE_PARTICLE_MASS		1.473e-12                           // (g^0.5)

#define ELECTRON_CHARGE				4.80e-10                            // (esu)
#define ELECTRON_CHARGE_POWER_4			5.308416e-38                        // (esu^4)

#define SOLAR_SURFACE_GRAVITY 			27400                               // (cm) (s^-2)
#define SOLAR_RADIUS          			6.96e10                             // (cm)
#define SOLAR_RADIUS_SQUARED			4.84e21                             // (cm^2)

#define LARGEST_DOUBLE				1e308                               // The largest number a double precision type can represent
#define SMALLEST_DOUBLE                         1e-300                              // The smallest number a double precision type can represent