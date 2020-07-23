// ****
// *
// * Include file for miscellaneous routines used by the hydrostatic code
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 07/20/2020
// *
// ****

void GetConfigurationParameters( PARAMETERS *pParams );

#ifdef USE_POLY_FIT_TO_GRAVITY
#else // USE_POLY_FIT_TO_GRAVITY
#define GRAVITY_POLY_ORDER		6
#define SUB_DOMAIN_STRUCTURE	"1\n0.0\t1.0\n"
void GenerateDefaultLoop( PARAMETERS Params );
#endif // USE_POLY_FIT_TO_GRAVITY

void WriteAMRFile( int iTotalSteps, double *s, double *T, double *nH, double *ne, PARAMETERS Params );
void WritePHYFile( int iTotalSteps, double *s, double *T, double *nH, double *ne, PARAMETERS Params );
void WriteSOLFile( double finalH0, PARAMETERS Params );
