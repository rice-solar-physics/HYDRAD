// ****
// *
// * Include file for miscellaneous routines used by the hydrostatic code
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 09/24/2019
// *
// ****

void GetConfigurationParameters( PARAMETERS *pParams );

#ifdef USE_POLY_FIT_TO_GRAVITY
#else // USE_POLY_FIT_TO_GRAVITY
void GenerateDefaultLoop( PARAMETERS Params );
#endif // USE_POLY_FIT_TO_GRAVITY

void WriteAMRFile( int iTotalSteps, double *s, double *T, double *nH, double *ne, PARAMETERS Params );
void WritePHYFile( int iTotalSteps, double *s, double *T, double *nH, double *ne, PARAMETERS Params );
void WriteSOLFile( double finalH0, PARAMETERS Params );
