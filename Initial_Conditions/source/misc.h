// ****
// *
// * Include file for miscellaneous routines used by the hydrostatic code
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 29/12/2012
// *
// ****

void GetConfigurationParameters( PARAMETERS *pParams );

#ifdef USE_TABULATED_GRAVITY
#else // USE_TABULATED_GRAVITY
void GenerateSemiCircularLoop( PARAMETERS Params );
#endif // USE_TABULATED_GRAVITY

void WriteAMRFile( int iTotalSteps, double *s, double *T, double *nH, double *ne, PARAMETERS Params );
void WritePHYFile( int iTotalSteps, double *s, double *T, double *nH, double *ne, PARAMETERS Params );
void WriteSOLFile( double finalH0, PARAMETERS Params );
