// ****
// *
// * A code to synthesise the emission detectable by a set
// * of observing instruments, given a set of response functions
// * and emitting ions
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 11/13/2015
// *
// ****


#include <stdio.h>

#include "forward.h"


int main( void )
{
PLOOP pLoop;
PFORWARD pForwardModel;
int i;

printf( "\n" );

pForwardModel = new CForward( "Forward_Model/config/forward_model.cfg" );
pForwardModel->WriteIonLineListToFile( "Results/Forward_Model/LINELIST.txt" );

// Forward model the first 600 seconds of evolution at a cadence of 5 seconds
for( i=0; i<596; i+=5 )
{
	printf( "Creating loop:\n\n" );
	pLoop = new CLoop( "Results", i, i+4, 1, 1E8 );

	pForwardModel->CreateVirtualDetector( pLoop );
	pForwardModel->ForwardModel( pLoop );

	delete pLoop;

	printf( "Writing detector files:\n\n" );
	pForwardModel->WriteDetectorToFile( "Results/Forward_Model/", i );

	pForwardModel->ResetVirtualDetector();
}

// printf( "\nWriting EM loci plots:\n\n" );
// pForwardModel->EM_Loci( 4.0, 8.0, 0.1, 10.0, "Results/Forward_Model/EM_Loci/" );

delete pForwardModel;

return 0;
}