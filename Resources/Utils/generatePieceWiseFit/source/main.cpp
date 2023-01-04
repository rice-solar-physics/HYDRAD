// ****
// *
// * A utility to calculate piece-wise polynomial fits to tabulated data
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 01/04/2023
// *
// ****

#include <stdio.h>
#include <stdlib.h>

#include "piecewisefit.h"

int main( int argc, char **argv )
{
	PPIECEWISEFIT pMagneticField;

	if( argc == 1 ) {
		printf( "\nEither:\n" );
		printf( "\t1. Input (data) and output (fit) files must be specified. E.g. generatePieceWiseFit fits/B(closed).txt fits/initial.amr.B\n" );
		printf( "\t2. An existing fit file must be specified. E.g. generatePieceWiseFit fits/initial.amr.B\n");
		exit( EXIT_SUCCESS );
	} else if( argc == 2 ) {
		pMagneticField = new CPieceWiseFit( argv[1] );
	} else {
		pMagneticField = new CPieceWiseFit( argv[1], argv[2] );
	}

	// Example useage

#ifdef VERBOSE
	pMagneticField->ShowPieceWiseFit();
#endif // VERBOSE
	pMagneticField->GetPieceWiseFit( 0.25 );
	pMagneticField->GetPieceWiseFit( 0.5 );
	pMagneticField->GetPieceWiseFit( 0.75 );

	pMagneticField->GetDerivative( 0.25, 1 );
	pMagneticField->GetDerivative( 0.5, 1 );
	pMagneticField->GetDerivative( 0.75, 1 );

	delete pMagneticField;

	return 0;
}