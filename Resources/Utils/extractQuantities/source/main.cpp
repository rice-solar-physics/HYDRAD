// ****
// *
// * A routine to extract particular sets of quantities from a series of data files in a
// * user-specified foirmat and write them into new files for ease of manipulation and plotting
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 01/04/2023
// *
// ****

#include <stdio.h>
#include <stdlib.h>

#include "../../../source/file.h"

int main( int argc, char **argv )
{
FILE *pCONFIGFile, *pINPUTFile, *pOUTPUTFile;
char **ppcDataType;
char szResultsDirectory[256], szRoot[256], szExtension[32], szINPUTFilename[256], szOUTPUTFilename[256];
char cBuffer;
double fBuffer;
int *piNumColumns, *piRow, *piColumn, iRange[2];
int iNumRows, iNumEntries, iEntry, iNumRecords;
int i, j, m, n;
int iBuffer;

if( argc == 1 ) {
	printf( "\nA configuration file must be specified. E.g. extractQuantities config.cfg\n");
	exit( EXIT_SUCCESS );
}

// Open the configuration file
pCONFIGFile = fopen( argv[1], "r" );
	// Get the directory containing the numerical results
	fscanf( pCONFIGFile, "%s", szResultsDirectory );
	// Get the filename structure of the files containing the data to be extracted
	fscanf(pCONFIGFile, "%s", szRoot );
	fscanf(pCONFIGFile, "%s", szExtension );
	// Get the number of rows per record in the files
	fscanf(pCONFIGFile, "%i", &iNumRows );
		// Allocate memory to store the number of columns per row and the data types of each entry
		piNumColumns = (int*)malloc( sizeof(int) * iNumRows );
		ppcDataType = (char**)malloc( sizeof(char*) * iNumRows );
		for( i=0; i<iNumRows; i++ ) {
			// Get the number of columns in the current row
			fscanf(pCONFIGFile, "%i", &(piNumColumns[i]) );
			// Allocate memory to store the data type of each entry in the current row
			ppcDataType[i] = (char*)malloc( sizeof(char) * piNumColumns[i] );
				for( j=0; j<piNumColumns[i]; j++ ) {
					// Get the data type of each entry in the current row
					fscanf( pCONFIGFile, "%s", &(ppcDataType[i][j]) );
				}
		}
	// Get the number of entries to extract from each record
	fscanf( pCONFIGFile, "%i", &iNumEntries );
		// Allocate memory to store the row and column number for each entry
		piRow = (int*)malloc( sizeof(int) * iNumEntries );
		piColumn = (int*)malloc( sizeof(int) * iNumEntries );	
		for( i=0; i<iNumEntries; i++ ) {
			// Get the row and column number for each entry
			fscanf( pCONFIGFile, "%i", &(piRow[i]) );
			fscanf( pCONFIGFile, "%c", &cBuffer );	// Get the comma separator
			fscanf( pCONFIGFile, "%i", &(piColumn[i]) );
		}
	fscanf( pCONFIGFile, "%i", &(iRange[0]) );
	fscanf( pCONFIGFile, "%i", &(iRange[1]) );
// Close the configuration file
fclose( pCONFIGFile );

printf( "\n%s/%s.%s\n", szResultsDirectory, szRoot, szExtension );
printf( "\niNumber of Rows = %i\n", iNumRows );
for( i=0; i<iNumRows; i++ ) {
	printf( "Number of Columns in Row %i = %i\n", i+1, piNumColumns[i] );
	printf( "\tColumn data types:" );
	for( j=0; j<piNumColumns[i]; j++ ) {
		printf( " %c", ppcDataType[i][j] );
	}
	printf( "\n" );
}
printf( "\nNumber of Entries = %i\n", iNumEntries );
printf( "\tEntries:" );
for( i=0; i<iNumEntries; i++ ) {
	printf( " [%i,%i]", piRow[i], piColumn[i] );
}
printf( "\n" );
printf( "\nFile Number Range: %i to %i\n", iRange[0], iRange[1] );

// Extract the entries from the range of files specified and write them into a new file
for( i=iRange[0]; i<=iRange[1]; i++ ) {

	// Get the number of records from the .amr file
	// Construct the .amr filename and open the file
	sprintf( szINPUTFilename, "%s/profile%i.amr", szResultsDirectory, i );
	pINPUTFile = fopen( szINPUTFilename, "r" );
		// Get the timestamp from the .amr file
		ReadDouble( pINPUTFile, &fBuffer );
		// Get the file number from the .amr file
		fscanf( pINPUTFile, "%i", &iBuffer );
		// Get the loop length from the .amr file
		ReadDouble( pINPUTFile, &fBuffer );
		// Get the number of records from the .amr file
		fscanf( pINPUTFile, "%i", &iNumRecords );
	fclose( pINPUTFile );

	// Construct the output filename and open the file
	sprintf( szOUTPUTFilename, "%s/profile%i.qts", szResultsDirectory, i );
	pOUTPUTFile = fopen( szOUTPUTFilename, "w" );
		// Construct the input filename and open the file
		sprintf( szINPUTFilename, "%s/%s%i.%s", szResultsDirectory, szRoot, i, szExtension );
		pINPUTFile = fopen( szINPUTFilename, "r" );
			// Extract the specified entries from each record and write them to a new file
			for( j=0; j<iNumRecords; j++ ) {
				// Reset the entry counter
				iEntry = 0;
				// Search the rows and columns of the current record for the particular entries
				for( m=1; m<=iNumRows; m++ ) {
					for( n=1; n<=piNumColumns[m-1]; n++ ) {

						// Switch on the data type (making it trivial to add types later on)
						switch ( ppcDataType[m-1][n-1] ) {
							
							case 'i' :
								// The entry is of integer type
								fscanf( pINPUTFile, "%i", &iBuffer );
							break;
								
							case 'f' :
								// The entry is of floating-point (double) type
								ReadDouble( pINPUTFile, &fBuffer );
							break;

							default :
							break;
						}
					
						if( piRow[iEntry] == m && piColumn[iEntry] == n ) {
							// We have located a particular entry
							// Insert a tab separator between columns in the new file
							if( iEntry > 0 )
								fprintf( pOUTPUTFile, "\t" );
							// Switch on the data type (making it trivial to add types later on)
							switch ( ppcDataType[m-1][n-1] ) {
							
								case 'i' :
									// The entry is of integer type
									fprintf( pOUTPUTFile, "%i", iBuffer );
								break;
								
								case 'f' :
									// The entry is of floating-point (double) type
									fprintf( pOUTPUTFile, "%.8e", fBuffer );
								break;

								default :
								break;
							}
							// Increment the entry counter
							iEntry++;
							// If all of the particular entries have been written to the new file then start a new row
							if( iEntry == iNumEntries )
								fprintf( pOUTPUTFile, "\n" );
						}
					}
				}
			}	
		// Close the input file	
		fclose( pINPUTFile );
	// Close the output file
	fclose( pOUTPUTFile );
}

// Free the allocated memory
free( piColumn );
free( piRow );
for( i=0; i<iNumRows; i++ )
	free( ppcDataType[i] );
free( ppcDataType );
free( piNumColumns );

return 0;
}