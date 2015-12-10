PRO CALC_EMISS_TABLES

szElementList = [ 'h', 'he', 'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne', $
                  'na', 'mg', 'al', 'si', 'p', 's', 'cl', 'ar', 'k', 'ca', $
                  'sc', 'ti', 'v', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn' ]

fTemperature = FINDGEN(41,start=4.0,increment=0.1)
fDensity = FINDGEN(17,start=7.0,increment=0.5)

; Loop over the first 30 elements
FOR iZ = 1, 30 DO BEGIN

  ; **** SPECIFY THE OUTPUT DIRECTORY FOR THE .em FILES ****
  szFileName = STRCOMPRESS( '/HYDRAD_v4/Radiation_Model/atomic_data/emissivities/chianti_v7/'+ $
    szElementList[iZ-1]+'.em', /remove_all )

  OPENW, 2, szFileName

  sZIonList = []
  iNumIons = 0
  FOR iSpec = 1, iZ DO BEGIN
    
    ; **** SUBSTITUTE THE FULL PATH TO YOUR ssw/packages/chianti/dbase/ DIRECTORY ****
    sDirectory = STRCOMPRESS( STRING( '/ssw/packages/chianti/dbase/', $
      szElementList[iZ-1], '/', szElementList[iZ-1], '_', iSpec, '/', szElementList[iZ-1], '_', $
      iSpec, '.wgfa' ), /remove_all )

    result = FILE_TEST( sDirectory )
    IF result EQ 1 THEN BEGIN
      iNumIons = iNumIons + 1
      szIonList = [szIonList, STRCOMPRESS( STRING( iSpec ), /remove_all )]
    ENDIF
    
  ENDFOR

  PRINTF, 2, 'The number of ions for which emissivity data is available.'
  PRINTF, 2, STRCOMPRESS( STRING( iNumIons ), /remove_all )
  PRINTF, 2
  PRINTF, 2, 'The particular ions in spectroscopic notation.'
  PRINTF, 2, szIonList

  ; Choose the range of spectroscopic numbers
  FOR iSpec = 1, iZ DO BEGIN

    sDirectory = STRCOMPRESS( STRING( '/ssw/packages/chianti/dbase/', $
      szElementList[iZ-1], '/', szElementList[iZ-1], '_', iSpec, '/', szElementList[iZ-1], '_', $
      iSpec, '.wgfa' ), /remove_all )

    result = FILE_TEST( sDirectory )
    IF result EQ 1 THEN BEGIN

      emiss = EMISS_CALC( iZ, iSpec, temp=fTemperature, dens=fDensity )
      emiss_total = TOTAL(emiss.em,3)
      szHeader = STRING( 'The total emissivity of ', szElementList[iZ-1], STRCOMPRESS( STRING( iSpec ), $
        /remove_all ), ' as a function of T and n.' )

      PRINTF, 2
      PRINTF, 2, szHeader
      PRINTF, 2
      PRINTF, 2, emiss_total
      
    ENDIF

  ENDFOR

  CLOSE, 2

ENDFOR

END