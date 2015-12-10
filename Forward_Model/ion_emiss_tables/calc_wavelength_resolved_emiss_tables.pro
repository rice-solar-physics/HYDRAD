PRO CALC_WAVELENGTH_RESOLVED_EMISS_TABLES

szElementList = [ 'h', 'he', 'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne', $
                  'na', 'mg', 'al', 'si', 'p', 's', 'cl', 'ar', 'k', 'ca', $
                  'sc', 'ti', 'v', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn' ]
  
fTemperature = FINDGEN(41,start=4.0,increment=0.1)
fDensity = FINDGEN(17,start=7.0,increment=0.5)

; Loop over the first 30 elements
FOR iZ = 1, 30 DO BEGIN     

  ; Loop over each ionization state, except the fully ionized state
  FOR iSpec=1,iZ DO BEGIN    

    ; **** SUBSTITUTE THE FULL PATH TO YOUR ssw/packages/chianti/dbase/ DIRECTORY ****
    sDirectory = STRCOMPRESS( STRING( '/ssw/packages/chianti/dbase/', $
      szElementList[iZ-1], '/', szElementList[iZ-1], '_', iSpec, '/', szElementList[iZ-1], '_', $
      iSpec,'.wgfa' ), /remove_all ) 

    result = FILE_TEST( sDirectory )  
    IF result EQ 1 THEN BEGIN    

      ; Calculate the emissivities for each line of the ionization state    
      emiss = EMISS_CALC( iZ, iSpec, temp=fTemperature, dens=fDensity )
    
      ; Open a file to write the data
      OPENW, 1, STRCOMPRESS( STRING( szElementList[iZ-1], '_', iSpec, '.em'), /remove_all )
    
      ; Write the ion label to the data file
      PRINTF, 1, STRCOMPRESS( STRING( szElementList[iZ-1], '_', iSpec ), /remove_all )
    
      ; Write the number of wavelengths to the data file
      iNumWvlns = SIZE( emiss.lambda, /n_elements )
      PRINTF, 1, iNumWvlns
    
      ; Write the wavelength range to the data file
      PRINTF, 1, emiss[0].lambda, emiss[iNumWvlns-1].lambda
    
      ; Write the list of wavelengths to the data file
      PRINTF, 1, emiss.lambda
    
      ; Write the number of densities to the data file
      PRINTF, 1, SIZE( fDensity, /n_elements )
    
      ; Write the list of densities to the data file
      PRINTF, 1, fDensity
    
      ; Write the number of temperatures to the data file
      PRINTF, 1, SIZE( fTemperature, /n_elements )
    
      ; Write the list of temperature to the data file
      PRINTF, 1, fTemperature
    
      ; For each wavelength, write the temperature and density-dependent 
      ; ion emissivities to the data file
      FOR w=1L, iNumWvlns DO BEGIN
        ; Write the wavelength to the data file
        PRINTF, 1, emiss[w-1].lambda
        ; Write the ion emissivities to the data file
        PRINTF, 1, emiss[w-1].em    
      ENDFOR
        
      CLOSE, 1    
    
    ENDIF

    ; Repeat the calculation for dielectronic recombination data
    ; **** SUBSTITUTE THE FULL PATH TO YOUR ssw/packages/chianti/dbase/ DIRECTORY ****
    sDirectory = STRCOMPRESS( STRING( '/ssw/packages/chianti/dbase/', $
      szElementList[iZ-1], '/', szElementList[iZ-1], '_', iSpec, 'd/', szElementList[iZ-1], '_', $
      iSpec, 'd.wgfa' ), /remove_all )

    result = FILE_TEST( sDirectory )
    IF result EQ 1 THEN BEGIN
      
      emiss_d = EMISS_CALC( iZ, iSpec, temp=fTemperature, dens=fDensity, /diel )
    
      OPENW, 1, STRCOMPRESS( STRING( szElementList[iZ-1], '_', iSpec, 'd.em' ), /remove_all )

      ; Write the ion label to the data file
      PRINTF, 1, STRCOMPRESS( STRING( szElementList[iZ-1], '_', iSpec ), /remove_all )
    
      ; Write the number of wavelengths to the data file
      iNumWvlns = SIZE( emiss_d.lambda, /n_elements )
      PRINTF, 1, iNumWvlns
    
      ; Write the wavelength range to the data file
      PRINTF, 1, emiss_d[0].lambda, emiss_d[iNumWvlns-1].lambda
    
      ; Write the list of wavelengths to the data file
      PRINTF, 1, emiss_d.lambda
    
      ; Write the number of densities to the data file
      PRINTF, 1, SIZE(fDensity,/n_elements)
    
      ; Write the list of densities to the data file
      PRINTF, 1, fDensity
    
      ; Write the number of temperatures to the data file
      PRINTF, 1, SIZE(fTemperature,/n_elements)
    
      ;Write the list of temperatures to the data file
      PRINTF, 1, fTemperature
    
      ; For each wavelength, write the temperature and density-dependent 
      ; ion emissivities to the data file
      FOR w=1L, iNumWvlns DO BEGIN
        ; Write the wavelength to the data file
        PRINTF, 1, emiss_d[w-1].lambda
        ; Write the ion emissivities to the data file
        PRINTF, 1, emiss_d[w-1].em    
      ENDFOR
        
      CLOSE, 1
    
    ENDIF
    
  ENDFOR

ENDFOR

END