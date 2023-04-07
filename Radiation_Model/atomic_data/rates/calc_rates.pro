PRO CALC_RATES
; Construct the ionization and recombination rate files for each element

; Create the element label array
szLabel = ['h','he','li','be','b','c','n','o','f','ne','na','mg','al','si','p','s','cl','ar','k','ca','sc','ti','v','cr','mn','fe','co','ni','cu','zn']

; Create the temperature arrays
flog10Te = [4.00,4.10,4.20,4.30,4.40,4.50,4.60,4.70,4.80,4.90,5.00,5.10,5.20,5.30,5.40,5.50,5.60,5.70,5.80,5.90,6.00,6.10,6.20,6.30,6.40,6.50,6.60,6.70,6.80,6.90,7.00,7.10,7.20,7.30,7.40,7.50,7.60,7.70,7.80,7.90,8.00]
fTe = 10.0^flog10Te

; Cycle through the elements and their ions to calculate the ionization and recombination rates as a function of temperature
FOR iA=1, 30 DO BEGIN
  ; Open a data file for the current element
  szFilename = STRCOMPRESS( STRING( szLabel[iA-1], '.rts' ), /REMOVE_ALL )
  OPENW, 1, szFilename
    ; Write the ionization rates to the data file
    FOR iZ=0, iA-1 DO BEGIN
      ; Write the element and charge state to the data file
      szIon = STRCOMPRESS( STRING( szLabel[iA-1], '+', iZ ), /REMOVE_ALL )
      PRINTF, 1, szIon
      ; Calculate the ionization rates
      fionRates = ioniz_rate( '', fTe, z=iA, ion=iZ+1 )
      ;fionRates[WHERE( fionRates LT 1.0E-30 )] = 1.0E-30
      FOR i=0, 40 DO BEGIN
        IF fionRates[i] LT 1.0E-30 THEN fionRates[i] = 1.0E-30
      ENDFOR
      ; Write the ionization rates to the data file
      PRINTF, 1, fionRates
    ENDFOR
    PRINTF, 1, ''
    ; Write the recombination rates to the data file
    FOR iZ=0, iA-1 DO BEGIN
      ; Write the element and charge state to the data file
      szIon = STRCOMPRESS( STRING( szLabel[iA-1], '+', iZ ), /REMOVE_ALL )
      PRINTF, 1, szIon
      ; Calculate the recombination rates
      frecRates = recomb_rate( '', fTe, z=iA, ion=iZ+2 )
      ;frecRates[WHERE( frecRates LT 1.0E-30 )] = 1.0E-30
      FOR i=0, 40 DO BEGIN
        IF frecRates[i] LT 1.0E-30 THEN frecRates[i] = 1.0E-30
      ENDFOR
      ; Write the recombination rates to the data file
      PRINTF, 1, frecRates
    ENDFOR
  CLOSE, 1
ENDFOR

END