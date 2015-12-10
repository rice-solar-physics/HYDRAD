PRO PLOT_SPECTRUM

; Double precision variables
fLambda0 = 0.0D
fLambda = 0.0D
fVelocity = 0.0D
fEquil = 0.0D
fNEquil = 0.0D

; floating point variables

; Integer variables
iPoints = 0
iCell = 0
i = 0

; String variables
szFilename = ' '

; Open the file containing the spectrum
szFilename = Dialog_Pickfile(/Read, Title='Choose .spec.det file', File='')
IF szFilename EQ '' THEN RETURN

PRINT, 'Rest wavelength: '
READ, fLambda0

PRINT, 'Show spectrum for detector cell (1 = foot-point): '
READ, iCell

OPENR, 1, szFilename

  READF, 1, iPoints

  ; Create the array to store the wavelengths
  fLambda = MAKE_ARRAY( iPoints, /DOUBLE, VALUE = 0.0D )
  ; Create the array to store the velocities
  fVelocity = MAKE_ARRAY( iPoints, /DOUBLE, VALUE = 0.0D )
  ; Create the array to store the equilibrium spectrum
  fEquil = MAKE_ARRAY( iPoints, /DOUBLE, VALUE = 0.0D )
  ; Create the array to store the non-equilibrium spectrum
  fNEquil = MAKE_ARRAY( iPoints, /DOUBLE, VALUE = 0.0D )
  
  ; Read the wavelengths from the data file
  READF, 1, fLambda
  
  ; Convert the wavelengths into velocities (km/s) - blue-shifts = negative velocities
  fVelocity = ( ( fLambda - fLambda0 ) / fLambda ) * (3e5)

  FOR i = 1, iCell, 1 DO BEGIN
    ; Read the location of the detector cell on the solar disk (in arcsec)
    READF, 1, fArcSec
    ; Read the equilibrium spectrum
    READF, 1, fEquil
    ; Read the non-equilibrium spectrum
    READF, 1, fNEquil
  ENDFOR

CLOSE, 1

; **** GRAPHICS ****
!P.MULTI = [0, 1, 2]
SET_PLOT, 'ps'
; ******************

DEVICE, FILE=szFilename+'.ps', /ENCAPSULATED, /INCHES, XSIZE=7.0, YSIZE=9.0

; Find the range for the y-axis scale
fYmax = fltarr(2)
fYmax[0] = MAX(fEquil)
fYmax[1] = MAX(fNEquil)

; Plot as a function of wavelength
PLOT, fLambda, fEquil, yrange=[0.0,MAX(fYmax)], linestyle=2, title='Intensity vs. Wavelength', xtitle='Wavelength (A)', ytitle='Log!l10!n Counts (DN pixel!u-1!n s!u-1!n)', xth=3, yth=3, th=3, charsize=1.5, charth=3
OPLOT, fLambda, fNEquil, linestyle=0, th=3

; Plot as a function of velocity
PLOT, fVelocity, fEquil, yrange=[0.0,MAX(fYmax)], linestyle=2, title='Intensity vs. Doppler-shift', xtitle='Velocity (km/s)', ytitle='Log!l10!n Counts (DN pixel!u-1!n s!u-1!n)', xth=3, yth=3, th=3, charsize=1.5, charth=3
OPLOT, fVelocity, fNEquil, linestyle=0, th=3

; **** GRAPHICS ****
DEVICE, /CLOSE
!P.MULTI = 0
SET_PLOT, 'win'
; ******************

END