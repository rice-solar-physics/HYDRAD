PRO PLOT_INTENSITY

; Double precision variables
fArcsec = 0.0D
fEquil = 0.0D
fNEquil = 0.0D

; floating point variables
temp1 = 0.0
temp2 = 0.0
temp3 = 0.0

; Integer variables
iSize = 1000
i = 0
iCount = 0

; String variables
szFilename = ' '

; Create the array to store the coordinates of the detector cell
fArcsec = MAKE_ARRAY( iSize, /DOUBLE, VALUE = 0.0D )
; Create the array to store the equilibrium intensities
fEquil = MAKE_ARRAY( iSize, /DOUBLE, VALUE = 0.0D )
; Create the array to store the non-equilibrium intensities
fNEquil = MAKE_ARRAY( iSize, /DOUBLE, VALUE = 0.0D )

; Open the file containing the intensities
szFilename = Dialog_Pickfile(/Read, Title='Choose .det file', File='')
IF szFilename EQ '' THEN RETURN

OPENR, 1, szFilename

  i = 0

  WHILE ( NOT EOF(1) ) DO BEGIN

    READF, 1, temp1, temp2, temp3
      
      fArcsec[i] = temp1  ; Coordinate (arcsec)
      fEquil[i] = temp2   ; Equilibrium intensities (DN pixel^-1 s^-1)
      fNEquil[i] = temp3  ; Non-equilibrium intensities (DN pixel^-1 s^-1)

      i = i + 1

  ENDWHILE

CLOSE, 1

iCount = i
FOR i = iCount, iSize - 1 DO BEGIN

  fArcsec[i] = fArcsec[iCount]
  fEquil[i] = fEquil[iCount]
  fNEquil[i] = fNEquil[iCount]

ENDFOR

; **** GRAPHICS ****
!P.MULTI = [0, 1, 2]
SET_PLOT, 'ps'
; ******************

DEVICE, FILE=szFilename+'.ps', /ENCAPSULATED, /INCHES, XSIZE=7.0, YSIZE=9.0

; Find the range for the y-axis scale
fYmax = fltarr(2)
fYmax[0] = MAX(fEquil)
fYmax[1] = MAX(fNEquil)

; Plot as a function of coordinate
PLOT, fArcsec, fEquil, yrange=[0.0,MAX(fYmax)], linestyle=2, title='Intensity vs. Coordinate', xtitle='Coordinate (arcsec)', ytitle='Log!l10!n Counts (DN pixel!u-1!n s!u-1!n)', xth=3, yth=3, th=3, charsize=1.5, charth=3
OPLOT, fArcsec, fNEquil, linestyle=0, th=3

; Plot as a function of detector cell number
PLOT, Indgen(iCount)+1, fEquil, yrange=[0.0,MAX(fYmax)], linestyle=2, psym=-2, title='Intensity vs. Detector Cell Number', xtitle='Detector Cell Number', ytitle='Log!l10!n Counts (DN pixel!u-1!n s!u-1!n)', xth=3, yth=3, th=3, charsize=1.5, charth=3
OPLOT, Indgen(iCount)+1, fNEquil, linestyle=0, psym=-4, th=3

; **** GRAPHICS ****
DEVICE, /CLOSE
!P.MULTI = 0
SET_PLOT, 'win'
; ******************

END