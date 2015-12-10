PRO PLOT_EM_LOCI

; Double precision variables
fTempVals = 0.0D
fCoord = 0.0D
fData = 0.0D
fPottaschData = 0.0D
fModelData = 0.0D
fLog10_Tmin = 0.0D
fLog10_Tmax = 0.0D
fLog10_EMmin = 0.0D
fLog10_EMmax = 0.0D

; Integer variables
iNumInstruments = 0
iNumVals = 0
iLines = 0
iItem = 0
iPixel = 0

; String variables
szFilename = ' '
szCorrectedFilename = ' '
szTitle = ' '

; Plot the emission measure loci for the 42nd pixel (apex of a 1E10 cm loop observed by EIS)
iPixel = 42

; If non-equilibrium ion population files exist then iLines = 2
iLines = 1

; Set the plot ranges
fLog10_Tmin = 5.5
fLog10_Tmax = 7.5
fLog10_EMmin = 20.0
fLog10_EMmax = 32.0

; **** GRAPHICS ****
szFilename = 'fig.eps'
szCorrectedFilename = 'fig_corrected.eps'
IF iLines EQ 1 THEN BEGIN
  !P.MULTI = 0
  SET_PLOT, 'ps'
  DEVICE, FILE=szFilename, /ENCAPSULATED, /PORTRAIT, /INCHES, XSIZE=7.0, YSIZE=4.5
ENDIF ELSE BEGIN
  !P.MULTI = [0, 1, 2, 0, 1]
  SET_PLOT, 'ps'
  DEVICE, FILE=szFilename, /ENCAPSULATED, /PORTRAIT, /INCHES, XSIZE=7.0, YSIZE=9.0
ENDELSE
; ******************

; Get the emission measure loci data
OPENR, 1, 'EM_Loci.dat'
  ; Get the number of instruments
  READF, 1, iNumInstruments
  ; Get the number of temperature values
  READF, 1, iNumVals
  ; Create an array to hold the temperature values
  fTempVals = make_array( iNumVals, /DOUBLE, VALUE = 0.0D )
  ; Read the temperature values
  READF, 1, fTempVals
  ; Create an array to hold the emission measures
  fData = make_array( iNumVals, (iNumInstruments*iLines), /DOUBLE, VALUE = 0.0D )
  FOR i=1, iPixel DO BEGIN
    ; Read the coordinate
    READF, 1, fCoord
    ; Read the emission measures
    READF, 1, fData
  ENDFOR
CLOSE, 1
fData = ALOG10(fData)

; Get the Pottasch emission measure
OPENR, 1, 'EM_Loci_Pottasch.dat'
; Get the number of instruments
READF, 1, iNumInstruments
fPottaschData = make_array( iNumInstruments, iLines+1, /DOUBLE, VALUE = 0.0D )
FOR i=1, iPixel DO BEGIN
  READF, 1, fCoord
  READF, 1, fPottaschData
ENDFOR
CLOSE, 1
fPottaschData[*,1] = ALOG10( fPottaschData[*,1] )
IF iLines GT 1 THEN BEGIN
fPottaschData[*,2] = ALOG10( fPottaschData[*,2] )
ENDIF

; Get the model emission measure
fModelData = make_array( 2, iNumVals, /DOUBLE, VALUE = 0.0D )
OPENR, 1, 'EM_Loci_model.dat'
  READF, 1, fModelData
CLOSE, 1
fModelData[1,*] = ALOG10( fModelData[1,*] )

; **** PANEL 1 ****
; Plot the equilibrium (true) emission measures
szTitle = STRING( 'EM Loci (Pixel = ', iPixel, ')' )
PLOT, fTempVals, fData[*,0], linestyle=1, xrange=[fLog10_Tmin,fLog10_Tmax], yrange=[fLog10_EMmin,fLog10_EMmax], title=szTitle, xtitle='Log!l10!n T (K)', ytitle='Log!l10!n EM (cm!u-5!n)', xth=3, yth=3, th=3, charsize=1.5, charth=3
iItem = iLines
FOR i=1, iNumInstruments-1 DO BEGIN
  OPLOT, fTempVals, fData[*,iItem], linestyle=1, th=3
  iItem = iItem + iLines
ENDFOR
; Over-plot the Pottasch emission measure
OPLOT, fPottaschData[*,0], fPottaschData[*,1], psym=1, th=3
; Over-plot the model emission measure
OPLOT, fModelData[0,*], fModelData[1,*], psym=4, th=3
; **** END OF PANEL 1 ****

; **** PANEL 2 ****
IF iLines GT 1 THEN BEGIN
; Plot the non-equilibrium emission measures
szTitle = STRING( 'EM Loci (Pixel = ', iPixel, ')' )
PLOT, fTempVals, fData[*,1], linestyle=1, xrange=[fLog10_Tmin,fLog10_Tmax], yrange=[fLog10_EMmin,fLog10_EMmax], title=szTitle, xtitle='Log!l10!n T (K)', ytitle='Log!l10!n EM (cm!u-5!n)', xth=3, yth=3, th=3, charsize=1.5, charth=3
iItem = iLines + 1
FOR i=1, iNumInstruments-1 DO BEGIN
  OPLOT, fTempVals, fData[*,iItem], linestyle=1, th=3
  iItem = iItem + iLines
ENDFOR
; Over-plot the Pottasch emission measure
OPLOT, fPottaschData[*,0], fPottaschData[*,2], psym=1, th=3
; Over-plot the model emission measure
OPLOT, fModelData[0,*], fModelData[1,*], psym=4, th=3
ENDIF
; **** END OF PANEL 2 ****

; **** GRAPHICS ****
DEVICE, /CLOSE
!P.MULTI = 0
SET_PLOT, 'win'
; ******************

; Rotate the file produced by IDL by 180 degrees if the /LANDSCAPE keyword is used
; FIXPS, szFilename, szCorrectedFilename

END