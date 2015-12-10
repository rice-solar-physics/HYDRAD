PRO ANIMATE_EM_LOCI

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
iNumPixels = 0
iFrameNumber = 0
iFPS = 0

; String variables
szFilename = ' '
szTitle = ' '

; Set the number of pixels in the emission measure loci file (for a 1E10 cm loop observed by EIS)
iNumPixels = 83

; If non-equilibrium ion population files exist then iLines = 2
iLines = 1

; Set the plot ranges
fLog10_Tmin = 5.5
fLog10_Tmax = 7.5
fLog10_EMmin = 20.0
fLog10_EMmax = 32.0

; **** GRAPHICS ****
IF iLines EQ 1 THEN BEGIN
  !P.MULTI = 0
  SET_PLOT, 'win'
  ;iXsize = 640
  ;iYsize = 480
  iXsize = 1024
  iYsize = 768
ENDIF ELSE BEGIN
  !P.MULTI = [0, 1, 2, 0, 1]
  SET_PLOT, 'win'
  ;iXsize = 640
  ;iYsize = 960
  iXsize = 1024
  iYsize = 1536
ENDELSE
; ******************

; Create the video object
szFilename = Dialog_Pickfile(/Write, Title='Save As', File='test.mp4')
IF szFilename EQ '' THEN RETURN
oVid = IDLffVideoWrite(szFilename)
vidStream = oVid.AddVideoStream(iXsize, iYsize, iFPS)

; Set up the plotting window
WINDOW, xsize = iXsize, ysize = iYsize

; Get the emission measure loci data
OPENR, 1, 'EM_Loci.dat'
; Get the Pottasch emission measure data
OPENR, 2, 'EM_Loci_Pottasch.dat'
  ; Get the number of instruments
  READF, 1, iNumInstruments
  ; Get the number of instruments
  READF, 2, iNumInstruments
  ; Create an array to hold the Pottasch emission measures
  fPottaschData = make_array( iNumInstruments, iLines+1, /DOUBLE, VALUE = 0.0D )
  ; Get the number of temperature values
  READF, 1, iNumVals
  ; Create an array to hold the temperature values
  fTempVals = make_array( iNumVals, /DOUBLE, VALUE = 0.0D )
  ; Read the temperature values
  READF, 1, fTempVals
  ; Create an array to hold the emission measures
  fData = make_array( iNumVals, (iNumInstruments*iLines), /DOUBLE, VALUE = 0.0D )

  ; Get the model emission measure
  fModelData = make_array( 2, iNumVals, /DOUBLE, VALUE = 0.0D )
  OPENR, 3, 'EM_Loci_model.dat'
    READF, 3, fModelData
  CLOSE, 3
  fModelData[1,*] = ALOG10( fModelData[1,*] )
  
  FOR iFrameNumber=1, iNumPixels DO BEGIN
    ; Read the coordinate
    READF, 1, fCoord
    ; Read the coordinate
    READF, 2, fCoord
    ; Read the emission measures
    READF, 1, fData
    ; Read the emission measures
    READF, 2, fPottaschData

    fData = ALOG10(fData)
    fPottaschData[*,1] = ALOG10( fPottaschData[*,1] )
    IF iLines GT 1 THEN BEGIN
      fPottaschData[*,2] = ALOG10( fPottaschData[*,2] )
    ENDIF

    ; **** PANEL 1 ****
    ; Plot the equilibrium (true) emission measures
    szTitle = STRING( 'EM Loci (Pixel = ', iFrameNumber, ')' )
    PLOT, fTempVals, fData[*,0], linestyle=0, xrange=[fLog10_Tmin,fLog10_Tmax], yrange=[fLog10_EMmin,fLog10_EMmax], title=szTitle, xtitle='Log!l10!n T (K)', ytitle='Log!l10!n EM (cm!u-5!n)', charsize=1.5
    iItem = iLines
    FOR i=1, iNumInstruments-1 DO BEGIN
      OPLOT, fTempVals, fData[*,iItem], linestyle=0
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
      szTitle = STRING( 'EM Loci (Pixel = ', iFrameNumber, ')' )
      PLOT, fTempVals, fData[*,1], linestyle=0, xrange=[fLog10_Tmin,fLog10_Tmax], yrange=[fLog10_EMmin,fLog10_EMmax], title=szTitle, xtitle='Log!l10!n T (K)', ytitle='Log!l10!n EM (cm!u-5!n)', charsize=1.5
      iItem = iLines + 1
      FOR i=1, iNumInstruments-1 DO BEGIN
        OPLOT, fTempVals, fData[*,iItem], linestyle=0
        iItem = iItem + iLines
      ENDFOR
      ; Over-plot the Pottasch emission measure
      OPLOT, fPottaschData[*,0], fPottaschData[*,2], psym=1, th=3
      ; Over-plot the model emission measure
      OPLOT, fModelData[0,*], fModelData[1,*], psym=4, th=3
    ENDIF
    ; **** END OF PANEL 2 ****

    ; Grab the screen image
    grab = TVRD(TRUE=1)
    ; Send the screen image to the MPEG object
    !NULL = oVid.Put(vidStream, grab)
  ENDFOR
CLOSE, 2
CLOSE, 1

; Close the video object
oVid.Cleanup

; **** GRAPHICS ****
!P.MULTI = 0
SET_PLOT, 'win'
; ******************

END