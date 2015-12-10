PRO ANIMATE_IB

; Double precision variables
Teq = 0.0D  ; Electron temperature
I1eq = 0.0D ; Equilibrium ionisation balance
I2eq = 0.0D
I3eq = 0.0D
T = 0.0D    ; Electron temperature
I1 = 0.0D   ; Ionisation balance
I2 = 0.0D
I3 = 0.0D

; floating point variables
Te = 0.0
temp1 = 0.0
temp2 = 0.0
temp3 = 0.0
temp4 = 0.0
temp5 = 0.0
temp6 = 0.0
temp7 = 0.0
temp8 = 0.0
temp9 = 0.0
temp10 = 0.0
temp11 = 0.0
temp12 = 0.0
temp13 = 0.0
temp14 = 0.0
temp15 = 0.0
temp16 = 0.0
temp17 = 0.0
temp18 = 0.0
temp19 = 0.0
temp20 = 0.0
temp21 = 0.0
temp22 = 0.0
temp23 = 0.0
temp24 = 0.0
temp25 = 0.0
temp26 = 0.0
temp27 = 0.0
temp28 = 0.0
temp29 = 0.0
x = 0.0
y = 0.0

; Integer variables
alength = 2000
i = 0
A = 0
; Set the window size
xsize = 640
ysize = 480
; Set the number of animation frames
num_frames = 121    ; i.e. 0 to 120 inclusive
frame_no = 0
fps = 0
count = 0
max_index = 0

; String variables
szFilename = ' '
szTitle = ' '

Teq = MAKE_ARRAY( 41, /DOUBLE, VALUE = 0.0D )
I1eq = MAKE_ARRAY( 41, /DOUBLE, VALUE = 0.0D )
I2eq = MAKE_ARRAY( 41, /DOUBLE, VALUE = 0.0D )
I3eq = MAKE_ARRAY( 41, /DOUBLE, VALUE = 0.0D )
T = MAKE_ARRAY( alength, /DOUBLE, VALUE = 0.0D )
I1 = MAKE_ARRAY( alength, /DOUBLE, VALUE = 0.0D )
I2 = MAKE_ARRAY( alength, /DOUBLE, VALUE = 0.0D )
I3 = MAKE_ARRAY( alength, /DOUBLE, VALUE = 0.0D )

; Open the equilibrium ionisation balance

OPENR, 1, 'he.bal'

FOR i = 0, 40 DO BEGIN

    ;READF, 1, Te, temp1, temp2 ;Hydrogen

    READF, 1, Te, temp1, temp2, temp3 ; Helium
    Teq[i] = Te
    I1eq[i] = temp1
    I2eq[i] = temp2
    I3eq[i] = temp3

    ;READF, 1, Te, temp1, temp2, temp3, temp4 ; Lithium
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5 ; Beryllium
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6 ;Boron
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7 ; Carbon
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8 ; Nitrogen
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9 ; Oxygen
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10 ; Fluorine
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11 ; Neon
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12 ; Sodium
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13 ; Magnesium
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14 ; Aluminium
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15 ; Silicon
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16 ; Phosphorus
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17 ; Sulphur
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18 ; Chlorine
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19 ; Argon
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20 ; Potassium
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21 ; Calcium
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22 ; Scandium
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23 ; Titanium
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24 ; Vanadium
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24, temp25 ; Chromium
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24, temp25, temp26 ; Manganese
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24, temp25, temp26, temp27 ; Iron
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24, temp25, temp26, temp27, temp28 ; Cobalt
    ;READF, 1, Te, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24, temp25, temp26, temp27, temp28, temp29 ; Nickel

ENDFOR

CLOSE, 1

; Create the video object
szFilename = Dialog_Pickfile(/Write, Title='Save As', File='test.mp4')
IF szFilename EQ '' THEN RETURN
oVid = IDLffVideoWrite(szFilename)
vidStream = oVid.AddVideoStream(xsize, ysize, fps)

; Set up the plotting window
WINDOW, xsize = xsize, ysize = ysize

; Compile the animation
FOR frame_no = 0, num_frames - 1 DO BEGIN

    ; Open the .phy file to get the electron temperature

    szFilename = STRCOMPRESS( STRING( 'profile', frame_no, '.phy' ), /REMOVE_ALL )

    OPENR, 1, szFilename

    i = 0

    WHILE ( NOT EOF(1) ) DO BEGIN

        READF, 1, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11

        T[i] = alog10(temp8)    ; Electron temperature

        i = i + 1

    ENDWHILE

    CLOSE, 1

    ; Open the ionisation balance file to get the ion populations

    szFilename = STRCOMPRESS( STRING( 'profile', frame_no, '.ine' ), /REMOVE_ALL )

    OPENR, 1, szFilename

    i = 0

    WHILE ( NOT EOF(1) ) DO BEGIN

        READF, 1, temp1

        ; Get the ion populations (choose those that are in the .ine files)

        READF, 1, A, temp1, temp2 ;Hydrogen

        READF, 1, A, temp1, temp2, temp3 ; Helium
        I1[i] = temp1
        I2[i] = temp2
        I3[i] = temp3

        ;READF, 1, A, temp1, temp2, temp3, temp4 ; Lithium
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5 ; Beryllium
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6 ;Boron
        READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7 ; Carbon
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8 ; Nitrogen
        READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9 ; Oxygen
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10 ; Fluorine
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11 ; Neon
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12 ; Sodium
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13 ; Magnesium
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14 ; Aluminium
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15 ; Silicon
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16 ; Phosphorus
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17 ; Sulphur
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18 ; Chlorine
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19 ; Argon
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20 ; Potassium
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21 ; Calcium
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22 ; Scandium
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23 ; Titanium
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24 ; Vanadium
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24, temp25 ; Chromium
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24, temp25, temp26 ; Manganese
        READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24, temp25, temp26, temp27 ; Iron
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24, temp25, temp26, temp27, temp28 ; Cobalt
        ;READF, 1, A, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24, temp25, temp26, temp27, temp28, temp29 ; Nickel

        i = i + 1

    ENDWHILE

    CLOSE, 1

    count = ( i - 1 ) / 2

    FOR i = count, alength - 1 DO BEGIN

        T[i] = T[count]
        I1[i] = I1[count]
        I2[i] = I2[count]
        I3[i] = I3[count]

    ENDFOR

    ; Plot the data
    szTitle = STRING( 'Time = ', LONG(10*frame_no), ' seconds' )
    PLOT, Teq, I1eq, linestyle=1, th=2, xrange=[4.0,7.0], yrange=[0.0,1.0], xth=2, yth=2, xtitle='Log!l10!n T (K)', ytitle='He Population Fraction', title=szTitle, charsize=2, charth=2
    OPLOT, Teq, I2eq, linestyle=1, th=2
    OPLOT, Teq, I3eq, linestyle=1, th=2
    OPLOT, T, I1, linestyle=0, th=2
    OPLOT, T, I2, linestyle=0, th=2
    OPLOT, T, I3, linestyle=0, th=2

    y = MAX( I1, max_index )
    x = T[max_index]
    XYOUTS, x, y, 'I', charsize=2, charth=2

    y = MAX( I2, max_index )
    x = T[max_index]
    XYOUTS, x, y, 'II', charsize=2, charth=2

    y = MAX( I3, max_index )
    x = T[max_index]
    XYOUTS, x, y, 'III', charsize=2, charth=2

    ; Grab the screen image
    grab = TVRD(TRUE=1)

    ; Send the screen image to the MPEG object
    !NULL = oVid.Put(vidStream, grab)

ENDFOR

; Close the video object
oVid.Cleanup

END