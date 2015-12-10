PRO ANIMATE
; Use .phy file as an example

; Double precision variables
x = 0.0D
v = 0.0D
w = 0.0D
y = 0.0D

; floating point variables
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

; Integer variables
alength = 2000
; Set the window size
xsize = 640
ysize = 480
; Set the number of animation frames
num_frames = 121    ; i.e. 0 to 120 inclusive
frame_no = 0
fps = 5
i = 0
count = 0

; String variables
szFilename = ' '
szTitle = ' '

x = MAKE_ARRAY( alength, /DOUBLE, VALUE = 0.0D )
v = MAKE_ARRAY( alength, /DOUBLE, VALUE = 0.0D )
w = MAKE_ARRAY( alength, /DOUBLE, VALUE = 0.0D )
y = MAKE_ARRAY( alength, /DOUBLE, VALUE = 0.0D )

; Create the video object
szFilename = Dialog_Pickfile(/Write, Title='Save As', File='test.mp4')
IF szFilename EQ '' THEN RETURN
oVid = IDLffVideoWrite(szFilename)
vidStream = oVid.AddVideoStream(xsize, ysize, fps)

; Set up the plotting window
WINDOW, xsize = xsize, ysize = ysize

; Compile the animation
FOR frame_no = 0, num_frames - 1 DO BEGIN

    szFilename = STRCOMPRESS( STRING( 'profile', frame_no, '.phy' ), /REMOVE_ALL )

    ; Open the profile to plot
    OPENR, 1, szFilename

    i = 0

    WHILE ( NOT EOF(1) ) DO BEGIN

        READF, 1, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11

        x[i] = temp1            ; Curvi-linear coordinate
        v[i] = alog10(temp4)    ; Electron number density
        w[i] = alog10(temp8)    ; Electron temperature
        y[i] = alog10(temp9)    ; Hydrogen temperature

        i = i + 1

    ENDWHILE

    CLOSE, 1

    count = i - 1

    FOR i = count, alength - 1 DO BEGIN

        x[i] = x[count]
        v[i] = v[count]
        w[i] = w[count]
        y[i] = y[count]

    ENDFOR

    ; Plot the data
    szTitle = STRING( 'Time = ', LONG(10*frame_no), ' seconds' )
    PLOT, x, v,  linestyle=0, th=2, xrange=[0.0,6e9], yrange=[4.0,12.0], xth=2, yth=2, xtitle='Spatial Coordinate (cm)', ytitle='Log!l10!n T (K) and n (cm!u-3!n)', title=szTitle, charsize=2, charth=2
    OPLOT, x, w, linestyle=0, th=2
    OPLOT, x, y, linestyle=4, th=2

    ; Grab the screen image
    grab = TVRD(TRUE=1)

    ; Send the screen image to the MPEG object
    !NULL = oVid.Put(vidStream, grab)

ENDFOR

; Close the video object
oVid.Cleanup

END