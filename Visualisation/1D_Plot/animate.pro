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
xsize = 1024
ysize = 1024
; Set the number of animation frames
num_frames = 501    ; i.e. 0 to 500 inclusive
frame_no = 0
fps = 25
i = 0
count = 0

; String variables
szFilename = ' '
szTitle = ' '

x = MAKE_ARRAY( alength, /DOUBLE, VALUE = 0.0D )
t = MAKE_ARRAY( alength, /DOUBLE, VALUE = 0.0D )
u = MAKE_ARRAY( alength, /DOUBLE, VALUE = 0.0D )
v = MAKE_ARRAY( alength, /DOUBLE, VALUE = 0.0D )
w = MAKE_ARRAY( alength, /DOUBLE, VALUE = 0.0D )
y = MAKE_ARRAY( alength, /DOUBLE, VALUE = 0.0D )
z = MAKE_ARRAY( alength, /DOUBLE, VALUE = 0.0D )
mi = MAKE_ARRAY( 2, /DOUBLE, VALUE = 0.0D )
ma = MAKE_ARRAY( 2, /DOUBLE, VALUE = 0.0D )

; Create the video object
szFilename = Dialog_Pickfile(/Write, Title='Save As', File='all.mp4')
IF szFilename EQ '' THEN RETURN
oVid = IDLffVideoWrite(szFilename)
vidStream = oVid.AddVideoStream(xsize, ysize, fps)

; Set up the plotting window
WINDOW, xsize = xsize, ysize = ysize
!P.MULTI = [0,2,2]

; Compile the animation
FOR frame_no = 0, num_frames - 1 DO BEGIN

    szFilename = STRCOMPRESS( STRING( 'profile', frame_no, '.phy' ), /REMOVE_ALL )

    ; Open the profile to plot
    OPENR, 1, szFilename

    i = 0

    WHILE ( NOT EOF(1) ) DO BEGIN

        READF, 1, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11

        x[i] = temp1            ; Curvi-linear coordinate
        t[i] = temp2            ; Bulk velocity
        u[i] = ALOG10(temp4)    ; Electron number density
        v[i] = ALOG10(temp6)    ; Electron pressure
        w[i] = ALOG10(temp7)    ; Hydrogen pressure
        y[i] = temp8            ; Electron temperature
        z[i] = temp9            ; Hydrogen temperature

        i = i + 1

    ENDWHILE

    CLOSE, 1

    count = i - 1

    FOR i = count, alength - 1 DO BEGIN

        x[i] = x[count]
        t[i] = t[count]
        u[i] = u[count]
        v[i] = v[count]
        w[i] = w[count]
        y[i] = y[count]
        z[i] = z[count]

    ENDFOR

    ; Scale the data
    x = x / 1E8   ; cm to Mm
    t = t / 1E5   ; cm/s to km/s
    ; u = u / 1.0 ; cm^-3 to cm^-3
    ; v = v / 1.0 ; dyne cm^-2 to dyne cm^-2
    ; w = w / 1.0 ; dyne cm^-2 to dyne cm^-2
    y = y / 1E6   ; K to MK
    z = z / 1E6   ; K to MK

    L = MAX(x)

    ; Plot the data
    szTitle = STRCOMPRESS( STRING( 'Time = ', (10.0/3600.0)*frame_no, ' hrs' ) )

    ; **** PANEL 1 ****
    PLOT, x, t,  linestyle=0, th=2, xrange=[MIN(x),MAX(x)], yrange=[MIN(t),MAX(t)], xth=2, yth=2, xtitle='Spatial Coordinate (Mm)', ytitle='Velocity (km s!u-1!n)', title=szTitle, charsize=2, charth=2
    ; **** PANEL 1 ****
    
    ; **** PANEL 2 ****
    PLOT, x, u,  linestyle=0, th=2, xrange=[MIN(x),MAX(x)], yrange=[MIN(u),MAX(u)], xth=2, yth=2, xtitle='Spatial Coordinate (Mm)', ytitle='Density (cm!u-3!n)', title=szTitle, charsize=2, charth=2
    ; **** PANEL 2 ****

    ; **** PANEL 3 ****
    mi[0] = MIN(v)
    mi[1] = MIN(w)
    ma[0] = MAX(v)
    ma[1] = MAX(w)
    ; PLOT, x, v,  linestyle=0, th=2, xrange=[MIN(x),MAX(x)], yrange=[MIN(mi),MAX(ma)], xth=2, yth=2, xtitle='Spatial Coordinate (Mm)', ytitle='Pressure (dyne cm!u-2!n)', title=szTitle, charsize=2, charth=2
    PLOT, x, v,  linestyle=0, th=2, xrange=[MIN(x),MAX(x)], yrange=[MIN(mi),0.0], xth=2, yth=2, xtitle='Spatial Coordinate (Mm)', ytitle='Pressure (dyne cm!u-2!n)', title=szTitle, charsize=2, charth=2
    OPLOT, x, w, linestyle=4, th=2
    ; **** PANEL 3 ****

    ; **** PANEL 4 ****
    ma[0] = MAX(y)
    ma[1] = MAX(z)
    PLOT, x, y,  linestyle=0, th=2, xrange=[MIN(x),MAX(x)], yrange=[0.0,MAX(ma)], xth=2, yth=2, xtitle='Spatial Coordinate (Mm)', ytitle='Temperature (MK)', title=szTitle, charsize=2, charth=2
    OPLOT, x, z, linestyle=4, th=2
    ; **** PANEL 4 ****

    ; Grab the screen image
    grab = TVRD(TRUE=1)

    ; Send the screen image to the MPEG object
    !NULL = oVid.Put(vidStream, grab)

ENDFOR

; Close the video object
oVid.Cleanup

END