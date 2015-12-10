PRO ANIMATE_DFN

; Floating point variables
time_elapsed = 0.0
length = 0.0
x = 0.0
y = 0.0
z = 0.0

; Integer variables
; Set the window size
xsize = 640
ysize = 480
num_frames = 121    ; i.e. 0 to 120 inclusive
fps = 5
profile_number = 0
grid_cells = 0
distribution_data_points = 1000
i = 0
j = 0

; String variables
szFilename = ' '
szTitle = ' '

; Create the video object
szFilename = Dialog_Pickfile(/Write, Title='Save As', File='test.mp4')
IF szFilename EQ '' THEN RETURN
oVid = IDLffVideoWrite(szFilename)
vidStream = oVid.AddVideoStream(xsize, ysize, fps)

; Set up the plotting window
WINDOW, xsize = xsize, ysize = ysize, PIXMAP=0

FOR frame_no = 0, num_frames-1 DO BEGIN

    szFilename = STRCOMPRESS( STRING( 'profile', frame_no, '.amr' ), /REMOVE_ALL )

    OPENR, 1, szFilename

        READF, 1, time_elapsed
        READF, 1, profile_number
        READF, 1, length

        READF, 1, grid_cells

    CLOSE, 1

    x = FltArr( grid_cells )
    y = FltArr( distribution_data_points )
    z = FltArr( distribution_data_points, grid_cells )

    szFilename = STRCOMPRESS( STRING( 'profile', frame_no, '.dfn' ), /REMOVE_ALL )

    OPENR, 1, szFilename

        READF, 1, x
        READF, 1, y
        READF, 1, z

    CLOSE, 1

    FOR i=0, distribution_data_points-1 DO BEGIN
        FOR j=0, grid_cells-1 DO BEGIN
            if( z(i,j) lt -50.0 ) THEN z(i,j) = -50.0
        ENDFOR
    ENDFOR

    szTitle = STRING( 'Time = ', LONG(10*frame_no), ' seconds' )
    SHADE_SURF, z, y, x, AX=45, AZ=30, xrange=[-2e10,2e10], yrange=[0.0,8e9], zrange=[-50.0,-10.0], xth=2, yth=2, zth=2, xtitle='Velocity (cm s!u-1!n)', ytitle='Spatial Coordinate (cm)', ztitle='Log!l10!n f(v)', title=szTitle, charsize=2, charth=2

    ; To plot a single distribution function
    ; PLOT, y, z[*,CELL_NUMBER], linestyle=0, th=2, xth=2, yth=2, xtitle='Velocity (cm s!u-1!n)', ytitle='Log!l10!n f(v)', title=szTitle, charsize=2, charth=2  ; Where CELL_NUMBER is the cell number (counting from 0)

    ; Grab the screen image
    grab = TVRD(TRUE=1)

    ; Send the screen image to the MPEG object
    !NULL = oVid.Put(vidStream, grab)

ENDFOR

; Close the video object
oVid.Cleanup

END