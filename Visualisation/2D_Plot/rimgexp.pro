Function	rImgExp, image, $
			xs, ys, $
			out_xs, out_ys, $
			x_ran, y_ran, $
			ASPECT=aspect, $
			INTERPOLATE=interp, $
			MASKVALUE=maskvalue, $
			PS_INTERP_SIZE=psis, $
			POSITION=p, $
			NO_EXPAND=no_expand, $
			HELP=help, $
                        _extra=extrapars, test=test, noerase=noerase, overplot=overplot

SccsId = '@(#)imgexp.pro 2.2 7/16/93 Fen Tamanaha'
;+
; NAME:
;	IMGEXP
;
; PURPOSE:
;	This function expands the array <Image> to fill the current plotting
;	window.  This routine works for both X and PostScript devices.  The
;	optional scales <XS> and <YS> are likewise transformed and returned
;	in option parameters <Out_XS> and <Out_YS>.
;
; CATEGORY:
;	Image expansion.
;
; CALLING SEQUENCE:
;	Result = IMGEXP(Image, XS, YS, Out_XS, Out_YS, X_Ran, Y_Ran)
;
; INPUTS:
;	Image:	Two-dimensional array to be expanded.
;
; OPTIONAL INPUTS:
;	XS:	Vector of x-axis values.  The length must equal the number of
;		rows in <Image>
;
;	YS:	Vector of y-axis values.  The length must equal the number of
;		columns in <Image>
;
; KEYWORD PARAMETERS:
;	ASPECT=	Set this keyword to the aspect ratio (width/height) of the
;		pixels.  /ASPECT is the same as ASPECT=1 and produces square
;		pixels.
;
;	/INTERPOLATE:
;		Set this switch to enable bilinear interpolation for pixels
;		in the expanded image.  See /PS_INTERP_SIZE for information
;		on using this switch on a PostScript device.
;
;	MASKVALUE=
;		Set this keyword to the value that uninterpolated pixels around
;		the border of the image should be given.  The default is 
;		-9999.0.  Interpolation is not performed beyond the centers of
;		the original pixels.
;
;	PS_INTERP_SIZE=
;		Since PostScript devices have scalable pixels it is necessary
;		to force expansion to at most this many pixels in either
;		dimension.  The default is 256.  (It's really more complicated
;		than this.  Read the code if you need to know.)
;
;	POSITION=
;		Set this keyword to the variable that is to hold the four-
;		element vector containing the device coordinates of the
;		plotting region that will contain the expanded image.  This
;		is designed to be used by subsequent TV and PLOT commands.
;
;	/NO_EXPAND:
;		Set this switch to prevent the image from being expanded
;		to fill the plotting window.  An aspect ration of 1:1 is
;		forced for PostScript output so that it conforms to the X
;		window view.
;
; OUTPUTS:
;	Result:	This function returns an expanded version of the input <Image>
;		possibly interpolated.
;
; OPTIONAL OUTPUTS:
;	Out_XS:	Vector of x-axis values corresponding the the expanded image.
;
;	Out_YS:	Vector of y-axis values corresponding the the expanded image.
;
;	X_Ran:	Two-element vector that contains the full x-axis range
;		including the width of the pixels.  It is designed to be used
;		as input to the PLOT command.
;	
;	Y_Ran:	Two-element vector that contains the full y-axis range
;		including the height of the pixels.  It is designed to be used
;		as input to the PLOT command.
;	
; RESTRICTIONS:
;	This routine may work for other devices, but it has only been tested
;	on 'X' and 'PS'.
;
; PROCEDURE:
;	Straight forward.  :-)
;
; EXAMPLE:
;	p = 0
;	big = IMGEXP(small, lon, lat, biglon, biglat, xr, yr, Position=p)
;	TVSCL, big, p(0), p(1), /Device
;	Plot, [0,1], /NoData, /NoErase, Position=p, /Device, $
;		XRange=xr, YRange=yr
;
;	junk = IMGSCL( )	;prints out a "Usage:" line
;
; MODIFICATION HISTORY:
; 	Written by:	Fen Tamanaha, July 9, 1993.  Release 2.1
;	July 16, 1993	Fen: (2.2) Added /No_Expand keyword
;       Jan 10, 2000    D. Finkbeiner - added _extra pass-through to Plot
;-

    On_Error, 2

;
; Go through the optional parameters.
;
    n_parms = N_Params()
    If ( Keyword_Set(help) ) Then n_parms = 0
    Case ( n_parms ) Of 
	0: Begin
	    Message, 'im = IMGEXP(image [,xs [,ys [,out_xs [,out_xs [,x_ran [,y_ran]]]]]]', /Info
	    Message, '            [,ASPECT=] [,/INTERPOLATE] [,MASKVALUE=]', /Info
	    Message, '            [,PS_INTERP_SIZE=] [,POSITION=] [,/NO_EXPAND]', /Info
	    Return, 0
	End
        1: Begin
            sz = Size(image)
            If ( sz(0) NE 2 ) Then Begin
                Message, '<image> must be an array.'
            EndIf
            xs = FIndGen(sz(1))
            ys = FIndGen(sz(2))
        End
        2: Begin
            sz = Size(image)
            If ( sz(0) NE 2 ) Then Begin
                Message, '<image> must be an array.'
            EndIf
            If ( N_Elements(xs) NE sz(1) ) Then Begin
                Message, '<xs> does not match <image> dimensions.'
            EndIf
            ys = FIndGen(sz(2))
        End
        3: Begin
            sz = Size(image)
            If ( sz(0) NE 2 ) Then Begin
                Message, '<image> must be an array.'
            EndIf
            If ( N_Elements(xs) NE sz(1) ) Then Begin
                Message, '<xs> does not match <image> dimensions.'
            EndIf
            If ( N_Elements(ys) NE sz(2) ) Then Begin
                Message, '<ys> does not match <image> dimensions.'
            EndIf
        End
	Else: Begin
            sz = Size(image)
            If ( sz(0) NE 2 ) Then Begin
                Message, '<image> must be an array.'
            EndIf
	End
    EndCase

;
; Establish image variables and determine the aspect ration.
;
    im_x_width = Float(sz(1))			;image width
    im_y_width = Float(sz(2))			;image height
    im_aspect = im_x_width / im_y_width		;image aspect (width/height)

;
; If MASKVALUE contains a value, then that value is assumed to be a flag
;	for data not to be used in scaling.  Is also the value that will
;	be used to "blank out" the border region which could not be 
;	interpolated.  If MASKVALUE does not contain a value then it will
;	be assigned to -9999.0 and used to "blank out" the border
;	region.  A warning is issued if interpolation and thus border
;	blanking will occur.
;
    If ( N_Elements(maskvalue) LE 0 ) Then Begin
	maskvalue = -9999.0
	If ( Keyword_Set(interp) ) Then Begin
	    msg = String(Format='("Warning: Uninterpolated border set to ", F7.0, ".")', maskvalue)
	    Message, msg, /Info
	EndIf
    EndIf

;
; No matter what keywords are set, the same axis ranges will be used.  To
;	account for pixel width, the ranges extend half a pixel width
;	beyond the specified centers.  Because this routine does not
;	interpolate beyond the pixel centers, interpolation will cause
;	this border region to be set to the background color.
;
    xs_delta = (xs(im_x_width-1) - xs(0)) / Float(im_x_width - 1.0)
    ys_delta = (ys(im_y_width-1) - ys(0)) / Float(im_y_width - 1.0)
    x_ran = [xs(0)-xs_delta/2.0,xs(im_x_width-1)+xs_delta/2.0]
    y_ran = [ys(0)-ys_delta/2.0,ys(im_y_width-1)+ys_delta/2.0]

;
; Use a dummy plot to determine the plot region, establish device variables,
;	and determine the aspect ratio.
;

; Modified 4 May, 2000 by Doug Finkbeiner
;  by removing /NoErase from the following plot command, display
;  behaves properly for !p.multi not [0,1,1]. 

;    Plot, [0,1], /NoData, XStyle=4, YStyle=4, /NoErase, _extra=extrapars
    if size(extrapars, /type) EQ 8 then begin
    if total(strmatch(tag_names(extrapars), 'xs*',/fold)) NE 0 then begin
		ind=where(strmatch(tag_names(extrapars), 'xst*', /fold))
		if (extrapars.(ind) AND 4) EQ 0 then extrapars.(ind)+=4
    endif
    if total(strmatch(tag_names(extrapars), 'xs*',/fold)) NE 0 then begin
                ind=where(strmatch(tag_names(extrapars), 'yst*', /fold))
                if (extrapars.(ind) AND 4) EQ 0 then extrapars.(ind)+=4
    endif
    endif
noerase1=keyword_set(noerase)
;if total(!p.multi NE 0) then begin;noerase1=0
    Plot, [0,1], /NoData, XStyle=4, YStyle=4, _extra=extrapars, noerase=noerase1

    dev_x_range = !X.Window * !D.X_VSize	;window range in device
    dev_y_range = !Y.Window * !D.Y_VSize	; coordinates
    dev_x_width = dev_x_range(1) - dev_x_range(0) + 1
    dev_y_width = dev_y_range(1) - dev_y_range(0) + 1
    dev_aspect = dev_x_width / dev_y_width	;device aspect (width/height)


if total(!p.multi NE 0) then begin;noerase1=0
   nimg=!p.multi[1]*!p.multi[2]
   if !P.multi[0] EQ 0 then nremain=nimg else nremain=!p.multi[0]
   ;if ~keyword_set(noerase1) then !p.multi[0]=nremain-1
 ;  if keyword_set(noerase) then STOP
   if ~keyword_set(overplot) then !p.multi[0]=nremain-1
endif

;
; If ASPECT has been set and is greater than zero, then it contains the 
;	aspect ratio of each pixel.  The pixel shape is maintained by 
;	altering the device coordinate widths of the plotting region.
; The aspect ratio is forced to 1:1 so that the behavior under PostScript
;	mimic the X display.
;
    If ( N_Elements(aspect) GT 0 ) Then Begin
	pix_aspect = aspect(0)
    EndIf Else Begin
	pix_aspect = 0.0
    EndElse
    If ( pix_aspect LT 0 ) Then Begin
	Message, 'Warning: ASPECT cannot be negative --- ignoring.', /Cont
	pix_aspect = 0.0
    EndIf

    If ( Keyword_Set(no_expand) ) Then Begin
	If ( pix_aspect NE 0.0 ) Then Begin
	    Message, 'Warning: ASPECT keyword ignored by /NO_EXPAND.', /Cont
	EndIf
	pix_aspect = 1.0			;force square pixels

	If ( Keyword_Set(interp) ) Then Begin
	    Message, 'Warning: INTERPOLATE keyword ignored by /NO_EXPAND.', $
									/Cont
	EndIf
    EndIf

    If ( pix_aspect GT 0.0 ) Then Begin
	aspect_ratio = im_aspect * pix_aspect / dev_aspect
	If ( aspect_ratio GT 1.0 ) Then Begin
	    dev_y_width = dev_y_width / aspect_ratio
        EndIf Else Begin
	    dev_x_width = dev_x_width * aspect_ratio
        EndElse
    EndIf
;
; Set the plotting window position.
;
    p = [dev_x_range(0),dev_y_range(0), $
		dev_x_range(0)+dev_x_width,dev_y_range(0)+dev_y_width]


;
; If the plotting device has scalable pixels, then manual expansion is not
;	necessary unless interpolation is desired.  If the the pixels are
;	not hardware scalable, then the expansion must be performed here.
; PostScript has scalable pixels while the windows do not.  Interpolation is
;	not allowed when /No_Expand is set.
;
    scalable = (!D.Flags And 1) NE 0
    If ( scalable ) Then Begin
	If ( Keyword_Set(interp) And Not Keyword_Set(no_expand) ) Then Begin
	;
	; Interpolation for scalable pixels is a bit tricky.  The approach
	;	taken here is to shrink the plotting region down to some
	;	reasonable size.  The image is then expanded over this
	;	intermediate region and when plotted, will expand to the
	;	appropriate size.  The keyword PS_INTERP_SIZE= can be used
	;	to override the default maximum size (256).
	;
	    If ( N_Elements(psis) GT 0 ) Then Begin
		width_limit = 1.0 * psis(0)	;max size of interp region
	    EndIf Else Begin
		width_limit = 256.0		;def max size of interp region
	    EndElse
            ws = dev_x_width / width_limit	;width reduction factor
            hs = dev_y_width / width_limit	;height reduction factor
            dev_fact = Max([hs,ws])		;select largest factor
            reg_x_width = dev_x_width / dev_fact    ;limit number of device
            reg_y_width = dev_y_width / dev_fact    ; pixels in interp region

	;
	; An index grid for interpolation is constructed.  INTERPOLATE is
	;	used to interpolate the image and the scales.
	;
            x_factor = reg_x_width / im_x_width
            y_factor = reg_y_width / im_y_width
            x_offset = (x_factor - 1.0) / x_factor / 2.0
            y_offset = (y_factor - 1.0) / y_factor / 2.0
            xi = FIndGen(reg_x_width) / x_factor - x_offset	;x interp index
            yi = FIndGen(reg_y_width) / y_factor - y_offset	;y interp index

            im = Interpolate(image, xi, yi, /Grid, Missing=maskvalue)

	EndIf Else Begin
	;
	; With scalable pixels the image can be put directly on the TV and
	;	the scales don't need any reworking.  But, an index grid 
	;	is still constructed so we can computer the scales.
	;
            xi = FIndGen(im_x_width) 		;x interp index
            yi = FIndGen(im_y_width) 		;y interp index

	    im = image

	EndElse

    EndIf Else Begin		;scalable pixels
    ;
    ; Pixels that are not scalable require that we actually increase the
    ;	size of the image to match the device widths.  Whether or not we
    ;	are interpolating the pixels, we do need to interpolate the
    ;	axis scales.  The following computes an interpolation grid.
    ;
      If ( Keyword_Set(no_expand) ) Then Begin
      ;
      ; If /No_Expand is set then return original image and scales.  The
      ;		ploting region is reduced to the image width.
      ;
        xi = FIndGen(im_x_width) 		;x interp index
        yi = FIndGen(im_y_width) 		;y interp index
	im = image
	p(2) = P(0) + im_x_width
	p(3) = p(1) + im_y_width

      EndIf Else Begin				;expand image
        x_factor = dev_x_width / im_x_width
        y_factor = dev_y_width / im_y_width
        x_offset = (x_factor - 1.0) / x_factor / 2.0
        y_offset = (y_factor - 1.0) / y_factor / 2.0
        xi = FIndGen(dev_x_width) / x_factor - x_offset	;x interp index
        yi = FIndGen(dev_y_width) / y_factor - y_offset	;y interp index

	If ( Keyword_Set(interp) ) Then Begin
	;
	; An index grid for interpolation is constructed.  INTERPOLATE is
	;	used to interpolate the image and the scales.
	;
            im = Interpolate(image, xi, yi, /Grid, Missing=maskvalue)

	EndIf Else Begin
	;
	; If interpolation is not required then we use a "simple" POLY_2D
	;	warping.  This works since without interpolation POLY_2D
	;	uses truncation (not nearest neighbor) to determine the
	;	value of the interpolated pixels.  Thus, the pixel center
	;	shifting logic is not necessary.
	    im = Poly_2D(image, [[0,0],[1.0/x_factor,0]], $
				[[0,1.0/y_factor],[0,0]], $
				0, dev_x_width, dev_y_width)
	EndElse
      EndElse

    EndElse	;not scalable

;
; Compute the expanded axis scales.
;
    out_xs = xi * xs_delta + xs(0)
    out_ys = yi * ys_delta + ys(0)
    Return, im
End
