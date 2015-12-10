Pro	rdisplay, image0, $
		xs0, ys0, $
		Title=t, XTitle=xt, YTitle=yt, $
		SUBTITLE=st, $
		MIN=mn, MAX=mx, $
		BOT=bot, TOP=top, $
		LOG=log_scaling, $
		LEVELS=l, $
		BREAKS=breaks, $
		WRAP=wrap, $
		ASPECT=aspect, $
		INTERPOLATE=interp, $
		MASKVALUE=maskvalue, $
		PSFINE=psfine, $
		NO_EXPAND=no_expand, $
		NOERASE=noerase, $
                HELP=help, $
                 XTICKFORMAT=xtickformat, $
                 YTICKFORMAT=ytickformat, $
                 reverse=reverse, $
                 _extra=extrapars, true=true, quiet=quiet, silent=silent, $
		range=range, position=position, rinterp=rinterp, $
		overplot=overplot, $
		enforce=enforce
if keyword_set(overplot) then begin
  noerase=1
endif
sz=size(image0, /dim)
image=image0
if keyword_set(xs0) then xs=xs0
if keyword_set(ys0) then ys=ys0

if keyword_set(enforce) AND keyword_set(xs) AND $
	keyword_set(ys) AND keyword_set(extrapars) then begin
  tgnames=tag_nameS(extrapars)
  m1=0
  m2=0
  for i=0, n_elementS(tgnames)-1 do m1+=strmatch(tgnames[i],'xr*',/fold)
  for i=0, n_elementS(tgnames)-1 do m2+=strmatch(tgnames[i],'yr*', /fold)
  if m1 then begin

     itag=where(strmatch(tgnames, 'xr*', /fold))
     wh=where(xs LT extrapars.(itag)[1] AND xs GT extrapars.(itag)[0])
     image=image[wh, *]
     xs=xs[wh]
  endif    
  if m2 then begin
     itag=where(strmatch(tgnames, 'yr*', /fold))
     wh=where(ys LT extrapars.(itag)[1] AND ys GT extrapars.(itag)[0])
     image=image[*,wh]
     ys=ys[wh]
  endif    
endif


if keyword_set(rinterp) then begin
  xs1= (findgen(sz[0]*rinterp+1-rinterp )/rinterp)#replicate(1, sz[0]*rinterp+1-rinterp)
  ys1= (findgen(sz[1]*rinterp+1-rinterp)/rinterp)##replicate(1, sz[1]*rinterp+1-rinterp)
  image=bilinear(image, xs1, ys1)
  image[where(finite(image) EQ 0)]=max(image)+1

  xs=interpol(xs, findgen(n_elements(xs)), xs1[*,0])
  ys=interpol(ys, findgen(n_elements(ys)), ys1[0,*])
endif 

;if keyword_set(xrange) EQ 0 then xrange=[0, sz[0]-1]
;if keyword_set(yrange) EQ 0 then xrange=[0, sz[1]-1]

pstore=!P
if keyword_set(quiet) then silent=1
if keyword_set(range) then begin
   mn=min(range)
   mx=max(range)
endif
if keyword_set(position) then !p.position=position
if total(!p.multi) NE 0 then noerase=1
SccsId = '@(#)display.pro 3.9 95/08/31 14:42:48 '
;+
; NAME:
;	RDISPLAY
;
; PURPOSE:
;	This procedure will display an image with the TV command that fills
;	the plotting window.  It handles scale, annotations, X and PostScript
;	devices, aspect ratios, logarithmic scaling, and interpolation.  
;	Masked values and values below <MIN> are mapped to !P.BACKGROUND.
;	Values above <MAX> are mapped to !P.COLOR.
;
; CATEGORY:
;	Image display.
;
; CALLING SEQUENCE:
;	RDISPLAY, Image, XS, YS
;
; INPUTS:
;	Image:	Two-dimensional array to be displayed.
;
; OPTIONAL INPUTS:
;	XS:	Vector of x-axis values.  The length must equal the number of
;		rows in <Image>
;
;	YS:	Vector of y-axis values.  The length must equal the number of
;		columns in <Image>
;
; KEYWORD PARAMETERS:
;	TITLE=	Set this keyword to a string containing the title annotation
;		to be used by PLOT.
;
;	XTITLE=	Set this keyword to a string containing the x-axis annotation
;		to be used by PLOT.
;
;	YTITLE=	Set this keyword to a string containing the y-axis annotation
;		to be used by PLOT.
;
;	SUBTITLE=
;		Set this keyword to a string containing the subtitle annotation
;		to be used by PLOT.
;
;	ASPECT=	Set this keyword to the aspect ratio (width/height) of the
;		pixels.  /ASPECT is the same as ASPECT=1 and produces square
;		pixels.
;
;	/NOAXIS: if set then no axes are plotted (added by rld)
;	
;	/INTERPOLATE:
;		Set this switch to enable bilinear interpolation for pixels
;		in the expanded image.  See /PS_FINE for information
;		on using this switch on a PostScript device.
;
;	RINTERP:set to inerpolate the image... set to be the factor
;		by which you interpolate
;
;	MASKVALUE=
;		Set this keyword to the value that pixels with bad data or
;		no data have been flagged with.  These will be mapped to 0B.
;
;	MIN=	The minimum value of <Image> to be considered.  If MIN is not
;		provided, <Image> is searched for its minimum value.  All
;		values less than MIN are set to !P.BACKGROUND in the Result.
;
;	MAX=	The maximum value of <Image> to be considered.  If MAX is not
;		provided, <Image> is searched for its maximum value.  All
;		values greater than MAX are set to !P.COLOR in the Result.
;
;	RANGE=	alternate method of setting MIN and MAX. Set as [MIN, MAX]
;
;	BOT=	The minimum value of the scaled result.  If TOP is not
;		specified, 0B is used.
;
;	TOP=	The maximum value of the scaled result.  If TOP is not
;		specified, 255B is used.
;
;	LEVELS=	Set this keyword to a vector of data value boundaries between
;		which all elements of <Image> have the same scaled byte
;		value.  e.g. LEVELS=[0,1,2,5] maps all values below 0 and
;		above 5 to 0B, map values between 0 and 1 to 1B, map values
;		between 1 and 2 to 128B, and map values between 2 and 5 to
;		255B.  This does not plot contours.
;
;	/LOG:	Set this switch to cause a logarithmic mapping.  This is
;		overridden by the LEVELS keyword.
;
;	PS_FINE=
;		Set to the approximate number of elements along an axis to 
;		use in interpolating on a PostScript device.  Interpolation
;		to the full indexable range of a PostScript device would be
;		unwise.  This is only useful with /INTERPOLATE and
;		will increase the size of the PostScript file.  The default
;		size is 256 if PS_FINE= is not set or <=0, 512 if 1 < PS_FINE=
;		< 512, and take from the keyword if PS_FINE > 512.
;
;	/NOERASE:
;		Set the switch to prevent output device from being erased
;		before the image, scales, and annotations are displayed.
;	/OVERPLOT:
;		Set this switch to allow overplotting of inlays if you are
;		using !P.multi
;	/NO_EXPAND:
;		Set this switch to prevent the image from being expanded
;		to fill the plotting window.  Scaling to byte type is still
;		performed.
;	/QUIET:	Set to suppress printing to command line
; SIDE EFFECTS:
;	TV display is altered.
;
; SUBORDINATE ROUTINES:
;	IMGSCL()
;	IMGEXP()
;
; RESTRICTIONS:
;	This routine may work for other devices, but it has only been tested
;	on 'X' and 'PS'.
;
; PROCEDURE:
;	Straight forward.  :-)
;
; EXAMPLE:
;	LoadCT, 3
;	image = SHIFT(DIST(20, 20), 10, 10)
;	scale = FINDGEN(20) - 10.0
;	DISPLAY, image, scale, scale, /INTERPOLATE, TITLE='!6Smooth Slope', $
;		/ASPECT
;	;Use CONTOUR with /OVERPLOT to overlay contours.
;	CONTOUR, image, scale, scale, LEVELS=1.0+FINDGEN(4)*2.0, /OVERPLOT
;
;	DISPLAY		;prints out a "Usage:" line
;
; MODIFICATION HISTORY:
; 	Written by:	Fen Tamanaha, July 10, 1993.  Release 3.1
;	July 13, 1993	Fen: (3.2) Fixed /No_Expand
;	July 16, 1993	Fen: (3.3) Really fixed /No_Expand
;       101210 - added QUIET and RANGE keywords RdS, changed name to rdisplay
;-
; 940425 - Added AUTOSTRETCH= to perform histgram equaliazation of the
;		color table.  This whole automated scaling ability needs
;		to be thought out much more carefully.  Right now, I think
;		equalizing the byte scaled image without setting the
;		clipping the scaling range causes as loss of dynamic
;		range in the displayed colors.
; 950717 - Modified to use the new IMGSCL which does the full byte range
;		and maps to !P.BACKGROUND and !P.COLOR for out of bounds
;		values.
; 950720 - Added SUBTITLE= keyword.
; 950801 - PS_FINE= value read from keyword.
; 950816 - Empty subtitle passed to PLOT if SUBTITLE= not defined.
; 950831 - Add pass-thru support for IMGSCL() keywords WRAP= and BREAKS=.
; 030125 - Autostretch feature removed - it depends on comp_dist which
;            cannot be found - DPF
; 101210 - added QUIET and RANGE keywords RdS, changed name to rdisplay
; 040611 - now actually works with position keyword, RdS

    On_Error, 0

;
; Validate arguments.
;
    nparms = N_Params()
    If ( Keyword_Set(help) ) Then nparms = 0	;force a "Usage:" line
    Case ( nparms ) Of
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
            Message, 'Usage: DISPLAY, image [,xs [,ys]] [,TITLE=] [,XTITLE=] [,YTITLE=]', /Info
	    Message, '           [,MIN=] [,MAX=] [,/LOG] [,LEVELS=]', /Info
            Message, '           [,ASPECT=] [,/INTERPOLATE] [MASKVALUE=]', /Info
	    Message, '           [,/NO_EXPAND] [,/NOERASE] [,/PSFINE]', /Info
            Return
        End
    EndCase

;
; Check parameter values and assign default values.
;
    If ( N_Elements(mn) GT 0 ) Then Begin
	minval = mn(0)
    EndIf Else Begin
	minval = Min(image, MAX=maxval)
    EndElse

    If ( N_Elements(mx) GT 0 ) Then Begin
	maxval = mx(0)
    EndIf

    If ( N_Elements(bot) LE 0 ) Then Begin
	bot = 0B
    EndIf

    If ( N_Elements(top) LE 0 ) Then Begin
	top = !D.Table_Size-1
    EndIf

;
; The plotting device must be erased to reset the system variables so that
;	IMGEXP will get the default values.  The /NOERASE keyword should
;	be used to prevent this.  One typical situation is when DISPLAY
;	is called after a !P.MULTI change.  An ERASE at this point would
;	destroy the previous plots.
;
    If ( Not Keyword_Set(noerase) ) Then Begin
	Erase
    EndIf

;
; If /PSFINE is set then up the intermediate interpolated image width.
;	This only has an effect on PostScript output.
;
    If (Keyword_Set(psfine) ) Then Begin
	If ( psfine GT 512 ) Then Begin
	    psis = psfine(0)
	EndIf Else If ( psfine GT 256 ) Then Begin
	    psis = 512.0
	EndIf
    EndIf
    if keyword_set(silent) EQ 0 then Print, "Expanding image ..."

    im =rImgExp(image, xs, ys, xscale, yscale, xrange, yrange, $
		Aspect=aspect, Interpolate=Keyword_Set(interp), $
		MaskValue=maskvalue, Position=dev_pos, PS_Interp_Size=psis, $
		No_Expand=Keyword_Set(no_expand), _extra=extrapars, $
		noerase=noerase, overplot=overplot)

    sz = Size(im)
    im_x_width = Float(sz(1))                   ;image width
    im_y_width = Float(sz(2))                   ;image height
;
; Determine the device coordinates of tahe plotting regions.
;

    dev_x_width = dev_pos[2] - dev_pos[0] + 1
    dev_y_width = dev_pos[3] - dev_pos[1] + 1
    If ( (im_x_width GT dev_x_width) Or (im_y_width GT dev_y_width) ) Then Begin
	Message, 'Error: Scaled image is larger than plotting window.'
    EndIf

;
; Convert a non-byte type image to byte with IMGSCL.  Masked values and
;	values below <minval> are mapped to !P.BACKGROUND.  Values above
;	<maxval> are mapped to !P.COLOR.
;

    If ( sz(sz(0)+1) GT 1 ) Then Begin
	if keyword_set(silent) EQ 0 then Print, "Scaling image ..."
; DPF 6 May 2000
        IF keyword_set(reverse) THEN BEGIN 
           byte_im = ImgScl(-im, Min=-maxval, Max=-minval, Top=top-bot, $
; DPF 10 Jan 2000
;			Breaks=breaks, Wrap=wrap, $
;			Wrap=wrap, $
                        Log=log_scaling, Levels=l, MaskValue=maskvalue) + $
                          Byte(bot)
        ENDIF ELSE BEGIN 
           byte_im = ImgScl(im, Min=minval, Max=maxval, Top=top-bot, $
                            Log=log_scaling, Levels=l, MaskValue=maskvalue) + $
             Byte(bot)
        ENDELSE  
    EndIf Else Begin
;	Message, '<Image> is already byte type. No scaling done.', /Info
	byte_im = im
    EndElse

;
; Put the image on the TV.
;
    TV, byte_im, /Device, dev_pos(0), dev_pos(1), $
		XSize=dev_pos(2)-dev_pos(0), YSize=dev_pos(3)-dev_pos(1)



;
; Manage the title and axis labels.
;
    If ( Keyword_Set(t) ) Then Begin
        title = String(t)
    EndIf Else Begin
        title = ' '
    EndElse

    If ( Keyword_Set(xt) ) Then Begin
        xtitle = String(xt)
    EndIf Else Begin
        xtitle = ' '
    EndElse

    If ( Keyword_Set(yt) ) Then Begin
        ytitle = String(yt)
    EndIf Else Begin
        ytitle = ' '
    EndElse

    If ( Keyword_Set(st) ) Then Begin
        subt = String(st)
    EndIf Else Begin
        subt = ' '
    EndElse

    IF n_elements(xtickformat) EQ 0 THEN xtickformat = ''
    IF n_elements(ytickformat) EQ 0 THEN ytickformat = ''

;
; Overplot annotations.
;
	if size(extrapars, /type) EQ 8 then begin
      if total(strmatch(tag_names(extrapars), 'xs*',/fold)) NE 0 then begin
                ind=where(strmatch(tag_names(extrapars), 'xst*', /fold))
                if (extrapars.(ind) AND 1) EQ 0 then extrapars.(ind)+=1
    endif
    if total(strmatch(tag_names(extrapars), 'xs*',/fold)) NE 0 then begin
                ind=where(strmatch(tag_names(extrapars), 'yst*', /fold))
                if (extrapars.(ind) AND 1) EQ 0 then extrapars.(ind)+=1
    endif
    endif

    Plot, [0,1], /NoErase, /NoData, XStyle=1, YStyle=1, $
                /Device, Position=dev_pos, $
                XRange=xrange, YRange=yrange, $
                Title=title, XTitle=xtitle, YTitle=ytitle, $
		SubTitle=subt, xtickformat=xtickformat, $
                ytickformat=ytickformat,font=0, _extra=extrapars

pmulti=!p.multi
!p=pstore
!p.multi=pmulti
    Return
End
