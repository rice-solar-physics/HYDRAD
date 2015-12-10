Function	ImgScl, array, $
			MIN=min_lvl, MAX=max_lvl, $
			TOP=top, $
			LEVELS=l, $
			LOG=log_scl, $
			MASKVALUE=maskvalue


SccsId = '@(#)imgscl.pro 2.1 7/12/93 Fen Tamanaha'
;+
; NAME:
;	IMGSCL
;
; PURPOSE:
;	This function mimics BYTSCL() except that it maps the input range
;	into a byte range from 1 through TOP.  A byte value of 0 is reserved
;	for elements containing MASKVALUE usually assigned for bad pixels or
;	those without data.  The function can also perform logarithmic scaling
;	of the data into byte values.  Use of the LEVELS keyword will scale
;	all value within a given level to a single byte value.
;
; CATEGORY:
;	Image scaling.
;
; CALLING SEQUENCE:
;	Result = IMGSCL(Array)
;
; INPUTS:
;	Array:	Two-dimensional array to be expanded.
;
; KEYWORD PARAMETERS:
;	MIN=	The minimum value of Array to be considered.  If MIN is not
;		provided, Array is searched for its minimum value.  All
;		values less than or equal to MIN are set to 1 in the Result.
;
;	MAX=	The maximum value of Array to be considered.  If MAX is not
;		provided, Array is searched for its maximum value.  All
;		values greater than or equal to MAX are set to TOP in the
;		Result.
;
;	TOP=	The maximum value of the scaled result.  If TOP is not
;		specified, 255 is used. Note that the minimum value of the
;		scaled result is always 1 (NOT 0 as in BYTSCL).
;
;	LEVELS=	Set this keyword to a vector of data value boundaries between
;		which all elements of the Array have the same scaled byte
;		value.  e.g. LEVELS=[0,1,2,5] maps all values below 0 and
;		above 5 to 0B, map values between 0 and 1 to 1B, map values
;		between 1 and 2 to 128B, and map values between 2 and 5 to
;		255B.
;
;	/LOG:	Set this switch to cause a logarithmic mapping.  This is
;		overridden by the LEVELS keyword.
;
;	MASKVALUE=
;		Set this keyword to the value that pixels with bad data or
;		no data have been flagged with.  These will be mapped to 0B.
;
; OUTPUTS:
;	Result:	This function returns a byte array between 1 and TOP for data
;		in Array between MIN and MAX.
;
; RESTRICTIONS:
;
; PROCEDURE:
;	Straight forward.  :-)
;
; EXAMPLE:
;	image = IMGSCL(array, Min=-1, Top=!D.Table_Size-1, /Log, Mask=-9999.0)
;	TV, image
;
;	image = IMGSCL(array, Levels=[0,1,2,4,8,16,32])
;	TV, image		;plot with 6 colors plus the background
;
; MODIFICATION HISTORY:
; 	Written by:	Fen Tamanaha, July 10, 1993.  Release 2.1
;-

    On_Error, 2

;
; Check parameters and keywords.
;
    If ( N_Params() NE 1 ) Then Begin
	Message, 'Usage: Image = ImgScl(array [,MIN=] [,MAX=] [,TOP=] [,LEVELS=]', /Info
	Message, '                      [,MASKVALUE=] [,/LOG])', /Info
	Return, 0B
    EndIf

    If ( (Size(array))(0) LT 1 ) Then Begin
	Message, 'Error: <array> must be an array.'
    EndIf

;
; Use MASKVALUE to determine the indices that contain valid data.
;
    If ( N_Elements(maskvalue) GT 0 )  Then Begin
	valid = Where(array NE maskvalue, valid_count)
	If ( valid_count LT 1 ) Then Begin
	    Message, 'Warning: <array> contains only masked values.', /Cont
	    Return, Byte(0.0 * array)
	EndIf
    EndIf Else Begin
	valid = LIndGen(N_Elements(array))
    EndElse

;
; If not passed, compute maximum and minimum data values for the non-masked
;	elements.
;
    If ( N_Elements(min_lvl) LE 0 ) Then Begin
	min_lvl = Min(array(valid), Max=tmp_max)
    EndIf
    If ( N_Elements(max_lvl) LE 0 ) Then Begin
	If ( N_Elements(tmp_max) EQ 0 ) Then Begin
	    max_lvl = Max(array(valid))
	EndIf Else Begin
	    max_lvl = tmp_max
	EndElse
    EndIf

    If ( N_Elements(top) LE 0 ) Then Begin
	top = 255B
    EndIf

    If ( Keyword_Set(log_scl) ) Then Begin
	scaling_type = 1
	low_clip = 0.01		;lowest value allowed with the ALog10
    EndIf Else Begin
	scaling_type = 0
    EndElse

    n_lvls = N_Elements(l)
    Case ( n_lvls ) Of
	0: Begin
	    ; Ignore this keyword.
	End
	1: Begin
	    Message, 'Error: LEVELS must contain at least 3 entries.'
	End
	2: Begin
	    Message, 'Error: LEVELS must contain at least 3 entries.'
	End
	Else: Begin
	    scaling_type = 2
	End
    EndCase

;
; Perform scaling from <MIN>:<MAX> to <1>:<TOP>.
;
    Case ( scaling_type ) Of
	0: Begin
	    image = BytScl(array, Min=min_lvl, Max=max_lvl, Top=top-1B) + 1B
	End
	1: Begin
	    tmp = 1.0 + (array <max_lvl >min_lvl) - min_lvl
	    tmp = ALog10(Temporary(tmp) >low_clip)
	    image = BytScl(tmp, Top=top-1B) + 1B
	End
	2: Begin
	;
	; Produce grayscale plateaus at the level breaks.  All level boundaries
	;	MUST be specified.  This requires a minimum of 3 elements
	;	in LEVELS (bottom of low bin, transition value, and top of
	;	high bin).  All values outside of these levels is mapped to
	;	the background color table entry of 0B.
	;
	    levels = l(Sort(l))
	    image = Byte(array * 0)
	    blk_size = Float(top + 1) / Float(n_lvls - 2)

	;
	; Assign each successive block the next color index.
	;
	    For step = 0, n_lvls-2 Do Begin
		index = Where( (array GE levels(step)) And $
			(array LE levels(step+1)), count)
		If ( count GT 0 ) Then Begin
		    image(index) = Byte((step)*blk_size) >1B <top
		EndIf Else Begin
		    msg = String(Format='("Warning: No array values between levels ", I2, " and ", I2, ".")', step, step+1)
		    Message, msg, /Cont
		EndElse
		color = Byte((step+1.0)*blk_size)
	    EndFor

	;
	; Set values that are out of range to 0B.
	;
	    index = Where( (array LT levels(0)) Or $
			(array GT levels(n_lvls-1)), count)
	    If ( count GT 0 ) Then Begin
		image(index) = 0B
	    EndIf

	End
	Else: Begin
	    Message, 'Program Bug: Scaling type in error.'
	End
    EndCase

;
; Set the pixels with a value of MASKVALUE to zero.
;
    If ( N_Elements(maskvalue) GT 0 ) Then Begin
	masked = Where(array EQ maskvalue(0), count)
	If ( count GT 0 ) Then Begin
	    image(masked) = 0B
	EndIf
    EndIf

    Return, image
End
