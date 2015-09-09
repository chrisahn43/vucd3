PRO INT_SPECTRA_ANNULUS,minradius,maxradius,xcen,ycen
;create an integrated spectrum with option to subtract off sky

IF (NOT KEYWORD_SET(minradius)) THEN minradius=0.
IF (NOT KEYWORD_SET(maxradius)) THEN maxradius=1.

infile='../vucd3_combine_best8.fits'
outfile='vucd3_combine_best8_int_fc'+STRTRIM(FIX(minradius),2)+'-'+STRTRIM(FIX(maxradius),2)+'.fits'

cube=DOUBLE(READFITS(infile,head,ext=1,/SILENT))
varcube=DOUBLE(READFITS(infile,head,ext=2,/SILENT))
cubesize=SIZE(cube,/dim)
lambda0=SXPAR(head,'CRVAL3')
dlambda=SXPAR(head,'CD3_3')
lambda=FINDGEN(cubesize[2])*dlambda+lambda0

outspec=DBLARR(cubesize[2])
outrawspec=DBLARR(cubesize[2])
outvar=DBLARR(cubesize[2])
outsky=DBLARR(cubesize[2])
FOR i=0,cubesize[2]-1 DO BEGIN
    skyval=MEDIAN(cube[5:20,67:81,i])
    outsky[i]=skyval*!PI*(maxradius^2-minradius^2)
    skyval=0

    APER,(cube[*,*,i]-skyval),xcen,ycen,maxflux,fluxerr,sky,skyerr,1.0,[maxradius],[20,25],[-32000,32000],/SILENT,SETSKYVAL=0.0,/flu,/EXACT,/NAN
    IF (minradius GT 0.5) THEN BEGIN
        APER,(cube[*,*,i]-skyval),xcen,ycen,minflux,fluxerr,sky,skyerr,1.0,[minradius],[20,25],[-32000,32000],/SILENT,SETSKYVAL=0.0,/flu,/EXACT,/NAN
        flux=maxflux-minflux
    ENDIF ELSE flux=maxflux

    APER,(cube[*,*,i]),xcen,ycen,maxrawflux,rawfluxerr,rawsky,rawskyerr,1.0,[maxradius],[20,25],[-32000,32000],/SILENT,SETSKYVAL=0.0,/fl,/EXACT,/NAN
    IF (minradius GT 0.5) THEN BEGIN
        APER,(cube[*,*,i]),xcen,ycen,minrawflux,rawfluxerr,rawsky,rawskyerr,1.0,[minradius],[20,25],[-32000,32000],/SILENT,SETSKYVAL=0.0,/flu,/EXACT,/NAN
        rawflux=maxrawflux-minrawflux
    ENDIF ELSE rawflux=maxrawflux

    APER,(varcube[*,*,i]),xcen,ycen,varflux,varfluxerr,varsky,varskyerr,1.0,[maxradius],[20,25],[-32000,32000],/SILENT,SETSKYVAL=0.0,/flux,/EXACT,/NAN 
    outspec[i]=flux
    outvar[i]=varflux   
    outrawspec[i]=rawflux
ENDFOR


ind=WHERE(lambda LT 2.01e4 or lambda GT 2.43e4)
outspec[ind]=0.0
outvar[ind]=0.0
outsky[ind]=0.0

outarr=[[outspec],[outvar],[outsky]]

mkhdr,headout,outarr
sxaddpar,headout,'CRVAL1',lambda0
sxaddpar,headout,'CRPIX1',1
sxaddpar,headout,'CD1_1',dlambda
SXADDPAR,head,'BUNIT','erg/cm2/s/A'
WRITEFITS,outfile,outarr,headout


;plot,lambda,MEDSMOOTH(outspec,5),xrange=[2.0e4,2.45e4],/xsty,ysty=16
plot,lambda,outspec/MEDIAN(outspec),xrange=[2.26e4,2.43e4],/xsty,ysty=16,yrange=[0.6,1.2]
redshift=1.0043
plots,[2.293529e4,2.293529e4]*redshift,!Y.CRANGE
plots,[2.322656e4,2.322656e4]*redshift,!Y.CRANGE
plots,[2.352458e4,2.352458e4]*redshift,!Y.CRANGE
plots,[2.382957e4,2.382957e4]*redshift,!Y.CRANGE
plots,[2.2078e4,2.2078e4]*redshift,!Y.CRANGE
;plots,[2.26435e4,2.26435e4]*redshift,!Y.CRANGE
;[FeIII] lines
;plots,[2.1551e4,2.1551e4]*redshift,!Y.CRANGE,color=220
;plots,[2.2184e4,2.2184e4]*redshift,!Y.CRANGE,color=220
;plots,[2.2427e4,2.2427e4]*redshift,!Y.CRANGE,color=220
;plots,[2.3485e4,2.3485e4]*redshift,!Y.CRANGE,color=220

;plots,[2.1895e4,2.1895e4]*redshift,!Y.CRANGE,color=220



;READCOL,'../h2_wavnum.dat',wavelength,FORMAT='F'
;nlines=N_ELEMENTS(wavelength)
;FOR i=0,nlines-1 DO plots,[wavelength[i],wavelength[i]]*redshift,!Y.CRANGE,color=100

ind=WHERE(lambda GT 2.28e4 and lambda LT 2.40e4)
print,MEDIAN(outspec[ind]/SQRT(outvar[ind]))
a=GET_KBRD()
;STOP
END
