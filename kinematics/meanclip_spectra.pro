FUNCTION MEANCLIP_SPECTRA,alllambda,specarr,SIGMA=sigma,MAXREJ=maxrej,MAXITER=maxiter,NPOINTS=npoints
;spectra should have different spectra specified by first dimension,
;flux as the second
dim=SIZE(specarr,/dim)
nspectra=LONG(dim[0])
nlambda=dim[1]

lambdaref=alllambda[0,*]
minlambda=MIN(lambdaref)
dlambda=lambdaref[2]-lambdaref[1]
highreslambda=FINDGEN((nlambda+1)*10)*dlambda/10.0+minlambda

IF NOT (KEYWORD_SET(sigma)) THEN sigma=2.0
IF NOT (KEYWORD_SET(maxrej)) THEN maxrej=2
IF NOT (KEYWORD_SET(maxiter)) THEN maxiter=2

ogspecarr=specarr
;first regrid normalize the spectra
medianarr=FLTARR(nspectra)
normspecarr=FLTARR(nspectra,nlambda)
FOR i=0,nspectra-1 DO BEGIN
    IF (ABS(lambdaref[0]-alllambda[i,0]) GT 0.01*(lambdaref[1]-lambdaref[0])) THEN BEGIN
        ;for lambdraef=alllambda[i,*] this can be wrong by 5% on cosmic
        ;rays, but its better than 1% for most features
        specarr[i,*]=resample_spectra(specarr[i,*],alllambda[i,*],lambdaref,oversample=10,npoints=npoints)
    ;spline interpolation gets everything better than 1% for spectra
    ;whose wavelengths match the reference, but is wrong in some ways
;    highres=INTERPOL(specarr[i,*],alllambda[i,*],highreslambda,/SPLINE)    
;    specarray[i,*]=INTERPOL(highres,highreslambda,lambdaref,/SPLINE)
;        plot,alllambda[i,*],(specarr[i,*]-ogspecarr[i,*])/ogspecarr[i,*],psym=4,xrange=[2.1e4,2.3e4],yrange=[-0.1,0.1]
    ENDIF
    medianarr[i]=MEDIAN(specarr[i,*])
    normspecarr[i,*]=specarr[i,*]/FLOAT(medianarr[i])
ENDFOR

stddevarr=FLTARR(nlambda)
outspec=FLTARR(nlambda)
npoints=FLTARR(nlambda)
FOR i=0,nlambda-1 DO BEGIN
    stddev=STDDEV(normspecarr[*,i])
    stddevarr[i]=stddev
    bad=WHERE(ABS(normspecarr[*,i]-MEDIAN(normspecarr[*,i])) GT stddev*sigma,nbad,comp=good)
    iter=0
    WHILE (nbad GT 0 AND nbad LT maxrej AND iter LT maxiter) DO BEGIN
        stddev=STDDEV(normspecarr[good,i])
;        print, stddevarr[i],stddev
        bad=WHERE(ABS(normspecarr[*,i]-MEDIAN(normspecarr[good,i])) GT stddev*sigma,nbad,comp=good)
        iter=iter+1
    ENDWHILE

;    IF (nbad GT 1) THEN print,i,nbad
    outspec[i]=MEAN(specarr[good,i])
    npoints[i]=nspectra-nbad
ENDFOR

RETURN,outspec

STOP
END
;Image     Xcen         Ycen          Lambda[1000]
;279       28.5906      31.0347       22133.300
;281       28.3702      29.2584       22133.300
;84       30.5430      27.8772       22134.970
;87       29.6298      25.3319       22134.970
;90       27.7812      24.3000       22134.970
;integer offsets for gemcombine
; 0   0  0
; 0  -2  0
; 2  -3 -1
; 1  -6 -1
;-1  -7 -1
