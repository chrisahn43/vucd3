@kband_indices.pro
@ppxf_wallace.pro
PRO LINEMAPS,min_sn=min_sn,CHECK=check,DOERROR=doerror

rv0=308. ;zeropoint in velocity
IF NOT (KEYWORD_SET(min_sn)) THEN min_sn=15.

infile='m60-ucd1_combine_feb20.fits'
ifu=READFITS(infile,head,ext=1)
var=READFITS(infile,head,ext=2)
tag='feb20'
size=SIZE(ifu,/DIM)
lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
lambda=FINDGEN(size[2])*dlambda+lambda0

caimage=FLTARR(size[0:1])
coimage=FLTARR(size[0:1])
mgimage=FLTARR(size[0:1])
mgaimage=FLTARR(size[0:1])
naimage=FLTARR(size[0:1])
brgimage=FLTARR(size[0:1])
caviiiimage=FLTARR(size[0:1])

IF KEYWORD_SET(DOERROR) THEN BEGIN
    caerrim=FLTARR(size[0:1])
    coerrim=FLTARR(size[0:1])
    mgerrim=FLTARR(size[0:1])
    mgaerrim=FLTARR(size[0:1])
    naerrim=FLTARR(size[0:1])
    brgerrim=FLTARR(size[0:1])
ENDIF

OPENW,1,'linemaps/allpix_'+tag+'.dat',width=1000.
snind=WHERE(lambda GT 22140. AND lambda LT 22866,nsnind)
FOR i=3,size[0]-4 DO BEGIN
FOR j=3,size[1]-4 DO BEGIN
signalnoise=MEDIAN(ifu[i,j,snind])/MEDIAN(SQRT(var[i,j,snind]))
;print,signalnoise
zero=WHERE(ifu[i,j,*] EQ 0. OR var[i,j,*] EQ 0.,nzero)
IF (nzero LT 100 AND signalnoise GT 3.) THEN BEGIN
    ppxf_wallace,lambda,REFORM(ifu[i,j,*]),REFORM(var[i,j,*]),4.1,outrv=vel,outdisp=disp,outsol=sol,initvel=rv0,/quiet
    IF (ABS(vel-rv0) GT 70.) THEN vel=rv0
    kind=kband_indices(lambda,REFORM(ifu[i,j,*]),velocity=vel,/PLOT)
    IF KEYWORD_SET(CHECK) THEN a=get_kbrd(1)
    mgimage[i,j]=kind[0]
    brgimage[i,j]=kind[1]
    naimage[i,j]=kind[2]
    caimage[i,j]=kind[3]
    mgaimage[i,j]=kind[4]    
    coimage[i,j]=kind[5]
    caviiiimage[i,j]=kind[5]
    IF KEYWORD_SET(DOERROR) THEN BEGIN
        nruns=100
        allind=FLTARR(N_ELEMENTS(kind),nruns)
        FOR k=0,nruns-1 DO BEGIN
            newspectra=ifu[i,j,*]+SQRT(var[i,j,*])*RANDOMN(seed,size[2])
            oneind=kband_indices(lambda,newspectra,velocity=vel)
            allind[*,k]=oneind
        ENDFOR 
        mgerrim[i,j]=STDDEV(allind[0,*])
        brgerrim[i,j]=STDDEV(allind[1,*])
        naerrim[i,j]=STDDEV(allind[2,*])
        caerrim[i,j]=STDDEV(allind[3,*])
        mgaerrim[i,j]=STDDEV(allind[4,*])
        coerrim[i,j]=STDDEV(allind[5,*])        
        print,i,j,kind[0]/STDDEV(allind[0,*]),kind[1]/STDDEV(allind[1,*]),kind[2]/STDDEV(allind[2,*]),kind[3]/STDDEV(allind[3,*]),kind[4]/STDDEV(allind[4,*]),(1.0-kind[5])/STDDEV(allind[5,*]),signalnoise,kind[3],STDDEV(allind[3,*]),FORMAT='(2I4,7F7.1,2F8.3)'
        printf,1,i,j,kind[0],STDDEV(allind[0,*]),kind[1],STDDEV(allind[1,*]),kind[2],STDDEV(allind[2,*]),kind[3],STDDEV(allind[3,*]),kind[4],STDDEV(allind[4,*]),kind[5],STDDEV(allind[5,*]),signalnoise,FORMAT='(2I4,12F9.4,F7.1)'
    ENDIF
ENDIF
ENDFOR
ENDFOR
CLOSE,1
FREE_LUN,1
;WRITEFITS,'linemaps/allpix_'+tag+'_mg.fits',mgimage
;WRITEFITS,'linemaps/allpix_'+tag+'_mga.fits',mgaimage
;WRITEFITS,'linemaps/allpix_'+tag+'_brg.fits',brgimage
;WRITEFITS,'linemaps/allpix_'+tag+'_na.fits',naimage
;WRITEFITS,'linemaps/allpix_'+tag+'_ca.fits',caimage
;WRITEFITS,'linemaps/allpix_'+tag+'_co.fits',coimage
WRITEFITS,'linemaps/allpix_'+tag+'_caviii.fits',caviiimage

IF KEYWORD_SET(DOERROR) THEN BEGIN
    WRITEFITS,'linemaps/allpix'+tag+'_mgerr.fits',mgerrim
    WRITEFITS,'linemaps/allpix'+tag+'_mgaerr.fits',mgaerrim
    WRITEFITS,'linemaps/allpix'+tag+'_brgerr.fits',brgerrim
    WRITEFITS,'linemaps/allpix'+tag+'_naerr.fits',naerrim
    WRITEFITS,'linemaps/allpix'+tag+'_caerr.fits',caerrim
    WRITEFITS,'linemaps/allpix'+tag+'_coerr.fits',coerrim
    WRITEFITS,'linemaps/allpix'+tag+'_cosn.fits',coimage/coerrim
    WRITEFITS,'linemaps/allpix'+tag+'_brgsn.fits',brgimage/brgerrim
ENDIF


END

