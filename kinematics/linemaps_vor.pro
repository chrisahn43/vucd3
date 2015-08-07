@kband_indices.pro
@ppxf_wallace.pro
PRO LINEMAPS_vor,min_sn=min_sn,CHECK=check,DOERROR=doerror

min_sn=25.

rv0=-61. ;zeropoint in velocity
IF NOT (KEYWORD_SET(min_sn)) THEN min_sn=15.

infile='ngc404_combine_near9.fits'
ifu=READFITS(infile,head,ext=1)
var=READFITS(infile,head,ext=2)
tag='sn25'
size=SIZE(ifu,/DIM)
lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
lambda=FINDGEN(size[2])*dlambda+lambda0
makex,ifu,x,y

minfit=2.28e4  ;part of the target spectra to extract
maxfit=2.395e4
snind=WHERE(lambda GT minfit AND lambda LT maxfit)
signal=FLTARR(size[0:1])
noise=FLTARR(size[0:1])
FOR i=0,size[0]-1 DO BEGIN
FOR j=0,size[1]-1 DO BEGIN
    uspec=ifu[i,j,snind]
    uvar=var[i,j,snind]
    signal[i,j]=MEDIAN(uspec)
    noise[i,j]=MEDIAN(SQRT(uvar))
ENDFOR
ENDFOR

signoise=signal/noise
logsn=ALOG10(signoise)
logsn=logsn-MIN(logsn[WHERE(FINITE(logsn) EQ 1,comp=bad)])
logsn[bad]=0.0

voronoi_2d_binning, x, y, signal, noise, min_sn, binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale, /QUIET,/NO_CVT ;-- S/N=10 failed without this

imsize=SIZE(signal,/dim)
maxval=MAX(signal,maxpos)
;xinit=maxpos MOD imsize[0]
;yinit=maxpos/imsize[0]
;GCNTRD,signal,xinit,yinit,xcen,ycen,4.5
;print,'Center',xcen,ycen
xcen=32.30 & ycen=28.43 ;corrected for dust emission
xoutbin=(xbar-xcen)*0.05
youtbin=(ybar-ycen)*0.05

;bin the spectra together
nbins=N_ELEMENTS(xnode)
bin_spectra=FLTARR(nbins,size[2])
bin_var=FLTARR(nbins,size[2])
measured_sn=FLTARR(nbins)
snind=WHERE(lambda GT 22140. AND lambda LT 22600.)
FOR i=0,nbins-1 DO BEGIN
    ind=WHERE(binnum EQ i,nind)
    xind=ind MOD size[0]
    yind=ind / size[0]
    tempspec=FLTARR(size[2])
    tempvar=FLTARR(size[2])
    FOR j=0,nind-1 DO BEGIN
        tempspec=tempspec+ifu[xind[j],yind[j],*]
        tempvar=tempvar+var[xind[j],yind[j],*]
    ENDFOR
    bin_spectra[i,*]=tempspec;/FLOAT(nind)
    bin_var[i,*]=tempvar
    MEANCLIP,tempspec[snind],meany,4.,subs=subs
    measured_sn[i]=MEAN(tempspec[snind[subs]])/STDDEV(tempspec[snind[subs]])
ENDFOR


caimage=FLTARR(size[0:1])
coimage=FLTARR(size[0:1])
mgimage=FLTARR(size[0:1])
mgaimage=FLTARR(size[0:1])
naimage=FLTARR(size[0:1])
brgimage=FLTARR(size[0:1])
;caviiiimage=FLTARR(size[0:1])
binimage=FLTARR(size[0:1])
IF KEYWORD_SET(DOERROR) THEN BEGIN
    caerrim=FLTARR(size[0:1])
    coerrim=FLTARR(size[0:1])
    mgerrim=FLTARR(size[0:1])
    mgaerrim=FLTARR(size[0:1])
    naerrim=FLTARR(size[0:1])
    brgerrim=FLTARR(size[0:1])
ENDIF

OPENW,1,'linemaps/vor_'+tag+'.dat',width=1000.
snind=WHERE(lambda GT 22140. AND lambda LT 22866,nsnind)
FOR i=0,nbins-1 DO BEGIN
inbin=WHERE(binnum EQ i)
signalnoise=MEDIAN(bin_spectra[i,snind])/MEDIAN(SQRT(bin_var[i,snind]))
;print,signalnoise
zero=WHERE(bin_spectra[i,*] EQ 0. OR bin_var[i,*] EQ 0.,nzero)
IF (nzero LT 100 AND signalnoise GT 3.) THEN BEGIN
    ppxf_wallace,lambda,REFORM(bin_spectra[i,*]),REFORM(bin_var[i,*]),4.1,outrv=vel,outdisp=disp,outsol=sol,initvel=rv0,/quiet
    IF (ABS(vel-rv0) GT 70.) THEN vel=rv0
    kind=kband_indices(lambda,REFORM(bin_spectra[i,*]),velocity=vel,/PLOT)
    IF KEYWORD_SET(CHECK) THEN a=get_kbrd(1)
    mgimage[inbin]=kind[0]
    brgimage[inbin]=kind[1]
    naimage[inbin]=kind[2]
    caimage[inbin]=kind[3]
    mgaimage[inbin]=kind[4]    
    coimage[inbin]=kind[5]
;    caviiiimage[inbin]=kind[6]
    IF KEYWORD_SET(DOERROR) THEN BEGIN
        nruns=100
        allind=FLTARR(N_ELEMENTS(kind),nruns)
        FOR k=0,nruns-1 DO BEGIN
            newspectra=bin_spectra[i,*]+SQRT(bin_var[i,*])*RANDOMN(seed,size[2])
            oneind=kband_indices(lambda,newspectra,velocity=vel)
            allind[*,k]=oneind
        ENDFOR 
        mgerrim[inbin]=STDDEV(allind[0,*])
        brgerrim[inbin]=STDDEV(allind[1,*])
        naerrim[inbin]=STDDEV(allind[2,*])
        caerrim[inbin]=STDDEV(allind[3,*])
        mgaerrim[inbin]=STDDEV(allind[4,*])
        coerrim[inbin]=STDDEV(allind[5,*])        
        print,i,kind[0]/STDDEV(allind[0,*]),kind[1]/STDDEV(allind[1,*]),kind[2]/STDDEV(allind[2,*]),kind[3]/STDDEV(allind[3,*]),kind[4]/STDDEV(allind[4,*]),(1.0-kind[5])/STDDEV(allind[5,*]),signalnoise,kind[3],STDDEV(allind[3,*]),FORMAT='(I4,7F7.1,2F8.3)'
        printf,1,i,kind[0],STDDEV(allind[0,*]),kind[1],STDDEV(allind[1,*]),kind[2],STDDEV(allind[2,*]),kind[3],STDDEV(allind[3,*]),kind[4],STDDEV(allind[4,*]),kind[5],STDDEV(allind[5,*]),signalnoise,FORMAT='(I4,12F9.4,F7.1)'
    ENDIF
ENDIF
ENDFOR
CLOSE,1
FREE_LUN,1
WRITEFITS,'linemaps/vor_'+tag+'_mg.fits',mgimage
WRITEFITS,'linemaps/vor_'+tag+'_mga.fits',mgaimage
WRITEFITS,'linemaps/vor_'+tag+'_brg.fits',brgimage
WRITEFITS,'linemaps/vor_'+tag+'_na.fits',naimage
WRITEFITS,'linemaps/vor_'+tag+'_ca.fits',caimage
WRITEFITS,'linemaps/vor_'+tag+'_co.fits',coimage
WRITEFITS,'linemaps/vor_'+tag+'_bin.fits',binnum

;WRITEFITS,'linemaps/vor_'+tag+'_caviii.fits',caviiimage

IF KEYWORD_SET(DOERROR) THEN BEGIN
    WRITEFITS,'linemaps/vor'+tag+'_mgerr.fits',mgerrim
    WRITEFITS,'linemaps/vor'+tag+'_mgaerr.fits',mgaerrim
    WRITEFITS,'linemaps/vor'+tag+'_brgerr.fits',brgerrim
    WRITEFITS,'linemaps/vor'+tag+'_naerr.fits',naerrim
    WRITEFITS,'linemaps/vor'+tag+'_caerr.fits',caerrim
    WRITEFITS,'linemaps/vor'+tag+'_coerr.fits',coerrim
    WRITEFITS,'linemaps/vor'+tag+'_cosn.fits',coimage/coerrim
    WRITEFITS,'linemaps/vor'+tag+'_brgsn.fits',brgimage/brgerrim
ENDIF


END

