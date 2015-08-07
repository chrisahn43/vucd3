@ppxf_wallace.pro
PRO LINEMAPS_SUB,min_sn=min_sn,CHECK=check,DOERROR=doerror

rv0=-61. ;zeropoint in velocity
IF NOT (KEYWORD_SET(min_sn)) THEN min_sn=15.

infile='ngc404_combine_near9_fluxcal.fits'
outfile='linemaps/ngc404_combine_near9_fluxcal_contsub.fits'
ifu=READFITS(infile,head,ext=1)
var=READFITS(infile,head2,ext=2)
tag='fc'
size=SIZE(ifu,/DIM)
lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
lambda=FINDGEN(size[2])*dlambda+lambda0

brgimage=FLTARR(size[0:1])
ca8image=FLTARR(size[0:1])
ca8image_off1=FLTARR(size[0:1])
ca8image_off2=FLTARR(size[0:1])

IF KEYWORD_SET(DOERROR) THEN BEGIN
    brgerrim=FLTARR(size[0:1])
    ca8errim=FLTARR(size[0:1])
ENDIF


snind=WHERE(lambda GT 22140. AND lambda LT 22866,nsnind)
fitind=WHERE(lambda GT 2.285e4 AND lambda LT 2.39e4,nfitind)
outifu=FLTARR(size[0],size[1],nfitind)
FOR i=3,size[0]-4 DO BEGIN
FOR j=3,size[1]-4 DO BEGIN
signalnoise=MEDIAN(ifu[i,j,snind])/MEDIAN(SQRT(var[i,j,snind]))
;print,signalnoise
zero=WHERE(ifu[i,j,*] EQ 0. OR var[i,j,*] EQ 0.,nzero)
IF (nzero LT 100 AND signalnoise GT 3.) THEN BEGIN
    ppxf_wallace,lambda,REFORM(ifu[i,j,*]),REFORM(var[i,j,*]),4.1,outrv=vel,outdisp=disp,outsol=sol,initvel=rv0,bestfit=bestfit,loglamout=loglamout,/quiet
    
    lambdaout=EXP(loglamout)
    dlambdaout=ABS(lambdaout-lambdaout[1:*])
    dlambdaout=[dlambdaout[0],dlambdaout]
    bestfit=bestfit*dlambda/dlambdaout
    bestfitlin=RESAMPLE_SPECTRA(bestfit,lambdaout,lambda)
    subspec=REFORM(ifu[i,j,*])-bestfitlin
    outifu[i,j,*]=subspec[fitind]
    plot,lambda,subspec,xrange=[2.31e4,2.33e4]
    caind=WHERE(lambda GT 23208. AND lambda LT 23222.)
    off1ind=WHERE(lambda GT 23180. AND lambda LT 23194.,noff1)
    off2ind=WHERE(lambda GT 23235. AND lambda LT 23249.,noff2)
    oplot,lambda[caind],subspec[caind],psym=4
    ca8image[i,j]=TOTAL(subspec[caind])*dlambda
    ca8image_off1[i,j]=TOTAL(subspec[off1ind])*dlambda
    ca8image_off2[i,j]=TOTAL(subspec[off2ind])*dlambda

    IF KEYWORD_SET(CHECK) THEN a=get_kbrd(1)
    IF KEYWORD_SET(DOERROR) THEN BEGIN
        nruns=100
        indarr=FLTARR(nruns)
        FOR k=0,nruns-1 DO BEGIN
            newspectra=subspec+SQRT(var[i,j,*])*RANDOMN(seed,size[2])
            indarr[k]=TOTAL(newspectra[caind])*dlambda
        ENDFOR 
        ca8errim[i,j]=STDDEV(indarr)
    ENDIF

ENDIF
ENDFOR
ENDFOR

WRITEFITS,'linemaps/allpix_'+tag+'_ca8.fits',ca8image
WRITEFITS,'linemaps/allpix_'+tag+'_ca8_off1.fits',ca8image_off1
WRITEFITS,'linemaps/allpix_'+tag+'_ca8_off2.fits',ca8image_off2


outhead=head
SXADDPAR,outhead,'CRVAL3',lambda[fitind[0]]
WRITEFITS,outfile,outifu,outhead


IF KEYWORD_SET(DOERROR) THEN BEGIN
    WRITEFITS,'linemaps/allpix_'+tag+'_ca8sn.fits',ca8image/ca8errim
;    WRITEFITS,'linemaps/allpix'+tag+'_cosn.fits',coimage/coerrim

ENDIF


END

