PRO EMLINE_VELMAP
;fit gaussian to determine the velocity and dispersion of emission
;lines


c=2.99d5
sigma_inst=4.1/2.35
centerwave=21218.334 ;wavelength from Black & van Dishoeck 1987, ApJ 322 412 
vlsr=13.25

ifu=READFITS('ngc404_combine_near9.fits',head,ext=1)
junk=READFITS('hst/ngc404_combine_int_astcor.fits',usehead)
fwhmfile='../sky/ngc404/fullcube_disp_med.fits'
fwhm=READFITS(fwhmfile)
sigma_inst=fwhm/2.35
imsize=SIZE(ifu,/dim)
lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
lambda=FINDGEN(imsize[2])*dlambda+lambda0

minlcont=-50.
maxlcont=-20.
minucont=20.
maxucont=50.
fitwidth=25.
lind=WHERE(lambda GE centerwave+minlcont AND lambda LE centerwave+maxlcont,nlind)
uind=WHERE(lambda GE centerwave+minucont AND lambda LE centerwave+maxucont,nuind)
fitind=WHERE(lambda GE centerwave-fitwidth AND lambda LE centerwave+fitwidth,nfitind)

intmap=FLTARR(imsize[0:1])
velmap=FLTARR(imsize[0:1])
dispmap=FLTARR(imsize[0:1])
loadct,5
inregion=WHERE(lambda GT centerwave+2*minlcont AND lambda LE centerwave+2*maxucont,ninregion)
outspec=FLTARR(ninregion)
FOR i=0,imsize[0]-1 DO BEGIN
FOR j=0,imsize[1]-1 DO BEGIN
    spec=ifu[i,j,*]
    background=MEDIAN([spec[lind],spec[uind]])
    noise=STDDEV([spec[lind],spec[uind]])
    outspec=[[outspec],[(spec-background)[inregion]]]
    IF (TOTAL(spec[fitind]-background) GT 10.*noise) THEN BEGIN        
        titstring=STRING(i)+'_'+STRING(j)
        plot,lambda,spec-background,title=titstring,xrange=[centerwave+2*minlcont,centerwave+2*maxucont],/xsty
        oplot,lambda[fitind],REPLICATE(0.,nfitind),thick=3
        fit=GAUSSFIT(lambda[fitind],spec[fitind]-background,par,nterms=3)
        oplot,lambda[fitind],fit,color=100
        vel=(par[1]-centerwave)/centerwave*c+vlsr
        disp=SQRT(par[2]^2-sigma_inst[i,j]^2)/centerwave*c
        IF (ABS(vel) LT 300. AND FINITE(disp) EQ 1) THEN BEGIN
            intmap[i,j]=par[0]
            velmap[i,j]=(par[1]-centerwave)/centerwave*c
            IF (par[2] GT sigma_inst[i,j]) THEN BEGIN
                dispmap[i,j]=SQRT(par[2]^2-sigma_inst[i,j]^2)/centerwave*c
            ENDIF ELSE print,par[2]-sigma_inst[i,j]

;                print,velmap[i,j],dispmap[i,j]
;        WAIT,0.05
        ENDIF
;        IF (i EQ 30) and (j EQ 27) THEN STOP
    ENDIF
ENDFOR
ENDFOR
SXDELPAR,head,'EXTNAME'
SXDELPAR,head,'EXTVER'


WRITEFITS,'h2_intmap.fits',intmap,usehead
WRITEFITS,'h2_velmap.fits',velmap,usehead
WRITEFITS,'h2_dispmap.fits',dispmap,usehead
WRITEFITS,'h2_allspec.fits',outspec

STOP
END
