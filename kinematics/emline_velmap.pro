PRO EMLINE_VELMAP,DOERROR=doerror
;fit gaussian to determine the velocity and dispersion of emission
;lines


c=2.99d5
usesigma=3.8/2.35
centerwave=21218.334 ;wavelength from Black & van Dishoeck 1987, ApJ 322 412 
tag='h2'
;centerwave=24065.914
;tag='h2q1'
;centerwave=21660.9
;tag='brg'
vned=1290.
centerwave=centerwave*(1.+vned/c)
vlsr=13.66


ifu=READFITS('m60-ucd1_combine_feb20.fits',head,ext=1)
var=READFITS('m60-ucd1_combine_feb20.fits',ext=2)
err=SQRT(var)
usehead=head
;junk=READFITS('hst/ngc4395_combine_near9_int_astcor.fits',usehead)
sky=0.;(MEDIAN(MEDIAN(ifu[54:59,2:13,*],dim=1),dim=1)+MEDIAN(MEDIAN(ifu[57:60,49:60,*],dim=1),dim=1)+MEDIAN(MEDIAN(ifu[61:64,45:56,*],dim=1),dim=1)+MEDIAN(MEDIAN(ifu[59:63,4:25,*],dim=1),dim=1))/4.


;sky=MEDIAN(MEDIAN(ifu[54:59,2:13,*],dim=1),dim=1)
;fwhmfile='../sky/ngc4395/fullcube_disp_med.fits'
;fwhm=READFITS(fwhmfile)
imsize=SIZE(ifu,/dim)
sigma_inst=FLTARR(imsize[0:1])
sigma_inst[*]=usesigma

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
interr=FLTARR(imsize[0:1])
velerr=FLTARR(imsize[0:1])
disperr=FLTARR(imsize[0:1])
loadct,5
inregion=WHERE(lambda GT centerwave+2*minlcont AND lambda LE centerwave+2*maxucont,ninregion)
outspec=FLTARR(ninregion)

;plot,lambda,sky,xrange=[centerwave+2*minlcont,centerwave+2*maxucont],/xsty
;oplot,lambda,sky2,color=100
FOR i=0,imsize[0]-1 DO BEGIN
FOR j=0,imsize[1]-1 DO BEGIN
    spec=ifu[i,j,*]-sky
    errspec=err[i,j,*]
    ospec=ifu[i,j,*]
    background=MEDIAN([spec[lind],spec[uind]])
    obackground=MEDIAN([ospec[lind],ospec[uind]])
    noise=STDDEV([spec[lind],spec[uind]])
    outspec=[[outspec],[(spec-background)[inregion]]]

;    IF (i EQ 38 AND j EQ 29) THEN BEGIN
;        titstring=STRING(i)+'_'+STRING(j)
;        plot,lambda,spec-background,title=titstring,xrange=[centerwave+2*minlcont,centerwave+2*maxucont],/xsty,yrange=[-50,300],/ysty
;        oplot,lambda,ospec-obackground,color=50        
;        a=GET_KBRD()
;    ENDIF
    IF (TOTAL(spec[fitind]-background) GT 10.*noise) THEN BEGIN        
;        ppxf_wallace_full,lambda,REFORM(ospec),REFORM(var[i,j,*]),fwhm[i,j],bestfit=bestfit,loglamout=loglamout,initvel=-70.
        titstring=STRING(i)+'_'+STRING(j)
        plot,lambda,spec-background,title=titstring,xrange=[centerwave+2*minlcont,centerwave+2*maxucont],/xsty,yrange=[-50,300],/ysty
        oplot,lambda,ospec-obackground,color=50
;        oplot,loglamout,bestfit-obackground,color=200
        oplot,lambda[fitind],REPLICATE(0.,nfitind),thick=3
        fit=GAUSSFIT(lambda[fitind],spec[fitind]-background,par,nterms=3)
        oplot,lambda[fitind],fit,color=100
        vel=(par[1]-centerwave)/centerwave*c+vlsr
        disp=SQRT(par[2]^2-sigma_inst[i,j]^2)/centerwave*c
;        IF (i EQ 38 AND j EQ 29) THEN a=GET_KBRD()
        IF (ABS(vel) LT 600. AND FINITE(disp) EQ 1) THEN BEGIN
            intmap[i,j]=par[0]
            velmap[i,j]=(par[1]-centerwave)/centerwave*c
            IF (par[2] GT sigma_inst[i,j]) THEN BEGIN
                dispmap[i,j]=SQRT(par[2]^2-sigma_inst[i,j]^2)/centerwave*c
            ENDIF ELSE print,par[2]-sigma_inst[i,j]

;                print,velmap[i,j],dispmap[i,j]
        WAIT,0.01
            IF KEYWORD_SET(DOERROR) THEN BEGIN
                nruns=50
                npixels=N_ELEMENTS(spec)
                int_arr=FLTARR(nruns)
                vel_arr=FLTARR(nruns)
                disp_arr=FLTARR(nruns)
                FOR k=0,nruns-1 DO BEGIN
                    newspec=spec[fitind]+err[i,j,fitind]*RANDOMN(seed,nfitind)
                    efit=GAUSSFIT(lambda[fitind],newspec-background,epar,nterms=3)
                    int_arr[k]=epar[0]
                    vel_arr[k]=(epar[1]-centerwave)/centerwave*c
                    IF (par[2] GT sigma_inst[i,j]) THEN $
                      disp_arr[k]=SQRT(epar[2]^2-sigma_inst[i,j]^2)/centerwave*c $
                      else disp_arr[k]=0.
                ENDFOR                
                interr[i,j]=STDDEV(int_arr)
                velerr[i,j]=STDDEV(vel_arr)
                disperr[i,j]=STDDEV(disp_arr[WHERE(FINITE(disp_arr) EQ 1)])
            ENDIF ;doerror           
        ENDIF ;(ABS(vel) and finite)
    ENDIF ;signal at all
ENDFOR
ENDFOR
SXDELPAR,usehead,'EXTNAME'
SXDELPAR,usehead,'EXTVER'

intfile='linemaps/'+tag+'_intmap.fits'
WRITEFITS,intfile,intmap,usehead
velfile='linemaps/'+tag+'_velmap.fits'
WRITEFITS,velfile,velmap,usehead
dispfile='linemaps/'+tag+'_dispmap.fits'
WRITEFITS,dispfile,dispmap,usehead
intfile='linemaps/'+tag+'_interr.fits'
WRITEFITS,intfile,interr,usehead
velfile='linemaps/'+tag+'_velerr.fits'
WRITEFITS,velfile,velerr,usehead
dispfile='linemaps/'+tag+'_disperr.fits'
WRITEFITS,dispfile,disperr,usehead

allfile='linemaps/'+tag+'_allspec.fits'
WRITEFITS,allfile,outspec

STOP
END
