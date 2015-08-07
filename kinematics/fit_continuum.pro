@kband_indices.pro
@ppxf_wallace.pro
PRO FIT_CONTINUUM


infile='m60-ucd1_combine_may18.fits'
ifu=READFITS(infile,head1,ext=1)
var=READFITS(infile,head,ext=2)
err=SQRT(var)
size=SIZE(ifu,/DIM)
lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
lambda=FINDGEN(size[2])*dlambda+lambda0


contimage=FLTARR(size[0:1])
conterrimage=FLTARR(size[0:1])
makex,contimage,xarr,yarr
xcenter=34.19
ycenter=26.23
rarr=SQRT((xarr-xcenter)^2+(yarr-ycenter)^2)

minlambda=2.05e4
maxlambda=2.28e4
lamind=WHERE(lambda GT 2.05e4 AND lambda LT 2.28e4)
lamcut=lambda[lamind]
FOR i=1,size[0]-2 DO BEGIN
FOR j=1,size[1]-2 DO BEGIN
    spec=ifu[i,j,lamind]
    specerr=err[i,j,lamind]
    MEANCLIP,spec,meany,4.,subs=subs
    meanspec=MEAN(spec[subs])
    spec=spec/meanspec   
    specerr=spec/meanspec
;    specerr[*]=1.

    fit=POLY_FIT(lamcut[subs],spec[subs],1,chisq=chisq,measure_errors=specerr[subs],sigma=sigma,yfit=yfit)
    contimage[i,j]=fit[1]*(maxlambda-minlambda)
    conterrimage[i,j]=sigma*(maxlambda-minlambda)
;    plot,lamcut,spec,yrange=[0.8,1.2],title=STRTRIM(chisq/N_ELEMENTS(subs),2)
;    oplot,lamcut[subs],yfit,thick=4
;    WAIT,0.02
ENDFOR
ENDFOR
WRITEFITS,'linemaps/contimage_may18.fits',contimage
WRITEFITS,'linemaps/conterrimage_may18.fits',conterrimage
STOP

coim=READFITS('linemaps/allpix_try1_co.fits')
myplot,file='linemaps/cont_co.ps'
ind=WHERE(rarr LT 25.)
plot,coim[ind],contimage[ind],psym=4,xrange=[0.15,0.35],yrange=[-0.05,0.25],/ysty,/xsty,xtitle='CO Linestrength',ytitle='Cont. Slope (2.28!9m!xm-2.05!9m!xm)',chars=2
ind=WHERE(rarr LT 4.)
loadct,13,/si
oplot,coim[ind],contimage[ind],psym=1,symsize=2,color=250,thick=5
legend,['All Data','R < 0.2"'],symsize=[1,2],color=[0,250],thick=[2,5],psym=[4,1],/right,chars=1.5
loadct,13,/si
device,/close
set_plot,'x'
STOP

END

