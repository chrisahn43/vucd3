@fitfunc_psf.pro
@make_resim.pro
PRO FIT_PSF2

COMMON PSFSHARE,scaling,hstsize,hsthead,psfsize

fitrad=12.
scaling=10
toppeak=10.

;prepare NIFS image
nifsfile='../vucd3_best8_cont_astcor.fits'
nifs=READFITS(nifsfile,nifshead)
SXADDPAR,nifshead,'EQUINOX','J2000'
nifssize=SIZE(nifs)
junk=MAX(nifs,maxpos)
xinit=maxpos MOD nifssize[1]
yinit=maxpos/nifssize[1]
GCNTRD,nifs,xinit,yinit,xcen,ycen,4.5
minx=xcen-fitrad & maxx=xcen+fitrad
miny=ycen-fitrad & maxy=ycen+fitrad
;minx=xcen & maxx=xcen+fitrad
;miny=ycen-fitrad/1.5+1 & maxy=ycen+fitrad/1.5
HEXTRACT,nifs,nifshead,nifs_fitim,nifs_fithead,minx,maxx,miny,maxy
WRITEFITS,'nifs_subim.fits',nifs_fitim,nifs_fithead
weight=nifs_fitim

;weight[*]=1

imsize=SIZE(nifs_fitim)
junk=MAX(nifs_fitim,maxpos)
xcenout=maxpos MOD imsize[1]
ycenout=maxpos/imsize[1]
makex,nifs_fitim,xarr,yarr,/ZERO
rarr=SQRT((xarr-xcenout)^2+(yarr-ycenout)^2)
ind=WHERE(rarr GT fitrad)
weight[ind]=0

;prepare HST image
hstfile='../HST_10137_03_ACS_HRC_F814W_drz.fits'
hst=READFITS(hstfile,hsthead,ext=1)
SXADDPAR,hsthead,'EQUINOX','J2000'
xhstinit=627. & yhstinit=659.
GCNTRD,hst,xhstinit,yhstinit,xhstcen,yhstcen,3.5
minx=xhstcen-fitrad*2. & maxx=xhstcen+fitrad*2.
miny=yhstcen-fitrad*2. & maxy=yhstcen+fitrad*2.
HASTROM,hst,hsthead,hst_noconvim,junk,nifs_fithead
;HEXTRACT,hst,hsthead,hst_fitim,hst_fithead,minx,maxx,miny,maxy
HASTROM,hst,hsthead,hst_fitim,hst_fithead,nifs_fithead
hst=hst_fitim
hsthead=hst_fithead

loadct,5,/silent
ind=WHERE(rarr GT 10. AND rarr LT 20.)
fluxscale=MEDIAN(nifs_fitim[ind])/MEDIAN(hst_noconvim[ind])

hstsize=SIZE(hst)
bighstim=REBIN(hst,hstsize[1]*scaling,hstsize[2]*scaling)

;prepare PSF iamge
psfsize=FIX(fitrad*2.*scaling+1.5)
;STOP
;single component PSF
nparam=4
initguess=[3.0,fluxscale,0.,0.0]
parinfo = replicate({value:0.D, fixed:0,limited:[0.D,0.D],limits:[0.D,0.D],mpmaxstep:0.D},nparam)
parinfo[*].value=initguess
parinfo[*].limited[*]=1
parinfo[0].limits=[0.5,10.0]
parinfo[1].limits=[fluxscale*0.1,fluxscale*10.0]
parinfo[2].limits=[-1.,1.]
parinfo[3].limits=[-1.,1.]
;singlepsf needs bighstim (x), nifs_fithead(y), nifs_fitim, psfsize,
;                hsthead, hstsize, scaling
fit=mpfit2dfun('singlepsf', bighstim, nifs_fithead, nifs_fitim,weight=weight, perror=err, xtol=1.D-4, maxiter=30, status=status, parinfo=parinfo, bestnorm=bestnorm)
fitim=SINGLEPSF(bighstim,nifs_fithead,fit)
resid=nifs_fitim-fitim
dof=FLOAT(N_ELEMENTS(WHERE(weight GT 0))-nparam)
redchisq=TOTAL(weight*resid^2)/dof
make_resim,'psf_residual2.fits',nifs_fitim,fitim,weight
WRITEFITS,'model2.fits',fitim,nifs_fithead
print,redchisq

rarr=rarr*0.05
myplot,file='fit_psf2.ps',ysize=9,yoff=1,/inches
!P.MULTI=[0,1,2]
plot,rarr,nifs_fitim,psym=1,symsize=0.5,/ylog,yrange=[1,2000],/ysty,xrange=[-0.05,1.8],/xsty,xtitle='Radius ["]',ytitle='Pixel Value',title='Single Gaussian PSF Model',xcharsize=0.00001,charsize=1.8,ymargin=[1,3]
oplot,rarr,hst_noconvim*fit[1],psym=4,symsize=0.5,color=50
oplot,rarr,fitim,psym=5,symsize=0.5,color=100
legend,['NIFS','Unconv HST','Conv HST'],psym=[1,4,5],color=[0,50,100],charsize=1.8,/right
plot,rarr,resid/nifs_fitim,psym=5,symsize=0.5,/nodata,/ysty,xrange=[-0.05,1.8],/xsty,xtitle='Radius ["]',ytitle='Resid',charsize=1.8,ymargin=[5,-1]
oplot,rarr,resid/nifs_fitim,psym=5,symsize=0.5,color=100
plots,!X.CRANGE,[0.0,0.0],linestyle=2,thick=3


plot,rarr,nifs_fitim,psym=1,symsize=0.5,/ylog,yrange=[1,2000],/ysty,xrange=[-0.05,0.4],/xsty,xtitle='Radius ["]',ytitle='Pixel Value',title='Single Gaussian PSF Model',xcharsize=0.00001,charsize=1.8,ymargin=[1,3]
oplot,rarr,hst_noconvim*fit[1],psym=4,symsize=1.0,color=50
oplot,rarr,fitim,psym=5,symsize=1.0,color=100
legend,['NIFS','Unconv HST','Conv HST'],psym=[1,4,5],color=[0,50,100],charsize=1.8,/right
plot,rarr,resid/nifs_fitim,psym=5,symsize=1.0,/nodata,/ysty,xrange=[-0.05,0.4],/xsty,xtitle='Radius ["]',ytitle='Resid',charsize=1.8,ymargin=[5,-1]
oplot,rarr,resid/nifs_fitim,psym=5,symsize=1.0,color=100
plots,!X.CRANGE,[0.0,0.0],linestyle=2,thick=3


nparam=6
initguess=[fit,10.0,0.1]
fitim=DOUBLEPSF(bighstim,nifs_fithead,initguess)
parinfo = replicate({value:0.D, fixed:0,limited:[0.D,0.D],limits:[0.D,0.D],mpmaxstep:0.D},nparam)
parinfo[*].value=initguess
parinfo[*].limited[*]=1
parinfo[0].limits=[0.5,10.0]
parinfo[1].limits=[fluxscale*0.1,fluxscale*10.0]
parinfo[2].limits=[-1.,1.]
parinfo[3].limits=[-1.,1.]
parinfo[4].limits=[2.0,2.*fitrad]
parinfo[5].limits=[0.00,1.0]
;singlepsf needs bighstim (x), nifs_fithead(y), nifs_fitim, psfsize,
;                hsthead, hstsize, scaling
fit=mpfit2dfun('doublepsf', bighstim, nifs_fithead, nifs_fitim,weight=weight, perror=err, xtol=1.D-4, maxiter=30, status=status, parinfo=parinfo, bestnorm=bestnorm)
fitim=DOUBLEPSF(bighstim,nifs_fithead,fit)
resid=nifs_fitim-fitim
dof=FLOAT(N_ELEMENTS(WHERE(weight GT 0))-nparam)
redchisq=TOTAL(weight*resid^2)/dof
make_resim,'psf_residual_double2.fits',nifs_fitim,fitim,weight
WRITEFITS,'model_double2.fits',fitim,nifs_fithead
print,redchisq

!P.MULTI=[0,1,2]
plot,rarr,nifs_fitim,psym=1,symsize=0.5,/ylog,yrange=[1,2000],/ysty,xrange=[-0.05,1.8],/xsty,xtitle='Radius ["]',ytitle='Pixel Value',title='Double Gaussian PSF Model',xcharsize=0.00001,charsize=1.8,ymargin=[1,3]
oplot,rarr,hst_noconvim*fit[1],psym=4,symsize=0.5,color=50
oplot,rarr,fitim,psym=5,symsize=0.5,color=100
legend,['NIFS','Unconv HST','Conv HST'],psym=[1,4,5],color=[0,50,100],charsize=1.8,/right
plot,rarr,resid/nifs_fitim,psym=5,symsize=0.5,/nodata,/ysty,xrange=[-0.05,1.8],/xsty,xtitle='Radius ["]',ytitle='Resid',charsize=1.8,ymargin=[5,-1]
oplot,rarr,resid/nifs_fitim,psym=5,symsize=0.5,color=100
plots,!X.CRANGE,[0.0,0.0],linestyle=2,thick=3


plot,rarr,nifs_fitim,psym=1,symsize=1.0,/ylog,yrange=[1,2000],/ysty,xrange=[-0.05,0.4],/xsty,xtitle='Radius ["]',ytitle='Pixel Value',title='Double Gaussian PSF Model',xcharsize=0.00001,charsize=1.8,ymargin=[1,3]
oplot,rarr,hst_noconvim*fit[1],psym=4,symsize=1.0,color=50
oplot,rarr,fitim,psym=5,symsize=1.0,color=100
legend,['NIFS','Unconv HST','Conv HST'],psym=[1,4,5],color=[0,50,100],charsize=1.8,/right
plot,rarr,resid/nifs_fitim,psym=5,symsize=1.0,/nodata,/ysty,xrange=[-0.05,0.4],/xsty,xtitle='Radius ["]',ytitle='Resid',charsize=1.8,ymargin=[5,-1]
oplot,rarr,resid/nifs_fitim,psym=5,symsize=1.0,color=100
plots,!X.CRANGE,[0.0,0.0],linestyle=2,thick=3

ind=WHERE(rarr LT 1.0)
print,STDDEV((resid/nifs_fitim)[ind])

device,/close
set_plot,'x'
STOP


END
