@fitfunc_psf.pro
@make_resim.pro
PRO FIT_PSF

COMMON PSFSHARE,scaling,hstsize,hsthead,psfsize

fitrad=20.;23
scaling=10
toppeak=10.

;prepare NIFS image
nifsfile='../m60-ucd1_cont_astcor.fits'
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
;miny=ycen-fitrad+1 & maxy=ycen+fitrad
HEXTRACT,nifs,nifshead,nifs_fitim,nifs_fithead,minx,maxx,miny,maxy
WRITEFITS,'nifs_subim.fits',nifs_fitim,nifs_fithead
weight=nifs_fitim
weight[*]=nifs_fitim

imsize=SIZE(nifs_fitim)
junk=MAX(nifs_fitim,maxpos)
xcenout=maxpos MOD imsize[1]
ycenout=maxpos/imsize[1]
makex,nifs_fitim,xarr,yarr,/ZERO
rarr=SQRT((xarr-xcenout)^2+(yarr-ycenout)^2)
rarr=rarr*0.05

;prepare HST image
hstfile='../g_lucy.fits'
hst=READFITS(hstfile,hsthead)
SXADDPAR,hsthead,'EQUINOX','J2000'
xhstinit=325. & yhstinit=329.

GCNTRD,hst,xhstinit,yhstinit,xhstcen,yhstcen,3.5
minx=xhstcen-fitrad*2. & maxx=xhstcen+fitrad*2.
miny=yhstcen-fitrad*2. & maxy=yhstcen+fitrad*2.
HASTROM,hst,hsthead,hst_noconvim,junk,nifs_fithead
HEXTRACT,hst,hsthead,hst_fitim,hst_fithead,minx,maxx,miny,maxy
hst=hst_fitim
hsthead=hst_fithead



loadct,5,/silent
;ind=WHERE(rarr GT 10. AND rarr LT 20.)
fluxscale=MEDIAN(nifs_fitim)/MEDIAN(hst_noconvim)

hstsize=SIZE(hst)
bighstim=REBIN(hst,hstsize[1]*scaling,hstsize[2]*scaling)

;prepare PSF iamge
psfsize=FIX(fitrad*2.*scaling+1.5)

;single component PSF
nparam=4
initguess=[3.0,fluxscale,0.5,0.5]
parinfo = replicate({value:0.D, fixed:0,limited:[0.D,0.D],limits:[0.D,0.D],mpmaxstep:0.D,relstep:0.0},nparam)
parinfo[*].value=initguess
parinfo[*].limited[*]=1
parinfo[0].limits=[0.5,10.0]
parinfo[1].limits=[fluxscale*0.1,fluxscale*10.0]
parinfo[2].limits=[-2.,2.]
parinfo[3].limits=[-2.,2.]
parinfo[*].relstep=0.2

;singlepsf needs bighstim (x), nifs_fithead(y), nifs_fitim, psfsize,
;                hsthead, hstsize, scaling
fit=mpfit2dfun('singlepsf', bighstim, nifs_fithead, nifs_fitim,weight=weight, perror=err, xtol=1.D-10, maxiter=30, status=status, parinfo=parinfo, bestnorm=bestnorm)
fitim=SINGLEPSF(bighstim,nifs_fithead,fit)
resid=nifs_fitim-fitim
dof=FLOAT(N_ELEMENTS(WHERE(weight GT 0))-nparam)
redchisq=TOTAL(weight*resid^2)/dof
make_resim,'psf_residual.fits',nifs_fitim,fitim,weight
WRITEFITS,'model.fits',fitim,nifs_fithead
print,redchisq


;make weight image masking out both HST and NIFS sources
jointim=nifs_fitim+fitim/2.
plot,rarr,jointim,psym=1,symsize=0.1,/ylog,yrange=[1,2000]
;used plot to find a fit-by-eye to the upper envelope of the nucleus
;this first (commented one) is not a very close shave, 2nd is more so
;testfunc=1.2*MAX(jointim)*EXP(-rarr/0.2)+0.1*MAX(jointim)*EXP(-rarr/0.7)
testfunc=1.3*MAX(jointim)*EXP(-rarr/0.18)+0.1*MAX(jointim)*EXP(-rarr/0.7)
oplot,rarr,testfunc,psym=4
sourceind=WHERE(jointim GT testfunc)
;weight[sourceind]=0

myplot,file='fit_psf.ps',ysize=9,yoff=1,/inches
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




;Gauss + Moffat
nparam=7
;take value from HIP116449 fit and subtract of HST psf
;initguess=[2.09,fit[1],fit[2],fit[3],20.06,0.922,4.765]
initguess=[2.0,fit[1],fit[2],fit[3],18.,1.1,4.765]
fitim=GAUSSMOFFAT(bighstim,nifs_fithead,initguess)
parinfo = replicate({value:0.D, fixed:0,limited:[0.D,0.D],limits:[0.D,0.D],mpmaxstep:0.D,step:0.0,relstep:0.0},nparam)
parinfo[*].value=initguess
parinfo[*].limited[*]=1
parinfo[0].limits=[0.5,10.0]
parinfo[1].limits=[fluxscale*0.1,fluxscale*10.0]
parinfo[2].limits=[-1.,2.]
parinfo[3].limits=[-1.,1.]
parinfo[4].limits=[2.0,5.*fitrad]
parinfo[5].limits=[0.00,5.]
parinfo[6].limits=[2.5,6]
;parinfo[0].fixed=1
parinfo[2].fixed=1
parinfo[3].fixed=1
;parinfo[4].fixed=1
;parinfo[5].fixed=1
parinfo[6].fixed=1

parinfo[*].relstep=0.5
;singlepsf needs bighstim (x), nifs_fithead(y), nifs_fitim, psfsize,
;                hsthead, hstsize, scaling
fit=mpfit2dfun('gaussmoffat', bighstim, nifs_fithead, nifs_fitim,weight=weight, perror=err, xtol=1.D-12, maxiter=30, status=status, parinfo=parinfo, bestnorm=bestnorm)
fitim=GAUSSMOFFAT(bighstim,nifs_fithead,fit,/WRITEPSF)
resid=nifs_fitim-fitim
dof=FLOAT(N_ELEMENTS(WHERE(weight GT 0))-nparam)
redchisq=TOTAL(weight*resid^2)/dof
make_resim,'psf_residual_moffat.fits',nifs_fitim,fitim,weight
WRITEFITS,'model_moffat.fits',fitim,nifs_fithead
print,redchisq

!P.MULTI=[0,1,2]
plot,rarr,nifs_fitim,psym=1,symsize=0.5,/ylog,yrange=[1,2000],/ysty,xrange=[-0.05,1.8],/xsty,xtitle='Radius ["]',ytitle='Pixel Value',title='Gauss + Moffat PSF Model',xcharsize=0.00001,charsize=1.8,ymargin=[1,3]
oplot,rarr,hst_noconvim*fit[1],psym=4,symsize=0.5,color=50
oplot,rarr,fitim,psym=5,symsize=0.5,color=100
legend,['NIFS','Unconv HST','Conv HST'],psym=[1,4,5],color=[0,50,100],charsize=1.8,/right
plot,rarr,resid/nifs_fitim,psym=5,symsize=0.5,/nodata,/ysty,xrange=[-0.05,1.8],/xsty,xtitle='Radius ["]',ytitle='Resid',charsize=1.8,ymargin=[5,-1]
oplot,rarr,resid/nifs_fitim,psym=5,symsize=0.5,color=100
plots,!X.CRANGE,[0.0,0.0],linestyle=2,thick=3


plot,rarr,nifs_fitim,psym=1,symsize=1.0,/ylog,yrange=[1,2000],/ysty,xrange=[-0.05,0.4],/xsty,xtitle='Radius ["]',ytitle='Pixel Value',title='Gauss + Moffat PSF Model',xcharsize=0.00001,charsize=1.8,ymargin=[1,3]
oplot,rarr,hst_noconvim*fit[1],psym=4,symsize=1.0,color=50
oplot,rarr,fitim,psym=5,symsize=1.0,color=100
legend,['NIFS','Unconv HST','Conv HST'],psym=[1,4,5],color=[0,50,100],charsize=1.8,/right
plot,rarr,resid/nifs_fitim,psym=5,symsize=1.0,/nodata,/ysty,xrange=[-0.05,0.4],/xsty,xtitle='Radius ["]',ytitle='Resid',charsize=1.8,ymargin=[5,-1]
oplot,rarr,resid/nifs_fitim,psym=5,symsize=1.0,color=100
plots,!X.CRANGE,[0.0,0.0],linestyle=2,thick=3

device,/close
set_plot,'x'




END
