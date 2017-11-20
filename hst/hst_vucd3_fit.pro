;######################################################################
;
; Copyright (C) 1999-2005, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; For details on the method see:
;   Cappellari M., 2002, MNRAS, 333, 400
;
; Updated versions of the software are available from my web page
; http://purl.org/cappellari/software
;
; If you have found this software useful for your
; research, I would appreciate an acknowledgment to
; `use of the MGE fitting software developed by Cappellari (2002)'.
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;######################################################################
;+
; NAME:
;       TEST_MGE_FIT
;
; AUTHOR:
;       Michele Cappellari, Astrophysics Sub-department, University of Oxford, UK
;
; PURPOSE:
;       Exercise all the routines in the MGE_FIT_SECTORS package.
;       This procedure is intended to be used as a template to be
;       customized for each particular MGE fitting problem.
;
; EXPLANATION:
;       Further information is available in
;       Cappellari M., 2002, MNRAS, 333, 400
;
; CALLING SEQUENCE:
;       TEST_MGE_FIT
;
; EXAMPLE:
;       The following command will take various minutes to complete,
;       while plotting intermediate results on the screen
;
;           test_mge_fit
;
; PROCEDURES USED:
;       The following procedures are contained in the TEST_MGE_FIT program.
;           FIT_M32     -- Perform an MGE fit of the galaxy M32
;           FIT_NGC4342 -- Perform an MGE fit of the galaxy NGC 4342
;           FIT_NGC4473 -- Perform an MGE fit of the galaxy NGC 4473
;           FIT_DOUBLE_POWERLAW_1D -- Fit a 1D MGE model
;           MGE_FIT_1D_HERNQUIST_MODEL -- Compute the circular velocity of an MGE 1D fit
;           FIT_NGC5831_TWIST -- Perform a twisted MGE fit of the galaxy NGC 5831
;
;       Other routines needed from the MGE_FIT_SECTORS package by
;       Michele Cappellari (http://purl.org/cappellari/idl)
;           SECTORS_PHOTOMETRY       -- perform photometry along sectors
;           MGE_FIT_SECTORS          -- do the actual MGE fit
;           MGE_PRINT_CONTOURS       -- plot the results
;           SECTORS_PHOTOMETRY_TWIST -- perform photometry along sectors with point symmetry
;           MGE_FIT_SECTORS_TWIST    -- do the actual MGE fit with possible isophote twist
;           MGE_PRINT_CONTOURS_TWIST -- plot the results with possible isophote twist
;           WFPC2_MGE_FIT            -- do an MGE fit of a WFPC2 image (PC1 + MOSAIC)
;           FIND_GALAXY              -- find galaxy center and position angle
;           MGE_FIT_1D               -- perform a 1D MGE fit
;           CAP_RANGE                -- returns a sequence of values
;
;       Other IDL routines needed:
;           BVLS  -- Michele Cappellari porting of Lawson & Hanson generalized NNLS
;                    http://purl.org/cappellari/idl
;           MPFIT -- Craig Markwardt porting of Levenberg-Marquardt MINPACK-1
;                    http://purl.com/net/mpfit
;           JAM package -- The Jeans Anisotropic MGE (JAM) package is required only
;                    if one wants to run the optional example MGE_FIT_1D_HERNQUIST_MODEL
;                    http://purl.org/cappellari/idl
;
;       Astronomy User's Library routines needed (http://idlastro.gsfc.nasa.gov):
;           FITS_READ
;           DIST_CIRCLE
;
; MODIFICATION HISTORY:
;       V1.0: Michele Cappellari, Leiden, January 2000
;       V2.0: Updated all examples, MC, Leiden July 2001
;       V2.1: Added 1D MGE fit example, MC, Leiden, 10 May 2003
;       V2.11: Removed the un-necessary parameter EPS from the routine
;           SECTORS_PHOTOMETRY_TWIST. MC, Leiden, 1 September 2004
;       V2.12: Replaced LOGRANGE keyword with the new MAGRANGE.
;           MC, Leiden, 1 May 2005
;       V2.2: Included a new example MGE_FIT_1D_HERNQUIST_MODEL.
;           MC, Oxford, 30 November 2008
;       V2.21: Uses CAP_RANGE. MC, Paranal, 8 November 2013
;-
;----------------------------------------------------------------------------
PRO fit_vucd3
;
  fits_read, '../data/HST_10137_03_ACS_HRC_F814W_drz.fits', img, h
  fits_read, 'sky_mask.fits',sky
  fits_read,'../data/HST_10137_03_ACS_HRC_F814W_drz.fits',wht,exten_no=2
  fits_read, '../data/HST_10137_03_ACS_HRC_F814W_drz.fits',mask,exten_no=3
;estimate sky level
mdrizz=[7.46042280032,7.9962758789,6.9261413824]
expt=350.0
mdrizzcount=mdrizz/expt
avgdrizz=mean(mdrizzcount)
img=img+avgdrizz
totcounts=img*expt
;stop
;totcounts=sqrt(totcounts)
imgsize=size(img,/dim)
for i=0,imgsize[0]-1 do begin
   for j=0,imgsize[1]-1 do begin
      if (totcounts[i,j] lt 0.) then totcounts[i,j]=mean(totcounts)
   endfor
endfor
err=sqrt(totcounts)


writefits,'vucd3_weightmap.fits',err
writefits,'vucd3_mask.fits',mask

img = (img - sky);*expt          ; subtract sky

writefits,'vucd3_skysubtract.fits',img
temp=img[431:831,462:862]
ngauss=20
minlevel=0.

find_galaxy,temp,majoraxis,eps,ang,xc,yc,FRACTION=.8

sectors_photometry,temp,eps,ang,xc,yc,radius,angle,counts,minlevel=minlevel

readcol,'tinytim_fits.dat',normpsf,sigmapsf,format='F,F'
MGE_fit_sectors,radius,angle,counts,eps,sol=sol,ngauss=ngauss,scale=scale,normpsf=normpsf,sigmapsf=sigmapsf
;stop


scale = 0.025 ; arcsec/pixel
extinct=0.034
zp=25.28697
msun=4.53
peak=sol[0,*]/(2*!PI*sol[1,*]^2*sol[2,*])
mu=zp+5*alog10(scale)-2.5*alog10(peak)-extinct
const=(64800/!PI)^2
intensity=const*(10^(0.4*(msun-mu)))
sigmaarc=sol[1,*]*scale
forprint,intensity,sigmaarc,sol[2,*],format='F,F,F',textout='vucd3_mge_output.dat'
;stop
;set_plot,'ps'
!P.Multi=[0,1,2]

fits_read,'./evstigneeva/galfit_model_fixed.fits',img,exten_no=1
fits_read,'./evstigneeva/galfit_model_fixed.fits',modelimg,exten_no=2
find_galaxy,img,m,e,a,xc,yc
xc=96 & yc=98
find_galaxy,modelimg,majoraxis,eps,ang,xc_mod,yc_mod
xc_mod=96 & yc_mod=98

fits_read,'./make_decon/serfixmodeli.fits',serfixi
find_galaxy,serfixi,m,e,a,xc_ser,yc_ser
radius=findgen(140)*0.4;(10^(findgen(85)*0.05))*scale
minrad=[0,findgen(139)*0.4];[0,(10^(findgen(84)*0.1))*scale]

;radius=radius/scale
;minrad=minrad/scale
photmag2=fltarr(n_elements(radius))
modelmag2=fltarr(n_elements(radius))
sersicmag=fltarr(n_elements(radius))
area2=fltarr(n_elements(radius))
for i=0,n_elements(radius)-1 do begin
   aper,img,xc,yc,maxflux_phot,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
   aper,img,xc,yc,minflux_phot,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
   aper,modelimg,xc_mod,yc_mod,maxflux_model,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
   aper,modelimg,xc_mod,yc_mod,minflux_model,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
   aper,serfixi,xc_ser,yc_ser,maxflux_ser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
   aper,serfixi,xc_ser,yc_ser,minflux_ser,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
   if (minrad[i] lt 0.025) then begin
      area2[i]=((!PI*(radius[i])^2)) ;*scale
      
      photflux=maxflux_phot
      modelflux=maxflux_model
      sersicflux=maxflux_ser
    endif else begin
      area2[i]=((!PI*(radius[i])^2)-(!PI*(minrad[i])^2)) ;*scale
      photflux=maxflux_phot-minflux_phot
      modelflux=maxflux_model-minflux_model
      sersicflux=maxflux_ser-minflux_ser
   endelse
   
    photmag2[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(photflux)-0.034
    modelmag2[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(modelflux)-0.034
    sersicmag[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(sersicflux)-extinct
 endfor


;set_plot,'x'

;stop
readcol,'./evstigneeva/vucd3_mge_outputsersic.dat',sersiclum,sersicsig,sersicq,sersicpa,format='D,D,D,D'
Msun=4.53d
fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/deconvolved_f814.fits',oldimg,head
hrotate,oldimg,head,img,newhead,1
find_galaxy,img,majoraxis,eps,ang,xci,yci


a=where(sersicq lt 0.8)
inintensity=sersiclum[a]
insigma=sersicsig[a]
inq=sersicq[a]
inpa=sersicpa[a]
mge2image,img,xci,yci,inintensity,insigma,inq,inpa,inmodel,zeropoint=zp,scale=scale,msun=msun

b=where(sersicq gt 0.8)
outintensity=sersiclum[b]
outsigma=sersicsig[b]
outq=sersicq[b]
outpa=sersicpa[b]
mge2image,img,xci,yci,outintensity,outsigma,outq,outpa,outmodel,zeropoint=zp,scale=scale,msun=msun

readcol,'vucd3_mge_output.dat',mgelum,mgesig,mgeq,format='F,F,F'
pa=fltarr(n_elements(mgeq))
pa[*]=0.
mge2image,img,xci,yci,mgelum,mgesig,mgeq,pa,mgemodel,zeropoint=zp,scale=scale,msun=msun
writefits,'./make_decon/mgedirectfit.fits',mgemodel
mgemag=fltarr(n_elements(radius))
sersicmag1=fltarr(n_elements(radius))
sersicmag2=fltarr(n_elements(radius))
for i=0,n_elements(radius)-1 do begin
   aper,inmodel,xci,yci,maxinflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
   aper,inmodel,xci,yci,mininflux,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
   aper,outmodel,xci,yci,maxoutflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
   aper,outmodel,xci,yci,minoutflux,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
   aper,mgemodel,xci,yci,maxmgeflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
   aper,mgemodel,xci,yci,minmgeflux,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
   if (minrad[i] lt 0.025) then begin
      area2[i]=((!PI*(radius[i])^2)) ;*scale
      sersicflux1=maxinflux
      sersicflux2=maxoutflux
      mgeflux=maxmgeflux
    endif else begin
      area2[i]=((!PI*(radius[i])^2)-(!PI*(minrad[i])^2)) ;*scale
      sersicflux1=maxinflux-mininflux
      sersicflux2=maxoutflux-minoutflux
      mgeflux=maxmgeflux-minmgeflux
   endelse
    sersicmag1[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(sersicflux1)
    sersicmag2[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(sersicflux2)
    mgemag[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(mgeflux)
 endfor

set_plot,'ps'
!P.Multi=[0,1,2]
radius=radius*scale
device,filename='surfbright_hstindivsersic_center.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7
;myplot,filename='surfbright_hstindivsersic_center.ps'
djs_plot,radius,photmag2,psym=2,xtitle='Radius ["]',ytitle='\mu_{F814W} [Mag/asec^2]',yran=[20.999,13],xran=[0,1.],charsize=1.5,charthick=4,xthick=3,ythick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],/ysty,ycharsize=0.95
djs_oplot,radius,sersicmag1,color='green',thick=3,linestyle=2 ;,psym=2
djs_oplot,radius,sersicmag2,color='blue',thick=3,linestyle=2
djs_oplot,radius,sersicmag,color='red',thick=3
djs_oplot,radius,modelmag2,color='cyan',thick=4
;djs_oplot,radius,mgemag,color='cyan',thick=4
items=['n=3.51','n=1.28']
lines=[2,2]
color=['green','blue']
al_legend,items,linestyle=lines,colors=color,/window,background_color='white',charthick=4,thick=3,/top,/right
xyouts,[.02],[20.5],['VUCD3'],charthick=3,charsize=1.5,/data
djs_plot,radius,photmag2-modelmag2,psym=2,ymargin=[4,0],position=[0.1,0.1,0.95,0.25],xtitle='Radius ["]',ytitle='\Delta \mu',charthick=4,xthick=3,ythick=2,charsize=1.5,ycharsize=0.65,xran=[0,1.],yran=[-0.055,0.055],/ysty
djs_oplot,radius,fltarr(n_elements(radius)),linestyle=2,thick=3
device,/close
!P.Multi=[0,1,1]
set_plot,'x'
stop

END
;----------------------------------------------------------------------------
;----------------------------------------------------------------------------
PRO hst_vucd3_fit
;
; This is the main routine to call in succession the MGE fits to
; M32, NGC4342, NGC4473, power-law and NGC5831, and measure the execution time.
; A run of this program takes: 
; - 691s on a Pentium III, 1.0GHz PC, with IDL 5.4.
; - 270s on a Pentium M, 1.6GHz PC, with IDL 6.1.
; - 150s on a Core2 Duo, 2.2GHz PC, with IDL 7.0.
; It was tested with IDL 5.6-8.1 under both Windows and Linux.
;
t = SYSTIME(1)
fit_vucd3
PRINT, 'Total computation time:', SYSTIME(1) - t, ' Seconds'

END
;----------------------------------------------------------------------------
