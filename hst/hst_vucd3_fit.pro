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
;img = img/expt
;stop
scale = 0.025 ; arcsec/pixel
ngauss = 15
minlevel = 0.0000;1;01 ; counts/pixel
imgsize=size(img,/dim)
readcol,'tinytim_fits.dat',normPSF,sigmaPSF,q,format='F,F,F'

; Here we use FIND_GALAXY directly inside the procedure. Usually you may want
; to experiment with different values of the FRACTION keyword, before adopting
; given values of Eps, Ang, Xc, Yc.

find_galaxy, img, majorAxis, eps, ang, xc, yc,FRACTION=0.13

; Perform galaxy photometry

sectors_photometry, img, eps, ang, xc, yc, radius, angle, counts, MINLEVEL=minlevel;,N_sectors=19,SECTOR_WIDTH=5

; Do the actual MGE fit 
;set_plot,'ps'
;device,filename='vucd3_photfits.ps'
MGE_fit_sectors, radius, angle, counts, eps, $
    NGAUSS=ngauss, SIGMAPSF=sigmaPSF, NORMPSF=normPSF, SOL=sol, SCALE=scale;,/BULGE_DISK
;device,/close
;set_plot,'x'
;stop

peakbright=(sol[0,*])/(2*!PI*((sol[1,*])^2)*sol[2,*])
extinct=0.034
surfbrightI=25.28697+5.*alog10(scale)-2.5*alog10(peakbright)-extinct
;http://www.stsci.edu/hst/acs/analysis/zeropoints/old_page/localZeropoints
Intensity=((64800./!PI)^2)*10^(0.4*(4.10-surfbrightI))
print,intensity
sigmaarc=(sol[1,*])*scale
q=sol[2,*]
forprint,intensity,sigmaarc,q,format='F,F,F',textout='vucd3_mge_output.dat'
modelweight=peakbright
modelsig=sol[1,*]
rad=10^((findgen(36)*0.073)+0.073);[findgen(19)*3+3,findgen(30)*10+60];,250,370]
minrad=10^(findgen(36)*0.073)-1;[findgen(20)*3,findgen(29)*10+60];,200,250]
N=n_elements(rad)
photflux=fltarr(N)
modelflux=fltarr(N)
area=fltarr(N)
photmag=fltarr(N)
modelmag=fltarr(N)
;Make surface brightness with errors on sky subtraction
fits_read, '../data/HST_10137_03_ACS_HRC_F814W_drz.fits', newimg, h
fits_read, 'sky1sigup_mask.fits',sky1sigup
fits_read, 'sky1sigdown_mask.fits',sky1sigdown
img1sigup=(newimg+avgdrizz)-sky1sigup
img1sigdown=(newimg+avgdrizz)-sky1sigdown
photmag1sigup=fltarr(N)
photmag1sigdown=fltarr(N)
;stop
zp=25.28697
modelimg=fltarr(imgsize[0],imgsize[1])
x=[reverse(findgen(xc))+1,findgen(imgsize[0]-xc)]
y=[reverse(findgen(yc))+1,findgen(imgsize[1]-yc)]
for i=0,n_elements(x)-1 do begin
   for j=0,n_elements(y)-1 do begin
      temp=0.
      for k=0,n_elements(modelweight)-1 do begin
         temp+= (modelweight[k]*exp(-(1./(2*modelsig[k]^2))*(x[i]^2+(y[j]^2/q[k]^2))))

      endfor
      modelimg[i,j]=temp
    endfor
endfor
stop
;writefits,'hst_mge_vucd3.fits',modelimg

for i=0,n_elements(rad)-1 do begin
   maxflux_phot=djs_phot(xc,yc,rad[i],0.,img,skyval=skyval)
   minflux_phot=djs_phot(xc,yc,minrad[i],0.,img,skyval=skyval)
   maxflux_model=djs_phot(xc,yc,rad[i],0.,modelimg,skyval=skyval)
   minflux_model=djs_phot(xc,yc,minrad[i],0.,modelimg,skyval=skyval)
   maxflux_phot1sigup=djs_phot(xc,yc,rad[i],0.,img1sigup,skyval=skyval)
   minflux_phot1sigup=djs_phot(xc,yc,minrad[i],0.,img1sigup,skyval=skyval)
   maxflux_phot1sigdown=djs_phot(xc,yc,rad[i],0.,img1sigdown,skyval=skyval)
   minflux_phot1sigdown=djs_phot(xc,yc,minrad[i],0.,img1sigdown,skyval=skyval)
   if (minrad[i] lt 1.1) then begin
      area[i]=((!PI*(rad[i])^2));*scale

      photflux[i]=maxflux_phot
      modelflux[i]=maxflux_model
      photflux1sigup=maxflux_phot1sigup
      photflux1sigdown=maxflux_phot1sigdown
   endif else begin
      area[i]=((!PI*(rad[i])^2)-(!PI*(minrad[i])^2));*scale
      photflux[i]=maxflux_phot-minflux_phot
      modelflux[i]=maxflux_model-minflux_model
      photflux1sigup=maxflux_phot1sigup-minflux_phot1sigup
      photflux1sigdown=maxflux_phot1sigdown-minflux_phot1sigdown
   endelse
   
   photmag[i]=zp+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(photflux[i])-0.034
   modelmag[i]=zp+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(modelflux[i])-0.034
   photmag1sigup[i]=zp+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(photflux1sigup)-0.034
   photmag1sigdown[i]=zp+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(photflux1sigdown)-0.034
   
endfor
stop
;mue1=18.19+2.5*alog10(2*!PI*(5.42*scale)^2)
;mue1=mue1+1.3
;mue2=18.5+2.5*alog10(2*!PI*(32.36*scale)^2)
;mue2=mue2+0.5
;bn1=1.999*4.5-0.327
;bn2=1.999*1.03-0.327
;mu1=mue1+((2.5*bn1)/(alog(10)))*((((rad)/(5.42))^(1./4.5))-1)
;mu2=mue2+((2.5*bn2)/(alog(10)))*((((rad)/(32.36))^(1./1.03))-1)
;int1=(10^(-0.4*(mu1-zp)))
;int2=(10^(-0.4*(mu2-zp)))
;inttot=int1+int2
;mutot=zp-2.5*alog10(inttot)
set_plot,'ps'
;device,filename='surfbright_hstvorigmge.ps',/color
!P.Multi=[0,1,2]
;djs_plot,rad*scale,photmag,psym=2,xtitle='Radius ["]',ytitle='\mu [Mag/sqare arcsecond]',yran=[23,13],xran=[0,11],charsize=1.5,charthick=4,xthick=3,ythick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95]
;djs_oplot,rad*scale,modelmag,color='blue',thick=3 ;,psym=2
;djs_oplot,rad*scale,modelmag,psym=4
;djs_plot,rad*scale,photmag-modelmag,psym=2,ymargin=[4,0],position=[0.1,0.1,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',charthick=4,xthick=3,ythick=2,charsize=1.5,ycharsize=0.5
;device,/close

;device,filename='surfbright_hstvorigmge_center.ps',/color
;djs_plot,rad*scale,photmag,psym=2,xtitle='Radius ["]',ytitle='\mu [Mag/sqare arcsecond]',yran=[21,13],xran=[0,1.],charsize=1.5,charthick=4,xthick=3,ythick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],/ysty
;djs_oplot,rad*scale,modelmag,color='blue',thick=3 ;,psym=2
;djs_oplot,rad*scale,modelmag,psym=4
;djs_plot,rad*scale,photmag-modelmag,psym=2,ymargin=[4,0],position=[0.1,0.1,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',charthick=4,xthick=3,ythick=2,charsize=1.5,ycharsize=0.5,xran=[0,1.],yran=[-0.3,0.75],/ysty
;device,/close
;stop

;device,filename='surfbright_hstvgalmge_center.ps',/color
;djs_plot,rad*scale,photmag,psym=2,xtitle='Radius ["]',ytitle='\mu [Mag/sqare arcsecond]',yran=[23,10],xran=[0,1.],charsize=1.5,charthick=4,xthick=3,ythick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95]
;djs_oplot,rad*scale,mutot,color='green',thick=3 ;,psym=2
;djs_oplot,rad*scale,mutot,psym=4
;djs_plot,rad*scale,photmag-mutot,psym=2,ymargin=[4,0],position=[0.1,0.1,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',charthick=4,xthick=3,ythick=2,charsize=1.5,ycharsize=0.5,xran=[0,1.],yran=[-0.75,0.75],/ysty
;device,/close
;stop

readcol,'vucd3_mge_outputsersic.dat',sersiclum,sersicsig,sersicq,sersicpa,format='D,D,D,D'
Msun=4.10d
const=(64800./!PI)^2
sersicsb=Msun-2.5*alog10(sersiclum/const)
sersicpeak=10^(-0.4*(sersicsb-zp-5*alog10(scale)+extinct))
radius=(10^(findgen(80)*0.025))*scale
minrad=[0,(10^(findgen(79)*0.025))*scale]
sersicflux=fltarr(n_elements(radius))
for i=0,n_elements(radius)-1 do begin
   tmp=0.
   for j=0,n_elements(sersicpeak)-1 do begin
      tmp+= (sersicpeak[j]*exp(-(1./(2*sersicsig[j]^2))*(radius[i])^2))
      sersicflux[i]=tmp
   endfor
endfor
readcol,'vucd3_mge_output.dat',lum,sig,q,format='F,F,F'
sb=Msun-2.5*alog10(lum/const)
peak=10^(-0.4*(sb-zp-5*alog10(scale)+extinct))
flux=fltarr(n_elements(radius))
for i=0,n_elements(radius)-1 do begin
   temp=0.
   for j=0,n_elements(peak)-1 do begin
      temp+= (peak[j]*exp(-(1./(2*sig[j]^2))*(radius[i])^2))
      flux[i]=temp
   endfor
endfor
mag=zp+5*alog10(scale)-2.5*alog10(flux)-extinct

radius=radius/scale
minrad=minrad/scale
photmag2=fltarr(n_elements(radius))
modelmag2=fltarr(n_elements(radius))
area2=fltarr(n_elements(radius))
for i=0,n_elements(radius)-1 do begin
   maxflux_phot=djs_phot(xc,yc,radius[i],0.,img,skyval=skyval)
   minflux_phot=djs_phot(xc,yc,minrad[i],0.,img,skyval=skyval)
   maxflux_model=djs_phot(xc,yc,radius[i],0.,modelimg,skyval=skyval)
   minflux_model=djs_phot(xc,yc,minrad[i],0.,modelimg,skyval=skyval)

   if (minrad[i] lt 0.025) then begin
      area2[i]=((!PI*(radius[i])^2)) ;*scale
      
      photflux=maxflux_phot
      modelflux=maxflux_model
    endif else begin
      area2[i]=((!PI*(radius[i])^2)-(!PI*(minrad[i])^2)) ;*scale
      photflux=maxflux_phot-minflux_phot
      modelflux=maxflux_model-minflux_model
   endelse
   
    photmag2[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(photflux)-0.034
    modelmag2[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(modelflux)-0.034
endfor
sersicmag=zp+5*alog10(scale)-2.5*alog10(sersicflux)-extinct
;set_plot,'x'
device,filename='surfbright_hstvgalmge_center.ps',/color
djs_plot,radius*scale,photmag2,psym=2,xtitle='Radius ["]',ytitle='\mu [Mag/sqare arcsecond]',yran=[21,13],xran=[0,1.],charsize=1.5,charthick=4,xthick=3,ythick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],/ysty
djs_oplot,radius*scale,sersicmag,color='green',thick=3 ;,psym=2
djs_oplot,radius*scale,sersicmag,psym=4
djs_plot,radius*scale,photmag2-sersicmag,psym=2,ymargin=[4,0],position=[0.1,0.1,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',charthick=4,xthick=3,ythick=2,charsize=1.5,ycharsize=0.5,xran=[0,1.],yran=[-0.07,0.8],/ysty
device,/close
device,filename='surfbright_hstvorigmge_center.ps',/color
djs_plot,radius*scale,photmag2,psym=2,xtitle='Radius ["]',ytitle='\mu [Mag/sqare arcsecond]',yran=[21,13],xran=[0,1.],charsize=1.5,charthick=4,xthick=3,ythick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],/ysty
djs_oplot,radius*scale,modelmag2,color='blue',thick=3 ;,psym=2
djs_oplot,radius*scale,modelmag2,psym=4
djs_plot,radius*scale,photmag2-modelmag2,psym=2,ymargin=[4,0],position=[0.1,0.1,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',charthick=4,xthick=3,ythick=2,charsize=1.5,ycharsize=0.5,xran=[0,1.],yran=[-0.4,0.8],/ysty
device,/close

device,filename='residual_comps.ps',/color
!P.Multi=[0,1,1]
djs_plot,radius*scale,photmag2-sersicmag,psym=2,xran=[0,1],yran=[-0.3,0.8],/ysty,ytitle='Residuals',xtitle='Radius ["]',charsize=1.5,charthick=4,xthick=3,ythick=3
djs_oplot,radius*scale,photmag2-modelmag2,psym=2,color='blue'
device,/close
stop
readcol,'vucd3_mge_outputsersic.dat',sersiclum,sersicsig,sersicq,sersicpa,format='D,D,D,D'
Msun=4.10d
const=(64800./!PI)^2
a=where(sersicq lt 0.8)
sersicsb1=Msun-2.5*alog10(sersiclum[a]/const)
sersicpeak1=10^(-0.4*(sersicsb1-zp-5*alog10(scale)+extinct))
b=where(sersicq gt 0.8)
sersicsb2=Msun-2.5*alog10(sersiclum[b]/const)
sersicpeak2=10^(-0.4*(sersicsb2-zp-5*alog10(scale)+extinct))

radius=(10^(findgen(80)*0.025))*scale
minrad=[0,(10^(findgen(79)*0.025))*scale]
sersicflux1=fltarr(n_elements(radius))
sersicflux2=fltarr(n_elements(radius))

for i=0,n_elements(radius)-1 do begin
   tmp=0.
   temp=0.
   for j=0,n_elements(sersicpeak1)-1 do begin
      tmp+= (sersicpeak1[j]*exp(-(1./(2*sersicsig[a[j]]^2))*(radius[i])^2))
      sersicflux1[i]=tmp
   endfor
   for j=0,n_elements(sersicpeak2)-1 do begin
      temp+=(sersicpeak2[j]*exp(-(1./(2*sersicsig[b[j]]^2))*(radius[i])^2))
      sersicflux2[i]=temp
   endfor 
endfor
sersicmag1=zp+5*alog10(scale)-2.5*alog10(sersicflux1)-extinct
sersicmag2=zp+5*alog10(scale)-2.5*alog10(sersicflux2)-extinct
!P.Multi=[0,1,2]
device,filename='surfbright_hstindivsersic_center.ps',/color
djs_plot,radius,photmag2,psym=2,xtitle='Radius ["]',ytitle='\mu [Mag/sqare arcsecond]',yran=[21,13],xran=[0,1.],charsize=1.5,charthick=4,xthick=3,ythick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],/ysty
djs_oplot,radius,sersicmag1,color='green',thick=3,linestyle=2 ;,psym=2
djs_oplot,radius,sersicmag2,color='blue',thick=3,linestyle=2
djs_oplot,radius,sersicmag,color='red',thick=3
djs_plot,radius,photmag2-sersicmag,psym=2,ymargin=[4,0],position=[0.1,0.1,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',charthick=4,xthick=3,ythick=2,charsize=1.5,ycharsize=0.5,xran=[0,1.],yran=[-0.07,0.8],/ysty
device,/close
!P.Multi=[0,1,1]
set_plot,'x'
stop

;set_plot,'ps'
;device,filename='fluxhst_fits.ps',/color
;djs_plot,rad*scale,modelflux,psym=2,xtitle='Radius ["]',ytitle='Flux',charsize=1.5,charthick=4,xthick=3,ythick=3
;djs_oplot,rad*scale,photflux,psym=2,color='blue'
;djs_oplot,rad,modelflux2,psym=2,color='green'
;device,/close
;set_plot,'x'
;stop
;djs_plot,rad,modelflux,psym=2
;djs_oplot,rad,photflux,psym=2,color='blue'
percentdiff=((photflux-modelflux)/photflux)*100
print,percentdiff
;stop

; Print the data-model contours comparison of the whole image

MGE_print_contours, img>minlevel, ang, xc, yc, sol, $
    FILE='hst_vucd3.ps', SCALE=scale, MAGRANGE=9, $
    SIGMAPSF=sigmaPSF, NORMPSF=normPSF, BINNING=7

; Print the data-model contours comparison of the central regions

s = SIZE(img)
img = img[xc-s[1]/9:xc+s[1]/9,yc-s[2]/9:yc+s[2]/9]
MGE_print_contours, img, ang, s[1]/9, s[2]/9, sol, $
    FILE='hst_vucd3_nuclear.ps', SCALE=scale, MAGRANGE=9, $
                    SIGMAPSF=sigmaPSF, NORMPSF=normPSF

;for i=1,40 do begin
;   modelrad=(0.01*findgen((300*i))+0.01)
;   temp=0.
;   for j=0,n_elements(modelrad)-2 do begin
;      for k=0,n_elements(modelweight)-1 do begin
;         temp+= (modelweight[k]*exp(-(1./(2*modelsig[k]^2))*((modelrad[j])^2)))*(2*!PI*modelrad[j])*(modelrad[j+1]-modelrad[j]);multiply by q
;      endfor
;   endfor
;   maxflux=temp
;   if i gt 1 then begin
;      minmodelrad=(0.01*findgen(300*(i-1))+0.01)
;      temp2=0.
;      for j=0,n_elements(minmodelrad)-2 do begin
;         for k=0,n_elements(modelweight)-1 do begin
;            temp2+= (modelweight[k]*exp(-(1./(2*modelsig[k]^2))*((minmodelrad[j])^2)))*(2*!PI*minmodelrad[j])*(minmodelrad[j+1]-minmodelrad[j]) ;multiply by q
;         endfor
;      endfor
;      minflux=temp2
;      modelflux2[i-1]=maxflux-minflux
;   endif else begin
;      modelflux2[i-1]=maxflux
;   endelse
   
;endfor
;stop
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
