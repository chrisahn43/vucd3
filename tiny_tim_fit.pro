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

;----------------------------------------------------------------------------
PRO fit_tinytim;ngc4342
;
; This procedure reproduces Figures 8-9 in Cappellari (2002)
;
; This example illustrates a simple MGE fit to one single HST/WFPC2/F814W image.
;
  fits_read, 'data/result00_psf.fits', img, h
  img=img*(1./total(img))
;skylev = 0.55 ; counts/pixel
;img = img - skylev ; subtract sky
scale = 0.0231 ; arcsec/pixel

ngauss = 5
minlevel = 5.e-6 ; counts/pixel

; Here we use FIND_GALAXY directly inside the procedure. Usually you may want
; to experiment with different values of the FRACTION keyword, before adopting
; given values of Eps, Ang, Xc, Yc.

find_galaxy, img, majorAxis, eps, ang, xc, yc, FRACTION=.9, /PLOT

; Perform galaxy photometry

sectors_photometry, img, eps, ang, xc, yc, radius, angle, counts, MINLEVEL=minlevel;,N_SECTORS=52,SECTOR_WIDTH=2.5

MGE_fit_sectors, radius, angle, counts, eps,$
                 SOL=sol, NGAUSS=ngauss, QBOUNDS=[0.9999999,1.],scale=scale;, /LINEAR;,/NEGATIVE
stop
modelnorm=1./total(sol[0,*]) ;make output sum 1
modelpeak=modelnorm*(sol[0,*]) ;make output sum 1
modelsig=sol[1,*] ;sigma
modelweight=modelpeak/(2*!PI*(modelsig)^2) ;in flux need A/2*pi*sigma^2
forprint,modelpeak,modelsig,sol[2,*],format='F,F,F',textout='tinytim_fits.dat' ;output for mge fits on hst image

photflux=fltarr(40)
modelflux=fltarr(40)
aperflux=fltarr(40)
rad=fltarr(40)
for i=1,40 do begin ;aperture photometry
   rad[i-1]=2*i
   asdf=djs_phot(xc,yc,rad[i-1],0.,img, skyval=skyval)
;   aper,img,xc,yc,flux,aperfluxerr,sky,skyerr,1.0,rad[i-1],[20,25],[-32000000,320000000],setskyval=0.0,/flux,/exact,/silent,/nan
   photflux[i-1]=asdf
;   aperflux[i-1]=flux
endfor


for i=1,40 do begin ;reconstruct gaussians
   modelrad=(0.01*findgen((200*i)+1))
   temp=0.
   for j=0,n_elements(modelrad)-2 do begin
      for k=0,n_elements(modelweight)-1 do begin
         temp+= (modelweight[k]*exp(-(1./(2*modelsig[k]^2))*((modelrad[j])^2)))*(2*!PI*modelrad[j])*(modelrad[j+1]-modelrad[j])
      endfor
   endfor
   modelflux[i-1]=temp
endfor

;set_plot,'ps'
;device,filename='tinytim_photfits.ps',/color
djs_plot,rad,modelflux,psym=2,yran=[0.2,1.1]
djs_oplot,rad,photflux,psym=2,color='blue'
percentdiff=((photflux-modelflux)/photflux)*100
print,percentdiff
;device,/close
;set_plot,'x'
stop
; Print the data-model contours comparison of the whole image

MGE_print_contours, img>minlevel, ang, xc, yc, sol, $
    FILE='tinytim.ps', SCALE=scale, MAGRANGE=9;, $
    ;SIGMAPSF=sigmaPSF, NORMPSF=normPSF, BINNING=7

; Print the data-model contours comparison of the central regions

s = SIZE(img)
img = img[xc-s[1]/9:xc+s[1]/9,yc-s[2]/9:yc+s[2]/9]
MGE_print_contours, img, ang, s[1]/9, s[2]/9, sol, $
    FILE='tinytim_nuclear.ps', SCALE=scale, MAGRANGE=9;, $
 ;   SIGMAPSF=sigmaPSF, NORMPSF=normPSF

END



;----------------------------------------------------------------------------
PRO tiny_tim_fit
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
;fit_m32
fit_tinytim;ngc4342
;fit_ngc4473
;fit_double_powerlaw_1d
;fit_ngc5831_twist
PRINT, 'Total computation time:', SYSTIME(1) - t, ' Seconds'

END

;----------------------------------------------------------------------------
