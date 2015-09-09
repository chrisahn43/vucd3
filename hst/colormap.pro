pro colormap

  fits_read,'../data/HST_10137_03_ACS_HRC_F606W_drz.fits',img_r,h_r
  fits_read,'sky_mask_f606.fits',sky_r
  mdrizzr=[7.49002591447,8.24492204734,7.3641947003]
  exptr=290.
  mdrizzcountr=mdrizzr/exptr
  avgdrizzr=mean(mdrizzcountr)
  img_r=img_r + avgdrizzr
  img_r=img_r - sky_r
  err_r=sqrt(img_r)

  fits_read,'../data/HST_10137_03_ACS_HRC_F814W_drz.fits',img_i,h_r
  fits_read, 'sky_mask.fits',sky_i
  mdrizz=[7.46042280032,7.9962758789,6.9261413824]
  expt=350.0
  mdrizzcount=mdrizz/expt
  avgdrizz=mean(mdrizzcount)
  img_i=img_i+avgdrizz
  img_i = img_i - sky_i               ; subtract sky
  err_i=sqrt(img_i)
  xc=627
  yc=660
;  exptime_r=870.
;  exptime_i=1050.
  scale=0.027
  zeropoint_r=25.9799
  zeropoint_i=25.28697
  rad=findgen(100)+1;10^((findgen(36)*0.073)+0.073)
  minrad=findgen(100);10^(findgen(36)*0.073)
  rflux=fltarr(n_elements(rad))
  iflux=fltarr(n_elements(rad))
  rmag=fltarr(n_elements(rad))
  imag=fltarr(n_elements(rad))
  rmagerr=fltarr(n_elements(rad))
  imagerr=fltarr(n_elements(rad))
  area=fltarr(n_elements(rad))
  for i=0,n_elements(rad)-1 do begin
     maxflux_r=djs_phot(xc,yc,rad[i],0.,img_r,skyval=skyval)
     minflux_r=djs_phot(xc,yc,minrad[i],0.,img_r,skyval=skyval)
     rflux[i]=maxflux_r-minflux_r
     maxfluxerr_r=djs_phot(xc,yc,rad[i],0.,err_r,skyval=skyval)
     minfluxerr_r=djs_phot(xc,yc,minrad[i],0.,err_r,skyval=skyval)
     rerrflux=maxfluxerr_r-minfluxerr_r
     maxflux_i=djs_phot(xc,yc,rad[i],0.,img_i,skyval=skyval)
     minflux_i=djs_phot(xc,yc,minrad[i],0.,img_i,skyval=skyval)  
     iflux[i]=maxflux_i-minflux_i
     maxfluxerr_i=djs_phot(xc,yc,rad[i],0.,err_i,skyval=skyval)
     minfluxerr_i=djs_phot(xc,yc,minrad[i],0.,err_i,skyval=skyval)
     ierrflux=maxfluxerr_i-minfluxerr_i
     area[i]=((!PI*(rad[i])^2)-(!PI*(minrad[i])^2))*scale

     rmag[i]=zeropoint_r + 5*alog10(scale) + 2.5*alog10(area[i]) - 2.5*alog10(rflux[i]);-0.034 ; + 2.5*alog10(exptime_r)
     rmagerr[i]=zeropoint_r + 5*alog10(scale) +((-2.5*alog10(rflux[i]+rerrflux))+(2.5*alog10(rflux[i]-rerrflux)))/2.
     imag[i]=zeropoint_i + 5*alog10(scale) + 2.5*alog10(area[i]) - 2.5*alog10(iflux[i]);-0.048 ; + 2.5*alog10(exptime_i)
     imagerr[i]=zeropoint_i + 5*alog10(scale) +((-2.5*alog10(iflux[i]+ierrflux))+(2.5*alog10(iflux[i]-ierrflux)))/2.

  endfor
  colormag=rmag-imag
  err=rmagerr-imagerr
  set_plot,'ps'
  device,filename='radialcolor.ps',/color
  djs_plot,rad*scale,colormag,psym=4,ytitle='(\mu_R - \mu_I) [mag/arcsec^2]',xtitle='Radius ["]',yran=[-0.2,0.8],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3
;  oploterr,rad*scale,colormag,err
  device,/close
  set_plot,'x'
  stop

END

