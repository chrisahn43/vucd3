pro colormap

  fits_read,'../data/HST_10137_03_ACS_HRC_F606W_drz.fits',img_r,h_r
  fits_read,'sky_mask_f606.fits',sky_r
  mdrizzr=[7.49002591447,8.24492204734,7.3641947003]
  exptr=290.
  mdrizzcountr=mdrizzr/exptr
  avgdrizzr=mean(mdrizzcountr)
  img_r=img_r + avgdrizzr
  
  totcountsr=img_r*exptr
  imgsizer=size(img_r,/dim)
  for i=0,imgsizer[0]-1 do begin
     for j=0,imgsizer[1]-1 do begin
        if (totcountsr[i,j] lt 0.) then totcountsr[i,j]=mean(totcountsr)
     endfor
  endfor
  err_r=sqrt(totcountsr)
  img_r=(img_r - sky_r);*exptr
  writefits,'vucd3_weightmap_r.fits',err_r
  writefits,'vucd3_skysubtract_r.fits',img_r
  stop
  fits_read,'../data/HST_10137_03_ACS_HRC_F814W_drz.fits',img_i,h_r
  fits_read, 'sky_mask.fits',sky_i
  mdrizz=[7.46042280032,7.9962758789,6.9261413824]
  expt=350.0
  mdrizzcount=mdrizz/expt
  avgdrizz=mean(mdrizzcount)
  img_i=img_i+avgdrizz
  totcountsi=img_i*expt
  imgsizei=size(img_i,/dim)
  for i=0,imgsizei[0]-1 do begin
     for j=0,imgsizei[1]-1 do begin
        if (totcountsi[i,j] lt 0.) then totcountsi[i,j]=mean(totcountsi)
     endfor
  endfor
  err_i=sqrt(totcountsi)
  img_i=(img_i - sky_i);*expt
  writefits,'vucd3_weightmap.fits',err_i
  writefits,'vucd3_skysubtract.fits',img_i
  stop
  xc=627
  yc=660
;  exptime_r=870.
;  exptime_i=1050.
  scale=0.025
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

pro colormap_models
  fits_read,'galfit_model.fits',img_i,h_i,exten_no=2
  fits_read,'galfit_model_r.fits',img_r,h_r,exten_no=2
  fits_read,'vucd3_weightmap.fits',err_i
  fits_read,'vucd3_weightmap_r.fits',err_r
  xc_i=98 & yc_i=99
  xc_r=98 & yc_r=98
  scale=0.025
  zeropoint_r=25.9799
  zeropoint_i=25.28697
  rad=findgen(100)+1
  minrad=findgen(100)
  rflux=fltarr(n_elements(rad))
  iflux=fltarr(n_elements(rad))
  rmag=fltarr(n_elements(rad))
  imag=fltarr(n_elements(rad))
  rmagerr=fltarr(n_elements(rad))
  imagerr=fltarr(n_elements(rad))
  area=fltarr(n_elements(rad))
  for i=0,n_elements(rad)-1 do begin
     maxflux_r=djs_phot(xc_r,yc_r,rad[i],0.,img_r,skyval=skyval)
     minflux_r=djs_phot(xc_r,yc_r,minrad[i],0.,img_r,skyval=skyval)
     rflux[i]=maxflux_r-minflux_r
     maxfluxerr_r=djs_phot(xc_r,yc_r,rad[i],0.,err_r,skyval=skyval)
     minfluxerr_r=djs_phot(xc_r,yc_r,minrad[i],0.,err_r,skyval=skyval)

     maxflux_i=djs_phot(xc_i,yc_i,rad[i],0.,img_i,skyval=skyval)
     minflux_i=djs_phot(xc_i,yc_i,minrad[i],0.,img_i,skyval=skyval)  
     iflux[i]=maxflux_i-minflux_i
     maxfluxerr_i=djs_phot(xc_i,yc_i,rad[i],0.,err_i,skyval=skyval)
     minfluxerr_i=djs_phot(xc_i,yc_i,minrad[i],0.,err_i,skyval=skyval)

     area[i]=((!PI*(rad[i])^2)-(!PI*(minrad[i])^2))
     ierrflux=(maxfluxerr_i-minfluxerr_i)/area[i]
     rerrflux=(maxfluxerr_r-minfluxerr_r)/area[i]

     rmag[i]=zeropoint_r + 5*alog10(scale) + 2.5*alog10(area[i]) - 2.5*alog10(rflux[i])-0.034
     rmagerr[i]=5*alog10(scale) +(((-2.5*alog10(rflux[i]+rerrflux))+(2.5*alog10(rflux[i]-rerrflux)))/2.)
     rmagerr[i]=5*alog10(scale) -2.5*(rerrflux/rflux[i])
     
     imag[i]=zeropoint_i + 5*alog10(scale) + 2.5*alog10(area[i]) - 2.5*alog10(iflux[i])-0.048
     imagerr[i]=5*alog10(scale) + (((-2.5*alog10(iflux[i]+ierrflux))+(2.5*alog10(iflux[i]-ierrflux)))/2.)
     imagerr[i]=5*alog10(scale) -2.5*(ierrflux/iflux[i])


  endfor
  colormag=rmag-imag
  err=rmagerr-imagerr

  set_plot,'ps'
  device,filename='radialcolor.ps',/color
  djs_plot,rad*scale,colormag,psym=4,ytitle='(\mu_R - \mu_I) [mag/arcsec^2]',xtitle='Radius ["]',yran=[-0.2,1.6],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3
  oploterr,rad*scale,colormag,err
  device,/close
  set_plot,'x'
  stop

END
