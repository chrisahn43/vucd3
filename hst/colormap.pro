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
;  stop
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
;  stop
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
     maxflux_i=djs_phot(xc,yc,rad[i],0.,img_i,skyval=skyval)
     minflux_i=djs_phot(xc,yc,minrad[i],0.,img_i,skyval=skyval)  
     iflux[i]=maxflux_i-minflux_i
     maxfluxerr_i=djs_phot(xc,yc,rad[i],0.,err_i,skyval=skyval)
     minfluxerr_i=djs_phot(xc,yc,minrad[i],0.,err_i,skyval=skyval)
     area[i]=((!PI*(rad[i])^2)-(!PI*(minrad[i])^2));*scale
     rerrflux=(maxfluxerr_r-minfluxerr_r)/area[i]
     ierrflux=(maxfluxerr_i-minfluxerr_i)/area[i]

     rmag[i]=zeropoint_r + 5*alog10(scale) + 2.5*alog10(area[i]) - 2.5*alog10(rflux[i])-0.034 ; + 2.5*alog10(exptime_r)
     rmagerr[i]=5*alog10(scale) -2.5*(rerrflux/rflux[i])

     imag[i]=zeropoint_i + 5*alog10(scale) + 2.5*alog10(area[i]) - 2.5*alog10(iflux[i])-0.048 ; + 2.5*alog10(exptime_i)
     imagerr[i]=5*alog10(scale) -2.5*(ierrflux/iflux[i])


  endfor
  colormag=rmag-imag
  err=rmagerr-imagerr
  set_plot,'ps'
  device,filename='radialcolor_orig.ps',/color
  djs_plot,rad*scale,colormag,psym=4,ytitle='(\mu_R - \mu_I) [mag/arcsec^2]',xtitle='Radius ["]',yran=[-0.2,1.0],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3
  oploterr,rad*scale,colormag,err
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
  device,filename='radialcolor_models.ps',/color
  djs_plot,rad*scale,colormag,psym=4,ytitle='(\mu_R - \mu_I) [mag/arcsec^2]',xtitle='Radius ["]',yran=[-0.2,1.6],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3
  oploterr,rad*scale,colormag,err
  device,/close
;  set_plot,'x'
  stop

                                  ;unconvolved mge sb comparison
                                ;F814W
  mue1=18.19+2.5*alog10(2*!PI*(5.42*scale)^2)
  mue1=mue1+1.3
  mue2=18.5+2.5*alog10(2*!PI*(32.36*scale)^2)
  mue2=mue2+0.5
  bn1=1.999*4.5-0.327
  bn2=1.999*1.03-0.327
  mu1=mue1+((2.5*bn1)/(alog(10)))*((((rad)/(5.42))^(1./4.5))-1)
  mu2=mue2+((2.5*bn2)/(alog(10)))*((((rad)/(32.36))^(1./1.03))-1)
  int1=(10^(-0.4*(mu1-zeropoint_i)))
  int2=(10^(-0.4*(mu2-zeropoint_i)))
  inttot=int1+int2
  mutot=zeropoint_i-2.5*alog10(inttot)

  mue1606=18.82+2.5*alog10(2*!PI*(4.66*scale)^2)
  mue1606=mue1606+1.3
  mue2606=19.2+2.5*alog10(2*!PI*(27.71*scale)^2)
  mue2606=mue2606+0.5
  bn1606=1.999*4.66-0.327
  bn2606=1.999*0.98-0.327
  mu1606=mue1606+((2.5*bn1606)/(alog(10)))*((((rad)/(4.66))^(1./4.46))-1)
  mu2606=mue2606+((2.5*bn2606)/(alog(10)))*((((rad)/(27.71))^(1./0.98))-1)
  int1606=(10^(-0.4*(mu1606-zeropoint_r)))
  int2606=(10^(-0.4*(mu2606-zeropoint_r)))
  inttot606=int1606+int2606
  mutot606=zeropoint_r-2.5*alog10(inttot606)
  stop
  colormu=mutot606-mutot
  device,filename='radialcolor_sersic.ps',/color
  djs_plot,rad*scale,colormu,psym=4,ytitle='(\mu_R - \mu_I) [mag/arcsec^2]',xtitle='Radius ["]',yran=[-0.2,1.6],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3
  oploterr,rad*scale,colormu,err
  device,/close


END
