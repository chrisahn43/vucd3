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
  fits_read,'galfit_model.fits',img_i,h_i,exten_no=1
  fits_read,'galfit_model_r.fits',img_r,h_r,exten_no=1
  fits_read,'galfit_model.fits',model_i,exten_no=2
  fits_read,'galfit_model_r.fits',model_r,exten_no=2
  fits_read,'vucd3_weightmap.fits',err_i
  fits_read,'vucd3_weightmap_r.fits',err_r
  xc_i=98 & yc_i=99
  xc_r=98 & yc_r=99;8
  scale=0.025
  zeropoint_r=25.9799
  zeropoint_i=25.28697
  rad=findgen(100)+1
  minrad=findgen(100)
  rflux=fltarr(n_elements(rad))
  rmodel=fltarr(n_elements(rad))
  iflux=fltarr(n_elements(rad))
  imodel=fltarr(n_elements(rad))
  rmag=fltarr(n_elements(rad))
  imag=fltarr(n_elements(rad))
  rmagmod=fltarr(n_elements(rad))
  imagmod=fltarr(n_elements(rad))
  rmagerr=fltarr(n_elements(rad))
  imagerr=fltarr(n_elements(rad))
  area=fltarr(n_elements(rad))
  for i=0,n_elements(rad)-1 do begin
     maxflux_r=djs_phot(xc_r,yc_r,rad[i],0.,img_r,skyval=skyval)
     minflux_r=djs_phot(xc_r,yc_r,minrad[i],0.,img_r,skyval=skyval)
     rflux[i]=maxflux_r-minflux_r
     maxfluxerr_r=djs_phot(xc_r,yc_r,rad[i],0.,err_r,skyval=skyval)
     minfluxerr_r=djs_phot(xc_r,yc_r,minrad[i],0.,err_r,skyval=skyval)
     maxfluxmod_r=djs_phot(xc_r,yc_r,rad[i],0.,model_r,skyval=skyval)
     minfluxmod_r=djs_phot(xc_r,yc_r,minrad[i],0.,model_r,skyval=skyval)
     rmodel[i]=maxfluxmod_r-minfluxmod_r
     
     maxflux_i=djs_phot(xc_i,yc_i,rad[i],0.,img_i,skyval=skyval)
     minflux_i=djs_phot(xc_i,yc_i,minrad[i],0.,img_i,skyval=skyval)  
     iflux[i]=maxflux_i-minflux_i
     maxfluxerr_i=djs_phot(xc_i,yc_i,rad[i],0.,err_i,skyval=skyval)
     minfluxerr_i=djs_phot(xc_i,yc_i,minrad[i],0.,err_i,skyval=skyval)
     maxfluxmod_i=djs_phot(xc_i,yc_i,rad[i],0.,model_i,skyval=skyval)
     minfluxmod_i=djs_phot(xc_i,yc_i,minrad[i],0.,model_i,skyval=skyval)
     imodel[i]=maxfluxmod_i-minfluxmod_i

     area[i]=((!PI*(rad[i])^2)-(!PI*(minrad[i])^2))
     ierrflux=(maxfluxerr_i-minfluxerr_i)/area[i]
     rerrflux=(maxfluxerr_r-minfluxerr_r)/area[i]

     rmag[i]=zeropoint_r + 5*alog10(scale) + 2.5*alog10(area[i]) - 2.5*alog10(rflux[i])-0.034
;     rmagerr[i]=5*alog10(scale) +(((-2.5*alog10(rflux[i]+rerrflux))+(2.5*alog10(rflux[i]-rerrflux)))/2.)
     rmagerr[i]=5*alog10(scale) -2.5*(rerrflux/rflux[i])
     rmagmod[i]=zeropoint_r+5*alog10(scale) + 2.5*alog10(area[i]) - 2.5*alog10(rmodel[i])-0.034
     
     imag[i]=zeropoint_i + 5*alog10(scale) + 2.5*alog10(area[i]) - 2.5*alog10(iflux[i])-0.048
;     imagerr[i]=5*alog10(scale) + (((-2.5*alog10(iflux[i]+ierrflux))+(2.5*alog10(iflux[i]-ierrflux)))/2.)
     imagerr[i]=5*alog10(scale) -2.5*(ierrflux/iflux[i])
     imagmod[i]=zeropoint_i + 5*alog10(scale) + 2.5*alog10(area[i]) - 2.5*alog10(imodel[i])-0.034

  endfor
  colormag=rmagmod-imagmod
  err=rmagerr-imagerr
  resid_i=imag-imagmod
  resid_r=rmag-rmagmod



                                  ;unconvolved mge sb comparison
                                ;F814W
  mue1=18.2+2.5*alog10(2*!PI*(5.37*scale)^2)
  mue1=mue1+1.3
  mue2=18.49+2.5*alog10(2*!PI*(32.17*scale)^2)
  mue2=mue2+0.5
  bn1=1.999*4.6-0.327
  bn2=1.999*1.05-0.327
  mu1=mue1+((2.5*bn1)/(alog(10)))*((((rad)/(5.37))^(1./4.6))-1)
  mu2=mue2+((2.5*bn2)/(alog(10)))*((((rad)/(32.17))^(1./1.05))-1)
  int1=(10^(-0.4*(mu1-zeropoint_i)))
  int2=(10^(-0.4*(mu2-zeropoint_i)))
  inttot=int1+int2
  mutot=zeropoint_i-2.5*alog10(inttot)
;  mue1=18.61+2.5*alog10(2*!PI*(2.99*scale)^2)
;  mue1=mue1+1.1
;  mue2=18.09+2.5*alog10(2*!PI*(26.56*scale)^2)
;  mue2=mue2+0.6
;  bn1=1.999*3.44-0.327
;  bn2=1.999*1.72-0.327
;  mu1=mue1+((2.5*bn1)/(alog(10)))*((((rad)/(2.99))^(1./3.44))-1)
;  mu2=mue2+((2.5*bn2)/(alog(10)))*((((rad)/(26.56))^(1./1.72))-1)
;  int1=(10^(-0.4*(mu1-zeropoint_i)))
;  int2=(10^(-0.4*(mu2-zeropoint_i)))
;  inttot=int1+int2
;  mutot=zeropoint_i-2.5*alog10(inttot)

;  mue1606=19.04+2.5*alog10(2*!PI*(3.36*scale)^2)
;  mue1606=mue1606+1.1
;  mue2606=18.95+2.5*alog10(2*!PI*(25.2*scale)^2)
;  mue2606=mue2606+0.6
;  bn1606=1.999*3.67-0.327
;  bn2606=1.999*1.25-0.327
;  mu1606=mue1606+((2.5*bn1606)/(alog(10)))*((((rad)/(3.36))^(1./3.67))-1)
;  mu2606=mue2606+((2.5*bn2606)/(alog(10)))*((((rad)/(25.2))^(1./1.25))-1)
;  int1606=(10^(-0.4*(mu1606-zeropoint_r)))
;  int2606=(10^(-0.4*(mu2606-zeropoint_r)))
;  inttot606=int1606+int2606
;  mutot606=zeropoint_r-2.5*alog10(inttot606)
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
  set_plot,'ps'
  device,filename='radialcolor_models.ps',/color
  djs_plot,rad*scale,colormag,psym=4,ytitle='(\mu_{F606W} - \mu_{F814W}) [mag/arcsec^2]',xtitle='Radius ["]',yran=[-0.2,1.6],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,xran=[0.,1.3]
  oploterr,rad*scale,colormag,err
  djs_oplot,rad*scale,colormu,color='blue',thick=3
  device,/close
  stop
  device,filename='residual_i.ps',/color
  djs_plot,rad*scale,resid_i,psym=4,ytitle='\Delta \mu',xtitle='Radius ["]',yran=[-0.5,0.5],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,thick=2,xran=[0.,1.3]
  djs_oplot,rad*scale,resid_r,psym=4,color='blue',thick=2
  device,/close
  set_plot,'x'
  
  stop

;  device,filename='radialcolor_sersic.ps',/color
;  djs_plot,rad*scale,colormu,psym=4,ytitle='(\mu_R - \mu_I) [mag/arcsec^2]',xtitle='Radius ["]',yran=[-0.2,1.6],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3
;  oploterr,rad*scale,colormu,err
;  device,/close


END

pro mlvscolor

  scale=0.025
  zeropoint_r_ab=25.9799
  zeropoint_i_ab=25.28697
  zeropoint_r_ve=25.90069
  zeropoint_i_ve=24.86147
  hstcolorin=0.61
  hstcolorout=0.74
  vsun=4.80
  isun=4.10
  readcol,'johnson_ssp.dat',johnz,johnage,johnmbol,u,b,v,r,i,j,h,k,format='F,F,F,F,F,F,F,F,F,F,F'
  readcol,'hstacs_ssp.dat',hstz,hstage,hstmbol,f220w,f250w,f330w,f334n,f435w,f475w,f550m,f555w,f606w,f625w,f658w,f660n,f775w,f814w,f850lp,f892n,format='F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  hstcolor=(f606w+zeropoint_r_ab-zeropoint_r_ve)-(f814w+zeropoint_i_ab-zeropoint_i_ve)
  johncolor=v-i
  djs_plot,hstcolor,johncolor,psym=4
  fit=linfit(hstcolor,johncolor)
  y=hstcolor*fit[1]+fit[0]
  djs_oplot,hstcolor,y,color='blue'

  colorin=(hstcolorin*fit[1]+fit[0])-(0.061-0.034)
  colorout=(hstcolorout*fit[1]+fit[0])-(0.061-0.034)
  stop
  
  readcol,'./bc03/chabrier/bc2003_hr_m62_chab_ssp.2color',logage2,vi,vj,vk,jh,/SILENT,FORMAT='F,X,X,X,X,F,F,F,X,F'
;  readcol,'~/Desktop/test.2color',logage2,vi,vj,vk,jh,/SILENT,FORMAT='F,X,X,X,X,F,F,F,X,F'

  readcol,'./bc03/chabrier/bc2003_hr_m62_chab_ssp.4color',logage4,bmag,vmag,mlb,mlv,/SILENT,FORMAT='F,X,F,F,F,F'
;  readcol,'~/Desktop/test.4color',logage4,bmag,vmag,mlb,mlv,/SILENT,FORMAT='F,X,F,F,F,F'

  c=where(vi gt colorin-0.003 and vi lt colorin+0.003) ; and logage2 gt 8.5)
;  d=where(vi gt colorout-0.083 and vi lt colorout+0.083); and logage2 gt 6.5)
  d=where(vi gt colorout-0.07 and vi lt colorout+0.07)
  mlvin=mean(mlv[c])
  mlvout=mean(mlv[d])
  mliin=mlvin*10^(-0.4*(colorin-(vsun-isun)))
  mliout=mlvout*10^(-0.4*(colorout-(vsun-isun)))
  print,mliin,mliout

  readcol,'./evstigneeva/vucd3_mge_outputsersic.dat',intensity,sigma,q,pa,format='D,D,D,D'
  a=where(q lt 0.8)
  intensity[a]=intensity[a]*mliin
  b=where(q gt 0.8)
  intensity[b]=intensity[b]*mliout
  forprint,intensity,sigma,q,pa,textout='./evstigneeva/vucd3_mge_outputsersic_mass.dat',format='D,D,D,D'

  stop
  colorin=0.61-(0.061-0.034)
  colorout=0.74-(0.061-0.034)
  vsun=4.80
  isun=4.53
  readcol,'~/bc03/bc03/hst_add.1ABmag',logage,u,g,r,i,z,f606w,f814w,f475w,f850lp,format='F,F,F,F,F,F,F,F,F,F'
  hstcolor=f606w-f814w
  readcol,'~/bc03/bc03/hst_add.4color',logage4,mbol,bmag,vmag,kmag,mliv,mrem,mret,mgal,sfr,mtot,mlbtot,mlvtot,format='F,F,F,F,F,F,F,F,F,F,F,F,F'
  inind=where(hstcolor gt colorin-0.001 and hstcolor lt colorin+0.001)
  outind=where(hstcolor gt colorout-0.05 and hstcolor lt colorout+0.05)
  iinlum=(10^(-0.4*(f814w[inind]-isun)))
  ioutlum=(10^(-0.4*(f814w[outind]-isun)))

  vinlum=(10^(-0.4*(vmag[inind]-vsun)))
  voutlum=(10^(-0.4*(vmag[outind]-vsun)))

  mlin=mlvtot[inind]*(vinlum/iinlum)
  mlin=mlin[0]
  mlout=mlvtot[outind]*(voutlum/ioutlum)
  mlout=mlout[0]
  print,mlin,mlout
  readcol,'./evstigneeva/vucd3_mge_outputsersic.dat',intensity,sigma,q,pa,format='D,D,D,D'
  a=where(q lt 0.8)
  intensity[a]=intensity[a]*mlin
  b=where(q gt 0.8)
  intensity[b]=intensity[b]*mlout
  forprint,intensity,sigma,q,pa,textout='./evstigneeva/vucd3_mge_outputsersic_nmass.dat',format='D,D,D,D'
  stop
;  colorin=(hstcolorin*fit[1]+fit[0])-(0.061-0.034)
;  colorout=(hstcolorout*fit[1]+fit[0])-(0.061-0.034)
;  f606w=f606w+zeropoint_r_ab-zeropoint_r_ve
;  f814w=f814w+zeropoint_i_ab-zeropoint_i_ve
;  hstcolor=f606w-f814w
;  inind=where(hstcolor gt colorin-0.001 and hstcolor lt colorin+0.001)
;  outind=where(hstcolor gt colorout-0.05 and hstcolor lt colorout+0.05)
;  isun=4.10
;  iinlum=(10^(-0.4*(f814w[inind]-isun)))
;  ioutlum=(10^(-0.4*(f814w[outind]-isun)))
;  mlin=mlvtot[inind]*(vinlum/iinlum)
;  mlin=mlin[0]
;  mlout=mlvtot[outind]*(voutlum/ioutlum)
;  mlout=mlout[0]
;  print,mlin,mlout
 
;  stop
;FIND COLORS AT RADII GIVEN BY MGE OUTPUTS FREE-FREE

  readcol,'./evstigneeva/vucd3_mge_outputsersic_free.dat',ilum,isig,iq,ipa,format='D,D,D,D'
  radius=isig/scale
  radius[24]=118.
  minradius=[0,radius[0:13],0,radius[15:23]]
  fits_read,'./make_decon/mgefreemodeli.fits',iimg
  fits_read,'./make_decon/mgefreemodelv.fits',vimg
  find_galaxy,iimg,m,e,a,xci,yci
  find_galaxy,vimg,m,e,a,xcv,ycv
  vmag=fltarr(n_elements(radius))
  imag=fltarr(n_elements(radius))

  for i=0,n_elements(radius)-1 do begin
     area=((!PI*radius[i]^2)-(!PI*minradius[i]^2))
     
     aper,iimg,xci,yci,imaxflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,iimg,xci,yci,iminflux,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,vimg,xcv,ycv,vmaxflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,vimg,xcv,ycv,vminflux,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     if (minradius[i] eq 0.) then begin
        iflux=imaxflux
        vflux=vmaxflux
     endif else begin
        iflux=imaxflux-iminflux
        vflux=vmaxflux-vminflux
;        stop
     endelse
     if (vflux lt 0.) then stop
     vmag[i]=zeropoint_r_ab+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(vflux)-0.061
     imag[i]=zeropoint_i_ab+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(iflux)-0.034
  endfor
  color=vmag-imag

  color=color*fit[1]+fit[0]
  readcol,'johnson_ssp.dat',johnz,johnage,johnmbol,u,b,v,r,i,j,h,k,format='F,F,F,F,F,F,F,F,F,F,F'

;  extinct=2.
;  vext=extinct*1.00600
;  v=v+vext
;  iext=extinct*0.60329
;  i=i+iext
  masstolight=fltarr(n_elements(color))
  for i=0,n_elements(color)-1 do begin
     cind=where(vi lt color[i]+0.001 and vi gt color[i]-0.001,count)
     if (count ge 1.) then begin
        masstolight[i]=mean(mlv[cind])
     endif else begin
        cind=where(vi lt color[i]+0.003 and vi gt color[i]-0.003,count)
        if (count ge 1.) then begin
           masstolight[i]=mean(mlv[cind])
        endif else begin
           cind=where(vi lt color[i]+0.01 and vi gt color[i]-0.01,count)
           if (count ge 1.) then begin
              masstolight[i]=mean(mlv[cind])
           endif else begin
              cind=where(vi lt color[i]+0.075 and vi gt color[i]-0.075,count)
              if (count ge 1.) then begin
                 masstolight[i]=mean(mlv[cind])
              endif else begin
                 cind=where(vi eq max(vi))
                 masstolight[i]=mean(mlv[cind])
              endelse 
           endelse
        endelse 
     endelse
     
  endfor
  masstolight=masstolight*10^(-0.4*(color-(vsun-isun)))
  readcol,'./evstigneeva/vucd3_mge_outputsersic_free.dat',intensity,sigma,q,pa,format='D,D,D,D'
  intensity=intensity*masstolight
  forprint,intensity,sigma,q,pa,textout='./evstigneeva/vucd3_mge_outputsersic_mass_free.dat',format='D,D,D,D'

  stop
  color=vmag-imag
  masstolight=fltarr(n_elements(color))
  readcol,'~/bc03/bc03/hst_add.1ABmag',logage,u,g,r,i,z,f606w,f814w,f475w,f850lp,format='F,F,F,F,F,F,F,F,F,F'
  hstcolor=f606w-f814w
  isun=4.53
  vsun=4.80
  readcol,'~/bc03/bc03/hst_add.4color',logage4,mbol,bmag,vmag,kmag,mliv,mrem,mret,mgal,sfr,mtot,mlbtot,mlvtot,format='F,F,F,F,F,F,F,F,F,F,F,F,F'

  for i=0,n_elements(color)-1 do begin
     cind=where(hstcolor lt color[i]+0.001 and hstcolor gt color[i]-0.001,count)
     if (count ge 1.) then begin
        ilum=mean(10^((-0.4*(f814w[cind]-isun))))
        vlum=mean(10^((-0.4*(vmag[cind]-vsun))))
        masstolight[i]=mean(mlvtot[cind]*(vlum/ilum))
     endif else begin
        cind=where(hstcolor lt color[i]+0.003 and hstcolor gt color[i]-0.003,count)
        if (count ge 1.) then begin
           ilum=mean(10^((-0.4*(f814w[cind]-isun))))
           vlum=mean(10^((-0.4*(vmag[cind]-vsun))))
           masstolight[i]=mean(mlvtot[cind]*(vlum/ilum))
        endif else begin
           cind=where(hstcolor lt color[i]+0.01 and hstcolor gt color[i]-0.01,count)
           if (count ge 1.) then begin
              ilum=mean(10^((-0.4*(f814w[cind]-isun))))
              vlum=mean(10^((-0.4*(vmag[cind]-vsun))))
              masstolight[i]=mean(mlvtot[cind]*(vlum/ilum))
           endif else begin
              cind=where(hstcolor lt color[i]+0.075 and hstcolor gt color[i]-0.075,count)
              if (count ge 1.) then begin
                 ilum=mean(10^((-0.4*(f814w[cind]-isun))))
                 vlum=mean(10^((-0.4*(vmag[cind]-vsun))))
                 masstolight[i]=mean(mlvtot[cind]*(vlum/ilum))
              endif else begin
                 cind=where(hstcolor eq max(hstcolor))
                 ilum=mean(10^((-0.4*(f814w[cind]-isun))))
                 vlum=mean(10^((-0.4*(vmag[cind]-vsun))))
                 masstolight[i]=mean(mlvtot[cind]*(vlum/ilum))
              endelse 
           endelse
        endelse 
     endelse
     
  endfor
  
  readcol,'./evstigneeva/vucd3_mge_outputsersic_free.dat',intensity,sigma,q,pa,format='D,D,D,D'
  intensity=intensity*masstolight
  forprint,intensity,sigma,q,pa,textout='./evstigneeva/vucd3_mge_outputsersic_nmass_free.dat',format='D,D,D,D'
  
  stop


; rad=findgen(66)+1
 ; mue1=18.2+2.5*alog10(2*!PI*(5.37*scale)^2)
 ; mue1=mue1+1.3
 ; mue2=18.49+2.5*alog10(2*!PI*(32.17*scale)^2)
 ; mue2=mue2+0.5
 ; bn1=1.999*4.6-0.327
 ; bn2=1.999*1.05-0.327
 ; mu1=mue1+((2.5*bn1)/(alog(10)))*((((rad)/(5.37))^(1./4.6))-1)
 ; mu2=mue2+((2.5*bn2)/(alog(10)))*((((rad)/(32.17))^(1./1.05))-1)
 ; int1=(10^(-0.4*(mu1-zeropoint_i)))
 ; int2=(10^(-0.4*(mu2-zeropoint_i)))
 ; inttot=int1+int2
 ; mutot=zeropoint_i-2.5*alog10(inttot)-0.02037518

 ; mue1606=18.82+2.5*alog10(2*!PI*(4.66*scale)^2)
 ; mue1606=mue1606+1.3
 ; mue2606=19.2+2.5*alog10(2*!PI*(27.71*scale)^2)
 ; mue2606=mue2606+0.5
 ; bn1606=1.999*4.66-0.327
 ; bn2606=1.999*0.98-0.327
 ; mu1606=mue1606+((2.5*bn1606)/(alog(10)))*((((rad)/(4.66))^(1./4.46))-1)
 ; mu2606=mue2606+((2.5*bn2606)/(alog(10)))*((((rad)/(27.71))^(1./0.98))-1)
 ; int1606=(10^(-0.4*(mu1606-zeropoint_r)))
 ; int2606=(10^(-0.4*(mu2606-zeropoint_r)))
 ; inttot606=int1606+int2606
 ; mutot606=zeropoint_r-2.5*alog10(inttot606)-0.0566446;;

;  colormu=mutot606-mutot
;  colorconv=((colormu*fit[1])+fit[0]) ;-(0.061-0.034)
;  djs_plot,rad,colormu
;  rad=rad*scale
;  re1=5.37*scale
;  re2=32.17*scale
;  a=where(rad lt re1+0.01 and rad gt re1-0.01)
;  b=where(rad lt re2+0.01 and rad gt re2-0.01)
;  a=where(rad lt 0.4)
;  b=where(rad ge 0.4)
;  agein=alog10(mean(johnage[a]))
;  ageout=alog10(mean(johnage[b]))
;  colorin=mean(colorconv[a])
;  colorout=mean(colorconv[b])
;  stop;

;  djs_plot,colormu,colorconv,psym=2
;  stop
;  mue1=18.6+2.5*alog10(2*!PI*(3.02*scale)^2)
;  mue1=mue1+1.1
;  mue2=18.1+2.5*alog10(2*!PI*(26.62*scale)^2)
;  mue2=mue2+0.6
;  bn1=1.999*3.46-0.327
;  bn2=1.999*1.71-0.327
;  mu1=mue1+((2.5*bn1)/(alog(10)))*((((rad)/(3.02))^(1./3.46))-1)
;  mu2=mue2+((2.5*bn2)/(alog(10)))*((((rad)/(26.62))^(1./1.71))-1)
;  int1=(10^(-0.4*(mu1-zeropoint_i)))
;  int2=(10^(-0.4*(mu2-zeropoint_i)))
;  inttot=int1+int2
;  mutot=zeropoint_i-2.5*alog10(inttot)-0.02037518

 ; mue1606=19.04+2.5*alog10(2*!PI*(3.36*scale)^2)
;  mue1606=mue1606+1.1
;  mue2606=18.95+2.5*alog10(2*!PI*(25.2*scale)^2)
;  mue2606=mue2606+0.6
;  bn1606=1.999*3.67-0.327
;  bn2606=1.999*1.25-0.327
;  mu1606=mue1606+((2.5*bn1606)/(alog(10)))*((((rad)/(3.36))^(1./3.67))-1)
;  mu2606=mue2606+((2.5*bn2606)/(alog(10)))*((((rad)/(25.2))^(1./1.25))-1)
;  int1606=(10^(-0.4*(mu1606-zeropoint_r)))
;  int2606=(10^(-0.4*(mu2606-zeropoint_r)))
;  inttot606=int1606+int2606
;  mutot606=zeropoint_r-2.5*alog10(inttot606)-0.0566446

;  colormu=mutot606-mutot
;  colorconv=((colormu*fit[1])+fit[0]) ;-(0.061-0.034)
;  re1=3.02*scale
;  re2=26.62*scale
;  a=where(rad lt re1+0.01 and rad gt re1-0.01)
;  b=where(rad lt re2+0.01 and rad gt re2-0.01)
;  agein=alog10(mean(johnage[a]))
;  ageout=alog10(mean(johnage[b]))
;  colorin=mean(colorconv[a])
;  colorout=mean(colorconv[b])

;  stop
;  readcol,'./bc03/chabrier/bc2003_hr_m62_chab_ssp.2color',logage2,vi,vj,vk,jh,/SILENT,FORMAT='F,X,X,X,X,F,F,F,X,F'
;  readcol,'./bc03/chabrier/bc2003_hr_m62_chab_ssp.4color',logage4,bmag,vmag,mlb,mlv,/SILENT,FORMAT='F,X,F,F,F,F'
;  c=where(logage2 gt agein-0.2 and logage2 lt agein+0.2 and vi gt colorin-0.2 and vi lt colorin+0.2)
;  d=where(logage2 gt ageout-0.1 and logage2 lt ageout+0.1 and vi gt colorout-0.2 and vi lt colorout+0.2)

;  print,mean(mlv[c]),mean(mlv[d])
;  stop
  
END

PRO colorprofile_all
  scale=0.025
  rad=findgen(100)+1
  zeropoint_r=25.9799
  zeropoint_i=25.28697
;ORIGINAL OUTPUTS
  mue1=18.2+2.5*alog10(2*!PI*(5.37*scale)^2)
  mue1=mue1+1.3
  mue2=18.49+2.5*alog10(2*!PI*(32.17*scale)^2)
  mue2=mue2+0.5
  bn1=1.999*4.6-0.327
  bn2=1.999*1.05-0.327
  mu1=mue1+((2.5*bn1)/(alog(10)))*((((rad)/(5.37))^(1./4.6))-1)
  mu2=mue2+((2.5*bn2)/(alog(10)))*((((rad)/(32.17))^(1./1.05))-1)
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

  colormu=mutot606-mutot
  set_plot,'ps'
  device,filename='color_allpsf.ps',/color
  djs_plot,rad,colormu,ytitle='(\mu_{F606W} - \mu_{F814W}) [mag/arcsec^2]',xtitle='Radius [pixels]',yran=[0.3,1.5],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,thick=3,/xlog;,xran=[0.9,100]

  mue1=18.6+2.5*alog10(2*!PI*(3.02*scale)^2)
  mue1=mue1+1.1
  mue2=18.1+2.5*alog10(2*!PI*(26.62*scale)^2)
  mue2=mue2+0.6
  bn1=1.999*3.46-0.327
  bn2=1.999*1.71-0.327
  mu1=mue1+((2.5*bn1)/(alog(10)))*((((rad)/(3.02))^(1./3.46))-1)
  mu2=mue2+((2.5*bn2)/(alog(10)))*((((rad)/(26.62))^(1./1.71))-1)
  int1=(10^(-0.4*(mu1-zeropoint_i)))
  int2=(10^(-0.4*(mu2-zeropoint_i)))
  inttot=int1+int2
  mutot=zeropoint_i-2.5*alog10(inttot)

  mue1606=19.04+2.5*alog10(2*!PI*(3.36*scale)^2)
  mue1606=mue1606+1.1
  mue2606=18.95+2.5*alog10(2*!PI*(25.2*scale)^2)
  mue2606=mue2606+0.6
  bn1606=1.999*3.67-0.327
  bn2606=1.999*1.25-0.327
  mu1606=mue1606+((2.5*bn1606)/(alog(10)))*((((rad)/(3.36))^(1./3.67))-1)
  mu2606=mue2606+((2.5*bn2606)/(alog(10)))*((((rad)/(25.2))^(1./1.25))-1)
  int1606=(10^(-0.4*(mu1606-zeropoint_r)))
  int2606=(10^(-0.4*(mu2606-zeropoint_r)))
  inttot606=int1606+int2606
  mutot606=zeropoint_r-2.5*alog10(inttot606)
  colormu_10=mutot606-mutot
  djs_oplot,rad,colormu_10,color='blue',thick=3,/xlog
  device,/close
  set_plot,'x'
stop
  mue1=18.2+2.5*alog10(2*!PI*(5.37*scale)^2)
  mue1=mue1+1.3
  mue2=18.49+2.5*alog10(2*!PI*(32.17*scale)^2)
  mue2=mue2+0.5
  bn1=1.999*4.6-0.327
  bn2=1.999*1.05-0.327
  mu1=mue1+((2.5*bn1)/(alog(10)))*((((rad)/(5.37))^(1./4.6))-1)
  mu2=mue2+((2.5*bn2)/(alog(10)))*((((rad)/(32.17))^(1./1.05))-1)
  int1=(10^(-0.4*(mu1-zeropoint_i)))
  int2=(10^(-0.4*(mu2-zeropoint_i)))
  inttot=int1+int2
  mutot=zeropoint_i-2.5*alog10(inttot)
  mue1606=18.73+2.5*alog10(2*!PI*(5.37*scale)^2)
  mue1606=mue1606+1.3
  mue2606=19.26+2.5*alog10(2*!PI*(32.17*scale)^2)
  mue2606=mue2606+0.5
  bn1606=1.999*4.6-0.327
  bn2606=1.999*1.05-0.327
  mu1606=mue1606+((2.5*bn1606)/(alog(10)))*((((rad)/(5.37))^(1./4.6))-1)
  mu2606=mue2606+((2.5*bn2606)/(alog(10)))*((((rad)/(32.17))^(1./1.05))-1)
  int1606=(10^(-0.4*(mu1606-zeropoint_r)))
  int2606=(10^(-0.4*(mu2606-zeropoint_r)))
  inttot606=int1606+int2606
  mutot606=zeropoint_r-2.5*alog10(inttot606)

  colormu=mutot606-mutot
  set_plot,'ps'
  device,filename='color_allpsf.ps',/color
  djs_plot,rad,colormu,ytitle='(\mu_{F606W} - \mu_{F814W}) [mag/arcsec^2]',xtitle='Radius [Pixels]',yran=[0.3,1.5],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,thick=3,/xlog;,xran=[0.9,100]

  mue1=18.6+2.5*alog10(2*!PI*(3.02*scale)^2)
  mue1=mue1+1.1
  mue2=18.1+2.5*alog10(2*!PI*(26.62*scale)^2)
  mue2=mue2+0.6
  bn1=1.999*3.46-0.327
  bn2=1.999*1.71-0.327
  mu1=mue1+((2.5*bn1)/(alog(10)))*((((rad)/(3.02))^(1./3.46))-1)
  mu2=mue2+((2.5*bn2)/(alog(10)))*((((rad)/(26.62))^(1./1.71))-1)
  int1=(10^(-0.4*(mu1-zeropoint_i)))
  int2=(10^(-0.4*(mu2-zeropoint_i)))
  inttot=int1+int2
  mutot=zeropoint_i-2.5*alog10(inttot)

  mue1606=19.14+2.5*alog10(2*!PI*(3.02*scale)^2)
  mue1606=mue1606+1.1
  mue2606=18.79+2.5*alog10(2*!PI*(26.62*scale)^2)
  mue2606=mue2606+0.6
  bn1606=1.999*3.46-0.327
  bn2606=1.999*1.71-0.327
  mu1606=mue1606+((2.5*bn1606)/(alog(10)))*((((rad)/(3.02))^(1./3.46))-1)
  mu2606=mue2606+((2.5*bn2606)/(alog(10)))*((((rad)/(26.62))^(1./1.71))-1)
  int1606=(10^(-0.4*(mu1606-zeropoint_r)))
  int2606=(10^(-0.4*(mu2606-zeropoint_r)))
  inttot606=int1606+int2606
  mutot606=zeropoint_r-2.5*alog10(inttot606)
  colormu_10=mutot606-mutot
  djs_oplot,rad,colormu_10,color='blue',thick=3,/xlog
  device,/close
  set_plot,'x'
stop

END
