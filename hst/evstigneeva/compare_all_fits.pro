pro compare_all_fits

  zeropt=25.28697
  extinct=0.034
  scale=0.025
  const=(64800./!PI)^2
  Msun=4.53d
  readcol,'vucd3_mge_outputsersic.dat',sersiclum,sersicsig,sersicq,sersicpa,format='D,D,D,D'
  sersicsb=Msun-2.5*alog10(sersiclum/const)
  sersicpeak=10^(-0.4*(sersicsb-zeropt-5*alog10(scale)+extinct))
  readcol,'vucd3_mge_outputsersic_free.dat',lum,sig,q,format='F,F,F'
  sb=Msun-2.5*alog10(lum/const)
  peak=10^(-0.4*(sb-zeropt-5*alog10(scale)+extinct))
  radius=[findgen(20)*0.05+0.05,(10^(findgen(70)*0.025+0.025))]*scale
  minradius=[findgen(21)*0.05,(10^(findgen(69)*0.025+0.025))]*scale
  radius=(radius+minradius)/2.
  ;Reconstruct MGE fits
  flux=fltarr(n_elements(radius))
  sersicflux=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     temp=0.
     for j=0,n_elements(peak)-1 do begin
        temp+= (peak[j]*exp(-(1./(2*sig[j]^2))*(radius[i])^2))

     endfor
     flux[i]=temp
  endfor
  for i=0,n_elements(radius)-1 do begin
     tmp=0.
     for j=0,n_elements(sersicpeak)-1 do begin
        tmp+= (sersicpeak[j]*exp(-(1./(2*sersicsig[j]^2))*(radius[i])^2))
 
     endfor
     sersicflux[i]=tmp
  endfor
  stop
  mgemag_free=zeropt+5*alog10(scale)-2.5*alog10(flux)-extinct
  mgemag_fixed=zeropt+5*alog10(scale)-2.5*alog10(sersicflux)-extinct
  diff=(mgemag_free-mgemag_fixed)

  fits_read,'galfit_model.fits',img,exten_no=1
  fits_read,'galfit_model.fits',freemod,exten_no=2
  fits_read,'galfit_model_fixed.fits',fixmod,exten_no=2
  xc=97.67 & yc=98.95
  radius=[findgen(20)*0.05+0.05,(10^(findgen(70)*0.025+0.025))]
  minradius=[findgen(21)*0.05,(10^(findgen(69)*0.025+0.025))]

  photmag_aper=fltarr(n_elements(radius))
  freemag=fltarr(n_elements(radius))
  fixmag=fltarr(n_elements(radius))
  ;DO PHOTOMETRY
  area=((!PI*(radius)^2)-(!PI*(minradius)^2))
  for i=0,n_elements(radius)-1 do begin
     aper,img,xc,yc,maxflux_aper,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.0,/flux,/exact
     aper,img,xc,yc,minflux_aper,mfluxerr,0.,mskyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.0,/flux,/exact

     aper,freemod,xc,yc,maxflux_free,freefluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.0,/flux,/exact
     aper,freemod,xc,yc,minflux_free,mfreefluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.0,/flux,/exact

     aper,fixmod,xc,yc,maxflux_fix,fixfluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.0,/flux,/exact
     aper,fixmod,xc,yc,minflux_fix,mfixfluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     if (minradius[i] lt 0.025) then begin
        photflux_aper=maxflux_aper
        freeflux=maxflux_free
        fixflux=maxflux_fix
     endif else begin
        photflux_aper=maxflux_aper-minflux_aper
        freeflux=maxflux_free-minflux_free
        fixflux=maxflux_fix-minflux_fix
     endelse
     

     photmag_aper[i]=zeropt+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(photflux_aper)-extinct
     freemag[i]=zeropt+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(freeflux)-extinct
     fixmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(fixflux)-extinct

;     stop
  endfor
  
  ;PHOTOMETRY ON MODELS

  !P.MULTI=[0,1,2]
  set_plot,'ps'
  device,filename='mge_sb_compare.ps',/color
  djs_plot,radius*scale,photmag_aper,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',charthick=4,xthick=3,ythick=3,thick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2
  djs_oplot,radius*scale,mgemag_fixed,color='blue',thick=3
  djs_oplot,radius*scale,mgemag_free,color='red',thick=3
  items=['Data','Free Unconvolved','Fixed Unconvolved']
  lines=[0,0,0]
  sym=[4,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,psym=sym,colors=color,/window,background_color='white',charthick=4,thick=3,/bottom,/left
  djs_plot,radius*scale,diff,psym=2,ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual (free-fixed)',charthick=4,xthick=3,ythick=3,charsize=1.5,ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.1,0.1]
  device,/close
  diff_conv=freemag-fixmag
  device,filename='sb_compare_convolve.ps',/color
  djs_plot,radius*scale,photmag_aper,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',charthick=4,xthick=3,ythick=3,thick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2
  djs_oplot,radius*scale,fixmag,color='blue',thick=3
  djs_oplot,radius*scale,freemag,color='red',thick=3
  items=['Data','Free Convolved','Fixed Convolved']
  lines=[0,0,0]
  sym=[4,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,psym=sym,colors=color,/window,background_color='white',charthick=4,thick=3,/bottom,/left
  djs_plot,radius*scale,diff_conv,psym=2,ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual (free-fixed)',charthick=4,xthick=3,ythick=3,charsize=1.5,ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.1,0.1]

  device,/close

  ratio=mgemag_free/mgemag_fixed

  device,filename='ratio_mge_sersic.ps',/color
  !P.MULTI=[0,1,1]
  djs_plot,radius*scale,ratio,psym=4,yran=[0.995,1.002],/ysty,xtitle='Radius ["]',ytitle='Free Sersic / Fixed Sersic',charthick=4,xthick=3,ythick=3,thick=3,xran=[0.,0.5],/xsty
  device,/close

  fits_read,'deconvolved_f814.fits',deconimg
  deconimg=deconimg/1050.
  find_galaxy,deconimg,majoraxis,eps,ang,xc_d,yc_d
  xc_d=121 & yc_d=121
  deconmag_aper=fltarr(n_elements(radius))

  for i=0,n_elements(radius)-1 do begin
     aper,deconimg,xc_d,yc_d,maxflux_aper,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,deconimg,xc_d,yc_d,minflux_aper,mfluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     if (minradius[i] lt 0.025) then begin
        photflux=maxflux_aper
     endif else begin
        photflux=maxflux_aper-minflux_aper
     endelse
     deconmag_aper[i]=zeropt+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(photflux)-extinct
  endfor
  
  device,filename='deconvolved_compare_mge.ps',/color
  !P.MULTI=[0,1,2]
  djs_plot,radius*scale,deconmag_aper,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',charthick=4,xthick=3,ythick=3,thick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,.15],/xsty,psym=2,title='Deconvolved Image vs. MGE fits'
  djs_oplot,radius*scale,mgemag_fixed,color='blue',thick=3
  djs_oplot,radius*scale,mgemag_free,color='red',thick=3
  items=['Data Deconvolved','Free Unconvolved','Fixed Unconvolved']
  lines=[0,0,0]
  sym=[4,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,psym=sym,colors=color,/window,background_color='white',charthick=4,thick=3,/bottom,/left
  djs_plot,radius*scale,deconmag_aper-mgemag_fixed,psym=2,ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',charthick=4,xthick=3,ythick=3,charsize=1.5,ycharsize=0.6,xran=[0,.15],/xsty,yran=[-0.35,0.6],color='blue'
  djs_oplot,radius*scale,deconmag_aper-mgemag_free,psym=2,color='red'

  device,/close

  zeropoint_i=25.28697
  bn1=1.999*3.25-0.327
  bn2=1.999*1.74-0.327
  radius=(radius+minradius)/2.
  mue1=18.58+5*alog10(2.75*scale)+2.5*alog10(2*!PI*3.25*(exp(bn1)/((bn1)^(2*3.25)))*GAMMA(2*3.25))
  mue2=17.99+5*alog10(24.71*scale)+2.5*alog10(2*!PI*1.74*(exp(bn2)/((bn2)^(2*1.74)))*GAMMA(2*1.74))
  mu1=mue1+((2.5*bn1)/(alog(10)))*((((radius*scale)/(2.75*scale))^(1./3.25))-1)
  mu2=mue2+((2.5*bn2)/(alog(10)))*((((radius*scale)/(24.71*scale))^(1./1.74))-1)
  int1=(10^(-0.4*(mu1-zeropoint_i)))
  int2=(10^(-0.4*(mu2-zeropoint_i)))
  inttot=int1+int2
  mutot=zeropoint_i-2.5*alog10(inttot)-extinct

  bn1=1.999*3.51-0.327
  bn2=1.999*1.28-0.327
  mue1=18.43+5*alog10(3.16*scale)+2.5*alog10(2*!PI*3.51*(exp(bn1)/((bn1)^(2*3.51)))*GAMMA(2*3.51))
  mue2=18.14+5*alog10(24.47*scale)+2.5*alog10(2*!PI*1.28*(exp(bn2)/((bn2)^(2*1.28)))*GAMMA(2*1.28))
  mu1=mue1+((2.5*bn1)/(alog(10)))*((((radius)/(3.16))^(1./3.51))-1)
  mu2=mue2+((2.5*bn2)/(alog(10)))*((((radius)/(24.47))^(1./1.28))-1)
  int1=(10^(-0.4*(mu1-zeropoint_i)))
  int2=(10^(-0.4*(mu2-zeropoint_i)))
  inttot=int1+int2
  mutot_fix=zeropoint_i-2.5*alog10(inttot)-extinct

  device,filename='deconvolved_compare_sersic.ps',/color
  djs_plot,radius*scale,deconmag_aper,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',charthick=4,xthick=3,ythick=3,thick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,0.15],/xsty,psym=2,title='Deconvolved Image vs. Sersic fits'
  djs_oplot,radius*scale,mutot,color='red',thick=3
  djs_oplot,radius*scale,mutot_fix,color='blue',thick=3
  items=['Data Deconvolved','Free Unconvolved','Fixed Unconvolved']
  lines=[0,0,0]
  sym=[4,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,psym=sym,colors=color,/window,background_color='white',charthick=4,thick=3,/bottom,/left
  djs_plot,radius*scale,deconmag_aper-mutot_fix,psym=2,ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',charthick=4,xthick=3,ythick=3,charsize=1.5,ycharsize=0.6,xran=[0,0.15],/xsty,yran=[-0.75,0.25],color='blue'
  djs_oplot,radius*scale,deconmag_aper-mutot,psym=2,color='red'

  device,/close
  set_plot,'x'
  !P.MULTI=[0,1,1]



  stop

END
