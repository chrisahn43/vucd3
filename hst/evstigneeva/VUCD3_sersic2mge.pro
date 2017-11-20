@sersic2mge.pro
pro VUCD3_sersic2mge
  scale=0.025
  Msun=4.53d
  A_B=0.034d
  zeropt=25.28697
  extinct=0.034

  re1=3.16d*scale
  re2=24.47d*scale
  sersic1={mag:18.43d,re:re1,n:3.51d,pa:19.04d,q:0.66d}
  sersic2={mag:18.14d,re:re2,n:1.28d,pa:18.43d,q:0.91d}
  sersics=[sersic1,sersic2]
  mge=sersics2mge(sersics,Msun,A_B=A_B)
  forprint,mge[*,0],mge[*,1],mge[*,2],mge[*,3],textout='vucd3_mge_outputsersic.dat',format='F,F,F,F'
  re1=2.75d*scale
  re2=24.71d*scale
  sersic1={mag:18.58d,re:re1,n:3.25d,pa:17.97d,q:0.62d}
  sersic2={mag:17.99d,re:re2,n:1.74d,pa:20.65d,q:0.89d}
  sersics=[sersic1,sersic2]
  mge=sersics2mge(sersics,Msun,A_B=A_B)
  forprint,mge[*,0],mge[*,1],mge[*,2],mge[*,3],textout='vucd3_mge_outputsersic_free.dat',format='F,F,F,F'

  Msun=4.73d
  A_B=0.061
  zeropt=25.9799
  re1=3.16d*scale
  re2=24.47d*scale
  sersic1={mag:19.05d,re:re1,n:3.51d,pa:19.04d,q:0.66d}
  sersic2={mag:18.89d,re:re2,n:1.28d,pa:18.43d,q:0.91d}
  sersics=[sersic1,sersic2]
  mge=sersics2mge(sersics,Msun,A_B=A_B)
  forprint,mge[*,0],mge[*,1],mge[*,2],mge[*,3],textout='vucd3_mge_outputsersic_free_606.dat',format='F,F,F,F'
  re1=2.75d*scale
  re2=24.71d*scale
  sersic1={mag:19.19d,re:re1,n:3.25d,pa:17.97d,q:0.62d}
  sersic2={mag:18.73d,re:re2,n:1.74d,pa:20.65d,q:0.89d}
  sersics=[sersic1,sersic2]
  mge=sersics2mge(sersics,Msun,A_B=A_B)
  forprint,mge[*,0],mge[*,1],mge[*,2],mge[*,3],textout='vucd3_mge_outputsersic_606.dat',format='F,F,F,F'

  
  
  
  stop
  const=(64800./!PI)^2
  readcol,'vucd3_mge_outputsersic.dat',sersiclum,sersicsig,sersicq,sersicpa,format='D,D,D,D'
  sersicsb=Msun-2.5*alog10(sersiclum/const)
  sersicpeak=10^(-0.4*(sersicsb-zeropt-5*alog10(scale)+extinct))
  readcol,'vucd3_mge_outputsersic_free.dat',lum,sig,q,format='F,F,F'
  sb=Msun-2.5*alog10(lum/const)
  peak=10^(-0.4*(sb-zeropt-5*alog10(scale)+extinct))
  radius=(10^(findgen(70)*0.025+0.025))*scale
  flux=fltarr(n_elements(radius))
  sersicflux=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     temp=0.
     for j=0,n_elements(peak)-1 do begin
        temp+= (peak[j]*exp(-(1./(2*sig[j]^2))*(radius[i])^2))
        flux[i]=temp
     endfor
  endfor
  for i=0,n_elements(radius)-1 do begin
     tmp=0.
     for j=0,n_elements(sersicpeak)-1 do begin
        tmp+= (sersicpeak[j]*exp(-(1./(2*sersicsig[j]^2))*(radius[i])^2))
        sersicflux[i]=tmp
     endfor
  endfor
  
  mag=zeropt+5*alog10(scale)-2.5*alog10(flux)-extinct
  sersicmag=zeropt+5*alog10(scale)-2.5*alog10(sersicflux)-extinct
  diff=(mag-sersicmag)

  
  fits_read,'galfit_model.fits',img,exten_no=1
  fits_read,'galfit_model.fits',freemod,exten_no=2
  fits_read,'galfit_model_fixed.fits',fixmod,exten_no=2
  xc=98 & yc=99
  radius=radius/scale
  minradius=[0,(10^(findgen(69)*0.025+0.025))]
  photmag=fltarr(n_elements(radius))
  photmag_aper=fltarr(n_elements(radius))
  freemag=fltarr(n_elements(radius))
  fixmag=fltarr(n_elements(radius))
  area=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     maxflux=djs_phot(xc,yc,radius[i],0.,img,skyval=skyval,/EXACT)
     aper,img,97.5,98.7,maxflux_aper,fluxerr,0.,skyerr,1,radius[i],-1,[-100000,100000],/silent,setskyval=0.0,/flux,/exact,/NAN
     aper,img,97.5,98.7,minflux_aper,mfluxerr,0.,mskyerr,1,minradius[i],-1,[-100000,100000],/silent,setskyval=0.0,/flux,/exact,/NAN
     minflux=djs_phot(xc,yc,minradius[i],0.,img,skyval=skyval,/EXACT)
     maxflux_fix=djs_phot(xc,yc,radius[i],0.,fixmod,skyval=skyval);,/EXACT)
     minflux_fix=djs_phot(xc,yc,minradius[i],0.,fixmod,skyval=skyval);,/EXACT)
     maxflux_free=djs_phot(xc,yc,radius[i],0.,freemod,skyval=skyval);,/EXACT)
     minflux_free=djs_phot(xc,yc,minradius[i],0.,freemod,skyval=skyval);,/EXACT)
     if (minradius[i] lt 0.025) then begin
        area[i]=((!PI*(radius[i])^2))
        photflux=maxflux
        photflux_aper=maxflux_aper
        modfluxfree=maxflux_free
        modfluxfix=maxflux_fix
     endif else begin
        area[i]=((!PI*(radius[i])^2)-(!PI*(minradius[i])^2))
        photflux=maxflux-minflux
        photflux_aper=maxflux_aper-minflux_aper
        modfluxfree=maxflux_free-minflux_free
        modfluxfix=maxflux_fix-minflux_fix
     endelse
     photmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(photflux)-extinct
     freemag[i]=zeropt+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(modfluxfree)-extinct
     fixmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(modfluxfix)-extinct
     photmag_aper[i]=zeropt+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(photflux_aper)-extinct

  endfor
  
     
  !P.MULTI=[0,1,2]
  set_plot,'ps'
  device,filename='mge_sb_compare.ps',/color
  djs_plot,radius*scale,photmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',charthick=4,xthick=3,ythick=3,thick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2
  djs_oplot,radius*scale,photmag_aper,color='green',thick=5
  djs_oplot,radius*scale,sersicmag,color='blue',thick=3
  djs_oplot,radius*scale,mag,color='red',thick=3
  items=['Data','Free Unconvolved','Fixed Unconvolved']
  lines=[0,0,0]
  sym=[4,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,psym=sym,colors=color,/window,background_color='white',charthick=4,thick=3,/bottom,/left
  djs_plot,radius*scale,diff,psym=2,ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual (free-fixed)',charthick=4,xthick=3,ythick=3,charsize=1.5,ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.1,0.1]

  device,/close
  diff=freemag-fixmag
  device,filename='sb_compare_convolve.ps',/color
  djs_plot,radius*scale,photmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',charthick=4,xthick=3,ythick=3,thick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2
  djs_oplot,radius*scale,fixmag,color='blue',thick=3
  djs_oplot,radius*scale,freemag,color='red',thick=3
  items=['Data','Free Convolved','Fixed Convolved']
  lines=[0,0,0]
  sym=[4,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,psym=sym,colors=color,/window,background_color='white',charthick=4,thick=3,/bottom,/left
  djs_plot,radius*scale,diff,psym=2,ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual (free-fixed)',charthick=4,xthick=3,ythick=3,charsize=1.5,ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.1,0.1]

  device,/close

  ratio=mag/sersicmag

  device,filename='ratio_mge_sersic.ps',/color
  !P.MULTI=[0,1,1]
  djs_plot,radius*scale,ratio,psym=4,yran=[0.995,1.002],/ysty,xtitle='Radius ["]',ytitle='Free Sersic / Fixed Sersic',charthick=4,xthick=3,ythick=3,thick=3,xran=[0.,0.5],/xsty
  device,/close

  fits_read,'deconvolved_f814.fits',deconimg
  deconimg=deconimg/1050.
  find_galaxy,deconimg,majoraxis,eps,ang,xc_d,yc_d
  xc_d=122 & yc_d=121
  deconmag=fltarr(n_elements(radius))
  deconmag_aper=fltarr(n_elements(radius))
  stop
  for i=0,n_elements(radius)-1 do begin
     maxflux_decon=djs_phot(xc_d,yc_d,radius[i],0.,deconimg,skyval=skyval,/EXACT)
     aper,deconimg,xc_d,yc_d,maxflux_aper,fluxerr,0.,skyerr,1,radius[i],-1,[-100000,100000],/silent,setskyval=0.0,/flux,/exact,/NAN
     aper,deconimg,xc_d,yc_d,minflux_aper,fluxerr,0.,skyerr,1,minradius[i],-1,[-100000,100000],/silent,setskyval=0.0,/flux,/exact,/NAN
     minflux_decon=djs_phot(xc_d,yc_d,minradius[i],0.,deconimg,skyval=skyval,/EXACT)
;     stop
     if (minradius[i] lt 0.025) then begin
        photflux_decon=maxflux_decon
        deconflux_aper=maxflux_aper
     endif else begin
        photflux_decon=maxflux_decon-minflux_decon
        deconflux_aper=maxflux_aper-minflux_aper
     endelse
     deconmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(photflux_decon)-extinct
     deconmag_aper[i]=zeropt+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(deconflux_aper)-extinct
  endfor
  
  device,filename='deconvolved_compare_mge.ps',/color
  !P.MULTI=[0,1,2]
  djs_plot,radius*scale,deconmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',charthick=4,xthick=3,ythick=3,thick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Deconvolved Image vs. MGE fits'
  djs_oplot,radius*scale,sersicmag,color='blue',thick=3
  djs_oplot,radius*scale,mag,color='red',thick=3
  items=['Data Deconvolved','Free Unconvolved','Fixed Unconvolved']
  lines=[0,0,0]
  sym=[4,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,psym=sym,colors=color,/window,background_color='white',charthick=4,thick=3,/bottom,/left
  djs_plot,radius*scale,deconmag-sersicmag,psym=2,ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',charthick=4,xthick=3,ythick=3,charsize=1.5,ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.35,0.6],color='blue'
  djs_oplot,radius*scale,deconmag-mag,psym=2,color='red'

  device,/close



  zeropoint_i=25.28697
  bn1=1.999*3.25-0.327
  bn2=1.999*1.74-0.327
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
  djs_plot,radius*scale,deconmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',charthick=4,xthick=3,ythick=3,thick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Deconvolved Image vs. Sersic fits'
  djs_oplot,radius*scale,mutot,color='red',thick=3
  djs_oplot,radius*scale,mutot_fix,color='blue',thick=3
  items=['Data Deconvolved','Free Unconvolved','Fixed Unconvolved']
  lines=[0,0,0]
  sym=[4,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,psym=sym,colors=color,/window,background_color='white',charthick=4,thick=3,/bottom,/left
  djs_plot,radius*scale,deconmag-mutot_fix,psym=2,ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',charthick=4,xthick=3,ythick=3,charsize=1.5,ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.75,0.2],color='blue'
  djs_oplot,radius*scale,deconmag-mutot,psym=2,color='red'

  device,/close
  set_plot,'x'
  !P.MULTI=[0,1,1]
  stop

;  device,filename='mgevsersic.ps',/color
;  djs_plot,radius,sersicmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',charthick=4,xthick=3,ythick=3,thick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95];,xran=[0,1.5],/xsty
;  djs_oplot,radius,mutot,color='blue',thick=3
;  djs_plot,radius,sersicmag-mutot,psym=2,ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',charthick=4,xthick=3,ythick=3,charsize=1.5,ycharsize=0.6;,xran=[0,1.5],/xsty
;  device,/close
;  set_plot,'x'
;  stop

END
