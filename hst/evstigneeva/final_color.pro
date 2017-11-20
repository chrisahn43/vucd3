pro final_color
  scale=0.025
  rad=findgen(60)+1
  minrad=findgen(60)
  zeropoint_r=25.9799
  zeropoint_i=25.28697
  fits_read,'galfit_model.fits',iimg,exten_no=1
  fits_read,'galfit_model_r.fits',rimg,exten_no=1
  fits_read,'galfit_model.fits',ifree,exten_no=2
  fits_read,'galfit_model_fixed.fits',ifixed,exten_no=2
  fits_read,'galfit_model_r.fits',rfree,exten_no=2
  fits_read,'galfit_model_rfixed.fits',rfixed,exten_no=2
  xc=98 & yc=99
  idata=fltarr(n_elements(rad))
  rdata=fltarr(n_elements(rad))
  imagfixed=fltarr(n_elements(rad))
  rmagfixed=fltarr(n_elements(rad))
  imagfree=fltarr(n_elements(rad))
  rmagfree=fltarr(n_elements(rad))
  for i=0,n_elements(rad)-1 do begin
     area=((!PI*(rad[i])^2)-(!PI*(minrad[i])^2))

     maxidata=djs_phot(xc,yc,rad[i],0.,iimg,skyval=skyval)
     minidata=djs_phot(xc,yc,minrad[i],0.,iimg,skyval=skyval)
     iflux=maxidata-minidata
     idata[i]=zeropoint_i + 5*alog10(scale) + 2.5*alog10(area) - 2.5*alog10(iflux)
     maxrdata=djs_phot(xc,yc,rad[i],0.,rimg,skyval=skyval)
     minrdata=djs_phot(xc,yc,minrad[i],0.,rimg,skyval=skyval)
     rflux=maxrdata-minrdata
     rdata[i]=zeropoint_r + 5*alog10(scale) + 2.5*alog10(area) - 2.5*alog10(rflux)
     
     maxifree=djs_phot(xc,yc,rad[i],0.,ifree,skyval=skyval)
     minifree=djs_phot(xc,yc,minrad[i],0.,ifree,skyval=skyval)
     ifluxfree=maxifree-minifree
     imagfree[i]=zeropoint_i + 5*alog10(scale) + 2.5*alog10(area) - 2.5*alog10(ifluxfree)

     maxrfree=djs_phot(xc,yc,rad[i],0.,rfree,skyval=skyval)
     minrfree=djs_phot(xc,yc,minrad[i],0.,rfree,skyval=skyval)
     rfluxfree=maxrfree-minrfree
     rmagfree[i]=zeropoint_r + 5*alog10(scale) + 2.5*alog10(area) - 2.5*alog10(rfluxfree)

     maxifixed=djs_phot(xc,yc,rad[i],0.,ifixed,skyval=skyval)
     minifixed=djs_phot(xc,yc,minrad[i],0.,ifixed,skyval=skyval)
     ifluxfixed=maxifixed-minifixed
     imagfixed[i]=zeropoint_i + 5*alog10(scale) + 2.5*alog10(area) - 2.5*alog10(ifluxfixed)

     maxrfixed=djs_phot(xc,yc,rad[i],0.,rfixed,skyval=skyval)
     minrfixed=djs_phot(xc,yc,minrad[i],0.,rfixed,skyval=skyval)
     rfluxfixed=maxrfixed-minrfixed
    rmagfixed[i]=zeropoint_r + 5*alog10(scale) + 2.5*alog10(area) - 2.5*alog10(rfluxfixed)
  
     

 endfor
  colordata=rdata-idata
  colorfree=rmagfree-imagfree
  rcolorfixed=rmagfixed-imagfree
  icolorfixed=rmagfree-imagfixed
  
;ORIGINAL OUTPUTS
;  sersic    : (  629.05,   661.21)   18.58      2.75    3.25    0.62    17.97
;  sersic    : ( {629.05}, {661.21})  17.99     24.71    1.74    0.89    20.65
  bn1=1.999*3.25-0.327
  bn2=1.999*1.74-0.327
  mue1=18.58+5*alog10(2.75*scale)+2.5*alog10(2*!PI*3.25*(exp(bn1)/((bn1)^(2*3.25)))*GAMMA(2*3.25))
;  mue1=18.58+2.5*alog10(2*!PI*(2.75*scale)^2)
;  mue1=mue1+1.25
  mue2=17.99+5*alog10(24.71*scale)+2.5*alog10(2*!PI*1.74*(exp(bn2)/((bn2)^(2*1.74)))*GAMMA(2*1.74))
;  mue2=17.99+2.5*alog10(2*!PI*(24.71*scale)^2)
;  mue2=mue2+1.;9
  mu1=mue1+((2.5*bn1)/(alog(10)))*((((rad)/(2.75))^(1./3.25))-1)
  mu2=mue2+((2.5*bn2)/(alog(10)))*((((rad)/(24.71))^(1./1.74))-1)
  int1=(10^(-0.4*(mu1-zeropoint_i)))
  int2=(10^(-0.4*(mu2-zeropoint_i)))
  inttot=int1+int2
  mutot=zeropoint_i-2.5*alog10(inttot)
;  sersic    : (  629.30,   661.02)   19.05      3.16    3.51    0.66    19.04
;  sersic    : ( {629.30}, {661.02})  18.89     24.47    1.28    0.91    18.43
  bn1606=1.999*3.51-0.327
  bn2606=1.999*1.28-0.327
  mue1606=19.05+5*alog10(3.16*scale)+2.5*alog10(2*!PI*3.51*(exp(bn1606)/((bn1606)^(2*3.51)))*GAMMA(2*3.51))
;  mue1606=19.05+2.5*alog10(2*!PI*(3.16*scale)^2)
;  mue1606=mue1606+1.3
  mue2606=18.89+5*alog10(24.47*scale)+2.5*alog10(2*!PI*1.28*(exp(bn2606)/((bn2606)^(2*1.28)))*GAMMA(2*1.28))
;  mue2606=18.89+2.5*alog10(2*!PI*(24.47*scale)^2)
;  mue2606=mue2606+0.8
  mu1606=mue1606+((2.5*bn1606)/(alog(10)))*((((rad)/(3.16))^(1./3.51))-1)
  mu2606=mue2606+((2.5*bn2606)/(alog(10)))*((((rad)/(24.47))^(1./1.28))-1)
  int1606=(10^(-0.4*(mu1606-zeropoint_r)))
  int2606=(10^(-0.4*(mu2606-zeropoint_r)))
  inttot606=int1606+int2606
  mutot606=zeropoint_r-2.5*alog10(inttot606)
  colormu_orig=mutot606-mutot
  set_plot,'ps'
  device,filename='color_all.ps',/color
  djs_plot,rad*scale,colormu_orig,ytitle='(\mu_{F606W} - \mu_{F814W}) [mag/arcsec^2]',xtitle='Radius ["]',yran=[0.3,1.],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,thick=3,xstyle=8,ymargin=[4,4],linestyle=2
  djs_oplot,rad*scale,colorfree,thick=3
  djs_oplot,rad*scale,colordata,thick=3,psym=4
;  stop
;Fixed F606 OUTPUTS:
;   sersic    : (  629.30,   661.02)   19.19     [2.75]  [3.25]  [0.62]  [17.97]
;   sersic    : ( {629.30}, {661.02})  18.73    [24.71]  [1.74]  [0.89]  [20.65]
  bn1606=1.999*3.25-0.327
  bn2606=1.999*1.74-0.327
  mue1606=19.19+5*alog10(2.75*scale)+2.5*alog10(2*!PI*3.25*(exp(bn1606)/((bn1606)^(2*3.25)))*GAMMA(2*3.25))
;  mue1606=19.19+2.5*alog10(2*!PI*(2.75*scale)^2)
;  mue1606=mue1606+1.25
  mue2606=18.73+5*alog10(24.71*scale)+2.5*alog10(2*!PI*1.74*(exp(bn2606)/((bn2606)^(2*1.74)))*GAMMA(2*1.74))
;  mue2606=18.73+2.5*alog10(2*!PI*(24.71*scale)^2)
;  mue2606=mue2606+1.
  mu1606=mue1606+((2.5*bn1606)/(alog(10)))*((((rad)/(2.75))^(1./3.25))-1)
  mu2606=mue2606+((2.5*bn2606)/(alog(10)))*((((rad)/(24.71))^(1./1.74))-1)
  int1606=(10^(-0.4*(mu1606-zeropoint_r)))
  int2606=(10^(-0.4*(mu2606-zeropoint_r)))
  inttot606=int1606+int2606
  mutot606_fixed=zeropoint_r-2.5*alog10(inttot606)
  colormu_fixed=mutot606_fixed-mutot
  djs_oplot,rad*scale,colormu_fixed,color='blue',thick=3,linestyle=2
  djs_oplot,rad*scale,rcolorfixed,thick=3,color='blue'
;Fixed F814 OUTPUTS:
; sersic    : (  629.05,   661.21)   18.43     [3.16]  [3.51]  [0.66]  [19.04]
; sersic    : ( {629.05}, {661.21})  18.14    [24.47]  [1.28]  [0.91]  [18.43]
  bn1=1.999*3.51-0.327
  bn2=1.999*1.28-0.327
  mue1=18.43+5*alog10(3.16*scale)+2.5*alog10(2*!PI*3.51*(exp(bn1)/((bn1)^(2*3.51)))*GAMMA(2*3.51))
;  mue1=18.43+2.5*alog10(2*!PI*(3.16*scale)^2)
;  mue1=mue1+1.3
  mue2=18.14+5*alog10(24.47*scale)+2.5*alog10(2*!PI*1.28*(exp(bn2)/((bn2)^(2*1.28)))*GAMMA(2*1.28))
;  mue2=18.14+2.5*alog10(2*!PI*(24.47*scale)^2)
;  mue2=mue2+0.8
  mu1=mue1+((2.5*bn1)/(alog(10)))*((((rad)/(3.16))^(1./3.51))-1)
  mu2=mue2+((2.5*bn2)/(alog(10)))*((((rad)/(24.47))^(1./1.28))-1)
  int1=(10^(-0.4*(mu1-zeropoint_i)))
  int2=(10^(-0.4*(mu2-zeropoint_i)))
  inttot=int1+int2
  mutot_fixed=zeropoint_i-2.5*alog10(inttot)
  colormu_fixed_i=mutot606-mutot_fixed
  djs_oplot,rad*scale,colormu_fixed_i,color='red',thick=3,linestyle=2
  djs_oplot,rad*scale,icolorfixed,thick=3,color='red'
  
  axis,xaxis=1,xtickv=rad,xcharsize=1.5,charthick=4,xthick=3,xtitle='Radius [Pixels]',/xsty,xran=[1,60]
  items=['Data','Convolved Model','Unconvolved Model']
  lines=[0,0,2]
  sym=[4,0,0]
  al_legend,items,linestyle=lines,psym=sym,/window,background_color='white',charthick=4,thick=3,/bottom,/right;,charsize=1.5
;  stop
                                ;No FOCUS outputs:
;  sersic    : (  629.05,   661.21)   18.59      2.73    3.24    0.62    18.01
;  sersic    : ( {629.05}, {661.21})  17.98     24.66    1.75    0.89    20.61
  bn1=1.999*3.24-0.327
  bn2=1.999*1.75-0.327
  mue1=18.59+5*alog10(2.73*scale)+2.5*alog10(2*!PI*3.24*(exp(bn1)/((bn1)^(2*3.24)))*GAMMA(2*3.24))
;  mue1=18.59+2.5*alog10(2*!PI*(2.73*scale)^2)
;  mue1=mue1+1.25
  mue2=17.98+5*alog10(24.66*scale)+2.5*alog10(2*!PI*1.75*(exp(bn2)/((bn2)^(2*1.75)))*GAMMA(2*1.75))
;  mue2=17.98+2.5*alog10(2*!PI*(24.66*scale)^2)
;  mue2=mue2+1.
  mu1=mue1+((2.5*bn1)/(alog(10)))*((((rad)/(2.73))^(1./3.24))-1)
  mu2=mue2+((2.5*bn2)/(alog(10)))*((((rad)/(24.66))^(1./1.75))-1)
  int1=(10^(-0.4*(mu1-zeropoint_i)))
  int2=(10^(-0.4*(mu2-zeropoint_i)))
  inttot=int1+int2
  mutot=zeropoint_i-2.5*alog10(inttot)
; sersic    : (  629.29,   661.02)   19.05      3.15    3.52    0.66    19.07
; sersic    : ( {629.29}, {661.02})  18.89     24.43    1.28    0.91    18.40
  bn1606=1.999*3.52-0.327
  bn2606=1.999*1.28-0.327
  mue1606=19.05+5*alog10(3.15*scale)+2.5*alog10(2*!PI*3.52*(exp(bn1606)/((bn1606)^(2*3.52)))*GAMMA(2*3.52))
;  mue1606=19.05+2.5*alog10(2*!PI*(3.15*scale)^2)
;  mue1606=mue1606+1.3
  mue2606=18.89+5*alog10(24.43*scale)+2.5*alog10(2*!PI*1.28*(exp(bn2606)/((bn2606)^(2*1.28)))*GAMMA(2*1.28))
;  mue2606=18.89+2.5*alog10(2*!PI*(24.43*scale)^2)
;  mue2606=mue2606+0.8
  mu1606=mue1606+((2.5*bn1606)/(alog(10)))*((((rad)/(3.15))^(1./3.52))-1)
  mu2606=mue2606+((2.5*bn2606)/(alog(10)))*((((rad)/(24.43))^(1./1.28))-1)
  int1606=(10^(-0.4*(mu1606-zeropoint_r)))
  int2606=(10^(-0.4*(mu2606-zeropoint_r)))
  inttot606=int1606+int2606
  mutot606=zeropoint_r-2.5*alog10(inttot606)
  colormu_orig_nofocus=mutot606-mutot

;  djs_oplot,rad*scale,colormu_orig_nofocus,linestyle=2,thick=3

  device,/close
  set_plot,'x'
  stop

END
