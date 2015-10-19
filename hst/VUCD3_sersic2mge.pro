@sersic2mge.pro
pro VUCD3_sersic2mge
  scale=0.025
  re1=5.42d*scale
  re2=32.36d*scale
  sersic1={mag:18.19d,re:re1,n:4.5d,pa:19.04d,q:0.71d}
  sersic2={mag:18.5d,re:re2,n:1.03d,pa:30.5d,q:0.97d}
  sersics=[sersic1,sersic2]
  Msun=4.10d
  A_B=0.034d
  zeropt=25.28697
  extinct=0.034
  mge=sersics2mge(sersics,Msun,A_B=A_B)
  forprint,mge[*,0],mge[*,1],mge[*,2],mge[*,3],textout='vucd3_mge_outputsersic.dat',format='F,F,F,F'
  stop
  const=(64800./!PI)^2
  readcol,'vucd3_mge_outputsersic.dat',sersiclum,sersicsig,sersicq,sersicpa,format='D,D,D,D'
  sersicsb=Msun-2.5*alog10(sersiclum/const)
  sersicpeak=10^(-0.4*(sersicsb-zeropt-5*alog10(scale)+extinct))
  readcol,'vucd3_mge_output.dat',lum,sig,q,format='F,F,F'
  sb=Msun-2.5*alog10(lum/const)
  peak=10^(-0.4*(sb-zeropt-5*alog10(scale)+extinct))
  radius=(10^(findgen(80)*0.025+0.025))*scale
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
  !P.MULTI=[0,1,2]
  set_plot,'ps'
  device,filename='mge_sb_compare.ps',/color
  djs_plot,radius,mag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',charthick=4,xthick=3,ythick=3,thick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty
  djs_oplot,radius,sersicmag,color='blue',thick=3
  djs_plot,radius,diff,psym=2,ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',charthick=4,xthick=3,ythick=3,charsize=1.5,ycharsize=0.6,xran=[0,1.5],/xsty
  device,/close
  mue1=18.19+2.5*alog10(2*!PI*(5.42*scale)^2)
  mue1=mue1+1.3
  mue2=18.5+2.5*alog10(2*!PI*(32.36*scale)^2)
  mue2=mue2+0.5
  bn1=1.999*4.5-0.327
  bn2=1.999*1.03-0.327
  mu1=mue1+((2.5*bn1)/(alog(10)))*((((radius)/(5.42*scale))^(1./4.5))-1)
  mu2=mue2+((2.5*bn2)/(alog(10)))*((((radius)/(32.36*scale))^(1./1.03))-1)
  int1=(10^(-0.4*(mu1-zeropt)))
  int2=(10^(-0.4*(mu2-zeropt)))
  inttot=int1+int2
  mutot=zeropt-2.5*alog10(inttot)
  
  device,filename='mgevsersic.ps',/color
  djs_plot,radius,sersicmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',charthick=4,xthick=3,ythick=3,thick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95];,xran=[0,1.5],/xsty
  djs_oplot,radius,mutot,color='blue',thick=3
  djs_plot,radius,sersicmag-mutot,psym=2,ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',charthick=4,xthick=3,ythick=3,charsize=1.5,ycharsize=0.6;,xran=[0,1.5],/xsty
  device,/close
  set_plot,'x'
  stop

END
