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
  radius=findgen(120)*0.025+0.025
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
  set_plot,'ps'
  device,filename='mge_sb_compare.ps',/color
  djs_plot,radius,mag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius,sersicmag,color='blue',thick=3
  device,/close
  set_plot,'x'
  stop
END
