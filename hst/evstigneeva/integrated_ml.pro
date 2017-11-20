pro integrated_ml
  scale=0.025
  readcol,'vucd3_mge_outputsersic.dat',lum,sig,q,pa,format='F,F,F,F'
  readcol,'vucd3_mge_outputsersic_mass.dat',m,format='F'
  radius=[0,(10^(findgen(80)*0.025))*scale]
  luminosity=fltarr(n_elements(radius))
  mass=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     temp=0.
     tmp=0.
     for j=0,n_elements(lum)-1 do begin
        temp+= (lum[j]*exp(-(1./(2*sig[j]^2))*(radius[i])^2))
     
        tmp+= (m[j]*exp(-(1./(2*sig[j]^2))*(radius[i])^2))
     
     endfor
     luminosity[i]=temp
     mass[i]=tmp
  endfor
  djs_plot,radius,luminosity,xran=[0,0.5]
  djs_oplot,radius,mass,color='blue'
  print,(total(mass))/(total(luminosity))
  stop
  
END
