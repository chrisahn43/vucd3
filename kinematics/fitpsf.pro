pro fitpsf
;this code is used to find the Kinematic PSF from the gaussian and moffat
;functions provided by ANIL
;  P(0) =              1.79948 ; width not taking into account HST PSF, including HST PSF, the width is 0.11” for the inner component
;  P(1) =              21.5393 overall scaling
;  P(2) =             0.463513 dx
;  P(3) =             0.442229 dy
;  P(4) =              21.6700 FWHM of the Moffat in pixels
;  P(5) =             0.810198 Light fraction Moffat (f_moffat=(1.245/(1+1.245)))
;  P(6) =              4.76500 Moffat power law index
;  rpsf=SQRT((xpsf-psfcenter-p[2]*scaling)^2+(ypsf-psfcenter-p[3]*scaling)^2)
;  psfwidth1=p[0]/2.35*scaling   ;in pixels
;  psfwidth2=p[4]/2.35*scaling   ;in pixels
;  psf1=exp(-(rpsf/psfwidth1)^2/2.)/(2.*!PI*psfwidth1^2)
;  psf1=psf1/TOTAL(psf1)
;  psf2=1./(1+rpsf/psfwidth2)^p[6]
;  psf2=p[5]*psf2/TOTAL(psf2)
;  psf=psf1+psf2
;print,TOTAL(psf)
;  psf=psf/TOTAL(psf)
;  psf=p[1]*psf
  scaling=0.05
  p=[1.79948,21.5393,0.463513,0.442229,21.670,0.810198,4.76500]
  r=range(0.01,2.5,100,/log)
  psfwidth1=0.11/2.35;p[0]/2.35*scaling
  psfwidth2=p[4]/2.35*scaling
  psf1=exp(-(r/psfwidth1)^2/2.)/(2.*!PI*psfwidth1^2)
  psf1=psf1/TOTAL(psf1)
  psf2=1./(1+r/psfwidth2)^p[6]
  psf2=p[5]*psf2/total(psf2)
  psf=psf1+psf2
  psf=psf/total(psf)
  psf=p[1]*psf

 ; mge_fit_1d,r,psf,ngauss=6,sol=sol
  mge_fit_1d,r,psf2,ngauss=10,sol=sol
  forprint,sol[0,*],sol[1,*],format='F,F',textout='temp.dat'
  readcol,'temp.dat',modelweight,modelsig,format='F,F'
  spawn, 'rm -rf temp.dat'
;  stop
  ;modelweight=(sol[0,*])        ;/(total(sol[0,*]))

;  modelweight=modelweight/(sqrt(2*!PI)*modelsig)
  modelsig=[0.0468085,modelsig]
  modelflux=fltarr(n_elements(r))
  ;modelpeak=(sol[0,*])/(total(sol[0,*]))
  modelpeak=((total(psf2)/(1+total(psf2))))*(modelweight)/(total(modelweight))
  modelpeak=[(1-total(modelpeak)),modelpeak]
  stop
  forprint,modelpeak,modelsig,format='F,F',textout='kinematic_psf.dat'
  modelsig=sol[1,*]
  for i=0,n_elements(r)-1 do begin
     temp=0.
     for j=0,n_elements(modelweight)-1 do begin
        temp+= (modelweight[j]*exp(-(1./(2*modelsig[j]^2))*(r[i])^2))
        modelflux[i]=temp
     endfor
  endfor

  djs_plot,r,psf2,psym=2
  djs_oplot,r,modelflux,color='blue'
  print,mean(((psf2-modelflux)/psf2)*100)
  stop
END
