PRO FITPSFNEW
                                ;THIS CODE IS USED TO FIND THE
                                ;KINEMATIC PSF FOR A DOUBLE GAUSSIAN
                                ;AND A GAUSSIAN AND MOFFAT FUNCTION
                                ;PROVIDED BY ANIL. SECOND TIME USING
                                ;KBAND IMAGE
;Gauss+gauss
;Iter      8   CHI-SQUARE =       972.35055          DOF = 523
;    P(0) =              3.23933  Inner gaussian width=0.162"
;    P(1) =              131.228 junk
;    P(2) =             0.507517 junk
;    P(3) =             0.425901 junk
;    P(4) =              13.9387 outer guassian width (0.697”)
;    P(5) =              1.23660 power (f_outer=1.236/(1+1.236))

;Gauss+Moffat
;Iter      4   CHI-SQUARE =       1032.9106          DOF = 524
;    P(0) =              2.75179  Inner Gaussian width=0.138"
;    P(1) =              129.819 normalization (junk)
;    P(2) =             0.510709 dx (junk)
;    P(3) =             0.418685 dy (junk)
;    P(4) =              21.6700 Moffat width
;    P(5) =              2.41817 Power (f_moffat=(2.418/(1+2.418)))
;    P(6) =              4.76500

;Gauss + Moffat function definition
;rpsf=SQRT((xpsf-psfcenter-p[2]*scaling)^2+(ypsf-psfcenter-p[3]*scaling)^2)
;psfwidth1=p[0]/2.35*scaling ;in pixels
;psfwidth2=p[4]/2.35*scaling ;in pixels
;psf1=exp(-(rpsf/psfwidth1)^2/2.)/(2.*!PI*psfwidth1^2)
;psf1=psf1/TOTAL(psf1)
;psf2=1./(1+rpsf/psfwidth2)^p[6]
;psf2=p[5]*psf2/TOTAL(psf2)
;psf=psf1+psf2
;print,TOTAL(psf)
;psf=psf/TOTAL(psf)
;psf=p[1]*psf  
  ;GAUSS + MOFFAT PSF
  scaling=0.05
  p=[2.75179,128.819,0.510709,0.418685,21.670,2.41817,4.76500]
  r=range(0.01,2.5,100,/log)
  psfwidth1=0.137589/2.35
  psfwidth2=p[4]/2.35*scaling
  psf1=exp(-(r/psfwidth1)^2/2.)/(2.*!PI*psfwidth1^2)
  psf1=psf1/TOTAL(psf1)
  psf2=1./(1+r/psfwidth2)^p[6]
  psf2=p[5]*psf2/total(psf2)
  psf=psf1+psf2
  psf=psf/total(psf)
  psf=p[1]*psf

  mge_fit_1d,r,psf2,ngauss=10,sol=sol
  forprint,sol[0,*],sol[1,*],format='F,F',textout='temp.dat'
  readcol,'temp.dat',modelweight,modelsig,format='F,F'
  spawn, 'rm -rf temp.dat'
  stop
;  modelweight=modelweight/(sqrt(2*!PI)*modelsig)
  modelflux=fltarr(n_elements(r))
  modelpeak=((total(psf2)/(1+total(psf2))))*(modelweight)/(total(modelweight))
  modelpeak=[(1-total(modelpeak)),modelpeak]
  modelsig=[psfwidth1,modelsig]
  forprint,modelpeak,modelsig,format='F,F',textout='newkinematic_psf_moffat.dat'
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

  p=[3.23933,131.228,0.507517,0.425901,13.9387,1.23660]
  psfwidth1=0.162/2.35
  psfwidth2=0.697/2.35
  psf1=exp(-(r/psfwidth1)^2/2.)/(2.*!PI*psfwidth1^2)
  psf1=psf1/TOTAL(psf1)
  psf2=exp(-(r/psfwidth2)^2/2.)/(2.*!PI*psfwidth2^2)
  psf2=p[5]*psf2/total(psf2)
  psf=psf1+psf2
  psf=psf/total(psf)
  lightfrac=p[5]/(1+p[5])
  modelpeak=[(1-lightfrac),lightfrac]
  modelsig=[psfwidth1,psfwidth2]
  forprint,modelpeak,modelsig,format='F,F',textout='newkinematic_psf_gauss.dat'

  stop
END
