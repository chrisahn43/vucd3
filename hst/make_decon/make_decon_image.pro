@sersic2image
@mge2image
PRO MAKE_DECON_IMAGE
  zeropoint=25.28697d
  zeropointr=25.9799d
  scale=0.025d
  Msun=4.53d
  Msunr=4.73d
  ;READ IN ORIGINAL IMAGES
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/galfit_model.fits',img,exten_no=1
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/galfit_model_r.fits',imgr,exten_no=1
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/deconvolved_f814.fits',oldimg,head
  hrotate,oldimg,head,img,newhead,1
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/deconvolved_f606.fits',oldimgr,headr
  hrotate,oldimgr,headr,imgr,newheadr,1
  imsizei=size(img,/dim)
  imsizer=size(imgr,/dim)
                                ;FIND THEIR CENTERS
  find_galaxy,img,majoraxis,eps,ang,xci,yci
  find_galaxy,imgr,majoraxis,eps,ang,xcr,ycr
                                ;CREATE X AND Y ARRAYS FOR USE IN
                                ;MODELS

                                ;CREATE IBAND FIXED MGE MODEL IMAGE
  readcol,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/vucd3_mge_outputsersic.dat',fixlum,fixsig,fixq,fixpa,format='D,D,D,D'
  mge2image,img,xci,yci,fixlum,fixsig,fixq,fixpa,model,zeropoint=zeropoint,scale=scale,msun=msun,fits='mgefixmodeli.fits'

                                ;CREATE IBAND FIXED SERSIC MODEL IMAGE
  mtot=[18.43,18.14]
  re=[3.16,24.47]
  n=[3.51,1.28]
  q=[0.66,0.91]
  pa=[19.04,18.43]
  sersic2image,img,xci,yci,mtot,re,n,q,pa,model,zeropoint=zeropoint,fits='serfixmodeli.fits'

                                ;CREATE IBAND FREE MGE MODEL IMAGE
  readcol,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/vucd3_mge_outputsersic_free.dat',freelum,freesig,freeq,freepa,format='D,D,D,D'
  mge2image,img,xci,yci,freelum,freesig,freeq,freepa,model,zeropoint=zeropoint,scale=scale,msun=msun,fits='mgefreemodeli.fits'

                                ;CREATE IBAND FREE SERSIC MODEL IMAGE
  mtot=[18.58,17.99]
  re=[2.75,24.71]
  n=[3.25,1.74]
  q=[0.62,0.89]
  pa=[17.97,20.65]
  sersic2image,img,xci,yci,mtot,re,n,q,pa,model,zeropoint=zeropoint,fits='serfreemodeli.fits'


                                ;CREATE VBAND FREE MGE MODEL
  
  readcol,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/vucd3_mge_outputsersic_free_606.dat',freelum,freesig,freeq,freepa,format='D,D,D,D'
  mge2image,imgr,xcr,ycr,freelum,freesig,freeq,freepa,model,zeropoint=zeropointr,scale=scale,msun=msunr,fits='mgefreemodelv.fits'

                                ;CREATE VBAND FREE SERSIC MODEL
  mtot=[19.05,18.89]
  re=[3.16,24.47]
  n=[3.51,1.28]
  q=[0.66,0.91]
  pa=[19.04,18.43]
  sersic2image,imgr,xcr,ycr,mtot,re,n,q,pa,model,zeropoint=zeropointr,fits='serfreemodelv.fits'

                                ;CREATE VBAND FIXED MGE MODEL
  
  readcol,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/vucd3_mge_outputsersic_606.dat',fixlum,fixsig,fixq,fixpa,format='D,D,D,D'
  mge2image,imgr,xcr,ycr,fixlum,fixsig,fixq,fixpa,model,zeropoint=zeropointr,scale=scale,msun=msunr,fits='mgefixmodelv.fits'

                                ;CREATE VBAND FIXED SERSIC MODEL
  
  mtot=[19.19,18.73]
  re=[2.75,24.71]
  n=[3.25,1.74]
  q=[0.62,0.89]
  pa=[17.97,20.65]
  sersic2image,imgr,xcr,ycr,mtot,re,n,q,pa,model,zeropoint=zeropointr,fits='serfixmodelv.fits'

  stop
  
;  makex,img,xarr,yarr,/zero
;  xi=xarr-xci
;  yi=yarr-yci
;  xi=[-1*(reverse(findgen(xci))+1),findgen(imsizei[0]-xci)]
;  yi=[-1*(reverse(findgen(yci))+1),findgen(imsizei[1]-yci)]
;  xr=[reverse(findgen(xcr))+1,findgen(imsizer[0]-xcr)]
;  yr=[reverse(findgen(ycr))+1,findgen(imsizer[1]-ycr)]
;  stop
  ;CREATE I BAND FIXED AND FREE MODEL IMAGES USING MGE FITS TO SERSIC
;  zeropt=25.28697d
;  extinct=0.034
;  scale=0.025d
;  const=(64800./!PI)^2
;  Msun=4.53d
;  readcol,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/vucd3_mge_outputsersic.dat',fixlum,fixsig,fixq,fixpa,format='D,D,D,D'
;  stop
;  fixsb=Msun-2.5*alog10(fixlum/const)
;  fixpeak=10^(-0.4*(fixsb-zeropt-5*alog10(scale)));+extinct))
;  fixsig=fixsig/scale
;  fixpa=fixpa+90.


;  fixpa=fixpa*!DtoR
;  fixmodeli=fltarr(imsizei[0],imsizei[1])

;  a=where(fixq lt 0.8)
;  inang=fixpa[a[0]]
;  xin=(xi*COS(inang))+(yi*SIN(inang))
;  yin=(-(xi*SIN(inang)))+(yi*COS(inang))
;  b=where(fixq gt 0.8)
;  outang=fixpa[b[0]]

;  xout=(xi*COS(outang))+(yi*SIN(outang))
;  yout=(-(xi*SIN(outang)))+(yi*COS(outang))

;  stop
  
;  for i=0,imsizei[0]-1 do begin
;     for j=0,imsizei[1]-1 do begin
;        in=0.
;        out=0.
;        for k=0,n_elements(fixpeak)-1 do begin
;           if (fixq[k] lt 0.8) then begin
;              in+= (fixpeak[k]*exp(-(1./(2*fixsig[k]^2))*(xin[i,j]^2+(yin[i,j]^2/fixq[k]^2))))

;           endif else begin
;              out+=(fixpeak[k]*exp(-(1./(2*fixsig[k]^2))*(xout[i,j]^2+(yout[i,j]^2/fixq[k]^2))))

 ;          endelse
 ;       endfor
 ;       fixmodeli[i,j]=in+out
;     endfor
;  endfor
;  writefits,'mgefixmodeli.fits',fixmodeli
  ;stop

                                ;MAKE SERSIC FIXED MODEL
;  fixmodelseri=fltarr(imsizei[0],imsizei[1])
;  bn1=1.999*3.51-0.327
;  bn2=1.999*1.28-0.327
;  mue1=18.43+5*alog10(3.16)+2.5*alog10(2*!PI*3.51*0.66*(exp(bn1)/((bn1)^(2*3.51)))*GAMMA(2*3.51))
;  mue2=18.14+5*alog10(24.47)+2.5*alog10(2*!PI*1.28*0.91*(exp(bn2)/((bn2)^(2*1.28)))*GAMMA(2*1.28))
;  for i=0,imsizei[0]-1 do begin
;     for j=0,imsizei[1]-1 do begin
;        inrad=sqrt(xin[i,j]^2+(yin[i,j]^2/0.66^2))
;        outrad=sqrt(xout[i,j]^2+(yout[i,j]^2/0.91^2))
 ;       mu1= mue1+((2.5*bn1)/(alog(10)))*((((inrad)/(3.16))^(1./3.51))-1)
;        mu2= mue2+((2.5*bn2)/(alog(10)))*((((outrad)/(24.47))^(1./1.28))-1)
;        int1=(10^(-0.4*(mu1-zeropt)));*(!PI*inrad^2)
;        int2=(10^(-0.4*(mu2-zeropt)));*(!PI*outrad^2)
;        fixmodelseri[i,j]=int1+int2
;     endfor
;  endfor

;  writefits,'serfixmodeli.fits',fixmodelseri
;  stop

;  readcol,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/vucd3_mge_outputsersic_free.dat',freelum,freesig,freeq,freepa,format='D,D,D,D'
;  freesb=Msun-2.5*alog10(freelum/const)
;  freepeak=10^(-0.4*(freesb-zeropt-5*alog10(scale)));+extinct))
;  freesig=freesig/scale
;  freepa=90.+freepa

;  freepa=freepa*!DtoR
;  freemodeli=fltarr(imsizei[0],imsizei[1])

;  a=where(freeq lt 0.8)
;  inang=freepa[a[0]]
;  xin=(xi*COS(inang))+(yi*SIN(inang))
;  yin=(-(xi*SIN(inang)))+(yi*COS(inang))
;  b=where(freeq gt 0.8)
;  outang=freepa[b[0]]
;  xout=xi*COS(outang)-yi*SIN(outang)
;  yout=xi*SIN(outang)+yi*COS(outang)

;  for i=0,imsizei[0]-1 do begin
;     for j=0,imsizei[1]-1 do begin
;        in=0.
;        out=0.
;        for k=0,n_elements(freepeak)-1 do begin
;           if (freeq[k] lt 0.8) then begin
;              in+= (freepeak[k]*exp(-(1./(2*freesig[k]^2))*(xin[i,j]^2+(yin[i,j]^2/freeq[k]^2))))

;           endif else begin
;              out+=(freepeak[k]*exp(-(1./(2*freesig[k]^2))*(xout[i,j]^2+(yout[i,j]^2/freeq[k]^2))))

;           endelse
;        endfor
;        freemodeli[i,j]=in+out
;     endfor
;  endfor
;  writefits,'mgefreemodeli.fits',freemodeli
  ;stop
                                ;MAKE SERSIC FREE MODEL
;  freemodelseri=fltarr(imsizei[0],imsizei[1])
;  bn1=1.999*3.25-0.327
;  bn2=1.999*1.74-0.327
;  mue1=18.58+5*alog10(2.75)+2.5*alog10(2*!PI*3.25*0.62*(exp(bn1)/((bn1)^(2*3.25)))*GAMMA(2*3.25))
;  mue2=17.99+5*alog10(24.71)+2.5*alog10(2*!PI*1.74*0.89*(exp(bn2)/((bn2)^(2*1.74)))*GAMMA(2*1.74))
;  for i=0,imsizei[0]-1 do begin
;     for j=0,imsizei[1]-1 do begin
;        inrad=sqrt(xin[i,j]^2+(yin[i,j]^2/0.62^2))
;        outrad=sqrt(xout[i,j]^2+(yout[i,j]^2/0.89^2))
;        mu1= mue1+((2.5*bn1)/(alog(10)))*((((inrad)/(2.75))^(1./3.25))-1)
;        mu2= mue2+((2.5*bn2)/(alog(10)))*((((outrad)/(24.71))^(1./1.74))-1)
;        int1=(10^(-0.4*(mu1-zeropt)));*(!PI*inrad^2)
;        int2=(10^(-0.4*(mu2-zeropt)));*(!PI*outrad^2)
;        freemodelseri[i,j]=int1+int2
;     endfor
;  endfor

 ; writefits,'serfreemodeli.fits',freemodelseri
  ;stop

END
