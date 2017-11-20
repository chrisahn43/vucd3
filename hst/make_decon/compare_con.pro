PRO COMPARE_CON
  scale=0.025
  zeropt=25.28697
  extinct=0.034
  fits_read,'mgefixmodeli.fits',fixmodeli
  fits_read,'serfixmodeli.fits',fixmodeliser
  fits_read,'mgefreemodeli.fits',freemodeli
  fits_read,'serfreemodeli.fits',freemodeliser

  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/galfit_model.fits',img,exten_no=1
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/galfit_model.fits',modelimg,exten_no=2

  find_galaxy,img,m,e,a,xc,yc
  find_galaxy,fixmodeli,m,e,a,xc_m,yc_m
;  xc_m=xc & yc_m=yc
  radius=(10^(findgen(39)*0.05+0.05))
  minradius=[0,(10^(findgen(38)*0.05+0.05))]
  
  mag=fltarr(n_elements(radius))
  modelmag=fltarr(n_elements(radius))
  fimodelmag=fltarr(n_elements(radius))
  fisersicmag=fltarr(n_elements(radius))
  frmodelmag=fltarr(n_elements(radius))
  frsersicmag=fltarr(n_elements(radius))
  imflux=fltarr(n_elements(radius))
  modflux=fltarr(n_elements(radius))
  fixserflux=fltarr(n_elements(radius))
  freeserflux=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     aper,img,xc,yc,maxflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,img,xc,yc,minflux,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     
     aper,modelimg,xc,yc,maxmodelflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,modelimg,xc,yc,minmodelflux,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     
     aper,fixmodeli,xc_m,yc_m,maxmodel,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,fixmodeli,xc_m,yc_m,minmodel,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,fixmodeliser,xc_m,yc_m,maxmodelser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,fixmodeliser,xc_m,yc_m,minmodelser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     
     aper,freemodeli,xc_m,yc_m,fmaxmodel,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,freemodeli,xc_m,yc_m,fminmodel,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,freemodeliser,xc_m,yc_m,fmaxmodelser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,freemodeliser,xc_m,yc_m,fminmodelser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     if (minradius[i] lt 0.025) then begin
        flux=maxflux
        model=maxmodelflux
        fimodel=maxmodel
        fisermodel=maxmodelser
        fmodel=fmaxmodel
        fsermodel=fmaxmodelser
        imflux[i]=maxflux
        modflux[i]=maxmodelflux
        fixserflux[i]=maxmodelser
        freeserflux[i]=fmaxmodelser

     endif else begin
        flux=maxflux-minflux
        model=maxmodelflux-minmodelflux
        fimodel=maxmodel-minmodel
        fisermodel=maxmodelser-minmodelser
        fmodel=fmaxmodel-fminmodel
        fsermodel=fmaxmodelser-fminmodelser
        imflux[i]=maxflux-minflux
        modflux[i]=maxmodelflux-minmodelflux
        fixserflux[i]=maxmodelser-minmodelser
        freeserflux[i]=fmaxmodelser-fminmodelser
        
     endelse

     area=((!PI*(radius[i])^2)-(!PI*(minradius[i])^2))

     mag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(flux)-extinct
     modelmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(model)-extinct
     fimodelmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(fimodel)
     fisersicmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(fisermodel) -extinct
     frmodelmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(fmodel)
     frsersicmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(fsermodel) -extinct

  endfor
  !P.MULTI=[0,1,2]
;  set_plot,'ps'
;  device,filename='convolved_vs_sersic_mge.ps',/color
  djs_plot,radius*scale,mag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,2.5],/xsty,psym=2,title='Convolved Image vs. Fixed Fits',charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,modelmag,color='blue',thick=3
  djs_oplot,radius*scale,fimodelmag,color='red',thick=3
  djs_oplot,radius*scale,fisersicmag,color='green',thick=3
  djs_plot,radius*scale,mag-modelmag,psym=2,ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,2.5],/xsty,yran=[-0.35,0.6] ,charthick=4,xthick=3,ythick=3,charsize=1.5
  djs_oplot,radius*scale,mag-fimodelmag,psym=2,color='red'
  djs_oplot,radius*scale,mag-fisersicmag,psym=2,color='green'
  stop
  djs_plot,radius*scale,mag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,2.5],/xsty,psym=2,title='Convolved Image vs. Free Fits',charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,modelmag,color='blue',thick=3
  djs_oplot,radius*scale,frmodelmag,color='red',thick=3
  djs_oplot,radius*scale,frsersicmag,color='green',thick=3
  djs_plot,radius*scale,mag-modelmag,psym=2,ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,2.5],/xsty,yran=[-0.35,0.6] ,charthick=4,xthick=3,ythick=3,charsize=1.5
  djs_oplot,radius*scale,mag-frmodelmag,psym=2,color='red'
  djs_oplot,radius*scale,mag-frsersicmag,psym=2,color='green'
;  device,/close
;  set_plot,'x'
  !P.MULTI=[0,1,1]



  stop
  
     



END
