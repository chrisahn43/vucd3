PRO COMPARE_DECON

  scale=0.025
  zeropt=25.28697
  extinct=0.034
  fits_read,'mgefixmodeli.fits',fixmodeli
  fits_read,'serfixmodeli.fits',fixmodeliser
  fits_read,'mgefreemodeli.fits',freemodeli
  fits_read,'serfreemodeli.fits',freemodeliser
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/deconvolved_f814.fits',oldimg,head

  hrotate,oldimg,head,deconimg,newhead,1

  find_galaxy,deconimg,m,e,a,xc,yc
  deconimg=deconimg/1050.

  find_galaxy,fixmodeli,m,e,a,xc_m,yc_m

  radius=(10^(findgen(40)*0.05+0.05))
  minradius=[0,(10^(findgen(39)*0.05+0.05))]
  deconmag=fltarr(n_elements(radius))
  modelmag=fltarr(n_elements(radius))
  sersicmag=fltarr(n_elements(radius))
  frmodelmag=fltarr(n_elements(radius))
  frsersicmag=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     aper,deconimg,xc,yc,maxflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,deconimg,xc,yc,minflux,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

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
        model=maxmodel
        sermodel=maxmodelser
        fmodel=fmaxmodel
        fsermodel=fmaxmodelser
     endif else begin
        flux=maxflux-minflux
        model=maxmodel-minmodel
        sermodel=maxmodelser-minmodelser
        fmodel=fmaxmodel-fminmodel
        fsermodel=fmaxmodelser-fminmodelser
     endelse

     area=((!PI*(radius[i])^2)-(!PI*(minradius[i])^2))

     deconmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(flux)-extinct
     modelmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(model);-extinct
     sersicmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(sermodel) -extinct
     frmodelmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(fmodel) ;-extinct
     frsersicmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(fsermodel) -extinct

  endfor
  !P.MULTI=[0,1,2]
  set_plot,'ps'
  device,filename='deconvolved_vs_sersic_mge_814.ps',/color
  djs_plot,radius*scale,deconmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Deconvolved Image vs. MGE fits',charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,modelmag,color='blue',thick=3
  djs_oplot,radius*scale,frmodelmag,color='red',thick=3
  djs_plot,radius*scale,deconmag-modelmag,psym=2,color='blue',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.35,0.6] ,charthick=4,xthick=3,ythick=3,charsize=1.5
  djs_oplot,radius*scale,deconmag-frmodelmag,psym=2,color='red'
;  stop
  djs_plot,radius*scale,deconmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Deconvolved Image vs. Sersic fits' ,charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,sersicmag,color='blue',thick=3
  djs_oplot,radius*scale,frsersicmag,color='red',thick=3
  djs_plot,radius*scale,deconmag-sersicmag,psym=2,color='blue',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.35,0.6] ,charthick=4,xthick=3,ythick=3,charsize=1.5
  djs_oplot,radius*scale,deconmag-frsersicmag,psym=2,color='red'
;  stop
  djs_plot,radius*scale,sersicmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Fixed Sersic Fits vs. MGE fits' ,charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,modelmag,color='blue',thick=3
  djs_plot,radius*scale,sersicmag-modelmag,psym=2,color='blue',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.05,0.05] ,charthick=4,xthick=3,ythick=3,charsize=1.5
;  stop
  djs_plot,radius*scale,frsersicmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Free Sersic Fits vs. MGE fits' ,charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,frmodelmag,color='red',thick=3
  djs_plot,radius*scale,frsersicmag-frmodelmag,psym=2,color='red',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.05,0.05] ,charthick=4,xthick=3,ythick=3,charsize=1.5
  device,/close
  set_plot,'x'
  stop
  zeropt=25.9799
  extinct=0.061
  fits_read,'mgefixmodelv.fits',fixmodelv
  fits_read,'serfixmodelv.fits',fixmodelvser
  fits_read,'mgefreemodelv.fits',freemodelv
  fits_read,'serfreemodelv.fits',freemodelvser
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/deconvolved_f606.fits',oldimg,head

  hrotate,oldimg,head,deconimg,newhead,1

  find_galaxy,deconimg,m,e,a,xc,yc
  deconimg=deconimg/870.

  find_galaxy,fixmodelv,m,e,a,xc_m,yc_m

  radius=(10^(findgen(40)*0.05+0.05))
  minradius=[0,(10^(findgen(39)*0.05+0.05))]
  deconmag=fltarr(n_elements(radius))
  modelmag=fltarr(n_elements(radius))
  sersicmag=fltarr(n_elements(radius))
  frmodelmag=fltarr(n_elements(radius))
  frsersicmag=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     aper,deconimg,xc,yc,maxflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,deconimg,xc,yc,minflux,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,fixmodelv,xc_m,yc_m,maxmodel,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,fixmodelv,xc_m,yc_m,minmodel,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,fixmodelvser,xc_m,yc_m,maxmodelser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,fixmodelvser,xc_m,yc_m,minmodelser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,freemodelv,xc_m,yc_m,fmaxmodel,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,freemodelv,xc_m,yc_m,fminmodel,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,freemodelvser,xc_m,yc_m,fmaxmodelser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,freemodelvser,xc_m,yc_m,fminmodelser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     if (minradius[i] lt 0.025) then begin
        flux=maxflux
        model=maxmodel
        sermodel=maxmodelser
        fmodel=fmaxmodel
        fsermodel=fmaxmodelser
     endif else begin
        flux=maxflux-minflux
        model=maxmodel-minmodel
        sermodel=maxmodelser-minmodelser
        fmodel=fmaxmodel-fminmodel
        fsermodel=fmaxmodelser-fminmodelser
     endelse

     area=((!PI*(radius[i])^2)-(!PI*(minradius[i])^2))

     deconmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(flux)-extinct
     modelmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(model);-extinct
     sersicmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(sermodel) -extinct
     frmodelmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(fmodel) ;-extinct
     frsersicmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(fsermodel) -extinct

  endfor
  !P.MULTI=[0,1,2]
  set_plot,'ps'
  device,filename='deconvolved_vs_sersic_mge_606.ps',/color
  djs_plot,radius*scale,deconmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Deconvolved Image vs. MGE fits',charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,modelmag,color='blue',thick=3
  djs_oplot,radius*scale,frmodelmag,color='red',thick=3
  djs_plot,radius*scale,deconmag-modelmag,psym=2,color='blue',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.35,0.6] ,charthick=4,xthick=3,ythick=3,charsize=1.5
  djs_oplot,radius*scale,deconmag-frmodelmag,psym=2,color='red'
;  stop
  djs_plot,radius*scale,deconmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Deconvolved Image vs. Sersic fits' ,charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,sersicmag,color='blue',thick=3
  djs_oplot,radius*scale,frsersicmag,color='red',thick=3
  djs_plot,radius*scale,deconmag-sersicmag,psym=2,color='blue',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.35,0.6] ,charthick=4,xthick=3,ythick=3,charsize=1.5
  djs_oplot,radius*scale,deconmag-frsersicmag,psym=2,color='red'
;  stop
  djs_plot,radius*scale,sersicmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Fixed Sersic Fits vs. MGE fits' ,charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,modelmag,color='blue',thick=3
  djs_plot,radius*scale,sersicmag-modelmag,psym=2,color='blue',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.05,0.05] ,charthick=4,xthick=3,ythick=3,charsize=1.5
;  stop
  djs_plot,radius*scale,frsersicmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Free Sersic Fits vs. MGE fits' ,charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,frmodelmag,color='red',thick=3
  djs_plot,radius*scale,frsersicmag-frmodelmag,psym=2,color='red',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.05,0.05] ,charthick=4,xthick=3,ythick=3,charsize=1.5
  device,/close
  set_plot,'x'
  !P.MULTI=[0,1,1]
  stop
  
END
