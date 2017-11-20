PRO COLORPROFILE
  scale=0.025
  radius=findgen(100)+1;(10^(findgen(40)*0.05+0.05))
  minradius=findgen(100);[0,(10^(findgen(39)*0.05+0.05))]
  zeropoint_v=25.9799
  zeropoint_i=25.28697
  ei=0.034
  ev=0.061
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/galfit_model.fits',iimg,exten_no=1
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/galfit_model_r.fits',vimg,exten_no=1
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/galfit_model.fits',ifree,exten_no=2
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/galfit_model_fixed.fits',ifixed,exten_no=2
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/galfit_model_r.fits',vfree,exten_no=2
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/galfit_model_rfixed.fits',vfixed,exten_no=2
  find_galaxy,iimg,m,e,a,xci,yci
  find_galaxy,vimg,m,e,a,xcv,ycv
  
  fits_read,'serfixmodeli.fits',serfixi
  fits_read,'serfreemodeli.fits',serfreei
  fits_read,'serfixmodelv.fits',serfixv
  fits_read,'serfreemodelv.fits',serfreev
  find_galaxy,serfixi,m,e,a,xci_m,yci_m
  find_galaxy,serfixv,m,e,a,xcv_m,ycv_m
  fits_read,'mgedirectfit.fits',mgei
  fits_read,'mgedirectfit_606.fits',mgev
  find_galaxy,mgei,m,e,a,xci_mge,yci_mge
  find_galaxy,mgev,m,e,a,xcv_mge,ycv_mge
  idata=fltarr(n_elements(radius))
  vdata=fltarr(n_elements(radius))
  imagfixed=fltarr(n_elements(radius))
  vmagfixed=fltarr(n_elements(radius))
  imagfree=fltarr(n_elements(radius))
  vmagfree=fltarr(n_elements(radius))
  imagserfixed=fltarr(n_elements(radius))
  imagserfree=fltarr(n_elements(radius))
  vmagserfixed=fltarr(n_elements(radius))
  vmagserfree=fltarr(n_elements(radius))
  imagmge=fltarr(n_elements(radius))
  vmagmge=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     area=((!PI*(radius[i])^2)-(!PI*(minradius[i])^2))

     aper,iimg,xci,yci,maxidata,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,iimg,xci,yci,minidata,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,vimg,xcv,ycv,maxvdata,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,vimg,xcv,ycv,minvdata,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,ifree,xci,yci,maxifree,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,ifree,xci,yci,minifree,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,ifixed,xci,yci,maxifixed,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,ifixed,xci,yci,minifixed,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,vfree,xcv,ycv,maxvfree,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,vfree,xcv,ycv,minvfree,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,vfixed,xcv,ycv,maxvfixed,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,vfixed,xcv,ycv,minvfixed,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,serfixi,xci_m,yci_m,maxifixser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,serfixi,xci_m,yci_m,minifixser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,serfreei,xci_m,yci_m,maxifreeser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,serfreei,xci_m,yci_m,minifreeser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,serfixv,xcv_m,ycv_m,maxvfixser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,serfixv,xcv_m,ycv_m,minvfixser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,serfreev,xcv_m,ycv_m,maxvfreeser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,serfreev,xcv_m,ycv_m,minvfreeser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,mgev,xcv_mge,ycv_mge,maxvmgeflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,mgev,xcv_mge,ycv_mge,minvmgeflux,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,mgei,xci_mge,yci_mge,maximgeflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,mgei,xci_mge,yci_mge,minimgeflux,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     if (minradius[i] lt 0.025) then begin
        iflux=maxidata
        vflux=maxvdata
        ifreeflux=maxifree
        ifixedflux=maxifixed
        vfreeflux=maxvfree
        vfixedflux=maxvfixed
        ifixserflux=maxifixser
        ifreeserflux=maxifreeser
        vfixserflux=maxvfixser
        vfreeserflux=maxvfreeser
        vmgeflux=maxvmgeflux
        imgeflux=maximgeflux
     endif else begin
        iflux=maxidata-minidata
        vflux=maxvdata-minvdata
        ifreeflux=maxifree-minifree
        ifixedflux=maxifixed-minifixed
        vfreeflux=maxvfree-minvfree
        vfixedflux=maxvfixed-minvfixed
        ifixserflux=maxifixser-minifixser
        ifreeserflux=maxifreeser-minifreeser
        vfixserflux=maxvfixser-minvfixser
        vfreeserflux=maxvfreeser-minvfreeser
        vmgeflux=maxvmgeflux-minvmgeflux
        imgeflux=maximgeflux-minimgeflux
     endelse
     idata[i]=zeropoint_i+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(iflux)-ei
     vdata[i]=zeropoint_v+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(vflux)-ev
     imagfixed[i]=zeropoint_i+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(ifixedflux)-ei
     vmagfixed[i]=zeropoint_v+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(vfixedflux)-ev
     imagfree[i]=zeropoint_i+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(ifreeflux)-ei
     vmagfree[i]=zeropoint_v+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(vfreeflux)-ev
     imagserfixed[i]=zeropoint_i+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(ifixserflux)-ei
     imagserfree[i]=zeropoint_i+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(ifreeserflux)-ei
     vmagserfixed[i]=zeropoint_v+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(vfixserflux)-ev
     vmagserfree[i]=zeropoint_v+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(vfreeserflux)-ev
     vmagmge[i]=zeropoint_v+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(vmgeflux)
     imagmge[i]=zeropoint_i+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(imgeflux)
  endfor

  colordata=vdata-idata
  colorfree=vmagfree-imagfree
  vcolorfixed=vmagfixed-imagfree
  icolorfixed=vmagfree-imagfixed
  sercolorfree=vmagserfree-imagserfree
  vsercolorfixed=vmagserfixed-imagserfree
  isercolorfixed=vmagserfree-imagserfixed
  mgecolor=vmagmge-imagmge
  set_plot,'ps'
  device,filename='color_all.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7
;  myplot,filename='color_all.ps'
  djs_plot,radius*scale,colordata,ytitle='(\mu_{F606W} - \mu_{F814W}) [Mag/asec^2]',xtitle='Radius ["]',yran=[0.3,1.],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,thick=5,xstyle=8,ymargin=[4,4],psym=4,xran=[0,1.5]
  djs_oplot,radius*scale,colorfree,thick=8
  djs_oplot,radius*scale,sercolorfree,thick=8,linestyle=2
  djs_oplot,radius*scale,vcolorfixed,thick=8,color='blue'
  djs_oplot,radius*scale,vsercolorfixed,thick=8,color='blue',linestyle=2
  djs_oplot,radius*scale,icolorfixed,thick=8,color='red'
  djs_oplot,radius*scale,isercolorfixed,thick=8,color='red',linestyle=2
;  djs_oplot,radius*scale,mgecolor,thick=4,color='green',linestyle=2
  rad=findgen(40)*2
  axis,xaxis=1,xtickv=rad,xcharsize=1.5,charthick=4,xthick=4,xtitle='Radius [Pixels]',/xsty,xran=[1,60]
  items=['Data','Convolved Model','Unconvolved Model'];,'Direct MGE Fit']
  lines=[0,0,2];,2]
  sym=[4,0,0];,0]
  colors=['black','black','black'];,'green']
  al_legend,items,linestyle=lines,psym=sym,color=colors,/window,background_color='white',charthick=4,thick=3,/bottom,/right ;,charsize=1.5
  xyouts,[0.02],[0.95],['VUCD3'],charthick=3,charsize=1.5,/data
  device,/close
  set_plot,'x'
  stop
  
END

PRO CROSS_CON
  scale=0.025
  radius=findgen(100)+1;(10^(findgen(40)*0.05+0.05))
  minradius=findgen(100);[0,(10^(findgen(39)*0.05+0.05))]
  zeropoint_v=25.9799
  zeropoint_i=25.28697
  ei=0.034
  ev=0.061
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/vucd3_skysubtract.fits',iimg
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/vucd3_skysubtract_r.fits',vimg
  find_galaxy,iimg,m,e,a,xci,yci
  find_galaxy,vimg,m,e,a,xcv,ycv
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/temp_psf.fits',ipsf
  ipsf=ipsf/TOTAL(ipsf)
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/temp_psf_r.fits',vpsf
  vpsf=vpsf/TOTAL(vpsf)
  iconv=convolve(iimg,vpsf)
  vconv=convolve(vimg,ipsf)
  find_galaxy,iconv,m,e,a,xci_c,yci_c
  find_galaxy,vconv,m,e,a,xcv_c,ycv_c
  idata=fltarr(n_elements(radius))
  vdata=fltarr(n_elements(radius))
  icross=fltarr(n_elements(radius))
  vcross=fltarr(n_elements(radius))

  for i=0,n_elements(radius)-1 do begin
     area=((!PI*(radius[i])^2)-(!PI*(minradius[i])^2))

     aper,iimg,xci,yci,maxidata,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,iimg,xci,yci,minidata,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,vimg,xcv,ycv,maxvdata,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,vimg,xcv,ycv,minvdata,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,iconv,xci_c,yci_c,maxiconv,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,iconv,xci_c,yci_c,miniconv,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,vconv,xcv_c,ycv_c,maxvconv,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,vconv,xcv_c,ycv_c,minvconv,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     if (minradius[i] lt 0.025) then begin
        iflux=maxidata
        vflux=maxvdata
        iconvflux=maxiconv
        vconvflux=maxvconv
     endif else begin
        iflux=maxidata-minidata
        vflux=maxvdata-minvdata
        iconvflux=maxiconv-miniconv
        vconvflux=maxvconv-minvconv
     endelse
     idata[i]=zeropoint_i+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(iflux)-ei
     vdata[i]=zeropoint_v+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(vflux)-ev
     icross[i]=zeropoint_i+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(iconvflux)-ei
     vcross[i]=zeropoint_v+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(vconvflux)-ev
     
  endfor

  colordata=vdata-idata
  colorconv=vcross-icross

  set_plot,'ps'
  device,filename='crosscolor.ps',/color
  djs_plot,radius*scale,colordata,ytitle='(\mu_{F606W} - \mu_{F814W}) [mag/arcsec^2]',xtitle='Radius ["]',yran=[0.3,1.],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,thick=3,xstyle=8,ymargin=[4,4],psym=4
  djs_oplot,radius*scale,colorconv,thick=3,psym=4,color='red',linestyle=2
  rad=findgen(50)*2+2
  axis,xaxis=1,xtickv=rad,xcharsize=1.5,charthick=4,xthick=4,xtitle='Radius [Pixels]',/xsty,xran=[1,100]
  items=['Data','Cross Convolved']
  lines=[0,0]
  sym=[4,4]
  color=['black','red']
  al_legend,items,linestyle=lines,psym=sym,colors=color,/window,background_color='white',charthick=4,thick=3,/bottom,/right;,charsize=1.5
  device,/close
  set_plot,'x'
  stop
 
END

