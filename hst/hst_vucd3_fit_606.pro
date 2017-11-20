PRO HST_VUCD3_FIT_606

  fits_read,'vucd3_skysubtract_r.fits',img,h
  temp=img[431:831,462:862]
  ngauss=20
  minlevel=0.

  find_galaxy,temp,majoraxis,eps,ang,xc,yc,FRACTION=.8

  sectors_photometry,temp,eps,ang,xc,yc,radius,angle,counts,minlevel=minlevel

  readcol,'tinytim_fits_606.dat',normpsf,sigmapsf,format='F,F'

  MGE_fit_sectors,radius,angle,counts,eps,sol=sol,ngauss=ngauss,scale=scale,normpsf=normpsf,sigmapsf=sigmapsf

  scale=0.025
  extinct=0.061
  zp=25.9799
  msun=4.73
  peak=sol[0,*]/(2*!PI*sol[1,*]^2*sol[2,*])
  mu=zp+5*alog10(scale)-2.5*alog10(peak)-extinct
  const=(64800/!PI)^2
  intensity=const*(10^(0.4*(msun-mu)))
  sigmaarc=sol[1,*]*scale
  forprint,intensity,sigmaarc,sol[2,*],format='F,F,F',textout='vucd3_mge_output_606.dat'

  set_plot,'ps'
  !P.Multi=[0,1,2]

  fits_read,'./evstigneeva/galfit_model_r.fits',img,exten_no=1
  fits_read,'./evstigneeva/galfit_model_r.fits',modelimg,exten_no=2
  find_galaxy,img,m,e,a,xc,yc
  xc=96 & yc=98
  find_galaxy,modelimg,m,e,a,xc_mod,yc_mod
  xc_mod=96 & yc_mod=98

  fits_read,'./make_decon/serfreemodelv.fits',serfreev
  find_galaxy,serfreev,m,e,a,xc_ser,yc_ser
  radius=(10^(findgen(85)*0.025))*scale
  minrad=[0,(10^(findgen(84)*0.025))*scale]
  radius=radius/scale
  minrad=minrad/scale
  photmag2=fltarr(n_elements(radius))
  modelmag2=fltarr(n_elements(radius))
  sersicmag=fltarr(n_elements(radius))
  area2=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     aper,img,xc,yc,maxflux_phot,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,img,xc,yc,minflux_phot,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,modelimg,xc_mod,yc_mod,maxflux_model,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,modelimg,xc_mod,yc_mod,minflux_model,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,serfreev,xc_ser,yc_ser,maxflux_ser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,serfreev,xc_ser,yc_ser,minflux_ser,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     if (minrad[i] lt 0.025) then begin
        area2[i]=((!PI*(radius[i])^2)) ;*scale
        
        photflux=maxflux_phot
        modelflux=maxflux_model
        sersicflux=maxflux_ser
     endif else begin
       area2[i]=((!PI*(radius[i])^2)-(!PI*(minrad[i])^2)) ;*scale
       photflux=maxflux_phot-minflux_phot
       modelflux=maxflux_model-minflux_model
       sersicflux=maxflux_ser-minflux_ser
   endelse
     
     photmag2[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(photflux)-0.034
     modelmag2[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(modelflux)-0.034
     sersicmag[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(sersicflux)-extinct
  endfor
  
  readcol,'./evstigneeva/vucd3_mge_outputsersic_free_606.dat',sersiclum,sersicsig,sersicq,sersicpa,format='D,D,D,D'
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/deconvolved_f606.fits',oldimg,head
  hrotate,oldimg,head,img,newhead,1
  find_galaxy,img,m,e,a,xcv,ycv

  a=where(sersicq lt 0.8)
  inintensity=sersiclum[a]
  insigma=sersicsig[a]
  inq=sersicq[a]
  inpa=sersicpa[a]
  mge2image,img,xcv,ycv,inintensity,insigma,inq,inpa,inmodel,zeropoint=zp,scale=scale,msun=msun

  b=where(sersicq gt 0.8)
  outintensity=sersiclum[b]
  outsigma=sersicsig[b]
  outq=sersicq[b]
  outpa=sersicpa[b]
  mge2image,img,xcv,ycv,outintensity,outsigma,outq,outpa,outmodel,zeropoint=zp,scale=scale,msun=msun

  readcol,'vucd3_mge_output_606.dat',mgelum,mgesig,mgeq,format='F,F,F'
  pa=fltarr(n_elements(mgeq))
  pa[*]=0.
  mge2image,img,xcv,ycv,mgelum,mgesig,mgeq,pa,mgemodel,zeropoint=zp,scale=scale,msun=msun
  writefits,'./make_decon/mgedirectfit_606.fits',mgemodel
  mgemag=fltarr(n_elements(radius))
  sersicmag1=fltarr(n_elements(radius))
  sersicmag2=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     aper,inmodel,xcv,ycv,maxinflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,inmodel,xcv,ycv,mininflux,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,outmodel,xcv,ycv,maxoutflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,outmodel,xcv,ycv,minoutflux,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,mgemodel,xcv,ycv,maxmgeflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,mgemodel,xcv,ycv,minmgeflux,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     if (minrad[i] lt 0.025) then begin
        area2[i]=((!PI*(radius[i])^2)) ;*scale
        sersicflux1=maxinflux
        sersicflux2=maxoutflux
        mgeflux=maxmgeflux
     endif else begin
        area2[i]=((!PI*(radius[i])^2)-(!PI*(minrad[i])^2)) ;*scale
        sersicflux1=maxinflux-mininflux
        sersicflux2=maxoutflux-minoutflux
        mgeflux=maxmgeflux-minmgeflux
     endelse
     sersicmag1[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(sersicflux1)
     sersicmag2[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(sersicflux2)
     mgemag[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(mgeflux)
  endfor
  !P.Multi=[0,1,2]
  radius=radius*scale
  device,filename='surfbright_hstindivsersic_center_606.ps',/color
  djs_plot,radius,photmag2,psym=2,xtitle='Radius ["]',ytitle='\mu_{F606W} [Mag/sqare arcsecond]',yran=[22,13.5],xran=[0,1.],charsize=1.5,charthick=4,xthick=3,ythick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],/ysty
  djs_oplot,radius,sersicmag1,color='green',thick=3,linestyle=2 ;,psym=2
  djs_oplot,radius,sersicmag2,color='blue',thick=3,linestyle=2
  djs_oplot,radius,sersicmag,color='red',thick=3
  djs_oplot,radius,modelmag2,color='purple',thick=4
;  djs_oplot,radius,mgemag,color='cyan',thick=4
  items=['n=3.51','n=1.28']
  lines=[2,2]
  color=['green','blue']
  al_legend,items,linestyle=lines,colors=color,/window,background_color='white',charthick=4,thick=3,/top,/right
  djs_plot,radius,photmag2-modelmag2,psym=2,ymargin=[4,0],position=[0.1,0.1,0.95,0.25],xtitle='Radius ["]',ytitle='\Delta \mu',charthick=4,xthick=3,ythick=2,charsize=1.5,ycharsize=0.5,xran=[0,1.],yran=[-0.1,0.1],/ysty
  device,/close
  !P.Multi=[0,1,1]
  set_plot,'x'
  stop


END
