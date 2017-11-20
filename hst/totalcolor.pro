PRO TOTALCOLOR
  scale=0.025
  radius=100.
  zeropoint_v=25.9799
  zeropoint_i=25.28697
  vsun=4.81
  isun=4.53
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/deconvolved_f814.fits',oldimg,head
  hrotate,oldimg,head,img,newhead,1
  find_galaxy,img,majoraxis,eps,ang,xci,yci
  readcol,'./evstigneeva/vucd3_mge_outputsersic.dat',iintensity,isigmaarc,iq,ipa,format='F,F,F,F'
  mge2image,img,xci,yci,iintensity,isigmaarc,iq,ipa,iimg,zeropoint=zeropoint_i,scale=scale,msun=isun

  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/deconvolved_f606.fits',oldimg,head
  hrotate,oldimg,head,img,newhead,1
  find_galaxy,img,majoraxis,eps,ang,xcv,ycv
  readcol,'./evstigneeva/vucd3_mge_outputsersic_free_606.dat',vintensity,vsigmaarc,vq,vpa,format='F,F,F,F'
  mge2image,img,xcv,ycv,vintensity,vsigmaarc,vq,vpa,vimg,zeropoint=zeropoint_v,scale=scale,msun=vsun

  aper,vimg,xcv,ycv,vflux,fluxerr,0.,skyerr,1,radius,-1,[1,1],/silent,setskyval=0.,/flux,/exact
  aper,iimg,xci,yci,iflux,fluxerr,0.,skyerr,1,radius,-1,[1,1],/silent,setskyval=0.,/flux,/exact

  area=!PI*(radius)^2

  vmag=zeropoint_v+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(vflux)
  imag=zeropoint_i+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(iflux)

  color=vmag-imag
  print,color
  stop
  
  
END
