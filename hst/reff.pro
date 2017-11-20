PRO REFF

  scale=0.025
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/make_decon/serfixmodeli.fits',iimg;,exten_no=1
  zeropoint=25.28697
  ;zeropoint=25.9799
  radius=findgen(100)+1
  find_galaxy,iimg,m,e,a,xc,yc
  iflux=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     aper,iimg,xc,yc,flux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     iflux[i]=flux
  endfor
  radius=radius*scale
  djs_plot,radius,iflux,psym=4,ytitle='counts',xtitle='Radius ["]',xran=[0,2.5],/xstyle
  stop
END
