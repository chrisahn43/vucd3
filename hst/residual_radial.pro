pro residual_radial
  fits_read,'galfit_model.fits',resid_i,exten_no=3
  fits_read,'galfit_model_r.fits',resid_r,exten_no=3
  xc_i=98 & yc_i=99
  xc_r=98 & yc_r=99
  scale=0.025
  rad=(findgen(150)+1)
  minrad=(findgen(150))
  iflux=fltarr(n_elements(rad))
  rflux=fltarr(n_elements(rad))
  for i=0,n_elements(rad)-1 do begin
     maxflux_i=djs_phot(xc_i,yc_i,rad[i],0.,resid_i,skyval=skyval)
     minflux_i=djs_phot(xc_i,yc_i,minrad[i],0.,resid_i,skyval=skyval)
     maxflux_r=djs_phot(xc_r,yc_r,rad[i],0.,resid_r,skyval=skyval)
     minflux_r=djs_phot(xc_r,yc_r,minrad[i],0.,resid_r,skyval=skyval)
     iflux[i]=maxflux_i-minflux_i
     rflux[i]=maxflux_r-minflux_r
  endfor

  djs_plot,rad,iflux,psym=4,ytitle='Residual in Model from Galfit',xtitle='Radius',charsize=1.5,charthick=4,xthick=3,ythick=3
  djs_oplot,rad,rflux,psym=4,color='blue'
  stop

END
