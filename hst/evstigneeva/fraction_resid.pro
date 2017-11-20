pro fraction_resid
  scale=0.025
  fits_read,'galfit_model.fits',img_i,exten_no=1,head
  find_galaxy,img_i,majoraxis,eps,ang,xc_i,yc_i
  fits_read,'galfit_model.fits',resid_i,exten_no=3
  fits_read,'galfit_model_r.fits',img_r,exten_no=1
  find_galaxy,img_r,majoraxis,eps,ang,xc_r,yc_r
  fits_read,'galfit_model_r.fits',resid_r,exten_no=3
  fits_read,'galfit_model_fixed.fits',img_i_fixed,exten_no=1
  find_galaxy,img_i_fixed,majoraxis,eps,ang,xc_i_fixed,yc_i_fixed
  fits_read,'galfit_model_fixed.fits',resid_i_fixed,exten_no=3
  fits_read,'galfit_model_rfixed.fits',img_r_fixed,exten_no=1
  find_galaxy,img_r_fixed,majoraxis,eps,ang,xc_r_fixed,yc_r_fixed
  fits_read,'galfit_model_rfixed.fits',resid_r_fixed,exten_no=3
  
  fraci=resid_i/img_i
;  fraci=abs(fraci)
  fracr=resid_r/img_r
;  fracr=abs(fracr)
  fraci_fixed=resid_i_fixed/img_i_fixed
;  fraci_fixed=abs(fraci_fixed)
  fracr_fixed=resid_r_fixed/img_r_fixed
;  fracr_fixed=abs(fracr_fixed)
  makex,fraci,x,y,/zero
  rarr=SQRT((x-xc_i)^2+(y-yc_i)^2)
                                ; we were looking at residuals divided
                                ; by the models. that is wrong we need
                                ; residuals divided by data.
  set_plot,'ps'
  device,filename='F606W_resid.ps',/color
  djs_plot,rarr,fracr,psym=5,yran=[-100,100],ytitle='\mu_{F606W} Fractional Residuals',xtitle='Radius [Pixels]',/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=0.5,color='red';,thick=3
  djs_oplot,rarr,fracr_fixed,psym=5,color='blue',symsize=0.5;,thick=3
;  stop
  device,/close
  device,filename='F606W_resid_center.ps',/color
  djs_plot,rarr,fracr,psym=1,yran=[-0.3,0.3],ytitle='\mu_{F606W} Fractional Residuals',xtitle='Radius [Pixels]',/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,xran=[0,50],color='red',symsize=0.4
  djs_oplot,rarr,fracr_fixed,psym=1,color='blue',symsize=0.4;,thick=3
;  stop
  device,/close
  device,filename='F814W_resid.ps',/color
  djs_plot,rarr,fraci,psym=4,yran=[-100,100],ytitle='\mu_{F814W} Fractional Residuals',xtitle='Radius [Pixels]',/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3;,thick=3
  djs_oplot,rarr,fraci_fixed,psym=4,color='blue';,thick=3
  device,/close
  device,filename='F814W_resid_center.ps',/color
  djs_plot,rarr,fraci,psym=1,yran=[-0.3,0.3],ytitle='\mu_{F814W} Fractional Residuals',xtitle='Radius [Pixels]',/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,xran=[0,50],symsize=0.4,color='red'
  djs_oplot,rarr,fraci_fixed,psym=1,color='blue',symsize=0.4;,thick=3
  device,/close
  set_plot,'x'
  stop
  
  rad=(findgen(60)+1)
  minrad=(findgen(60))
  residi=fltarr(n_elements(rad))
  residr=fltarr(n_elements(rad))
  residi_fixed=fltarr(n_elements(rad))
  residr_fixed=fltarr(n_elements(rad))
  for i=0,n_elements(rad)-1 do begin
     maxresid_i=djs_phot(xc_i,yc_i,rad[i],0.,fraci,skyval=skyval)
     minresid_i=djs_phot(xc_i,yc_i,minrad[i],0.,fraci,skyval=skyval)
     residi[i]=maxresid_i-minresid_i
     maxresid_r=djs_phot(xc_r,yc_r,rad[i],0.,fracr,skyval=skyval)
     minresid_r=djs_phot(xc_r,yc_r,minrad[i],0.,fracr,skyval=skyval)
     residr[i]=maxresid_r-minresid_r
     maxresid_i_fixed=djs_phot(xc_i_fixed,yc_i_fixed,rad[i],0.,fraci_fixed,skyval=skyval)
     minresid_i_fixed=djs_phot(xc_i_fixed,yc_i_fixed,minrad[i],0.,fraci_fixed,skyval=skyval)
     residi_fixed[i]=maxresid_i_fixed-minresid_i_fixed
     maxresid_r_fixed=djs_phot(xc_r_fixed,yc_r_fixed,rad[i],0.,fracr_fixed,skyval=skyval)
     minresid_r_fixed=djs_phot(xc_r_fixed,yc_r_fixed,minrad[i],0.,fracr_fixed,skyval=skyval)
     residr_fixed[i]=maxresid_r_fixed-minresid_r_fixed     
  endfor
  set_plot,'ps'
  device,filename='fractionresid.ps',/color
  djs_plot,rad*scale,residi,yran=[-100,100],xtitle='Radius ["]',ytitle='Fractional Residuals',charsize=1.5,charthick=4,xthick=3,ythick=3,/ysty,thick=3;,linestyle=1
  djs_oplot,rad*scale,residi_fixed,color='blue',thick=3;,linestyle=2
  djs_oplot,rad*scale,residr,color='green',thick=3;,linestyle=3
  djs_oplot,rad*scale,residr_fixed,color='red',thick=3;,linestyle=4
  items=['F814W out','F606W out','F814W fixed','F606W fixed']
  color=['black','green','blue','red']
  lines=[0,0,0,0]
  al_legend,items,colors=color,linestyle=lines,/window,background_color='white'
  device,/close
  
                                ;central arcsecond best fit
  device,filename='fracresid_best.ps',/color
  djs_plot,rad*scale,residr,color='green',xtitle='Radius ["]',ytitle='Fractional Residuals',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=3,yran=[-25,25],/ysty
  djs_oplot,rad*scale,residi_fixed,color='blue',thick=3
  items=['F606W','F814W fixed']
  color=['green','blue']
  lines=[0,0]
  al_legend,items,colors=color,linestyle=lines,/window,background_color='white'
  device,/close
  set_plot,'x'
  stop
END
