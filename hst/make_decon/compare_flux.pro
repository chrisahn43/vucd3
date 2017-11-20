PRO COMPARE_FLUX

  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/galfit_model.fits',img,exten_no=1
  find_galaxy,img,m,e,a,xc,yc
  makex,img,ximg,yimg,/zero
  ximg=ximg-xc
  yimg=yimg-yc
  

  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/deconvolved_f814.fits',oldimg,head
  hrotate,oldimg,head,deconimg,newhead,1
  find_galaxy,deconimg,m,e,a,xc_d,yc_d
  deconimg=deconimg/1050.
  makex,deconimg,deconx,decony,/zero
  deconx=deconx-xc_d
  decony=decony-yc_d

  fits_read,'mgefixmodeli.fits',fixmodeli
  find_galaxy,fixmodeli,m,e,a,xc_m,yc_m
  makex,fixmodeli,xmod,ymod,/zero
  xmod=xmod-xc_m
  ymod=ymod-yc_m

                                ;X-AXIS PLOT COMPARING THE THREE
  djs_plot,ximg[*,yc],img[*,yc],xran=[-10,10],yran=[0,50]
  djs_oplot,deconx[*,yc_d],deconimg[*,yc_d],color='blue'
  djs_oplot,xmod[*,yc_m],fixmodeli[*,yc_m],color='red'
  stop
  djs_plot,yimg[xc,*],img[xc,*],linestyle=2,xran=[-10,10],yran=[0,50]
  djs_oplot,decony[xc_d,*],deconimg[xc_d,*],color='blue',linestyle=2
  djs_oplot,ymod[xc_m,*],fixmodeli[xc_m,*],color='red',linestyle=2
  stop
  djs_plot,ximg[*,yc],img[*,yc],xran=[35,45],psym=2;,yran=[0,50]
  djs_oplot,deconx[*,yc_d],deconimg[*,yc_d],color='blue',psym=2
  djs_oplot,xmod[*,yc_m],fixmodeli[*,yc_m],color='red',psym=2

  djs_plot,yimg[xc,*],img[xc,*],psym=2,xran=[35,45],yran=[0,0.2]
  djs_oplot,decony[xc_d,*],deconimg[xc_d,*],color='blue',psym=2
  djs_oplot,ymod[xc_m,*],fixmodeli[xc_m,*],color='red',psym=2

END
