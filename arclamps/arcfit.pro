pro arcfit

  arc=readfits('arc_combine_best8.fits',header,ext=1)
  errarc=readfits('arc_combine_best8.fits',ext=2)
  arcsize=size(arc,/dim)
  crval=sxpar(header,'CRVAL3')
  cdelt=sxpar(header,'CD3_3')
  lambda=findgen(arcsize[2])*cdelt+crval
  onedarc=median(median(arc[*,*,*],dim=1),dim=1)
  onederr=median(median(errarc[*,*,*],dim=1),dim=1)
  djs_plot,lambda,onedarc,xran=[20050,24300],/xsty,yran=[0,1000],/ysty
  stop

  readcol,'anil_ArXe_K.dat',arclamb
  arcfwhm=fltarr(n_elements(arclamb))

  for i=0,n_elements(arclamb)-1 do begin
     tmp=arclamb[i]
     a=where(lambda ge tmp-6. and lambda le tmp+6.,c)
     if c gt 0. then begin
        x=lambda[a]
        y=onedarc[a]
        err=onederr[a]
        yfit=gaussfit(x,y,r,nterms=4,chisq=chisq,measure_errors=err)
        djs_plot,lambda,onedarc,xran=[tmp-50,tmp+50],yran=[0,max(y)]
        djs_oplot,x,yfit,color='blue'
        sig=2.35482*r[2]
        arcfwhm[i]=sig
        print,i,tmp,r[1],sig
;        stop
        
     endif
  endfor
  stop
  forprint,arclamb,arcfwhm,format='F,F',textout='linefwhm.dat'

  r=fltarr(n_elements(arclamb))
  r=arclamb/arcfwhm
  djs_plot,arclamb,r,psym=2
  fit=linfit(arclamb,r)
  y=fit[0]+fit[1]*arclamb
  djs_oplot,arclamb,y,color='blue'
  stop
END
