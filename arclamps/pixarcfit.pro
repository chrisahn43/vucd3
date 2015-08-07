pro pixarcfit
  arc=readfits('arc_combine_best8.fits',header,ext=1)
  errarc=readfits('arc_combine_best8.fits',ext=2)
  arcsize=size(arc,/dim)
  crval=sxpar(header,'CRVAL3')
  cdelt=sxpar(header,'CD3_3')
  wave=findgen(arcsize[2])*cdelt+crval

  readcol,'linefwhm.dat',fitwave,fw,format='F,F'

  fwhm=[-99.]
  lamb=[-99.]
  image_f=fltarr(arcsize[0],arcsize[1])
  spec_f=fltarr(arcsize[0],arcsize[1])
  for i=0,arcsize[0]-1 do begin
     for j=0,arcsize[1]-1 do begin
        djs_plot,wave,arc[i,j,*]
        for k=0,n_elements(fitwave)-1 do begin
           tmp=fitwave[k]
           a=where(wave ge tmp-6 and wave le tmp+6,c)
           x=wave[a]
           temp=arc[i,j,*]
           y=temp[a]
           errtemp=errarc[i,j,*]
           err=errtemp[a]
           if max(y) gt 0. then begin
              yfit=gaussfit(x,y,r,nterms=4,chisq=chisq);,measure_error=err)
              djs_oplot,x,yfit,color='blue'
              sig=2.35482*r[2]
              fwhm=[fwhm,sig]
              lamb=[lamb,tmp]
           endif
        endfor
        
        if n_elements(fwhm) gt 2. then begin
           fwhm=fwhm[1:n_elements(fwhm)-1]
           lamb=lamb[1:n_elements(lamb)-1]
           spec_res=lamb/fwhm
           image_f[i,j]=median(fwhm)
           spec_f[i,j]=median(spec_res)
           fwhm=[0.]
           lamb=[0.]
        endif else begin
           spec_f[i,j]=0.
           image_f[i,j]=0.
        endelse
       endfor
  endfor
  
  stop
  name='fwhm_map.fits'
  writefits,name,image_f,hdr
  stop
           
END

