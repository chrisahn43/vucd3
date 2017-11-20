PRO MAKE_KBAND_IMAGE
                                ;READ IN SSPs AND CONVERT HST FILTERS TO AB MAG
  zeropoint_v_ab=25.9799
  zeropoint_i_ab=25.28697
  zeropoint_v_ve=25.90069
  zeropoint_i_ve=24.86147
  extinct=0.15
  readcol,'../johnson_ssp.dat',johnz,johnage,johnmbol,u,b,v,r,i,j,h,k,format='F,F,F,F,F,F,F,F,F,F,F'
  readcol,'../hstacs_ssp.dat',hstz,hstage,hstmbol,f220w,f250w,f330w,f334n,f435w,f475w,f550m,f555w,f606w,f625w,f658w,f660n,f775w,f814w,f850lp,f892n,format='F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  vext=extinct*1.00600
  v=v+vext
  kext=extinct*0.11471
  k=k+kext
  f606ext=extinct*0.92860
  f606w=f606w+f606ext
  f814ext=extinct*0.59927
  f814w=f814w+f814ext
  ind=where(hstage gt 1.e9)
  f606w=f606w[ind]
  f814w=f814w[ind]
  k=k[ind]
  f606wab=(f606w+zeropoint_v_ab-zeropoint_v_ve)
  f814wab=(f814w+zeropoint_i_ab-zeropoint_i_ve)
  color=f606wab-f814wab
  kcolor=f814wab-k
  stop
  fits_read,'serfixmodeli.fits',iimg
  fits_read,'serfreemodelv.fits',vimg
  vimg=shift(vimg,0,1)
  iimg=iimg[1:*,1:*]
  vimg=vimg[1:*,1:*]
  imgsize=size(iimg,/dim)
  vmag=fltarr(imgsize[0],imgsize[1])
  imag=fltarr(imgsize[0],imgsize[1])


  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        vmag[i,j]=zeropoint_v_ab-2.5*alog10(vimg[i,j]);-0.061
        imag[i,j]=zeropoint_i_ab-2.5*alog10(iimg[i,j]);-0.034
     endfor
  endfor
  diff=vmag-imag
  kimg=fltarr(imgsize[0],imgsize[1])

  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        cind=where(diff[i,j] gt color-0.015 and diff[i,j] lt color+0.015,count)
        if (count eq 1) then begin
           i_k=kcolor[cind]
           kimg[i,j]=iimg[i,j]*(10^(-0.4*i_k))
        endif
        if (count gt 1) then begin
           i_k=mean(kcolor[cind])
           kimg[i,j]=iimg[i,j]*(10^(-0.4*i_k))

        endif
        if (count eq 0) then stop
     endfor
  endfor
  
  writefits,'kbandimage.fits',kimg
  stop
END
