pro psfdrizz_r
  fits_read,'f606_flt.fits',zero,header,exten_no=0
  fits_read,'f606_flt.fits',img,header_sci,exten_no=1
  fits_read,'f606_flt.fits',err,header_err,exten_no=2
  fits_read,'f606_flt.fits',dq,header_dq,exten_no=3
  outhead=header
  imgsize=size(img,/dim)
  fits_read,'psf_606_0focus.fits',psf
  find_galaxy,img,majorAxis,eps,ang,xc,yc
  find_galaxy,psf,majoraxis1,eps1,ang1,xcpsf,ycpsf
  img[*,*]=0.0
  img[xc-(xcpsf-1):xc+xcpsf,yc-(ycpsf-1):yc+ycpsf]=psf
  outim='psfonimg_r.fits'
  fits_write,outim,zero,outhead
  fits_open,outim,fcbout,/update
  fits_write,fcbout,img,header_sci,extname='SCI'
  fits_write,fcbout,err,header_err,extname='ERR'
  fits_write,fcbout,dq,header_dq,extname='DQ'
  fits_close,fcbout
    
  stop

END
