@myfunct.pro
pro mask_skysubtract
  infile='../data/HST_10137_03_ACS_HRC_F814W_drz.fits'
  fits_read,infile, img, h
  imgsize=size(img,/dim)
  fits_read,infile,mask,exten_no=3
  mdrizz=[7.46042280032,7.9962758789,6.9261413824]
  expt=350.0
  mdrizzcount=mdrizz/expt
  avgdrizz=mean(mdrizzcount)
  img=img+avgdrizz
  xcen=627 & ycen=660
  r=radial_dist(imgsize[0],imgsize[1],xcen,ycen)
  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        if r[i,j] le 100. then mask[i,j]=0
     endfor
  endfor
  xcenstar1=702 & ycenstar1=839
  rstar1=radial_dist(imgsize[0],imgsize[1],xcenstar1,ycenstar1)
  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        if rstar1[i,j] le 10. then mask[i,j]=0
     endfor
  endfor
  xcenstar2=539 & ycenstar2=471
  rstar2=radial_dist(imgsize[0],imgsize[1],xcenstar2,ycenstar2)
  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        if rstar2[i,j] le 10. then mask[i,j]=0
     endfor
  endfor
  xcenstar3=1168 & ycenstar3=265
  rstar3=radial_dist(imgsize[0],imgsize[1],xcenstar3,ycenstar3)
  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        if rstar3[i,j] le 10. then mask[i,j]=0
     endfor
  endfor
  newmask=fltarr(imgsize[0],imgsize[1])
  temperr=img*1050.
  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        if (mask[i,j] eq 7 and temperr[i,j] gt 0.) then newmask[i,j]=1./sqrt(temperr[i,j]) else newmask[i,j]=0.;1.;0.0
     endfor
  endfor
;  stop
;  err=sqrt(img)
  x = (findgen(1300)) # (fltarr(1300) +1)
  y = (fltarr(1300)+1) # (findgen(1300))
  z=img
  p0=[0.02,1.e-6,-4.e-6];,0.,1.e-7,-1.e-7,0.,0.,0.,0.]
  aa=mpfit2dfun('myfunct', x, y, z, err, p0, weights=newmask,covar=covar,perror=perror)
  stop
  plane=aa[0]+aa[1]*x+aa[2]*y
  plane1sigup=(aa[0]+perror[0])+(aa[1])*x+(aa[2])*y
  plane1sigdown=(aa[0]-perror[0])+(aa[1])*x+(aa[2])*y
  writefits,'sky_mask.fits',plane
  writefits,'sky1sigup_mask.fits',plane1sigup
  writefits,'sky1sigdown_mask.fits',plane1sigdown
  stop

END
