PRO TINY_TIM_FIT_606
  fits_read,'./evstigneeva/temp_psf_r.fits',img,h
  scale=0.025
  ngauss=10
  minlevel=1.e-40

  find_galaxy,img,majoraxis,eps,ang,xc,yc,FRACTION=0.3,/PLOT
  sectors_photometry, img, eps, ang, xc, yc, radius, angle, counts, MINLEVEL=minlevel
  MGE_fit_sectors, radius, angle, counts, eps,$
                   SOL=sol, NGAUSS=ngauss, QBOUNDS=[0.9999999,1.],scale=scale
  modelnorm=1./TOTAL(sol[0,*])
  modelpeak=modelnorm*(sol[0,*])
  modelsig=sol[1,*]
  forprint,modelpeak,modelsig,sol[2,*],format='F,F,F',textout='tinytim_fits_606.dat'
  stop

END
