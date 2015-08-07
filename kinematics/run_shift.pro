@shift_cube.pro
PRO RUN_SHIFT

reffile='CatfbrgnN20150502S0117.fits'

cube = mrdfits(reffile,1,h1)

imsize=SIZE(cube,/dim)
nlambda=imsize[2]
lambda0=SXPAR(h1,'CRVAL3')
dlambda=SXPAR(h1,'CD3_3')
lambda=FINDGEN(imsize[2])*dlambda+lambda0
;Feb 1 23.681
;May 3 -16.551
;correct everything to May 18 data
corfiles=['CatfbrgnN20150201S0190','CatfbrgnN20150201S0191','CatfbrgnN20150201S0192','CatfbrgnN20150201S0194']


;veldiff=vbarycor_Object-vbarycor_Reference 
veldiff=23.681+16.551
FOR i=0,N_ELEMENTS(corfiles)-1 DO SHIFT_CUBE,corfiles[i],lambda,VELDIFF=veldiff
END
