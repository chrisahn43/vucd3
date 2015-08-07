PRO SHIFT_CUBE, incubestem, lambdaref,VELDIFF=veldiff

incubefile=incubestem+'.fits'
outcubefile=incubestem+'_shift.fits'
ext0 = mrdfits(incubefile,0,h0)
cube = mrdfits(incubefile,1,h1)
var = mrdfits(incubefile,2,h2)
dq = mrdfits(incubefile,3,h3)

imsize=SIZE(cube,/dim)
nlambda=imsize[2]
lambda0=SXPAR(h1,'CRVAL3')
dlambda=SXPAR(h1,'CD3_3')
lambda=FINDGEN(imsize[2])*dlambda+lambda0
IF (KEYWORD_SET(veldiff)) THEN BEGIN
    lold=lambda
    lambda=lambda*(1.+veldiff/2.99792458d5)
    print, veldiff, lambda[0],lold[0]
ENDIF
oversample=10
highreslambda=FINDGEN((nlambda+1)*oversample)*dlambda/FLOAT(oversample)+lambda0

outcube=cube
outvar=var
outdq=dq
FOR i=0,imsize[0]-1 DO BEGIN
print,'Column ',i
FOR j=0,imsize[1]-1 DO BEGIN
    ;spline interpolation gets everything better than 1% for spectra
    ;whose wavelengths match the reference
;    highres=INTERPOL(cube[i,j,*],lambda,highreslambda,/SPLINE)    
;    outcube[i,j,*]=INTERPOL(highres,highreslambda,lambdaref,/SPLINE)
    outcube[i,j,*]=resample_spectra(cube[i,j,*],lambda,lambdaref,oversample=10,npoints=npoints)
    outvar[i,j,*]=resample_spectra(var[i,j,*],lambda,lambdaref,oversample=10,npoints=npoints)
    outdq[i,j,*]=resample_spectra(dq[i,j,*],lambda,lambdaref,oversample=10,npoints=npoints)
    
ENDFOR
ENDFOR

SXADDPAR,h1,'CRVAL3',lambdaref[0]
mwrfits,ext0,outcubefile,h0,/create
mwrfits,outcube,outcubefile,h1
mwrfits,var,outcubefile,h2
mwrfits,dq,outcubefile,h3

loadct,12,/silent
plot,lambda,cube[28,28,*],xrange=[2.29e4,2.30e4],ysty=16
oplot,lambdaref,outcube[28,28,*],color=100
loadct,0,/silent



END
