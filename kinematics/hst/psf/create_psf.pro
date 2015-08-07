PRO CREATE_PSF

psfradius=1.5 ; in arcseconds
pixelsize=0.01 ;in arcsec, NIFS pixel size is 0.05"

FWHM1=0.086 ;in arcseconds
FWHM2=0.95 ;in arcseconds
FRAC2=0.97
n=4.765

;set up pixel image
npsfpix=(psfradius/pixelsize)*2+1
psfcenter=psfradius/pixelsize
psfim=FLTARR(npsfpix,npsfpix)
pixarr=FINDGEN(npsfpix)
xarr=pixarr#(pixarr*0+1)
yarr=(pixarr*0+1)#pixarr
rarr=SQRT((xarr-psfcenter)^2+(yarr-psfcenter)^2)

;make Gauss+Moffat PSF
psfwidth1=FWHM1 / (2.35 * pixelsize)
psf1=exp(-(rarr/psfwidth1)^2/2.)/(2.*!PI*psfwidth1^2)
psfwidth2=FWHM2 / (2.35 * pixelsize)
psf2=1./(1+rarr/psfwidth2)^n
psf2=FRAC2*psf2/TOTAL(psf2)
psf=psf1+psf2
psf=psf/TOTAL(psf)

WRITEFITS,'psf_gauss_moffat_over5.fits',psf

END
