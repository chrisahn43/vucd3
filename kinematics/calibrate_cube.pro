PRO CALIBRATE_CUBE
;apply spectrophotometric calibration to the data cube.  

infile='ngc4395_combine_near9.fits'
outfile='ngc4395_combine_near9_fluxcal.fits'

cal=readplainfits('../flux_cal/fluxcal_hip1123_sep21.fits',/noerr)
callambda=cal[0,*]
calflux=cal[1,*]


zeroext=READFITS(infile,zerohead,ext=0)
im=READFITS(infile,head,ext=1)
varim=READFITS(infile,head2,ext=2)
dqim=READFITS(infile,head3,ext=3)
npixim=READFITS(infile,head4,ext=4)
sigim=READFITS(infile,head5,ext=5)
flagim=READFITS(infile,head6,ext=6)

imsize=SIZE(im,/DIM)
exptime=SXPAR(zerohead,'EXPTIME')
im=im/exptime
varim=varim/exptime^2
lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
lambda=FINDGEN(imsize[2])*dlambda+lambda0

IF (callambda[0] NE lambda[0]) THEN intcalflux=INTERPOL(calflux,callambda,lambda) ELSE intcalflux=calflux


outim=im
outvarim=varim
FOR i=0,imsize[0]-1 DO BEGIN
FOR j=0,imsize[1]-1 DO BEGIN
    outim[i,j,*]=im[i,j,*]*intcalflux
    outvarim[i,j,*]=varim[i,j,*]*intcalflux^2
ENDFOR
ENDFOR

SXADDPAR,head,'BUNIT','erg/cm2/s/A'

fits_write,outfile,zeroext,zerohead
fits_open,outfile,fcbout,/update
fits_write,fcbout,outim,head,extname='SCI'
fits_write,fcbout,outvarim,head2,extname='VAR'
fits_write,fcbout,dqim,head3,extname='DQ'
fits_write,fcbout,npixim,head4,extname='NPIX'
fits_write,fcbout,sigmaim,head5,extname='SIG'
fits_write,fcbout,flagim,head6,extname='FLAG'

fits_close,fcbout
STOP
END
