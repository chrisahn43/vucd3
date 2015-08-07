FUNCTION SINGLEPSF, x, y, p,NOCONV=NOCONV

;+
; NAME:  singlepsf
;
; INPUTS:
;        x:  scaled up HST image
;        y:  header for NIFS image
;        p: 4 element array. 
;        p[0]=FWHM in pixels
;        p[1]=scaling between images
;        p[2]=X offset in pixels
;        p[3]=Y offset in pixels

COMMON PSFSHARE,scaling,hstsize,hsthead,psfsize
;STOP
;prepare PSF
psfcenter=FIX(psfsize/2.)
psf=DBLARR(psfsize,psfsize)
makex,psf,xpsf,ypsf
psfwidth=p[0]/2.35*scaling ;in pixels
rpsf=SQRT((xpsf-psfcenter-p[2]*scaling)^2+(ypsf-psfcenter-p[3]*scaling)^2)
psf=exp(-(rpsf/psfwidth)^2/2.)/(2.*!PI*psfwidth^2)
psf=p[1]*psf/TOTAL(psf)


convresult=CONVOLVE(x,psf)
fullresult=REBIN(convresult,hstsize[1],hstsize[2])
HASTROM,fullresult,hsthead,result,outhead,y
;STOP  

RETURN, result

END

FUNCTION DOUBLEPSF, x, y, p,NOCONV=NOCONV

;+
; NAME:  doublepsf
;
; INPUTS:
;        x:  scaled up HST image
;        y:  header for NIFS image
;        p: 4 element array. 
;        p[0]=FWHM of first component
;        p[1]=overall scaling
;        p[2]=X offset in pixels
;        p[3]=Y offset in pixels
;        p[4]=FWHM of second component
;        p[5]=ratio of scaling of second component
COMMON PSFSHARE,scaling,hstsize,hsthead,psfsize
;STOP
;prepare PSF
psfcenter=FIX(psfsize/2.)+0.5
psf=DBLARR(psfsize,psfsize)
makex,psf,xpsf,ypsf
rpsf=SQRT((xpsf-psfcenter-p[2]*scaling)^2+(ypsf-psfcenter-p[3]*scaling)^2)
psfwidth1=p[0]/2.35*scaling ;in pixels
psfwidth2=p[4]/2.35*scaling ;in pixels
psf=(exp(-(rpsf/psfwidth1)^2/2.)/(2.*!PI*psfwidth1^2)+p[5]*exp(-(rpsf/psfwidth2)^2/2.)/(2.*!PI*psfwidth2^2))/(1.+p[5])
psf=p[1]*psf/TOTAL(psf)


convresult=CONVOLVE(x,psf)
fullresult=REBIN(convresult,hstsize[1],hstsize[2])
HASTROM,fullresult,hsthead,result,outhead,y
; writefits,'psf.fits',psf 
;STOP  
RETURN, result

END


FUNCTION GAUSSMOFFAT, x, y, p,NOCONV=NOCONV,WRITEPSF=WRITEPSF

;+
; NAME:  doublepsf
;
; INPUTS:
;        x:  scaled up HST image
;        y:  header for NIFS image
;        p: 4 element array. 
;        p[0]=FWHM of first component
;        p[1]=overall scaling
;        p[2]=X offset in pixels
;        p[3]=Y offset in pixels
;        p[4]=FWHM of second component
;        p[5]=ratio of scaling of second component
;        p[6]=Moffat n
;Moffat from, 
;http://www.mpa-garching.mpg.de/~dimitri/budda/musings.html

COMMON PSFSHARE,scaling,hstsize,hsthead,psfsize
;STOP
;prepare PSF, he claims default is 4.765
psfcenter=FIX(psfsize/2.)+0.5
psf=DBLARR(psfsize,psfsize)
makex,psf,xpsf,ypsf
rpsf=SQRT((xpsf-psfcenter-p[2]*scaling)^2+(ypsf-psfcenter-p[3]*scaling)^2)
psfwidth1=p[0]/2.35*scaling ;in pixels
psfwidth2=p[4]/2.35*scaling ;in pixels
psf1=exp(-(rpsf/psfwidth1)^2/2.)/(2.*!PI*psfwidth1^2)
psf2=1./(1+rpsf/psfwidth2)^p[6]
psf2=p[5]*psf2/TOTAL(psf2)
psf=psf1+psf2
;print,TOTAL(psf)
psf=psf/TOTAL(psf)
psf=p[1]*psf

convresult=CONVOLVE(x,psf)
fullresult=REBIN(convresult,hstsize[1],hstsize[2])
HASTROM,fullresult,hsthead,result,outhead,y
IF (KEYWORD_SET(writepsf)) THEN writefits,'psf_moffat.fits',psf 
;STOP  
RETURN, result

END

