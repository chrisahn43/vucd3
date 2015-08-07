@meanclip_spectra
PRO ANALYZE_CUBE

READCOL,'Catlist',imfiles,FORMAT='A'
tags=STRTRIM(INDGEN(N_ELEMENTS(imfiles)),2)
nfiles=N_ELEMENTS(imfiles)
veldiffs=REPLICATE(0.0,nfiles);[0.0,0.0,0.0,0.0,0.0,0.0,0.0]

;define region for sky subtraction
xskymin=44 & xskymax=58
yskymin=44 & yskymax=58
allspec=FLTARR(nfiles,2040)
alllambda=FLTARR(nfiles,2040)
loadct,4
FOR i=0,nfiles-1 DO BEGIN
    im=READFITS(imfiles[i],head,ext=1,/SILENT)
    junk=READFITS(imfiles[i],head0,ext=0,/SILENT)
    imsize=SIZE(im,/dim)
    lambda0=SXPAR(head,'CRVAL3')
    dlambda=SXPAR(head,'CD3_3')
;    print,SXPAR(head0,'PA')
    IF (i eq 0) THEN BEGIN
        reflambda0=lambda0
        refdlambda=dlambda
    ENDIF
    lambda=FINDGEN(imsize[2])*dlambda+lambda0
    lambda=lambda*(1.-veldiffs[i]/2.99792458d5)
    alllambda[i,*]=lambda

    skyintermediate=MEDIAN(im[xskymin:xskymax,yskymin:yskymax,*],DIM=1)
    skyspec=MEDIAN(skyintermediate,DIM=1)
    allskyspec=REBIN(REFORM(skyspec,1,1,imsize[2]),imsize[0],imsize[1],imsize[2])
    skysubspec=im-allskyspec
    outim=TOTAL(skysubspec,3)
    outfile='skysub_'+tags[i]+'.fits'
;    writefits,outfile,outim/FLOAT(imsize[2])

    maxval=MAX(outim,maxpos)
    xinit=maxpos MOD imsize[0]
    yinit=maxpos/imsize[0]
    GCNTRD,outim,xinit,yinit,xcen,ycen,4.5

    ;create spectra using aperture photometry at each wavelength
    outspec=FLTARR(imsize[2])
    FOR j=0,imsize[2]-1 DO BEGIN
        APER,skysubspec[*,*,j],xcen,ycen,flux,fluxerr,sky,skyerr,1.0,[2.],[15,20],[-32000,32000],/SILENT,SETSKYVAL=0.0,/flux
        outspec[j]=flux[0]
    ENDFOR
    ind=WHERE(lambda GT 22140. AND lambda LT 22866,nind)
    norm=MEDIAN(outspec[ind])/4500.
    outnormspec=outspec/norm
    print,tags[i],' ',xcen+1,ycen+1,lambda[1000],TOTAL(im[10:50,10:50,ind]),norm,FORMAT='(A7,A1,2F8.2,F10.3,E10.2,F8.3)'

loadct,4,/silent
    IF i EQ 0 THEN $
      plot,lambda,outnormspec,yrange=[2e3,6e3],/ysty,xrange=[2.28e4,2.35e4] ELSE $
      oplot,lambda,outnormspec+i*50.;,color=i*20

    allspec[i,*]=outnormspec
ENDFOR

centerspectra=MEANCLIP_SPECTRA(alllambda,allspec,npoints=npoints)
;plot,lambda,centerspectra,yrange=[2e3,7e3],/ysty,xrange=[2.25e4,2.4e4],thick=
;oplot,alllambda[0,*],centerspectra-1000,thick=2

mkhdr,headout,centerspectra
sxaddpar,headout,'CRVAL1',reflambda0
sxaddpar,headout,'CRPIX1',1
sxaddpar,headout,'CD1_1',refdlambda
WRITEFITS,'center_spec.fits',centerspectra,headout

;plot,alllambda[0,*],centerspectra,xrange=[2.29e4,2.31e4],ysty=16
STOP
;this plot shows that the  data is have gridding problems because of
;the jagged profile down the right side.

plot,outim[*,27]  
oplot,outim[*,27]

totpix=FLOAT(xcentermax-xcentermin+1)*(ycentermax-ycentermin+1)
centerspec=TOTAL(TOTAL(outim[xcentermin:xcentermax,ycentermin:ycentermax,*],1),1)/totpix
mkhdr,headout,centerspec
sxaddpar,headout,'CRVAL1',lambda0
sxaddpar,headout,'CRPIX1',1
;this doesn't work, need to use CD1_1
sxaddpar,headout,'CRDELT1',dlambda
ind=WHERE(lambda GT 24200.)
centerspec[ind]=0.
plot,lambda,centerspec,yrange=[0,300.]
writefits,'centerspec.fits',centerspec,headout
;STOP
xleftmin=xcen-5
xleftmax=xcen-4
yleftmin=ycen-3
yleftmax=ycen+3
totpix=FLOAT(xleftmax-xleftmin+1)*(yleftmax-yleftmin+1)
leftspec=TOTAL(TOTAL(outim[xleftmin:xleftmax,yleftmin:yleftmax,*],1),1)/totpix
writefits,'leftspec.fits',leftspec,headout

xrightmin=xcen+4
xrightmax=xcen+5
yrightmin=ycen-3
yrightmax=ycen+3
totpix=FLOAT(xrightmax-xrightmin+1)*(yrightmax-yrightmin+1)
rightspec=TOTAL(TOTAL(outim[xrightmin:xrightmax,yrightmin:yrightmax,*],1),1)/totpix
writefits,'rightspec.fits',rightspec,headout

plot,lambda,leftspec,yrange=[0,140],xrange=[2.29e4,2.40e4],/xsty,/nodata,/ysty
;oplot,lambda,SMOOTH(leftspec,2),color=100
;oplot,lambda,SMOOTH(rightspec-MEAN(rightspec-leftspec),2),color=200
oplot,lambda,leftspec,color=100
oplot,lambda,(rightspec-MEAN(rightspec-leftspec)),color=200
loadct,0
oplot,lambda,(centerspec-MEAN(centerspec-leftspec))

myplot,file='analyze_cube.ps',xsize=8,ysize=3,/inches
plot,lambda/1.e4,centerspec,xrange=[2.29,2.40],/xsty,xtitle='Wavelength [!9m!xm]',ytitle='Flux',charsize=1.8,yrange=[100,180]
device,/close
set_plot,'x'

STOP
END
