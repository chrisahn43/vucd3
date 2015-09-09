pro int_spectra_own,minradius,maxradius,xcen,ycen
  IF (NOT KEYWORD_SET(minradius)) then minradius=0.
  IF (NOT KEYWORD_SET(maxradius)) then maxradius=1.
  infile='../vucd3_combine_best8.fits'
  outfile='vucd3_combine_best8_int_fc'+STRTRIM(FIX(minradius),2)+'-'+STRTRIM(FIX(maxradius),2)+'.fits'

  cube=DOUBLE(READFITS(infile,head,ext=1,/SILENT))
  varcube=DOUBLE(READFITS(infile,head,ext=2,/SILENT))
  cubesize=SIZE(cube,/dim)
  lambda0=SXPAR(head,'CRVAL3')
  dlambda=SXPAR(head,'CD3_3')
  lambda=FINDGEN(cubesize[2])*dlambda+lambda0
  
  outspec=DBLARR(cubesize[2])
  outrawspec=DBLARR(cubesize[2])
  outvar=DBLARR(cubesize[2])
  outsky=DBLARR(cubesize[2])
  for i=0,cubesize[2]-1 do begin
     skyval=MEDIAN(cube[5:20,67:81,i])
     outsky[i]=skyval*!PI*(maxradius^2-minradius^2)
     skyval=0
     

     maxflux=djs_phot(xcen,ycen,maxradius,0.,cube[*,*,i]-skyval)
     if (minradius gt 0.5) then begin
        minflux=djs_phot(xcen,ycen,minradius,0.,cube[*,*,i]-skyval)
        flux=maxflux-minflux
     endif else flux=maxflux
     stop
     maxrawflux=djs_phot(xcen,ycen,maxradius,0.,cube[*,*,i])
     if (minradius gt 0.5) then begin
        minrawflux=djs_phot(xcen,ycen,minradius,0.,cube[*,*,i])
        rawflux=maxrawflux-minrawflux
     endif else rawflux=maxrawflux

     maxvar=djs_phot(xcen,ycen,maxradius,0.,varcube[*,*,i])
     outspec[i]=flux
     outvar[i]=maxvar
     outrawspec[i]=rawflux
  ENDFOR
  ind=WHERE(lambda LT 2.01e4 or lambda GT 2.43e4)
outspec[ind]=0.0
outvar[ind]=0.0
outsky[ind]=0.0

outarr=[[outspec],[outvar],[outsky]]

mkhdr,headout,outarr
sxaddpar,headout,'CRVAL1',lambda0
sxaddpar,headout,'CRPIX1',1
sxaddpar,headout,'CD1_1',dlambda
SXADDPAR,head,'BUNIT','erg/cm2/s/A'
WRITEFITS,outfile,outarr,headout


;plot,lambda,MEDSMOOTH(outspec,5),xrange=[2.0e4,2.45e4],/xsty,ysty=16
plot,lambda,outspec/MEDIAN(outspec),xrange=[2.26e4,2.43e4],/xsty,ysty=16,yrange=[0.6,1.2]
redshift=1.0043
plots,[2.293529e4,2.293529e4]*redshift,!Y.CRANGE
plots,[2.322656e4,2.322656e4]*redshift,!Y.CRANGE
plots,[2.352458e4,2.352458e4]*redshift,!Y.CRANGE
plots,[2.382957e4,2.382957e4]*redshift,!Y.CRANGE
plots,[2.2078e4,2.2078e4]*redshift,!Y.CRANGE
;plots,[2.26435e4,2.26435e4]*redshift,!Y.CRANGE
;[FeIII] lines
;plots,[2.1551e4,2.1551e4]*redshift,!Y.CRANGE,color=220
;plots,[2.2184e4,2.2184e4]*redshift,!Y.CRANGE,color=220
;plots,[2.2427e4,2.2427e4]*redshift,!Y.CRANGE,color=220
;plots,[2.3485e4,2.3485e4]*redshift,!Y.CRANGE,color=220

;plots,[2.1895e4,2.1895e4]*redshift,!Y.CRANGE,color=220



;READCOL,'../h2_wavnum.dat',wavelength,FORMAT='F'
;nlines=N_ELEMENTS(wavelength)
;FOR i=0,nlines-1 DO plots,[wavelength[i],wavelength[i]]*redshift,!Y.CRANGE,color=100

ind=WHERE(lambda GT 2.28e4 and lambda LT 2.40e4)
print,MEDIAN(outspec[ind]/SQRT(outvar[ind]))
;a=GET_KBRD()
;STOP
END

