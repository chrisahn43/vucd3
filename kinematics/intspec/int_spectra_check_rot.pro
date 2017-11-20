PRO INT_SPECTRA_CHECK_ROT,minradius,maxradius,xcen,ycen
  IF (NOT KEYWORD_SET(minradius)) THEN minradius=0.
  IF (NOT KEYWORD_SET(maxradius)) THEN maxradius=1.
  
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

  shift=[0,0.875,2,5,10.875]
  temprad=maxradius-0.3
  ind=where(shift lt temprad)
  shift=shift[ind]
;  stop
  outspec=DBLARR(cubesize[2])
  outrawspec=DBLARR(cubesize[2])
  outvar=DBLARR(cubesize[2])
  outsky=DBLARR(cubesize[2])
  for i=0,cubesize[2]-1 do begin
     skyval=MEDIAN(cube[5:20,67:81,i])
     outsky[i]=skyval*!PI*(maxradius^2-minradius^2)
     skyval=0
     
     aper,(cube[*,*,i]-skyval),xcen,ycen-max(shift),maxflux,fluxerr,sky,skyerr,1.0,[maxradius],[20,25],[-32000,32000],/silent,setskyval=0.,/flux,/exact,/NAN
     if (minradius gt 0.5) then begin
        aper,(cube[*,*,i]-skyval),xcen,ycen-shift[n_elements(shift)-2],minflux,fluxerr,sky,skyerr,1.0,[minradius],[20,25],[-32000,32000],/silent,setskyval=0.,/flux,/exact,/NAN
        flux=maxflux-minflux
     endif else flux=maxflux
     aper,(cube[*,*,i]),xcen,ycen-max(shift),maxrawflux,rawfluxerr,rawsky,rawskyerr,1.0,[maxradius],[20,25],[-32000,32000],/silent,setskyval=0.,/flux,/exact,/NAN
     if (minradius gt 0.5) then begin
        aper,(cube[*,*,i]),xcen,ycen-shift[n_elements(shift)-2],minrawflux,rawfluxerr,rawsky,rawskyerr,1.0,[minradius],[20,25],[-32000,32000],/silent,setskyval=0.,/flux,/exact,/NAN
        rawflux=maxrawflux-minrawflux
     endif else rawflux=maxrawflux

     aper,(varcube[*,*,i]),xcen,ycen-max(shift),varflux,varfluxerr,varsky,varskyerr,1.0,[maxradius],[20,25],[-32000,32000],/silent,setskyval=0.,/flux,/exact,/nan
     outspec[i]=flux
     outvar[i]=varflux
     outrawspec[i]=rawflux
  endfor
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
  
  
  plot,lambda,outspec/MEDIAN(outspec),xrange=[2.26e4,2.43e4],/xsty,ysty=16,yrange=[0.6,1.2]
  redshift=1.0043
  plots,[2.293529e4,2.293529e4]*redshift,!Y.CRANGE
  plots,[2.322656e4,2.322656e4]*redshift,!Y.CRANGE
  plots,[2.352458e4,2.352458e4]*redshift,!Y.CRANGE
  plots,[2.382957e4,2.382957e4]*redshift,!Y.CRANGE
  plots,[2.2078e4,2.2078e4]*redshift,!Y.CRANGE
  ind=WHERE(lambda GT 2.28e4 and lambda LT 2.40e4)
  print,MEDIAN(outspec[ind]/SQRT(outvar[ind]))

END
