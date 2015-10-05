pro int_spectra_pixel,minradius,maxradius,xcen,ycen
;  IF (NOT KEYWORD_SET(minradius)) then minradius=0.
;  IF (NOT KEYWORD_SET(maxradius)) then maxradius=1.
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
  r=radial_dist(cubesize[0],cubesize[1],xcen,ycen)
  for k=0,cubesize[2]-1 do begin
     for i=0,cubesize[0]-1 do begin
        for j=0,cubesize[1]-1 do begin
        
           skyval=MEDIAN(cube[5:20,67:81,i])
           outsky[k]=skyval*!PI*(maxradius^2-minradius^2)
           skyval=0
           if (r[i,j] lt maxradius) and (r[i,j] gt minradius) then begin
              outspec[k]+=cube[i,j,k]
              outvar[k]+=varcube[i,j,k]
           endif
        endfor
     endfor
  endfor
;  stop
  

;     maxflux=squarephot(xcen,ycen,maxradius,0.,cube[*,*,i]-skyval)

;     if (minradius gt 0.5) then begin
;        minflux=squarephot(xcen,ycen,minradius,0.,cube[*,*,i]-skyval)
;        flux=maxflux-minflux
;     endif else flux=maxflux
 ;    if (flux gt 0) then stop
;     maxrawflux=squarephot(xcen,ycen,maxradius,0.,cube[*,*,i])
;     if (minradius gt 0.5) then begin
;        minrawflux=squarephot(xcen,ycen,minradius,0.,cube[*,*,i])
;        rawflux=maxrawflux-minrawflux
;     endif else rawflux=maxrawflux

;     maxvar=squarephot(xcen,ycen,maxradius,0.,varcube[*,*,i])
;     if (maxvar lt 0.) then maxvar=maxvar*(-1.)
 ;    if (minradius gt 0.5) then begin
 ;       minvar=squarephot(xcen,ycen,minradius,0.,varcube[*,*,i])
 ;       varflux=maxvar-minvar
 ;    endif else varflux=maxvar
;     outspec[i]=flux
;     outvar[i]=maxvar
;     outrawspec[i]=rawflux
;  ENDFOR
 

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
;  stop
END
