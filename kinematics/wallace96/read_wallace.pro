PRO READ_WALLACE,star
!P.MULTI=[0,1,1]
files=file_search(star,'*.dat')
nfiles=N_ELEMENTS(files)
thresh=1.0

alllambda=[0.]
allspec=[0.]
FOR i=1,nfiles-1 DO BEGIN
    OPENR,1,files[i]
    READF,1,cont
    READF,1,doppler
    wavearr=FLTARR(1551)
    spec1arr=FLTARR(1551)
    spec2arr=FLTARR(1551)
    FOR j=0,1550 DO BEGIN
        READF,1,wavenumber,spec1;,spec2
        wavearr[j]=wavenumber+doppler
        spec1arr[j]=spec1
;        spec2arr[j]=spec2
    ENDFOR    
    CLOSE,1
    FREE_LUN,1
    lambda=1/wavearr*1.e8
    ind=WHERE(spec1arr GT -0.5)
;    plot,lambda[ind],spec1arr[ind]
    alllambda=[alllambda,lambda[ind]]
    allspec=[allspec,spec1arr[ind]]
;    print,cont,doppler
ENDFOR
alllambda=alllambda[1:*]
allspec=allspec[1:*]
uniqind=UNIQ(alllambda,SORT(alllambda))
outlambda=alllambda[uniqind]
outspec=allspec[uniqind]
plot,outlambda,outspec,xrange=[2.34e4,2.38e4]
difflambda=outlambda-outlambda[1:*]
ind=WHERE(outlambda GT 2.28e4 AND outlambda LT 2.369e4,nind)
mediandiff=0.0957031;ABS(MEDIAN(difflambda))
print,mediandiff,MAX(ABS(difflambda[ind]),maxpos)
;oplot,[outlambda[ind[maxpos]]],[outspec[ind[maxpos]]],psym=1,symsize=3
print,outlambda[ind[maxpos]]
;STOP

oversample=10
minlambda=20195.0;MIN(outlambda)
maxlambda=23961.0;MAX(outlambda)
nlambda=LONG((maxlambda-minlambda)/mediandiff)
noversample=nlambda*oversample
intlambda=minlambda+FINDGEN(noversample)*(mediandiff/oversample)
intspec=INTERPOL(outspec,outlambda,intlambda)
;intspec=INTERPOL(outspec,outlambda,intlambda,/SPLINE)
finallambda=minlambda+(FINDGEN(nlambda)+0.5)*(mediandiff)
reform=REFORM(intspec,oversample,nlambda)
finalspectra=TOTAL(reform,1)/FLOAT(oversample)
;oplot,intlambda,intspec;,xrange=[2.29e4,2.295e4]
errspec=INTARR(nlambda)
FOR i=0l,nlambda-1 DO BEGIN
    mindist=MIN(ABS(finallambda[i]-outlambda))
    IF (mindist GT thresh) THEN errspec[i]=1
ENDFOR
oplot,finallambda,errspec,thick=2,color=200
loadct,12,/silent
oplot,finallambda,finalspectra,color=100
nsmooth=FIX(5.5/mediandiff)
pixsize=0.8
simlambda=minlambda+FINDGEN(LONG(nlambda*(mediandiff/pixsize)))*pixsize
simspectra=RESAMPLE_SPECTRA(SMOOTH(finalspectra,nsmooth),finallambda,simlambda,oversample=100)
plot,simlambda,simspectra,color=200,xrange=[2.34e4,2.39e4]
oplot,finallambda,errspec,thick=2,color=200
loadct,0,/silent

outfile=star+'_dop.fits'
mkhdr,headout,finalspectra
sxaddpar,headout,'CRVAL1',finallambda[0]
sxaddpar,headout,'CRPIX1',1
sxaddpar,headout,'CD1_1',mediandiff
WRITEFITS,outfile,finalspectra,headout

outfile=star+'_err.fits'
mkhdr,headout,errspec
sxaddpar,headout,'CRVAL1',finallambda[0]
sxaddpar,headout,'CRPIX1',1
sxaddpar,headout,'CD1_1',mediandiff
WRITEFITS,outfile,errspec,headout

END
