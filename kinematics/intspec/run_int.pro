@int_spectra_annulus.pro
PRO RUN_INT


fwhmfile='../../sky/vucd3/fwhm_map.fits'
contfile='../lum_model/vucd3_combine_best8_cont.fits'
fwhm=READFITS(fwhmfile)
cont=READFITS(contfile)

radii=[0,1,1.75,2.75,4,8,12]
nradii=N_ELEMENTS(radii)

xcen=43.95
ycen=41.71
makex,fwhm,xarr,yarr,/zero
rarr=SQRT((xarr-xcen)^2+(yarr-ycen)^2)
OPENW,1,'spectra_fwhm.dat'
FOR i=0,nradii-2 DO BEGIN
   ind=WHERE(rarr GT radii[i] AND rarr LE radii[i+1])   
   avfwhm=TOTAL(fwhm[ind]*cont[ind])/TOTAL(cont[ind])
   avrad=TOTAL(rarr[ind]*cont[ind])/TOTAL(cont[ind])
   outfile='vucd3_combine_best8_int_fc'+STRTRIM(FIX(radii[i]),2)+'-'+STRTRIM(FIX(radii[i+1]),2)+'.fits'

   printf,1,radii[i],radii[i+1],avrad,avfwhm,' ',outfile,FORMAT='(2F5.1,2F7.3,2A)'
ENDFOR
close,1



FOR i=0,nradii-2 DO int_spectra_annulus,radii[i],radii[i+1],xcen,ycen

STOP
END
