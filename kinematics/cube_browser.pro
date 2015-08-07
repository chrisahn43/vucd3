PRO CUBE_BROWSER


minlambda=2.2e4
maxlambda=2.4e4

infile='CatfbrgnN20140518S0045.fits';
cube=READFITS(infile,head,ext=1,/SILENT)
varcube=READFITS(infile,head,ext=2,/SILENT)
cubesize=SIZE(cube,/dim)
lambda0=SXPAR(head,'CRVAL3')
dlambda=SXPAR(head,'CD3_3')
lambda=FINDGEN(cubesize[2])*dlambda+lambda0
lind=WHERE(lambda GT minlambda AND lambda LT maxlambda,nlind)

window,0,xsize=cubesize[0]*5.,ysize=cubesize[1]*5.,title='Integrated Light'
window,1,xsize=1500,ysize=500,title='Spectrum'

;make total intensity image
outim=TOTAL(cube[*,*,lind],3)/FLOAT(nlind)
wset,0
plot,FINDGEN(cubesize[0]),FINDGEN(cubesize[1]),charsize=0.0001,/nodata,xmargin=[0,0],ymargin=[0,0],/xsty,/ysty,/iso
imgunder,SQRT(outim),/scale
key=''
TVCRS,30.,30.,/data
; to compare to Maraston
binsize=100.
nbin=(maxlambda-minlambda)/binsize+1
binlam=FINDGEN(nbin)*binsize+minlambda
WHILE (key NE 'q') DO BEGIN
    CURSOR,x,y,3,/data
    x=FIX(x+0.5) & y=FIX(y+0.5)
    wset,1
    titstring='X='+STRTRIM(x,2)+' Y='+STRTRIM(y,2)
    
    plot,lambda[lind],(cube[x,y,lind])/MEDIAN(cube[x,y,lind]),/xst,ysty=16,xtitle='Wavelength',ytitle='Flux',yrange=[0.5,1.5],title=titstring
; to compare to Maraston
;    binspec=resample_spectra(cube[x,y,lind]/MEDIAN(cube[x,y,lind]),lambda[lind],binlam)
;    plot,binlam,binspec,/xst,ysty=16,xtitle='Wavelength',ytitle='Flux',yrange=[0.5,1.5],title=titstring
;    oplot,lambda[lind],(varcube[x,y,lind])/MEDIAN(varcube[x,y,lind])*0.5,color=150
    wset,0
    plot,FINDGEN(cubesize[0]),FINDGEN(cubesize[1]),charsize=0.0001,/nodata,xmargin=[0,0],ymargin=[0,0],/xsty,/ysty,/iso
    imgunder,SQRT(outim),/scale

    key=GET_KBRD(0)
ENDWHILE
STOP
END
