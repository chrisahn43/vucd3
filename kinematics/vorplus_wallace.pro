@ppxf_wallace.pro
PRO VORPLUS_WALLACE,MIN_SN=min_sn

;binning in elliptical sections beyond maxvor (~0.5")
maxvor=10. ;in pixels
axrat=0.85 ; in 1D HST model, ellip~0.15 from 5-14 pixels
radbin=[10,12,15,20]
nradbin=N_ELEMENTS(radbin)-1
pa=-49.5/!RADEG
angbinsize=!PI/2.
angbin=FINDGEN(5)*angbinsize
nangbin=N_ELEMENTS(angbin)-1

IF NOT (KEYWORD_SET(min_sn)) THEN min_sn=25.


infile='m60-ucd1_combine_best9.fits'
fwhmfile='../sky/m60-ucd1/fullcube_disp_med.fits'
contfile='lum_model/m60-ucd1_combine_best9_cont.fits'
initvel=1290.
minx=19
miny=19
;IRAF center measurement _cont.fits is 40.26 42.90
xcen=39.26-minx
ycen=41.90-miny

ifu=READFITS(infile,head,ext=1)
var=READFITS(infile,ext=2)
fwhm=READFITS(fwhmfile)
cont=READFITS(contfile)
ifu=ifu[minx:65,miny:65,*]
var=var[minx:65,miny:65,*]
fwhm=fwhm[minx:65,miny:65]
cont=cont[minx:65,miny:65]
size=SIZE(ifu,/DIM)
tag='best9'
lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
lambda=FINDGEN(size[2])*dlambda+lambda0
makex,ifu,x,y,/zero
x=x-xcen
y=y-ycen
xrot=x*COS(pa)-y*SIN(pa)
yrot=x*SIN(pa)+y*COS(pa)
rellip=SQRT(xrot^2+(yrot/axrat)^2)
;create angle offset from major axis by angbinsize/2.
pabin=pa-angbinsize/2.
xang=x*COS(pabin)-y*SIN(pabin)
yang=x*SIN(pabin)+y*COS(pabin)
ang_rot=ATAN(yang,xang)+!PI
;WRITEFITS,'test.fits',ang_rot
;from -!PI to +!PI, major axis in +x +y is 0
;STOP

outplotfile='vor_out/newbin_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'.ps'
!P.THICK=2
!P.CHARTHICK=2
!X.THICK=2
!Y.THICK=2
;set_plot,'ps'
;device,file=outplotfile,ysize=9,yoff=1,/inches,/color

;myplot_nott,file=outplotfile,ysize=9,yoff=1,/inches
!P.MULTI=[0,1,2]

;calculate signal to noise
minfit=2.28e4  ;part of the target spectra to extract
maxfit=2.395e4
fitind=WHERE(lambda GT minfit AND lambda LT maxfit)
signal=FLTARR(size[0:1])
noise=FLTARR(size[0:1])
FOR i=0,size[0]-1 DO BEGIN
FOR j=0,size[1]-1 DO BEGIN
    uspec=ifu[i,j,fitind]
    uvar=var[i,j,fitind]
    signal[i,j]=MEDIAN(uspec)
    noise[i,j]=MEDIAN(SQRT(uvar))
ENDFOR
ENDFOR
;imsize=SIZE(signal,/dim)
;maxval=MAX(signal,maxpos)
;xinit=maxpos MOD imsize[0]
;yinit=maxpos/imsize[0]
;GCNTRD,signal,xinit,yinit,xcen,ycen,4.5
;print,'Center',xcen,ycen
;Clipped Center      20.2577      22.8544



signoise=signal/noise
;WRITEFITS,'vor_out/sn.fits',signoise
logsn=ALOG10(signoise)
logsn=logsn-MIN(logsn[WHERE(FINITE(logsn) EQ 1,comp=bad)])
IF (bad[0] GE 0) THEN logsn[bad]=0.0
plot,FINDGEN(10),/nodata,/iso
loadct,13,/silent
imgunder,logsn*255./MAX(logsn)
loadct,0,/silent

binind=WHERE(rellip LE maxvor)

voronoi_2d_binning, x[binind], y[binind], signal[binind], noise[binind], min_sn, binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale, /PLOT, /QUIET,/NO_CVT ;-- S/N=10 failed without this
;reassign S/N<20 pixels to nearest bins
allbinnum=INTARR(size[0],size[1])
allbinnum[*,*]=999;!VALUES.F_NAN
allbinnum[binind]=binNum
WRITEFITS,'test2.fits',allbinnum
STOP
ibin=INDGEN(N_ELEMENTS(xbar))
badind=WHERE(sn LT 20.,nbad)
FOR i=0,nbad-1 DO BEGIN
   notbin=WHERE(ibin NE badind[i])
   inbin=WHERE(binnum EQ badind[i],ninbin)
   FOR j=0,ninbin-1 DO BEGIN
      xsub=binind[inbin[j]]
      dist=SQRT((x[xsub]-xbar[notbin])^2+(y[xsub]-ybar[notbin])^2)
      mindist=MIN(dist,closest)
      print,mindist,x[xsub],xbar[notbin[closest]],y[xsub],ybar[notbin[closest]]

      binnum[inbin[j]]=ibin[notbin[closest]]
      innewbin=binind[WHERE(binnum EQ ibin[notbin[closest]],npixnew)]      
      xbar[notbin[closest]]=TOTAL(signal[innewbin]*x[innewbin])/TOTAL(signal[innewbin])
      ybar[notbin[closest]]=TOTAL(signal[innewbin]*y[innewbin])/TOTAL(signal[innewbin])
      npixels[notbin[closest]]=npixnew
   ENDFOR
   xbar=xbar[notbin]
   ybar=ybar[notbin]
   ibin=ibin[notbin]
   npixels=npixels[notbin]
ENDFOR
   
xyouts,0,0,'!6'

allbinnum=INTARR(size[0],size[1])
allbinnum[*,*]=999;!VALUES.F_NAN
allbinnum[binind]=binNum
curbin=MAX(binnum)+1


;now do radial section bins
FOR i=0,nradbin-1 DO BEGIN
FOR j=0,nangbin-1 DO BEGIN
   inbin=WHERE(rellip GT radbin[i] AND rellip LE radbin[i+1] AND ang_rot GE angbin[j] AND ang_rot LT angbin[j+1],ninbin)
   print,radbin[i],angbin[j],curbin,ninbin
   xbar=[xbar,TOTAL(cont[inbin]*x[inbin])/TOTAL(cont[inbin])]
   ybar=[ybar,TOTAL(cont[inbin]*y[inbin])/TOTAL(cont[inbin])]
   npixels=[npixels,ninbin]
   ibin=[ibin,curbin]
   sn=[sn,TOTAL(signal[inbin])/SQRT(TOTAL(noise[inbin]^2))]
   allbinnum[inbin]=curbin
   curbin=curbin+1
ENDFOR
ENDFOR


xoutbin=xbar;(xbar-xcen)*0.05
youtbin=ybar;(ybar-ycen)*0.05
binnum=allbinnum
WRITEFITS,'test.fits',binnum


;bin the spectra together
nbins=N_ELEMENTS(xbar)
bin_spectra=FLTARR(nbins,size[2])
bin_var=FLTARR(nbins,size[2])
bin_fwhm=FLTARR(nbins)
measured_sn=FLTARR(nbins)
tsn=FLTARR(nbins)
snind=WHERE(lambda GT 22140. AND lambda LT 22600.)
tsnind=WHERE(lambda GT 22950. AND lambda LT 23950.)
ind=WHERE(fwhm GT 0.0)
medianfwhm=MEDIAN(fwhm[ind])
print,'MEDIAN FWHM',MEDIAN(fwhm)
FOR i=0,nbins-1 DO BEGIN
    ind=WHERE(binnum EQ ibin[i],nind)
    xind=ind MOD size[0]
    yind=ind / size[0]
    tempspec=FLTARR(size[2])
    tempvar=FLTARR(size[2])
    tempfwhm=0.
    FOR j=0,nind-1 DO BEGIN
        tempspec=tempspec+ifu[xind[j],yind[j],*]
        tempvar=tempvar+var[xind[j],yind[j],*]
        IF (tempfwhm GT 0.0) THEN tempfwhm=tempfwhm+fwhm[xind[j],yind[j],*] ELSE tempfwhm=tempfwhm+medianfwhm
                                   
    ENDFOR
    bin_spectra[i,*]=tempspec;/FLOAT(nind)
    bin_var[i,*]=tempvar
    bin_fwhm[i]=tempfwhm/FLOAT(nind)
    MEANCLIP,tempspec[snind],meany,4.,subs=subs
    measured_sn[i]=MEAN(tempspec[snind[subs]])/STDDEV(tempspec[snind[subs]])
    tsn[i]=MEDIAN(tempspec[tsnind]/SQRT(tempvar[tsnind]))
    print,ibin[i],npixels[i],nind,xbar[i],ybar[i],measured_sn[i],tsn[i]
    plot,lambda,tempspec,xrange=[2.25e4,2.4e4],title=STRTRIM(i,2)
;    print,xoutbin[i],youtbin[i],tsignoise/sn[i],nPixels[i],FORMAT='(3F7.2,I4)'
;    titstring=STRING(xoutbin[i],FORMAT='(F5.2)')+','+STRING(youtbin[i],FORMAT='(F5.2)')
;    plot,lambda,SMOOTH(bin_spectra[i,*],4),title=titstring,xrange=[2.29e4,2.4e4]
;    junk=GET_KBRD()
ENDFOR

STOP

binimage=FLTARR(size[0:1])
velbin=FLTARR(nbins)
velimage=FLTARR(size[0:1])
dispbin=FLTARR(nbins)
dispimage=FLTARR(size[0:1])
coimage=FLTARR(size[0:1])
caimage=FLTARR(size[0:1])
mgimage=FLTARR(size[0:1])
naimage=FLTARR(size[0:1])
brgimage=FLTARR(size[0:1])
snimage=FLTARR(size[0:1])
spimage=FLTARR(size[0:1])
spbin=FLTARR(nbins)
spdispbin=FLTARR(nbins)
lumimage=FLTARR(size[0:1])
lumbin=FLTARR(nbins)
chiimage=FLTARR(size[0:1])
chibin=FLTARR(nbins)
outfile='vor_out/newbin_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'.dat' 
!P.MULTI=[0,1,2]
OPENW,1,outfile,WIDTH=1000
FOR i=0,nbins-1 DO BEGIN
    xyouts,0,0,'!6'
    medspec=MEDIAN(bin_spectra[i,*])
    badspec=WHERE(REFORM(bin_spectra[i,*]) LT medspec/2. OR REFORM(bin_spectra[i,*]) GT medspec*2. OR REFORM(bin_spectra[i,*]) EQ 0.0,nbadspec)
    ind=WHERE(binnum EQ ibin[i])    
    IF (nbadspec GT 100) THEN BEGIN
        binimage[ind]=i
        velbin[i]=-999.
        dispbin[i]=-999.
        mgimage[ind]=-999.
        brgimage[ind]=-999.
        naimage[ind]=-999.
        caimage[ind]=-999.
        coimage[ind]=-999
        plot,lambda,REFORM(bin_spectra[i,*])
        plot,lambda,REFORM(bin_spectra[i,*])
        printf,1,ibin[i],xoutbin[i],youtbin[i],npixels[i],measured_sn[i],sn[i],REPLICATE(-999.,20),FORMAT='(I4,2F7.3,I5,22F9.3)'
        GOTO,SKIP
    ENDIF ELSE BEGIN
    ppxf_wallace,lambda,REFORM(bin_spectra[i,*]),REFORM(bin_var[i,*]),bin_fwhm[i],outrv=vel,outdisp=disp,outsol=sol,binnum=ibin[i],initvel=initvel,nodispchi=nodispchi,outerror=err,/doerror
    spimage[ind]=sol[7]
    spbin[i]=sol[7]
    spdispbin[i]=sol[8]
    lumimage[ind]=sol[9]
    lumbin[i]=sol[9]
    chiimage[ind]=sol[6]
    chibin[i]=sol[6]
        printf,1,ibin[i],xoutbin[i],youtbin[i],npixels[i],measured_sn[i],sn[i],sol[6],sol[0],err[0],err[1],sol[1],err[2],err[3],sol[2],err[4],err[5],sol[3],err[6],err[7],sol[7],err[8],err[9],sol[9],err[10],err[11],nodispchi,FORMAT='(I4,2F7.3,I5,22F9.3)'
;       printf,1,i,xoutbin[i],youtbin[i],npixels[i],measured_sn[i],sn[i],sol[6],sol[0],sol[1],sol[2],sol[3],sol[7],sol[9],nodispchi,FORMAT='(I4,2F7.3,I5,10F9.3)'
    binimage[ind]=ibin[i]
    velbin[i]=vel
    dispbin[i]=disp
    snimage[ind]=measured_sn[i]
    kind=kband_indices(lambda,bin_spectra[i,*],velocity=vel)
;    IF (sn[i] GE 0.8*min_sn AND npixels[i] LT 90.) THEN BEGIN
    IF (sn[i] GE 0.1*min_sn) THEN BEGIN
        velimage[ind]=velbin[i]
        dispimage[ind]=dispbin[i]
;        indices=['MgI','BrG','NaI','CaI','CO20']
        mgimage[ind]=kind[0]
        brgimage[ind]=kind[1]
        naimage[ind]=kind[2]
        caimage[ind]=kind[3]
        coimage[ind]=kind[4]
    ENDIF ELSE BEGIN
        velimage[ind]=-999.
        dispimage[ind]=-999.
        mgimage[ind]=-999.
        brgimage[ind]=-999.
        naimage[ind]=-999.
        caimage[ind]=-999.
        coimage[ind]=-999.
    ENDELSE
    print,xoutbin[i],youtbin[i],nPixels[i],spbin[i],spdispbin[i],lumbin[i],FORMAT='(2F7.2,I4,3F7.1)'
    print,''
    IF KEYWORD_SET(CHECK) THEN a=get_kbrd(1)
ENDELSE
    SKIP:
ENDFOR
CLOSE,1

WRITEFITS,'vor_out/newbin_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'_mg.fits',mgimage
WRITEFITS,'vor_out/newbin_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'_brg.fits',brgimage
WRITEFITS,'vor_out/newbin_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'_na.fits',naimage
WRITEFITS,'vor_out/newbin_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'_ca.fits',caimage
WRITEFITS,'vor_out/newbin_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'_co.fits',coimage
WRITEFITS,'vor_out/newbin_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'_sn.fits',snimage
WRITEFITS,'vor_out/newbin_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'_sp.fits',spimage
WRITEFITS,'vor_out/newbin_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'_lum.fits',lumimage

binfile='vor_out/newbin_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'_bin.fits'
WRITEFITS,binfile,binimage    
velfile='vor_out/vel_newbin_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'.fits'
WRITEFITS,velfile,velimage
dispfile='vor_out/disp_newbin_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'.fits'
WRITEFITS,dispfile,dispimage


vsys=velimage[FIX(xcen+0.5),FIX(ycen+0.5)]
velbin=velbin-vsys
velimage=velimage-vsys
;loadct,13
;plot_velfield,xbar,ybar,velbin,range=[-50,50]

cont=signal

DEVICE,/close
!P.MULTI=[0,1,1]

STOP

outpsfile='vor_out/voronoi_wallace_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'.ps'
myplot,file=outpsfile,xsize=8,ysize=8,/inches,bits=8,xoffset=0.25,yoffset=1.5
loadct,0
xcoord=(FINDGEN(size[0])-xcen)*0.05
ycoord=(FINDGEN(size[1])-ycen)*0.05

;first do velocities
minlev=-40 & maxlev=40
nlev=100.
dellev=(maxlev-minlev)/nlev
vellevels=findgen(nlev)*dellev+minlev
velcontimage=velimage
;bad=WHERE(velimage LT -200.)
;velcontimage[bad]=!VALUE.F_NAN
low=WHERE(velimage GT -200. AND velimage LE minlev,nlow)
IF (nlow GT 0) THEN velcontimage[low]=minlev
high=WHERE(velimage GT maxlev,nhigh)
IF (nhigh GT 0) THEN velcontimage[high]=maxlev
CONTOUR, signal,xcoord,ycoord, levels=levels,/fill,/xs,/ys,/nodata,charsize=1.8,xtitle='X offset ["]',ytitle='Y offset ["]',xminor=2,/iso,xmargin=[6,8]
LOADCT, 13      
CONTOUR, velcontimage,xcoord,ycoord,levels=vellevels,/fill,/overplot,c_colors=findgen(n_elements(vellevels))*256./n_elements(vellevels)
logcont=alog10(cont)
levels=ALOG10(MAX(cont)*[0.05,0.1,0.5])
contour,logcont,xcoord,ycoord,/overplot,levels=levels,thick=4
COLORBAR,range=[min(vellevels),max(vellevels)+dellev],format='(I4.0)',/vertical,right=1,divisions=5,charsize=1.8,title='km/sec',POSITION=[!X.WINDOW[1],!Y.WINDOW[0],!X.WINDOW[1]+0.05,!Y.WINDOW[1]]
;tvcircle,0.1,-1.3,1.65,0,/data,/fill
;xyouts,-1.15,1.6,'Approx. FWHM (0.2")',charsize=1.8

;now do dispersion
minlev=0 & maxlev=50
nlev=50.
dellev=(maxlev-minlev)/nlev
vellevels=findgen(nlev)*dellev+minlev
CONTOUR, signal,xcoord,ycoord, levels=levels,/fill,/xs,/ys,/nodata,charsize=1.8,xtitle='X offset ["]',ytitle='Y offset ["]',xminor=2,/iso,xmargin=[6,8]
LOADCT, 13      
CONTOUR, dispimage,xcoord,ycoord,levels=vellevels,/fill,/overplot,c_colors=findgen(n_elements(vellevels))*256./n_elements(vellevels)
logcont=alog10(cont)
levels=ALOG10(MAX(cont)*[0.05,0.1,0.5])
contour,logcont,xcoord,ycoord,/overplot,levels=levels,thick=3
COLORBAR,range=[min(vellevels),max(vellevels)+dellev],format='(I4.0)',/vertical,right=1,divisions=5,charsize=1.8,title='km/sec',POSITION=[!X.WINDOW[1],!Y.WINDOW[0],!X.WINDOW[1]+0.05,!Y.WINDOW[1]]
tvcircle,0.1,-1.3,1.65,0,/data,/fill
xyouts,-1.15,1.6,'Approx. FWHM (0.2")',charsize=1.8
loadct,0
device,/close
set_plot,'x'

STOP
END


