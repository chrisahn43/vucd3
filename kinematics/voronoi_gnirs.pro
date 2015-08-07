@ppxf_gnirs.pro
PRO VORONOI_GNIRS,MIN_SN=min_sn

IF NOT (KEYWORD_SET(min_sn)) THEN min_sn=25.

infile='ngc404_combine_near9.fits'
fwhmfile='../sky/ngc404/fullcube_disp_med.fits'
initvel=-61 ;NED value

ifu=READFITS(infile,head,ext=1)
var=READFITS(infile,ext=2)
fwhm=READFITS(fwhmfile)
tag='newfwhm'
size=SIZE(ifu,/DIM)
lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
lambda=FINDGEN(size[2])*dlambda+lambda0
makex,ifu,x,y


outplotfile='vor_out/gnirs_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'.ps'
myplot_nott,file=outplotfile,ysize=9,yoff=1,/inches
!P.MULTI=[0,1,2]

;calculate signal to noise
minfit=2.28e4  ;part of the target spectra to extract
maxfit=2.395e4
snind=WHERE(lambda GT minfit AND lambda LT maxfit)
signal=FLTARR(size[0:1])
noise=FLTARR(size[0:1])
FOR i=0,size[0]-1 DO BEGIN
FOR j=0,size[1]-1 DO BEGIN
    uspec=ifu[i,j,snind]
    uvar=var[i,j,snind]
    signal[i,j]=MEDIAN(uspec)
    noise[i,j]=MEDIAN(SQRT(uvar))
ENDFOR
ENDFOR

signoise=signal/noise
WRITEFITS,'vor_out/sn.fits',signoise
logsn=ALOG10(signoise)
logsn=logsn-MIN(logsn[WHERE(FINITE(logsn) EQ 1,comp=bad)])
logsn[bad]=0.0
plot,FINDGEN(10),/nodata,/iso
loadct,13,/silent
imgunder,logsn*255./MAX(logsn)
loadct,0,/silent

voronoi_2d_binning, x, y, signal, noise, min_sn, binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale, /PLOT, /QUIET,/NO_CVT ;-- S/N=10 failed without this

imsize=SIZE(signal,/dim)
maxval=MAX(signal,maxpos)
xinit=maxpos MOD imsize[0]
yinit=maxpos/imsize[0]
GCNTRD,signal,xinit,yinit,xcen,ycen,4.5
print,'Center',xcen,ycen
xoutbin=(xbar-xcen)*0.05
youtbin=(ybar-ycen)*0.05

;bin the spectra together
nbins=N_ELEMENTS(xnode)
bin_spectra=FLTARR(nbins,size[2])
bin_var=FLTARR(nbins,size[2])
bin_fwhm=FLTARR(nbins)
measured_sn=FLTARR(nbins)
snind=WHERE(lambda GT 22140. AND lambda LT 22600.)
FOR i=0,nbins-1 DO BEGIN
    ind=WHERE(binnum EQ i,nind)
    xind=ind MOD size[0]
    yind=ind / size[0]
    tempspec=FLTARR(size[2])
    tempvar=FLTARR(size[2])
    tempfwhm=0.
    FOR j=0,nind-1 DO BEGIN
        tempspec=tempspec+ifu[xind[j],yind[j],*]
        tempvar=tempvar+var[xind[j],yind[j],*]
        tempfwhm=tempfwhm+fwhm[xind[j],yind[j],*]
    ENDFOR
    bin_spectra[i,*]=tempspec;/FLOAT(nind)
    bin_var[i,*]=tempvar
    bin_fwhm[i]=tempfwhm/FLOAT(nind)
    MEANCLIP,tempspec[snind],meany,4.,subs=subs
    measured_sn[i]=MEAN(tempspec[snind[subs]])/STDDEV(tempspec[snind[subs]])
;    print,xoutbin[i],youtbin[i],tsignoise/sn[i],nPixels[i],FORMAT='(3F7.2,I4)'
;    titstring=STRING(xoutbin[i],FORMAT='(F5.2)')+','+STRING(youtbin[i],FORMAT='(F5.2)')
;    plot,lambda,SMOOTH(bin_spectra[i,*],4),title=titstring,xrange=[2.29e4,2.4e4]
;    junk=GET_KBRD()
ENDFOR

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
outfile='vor_out/gnirs_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'.dat' 

!P.MULTI=[0,1,2]
OPENW,1,outfile,WIDTH=1000
FOR i=0,nbins-1 DO BEGIN
    ind=WHERE(binnum EQ i)
    ppxf_gnirs,lambda,REFORM(bin_spectra[i,*]),REFORM(bin_var[i,*]),bin_fwhm[i],outrv=vel,outdisp=disp,outsol=sol,binnum=i,initvel=initvel,nodispchi=nodispchi,outerror=err,/doerror
    spimage[ind]=sol[7]
    spbin[i]=sol[7]
    spdispbin[i]=sol[8]
    lumimage[ind]=sol[9]
    lumbin[i]=sol[9]
    chiimage[ind]=sol[6]
    chibin[i]=sol[6]
        printf,1,i,xoutbin[i],youtbin[i],npixels[i],measured_sn[i],sn[i],sol[6],sol[0],err[0],err[1],sol[1],err[2],err[3],sol[2],err[4],err[5],sol[3],err[6],err[7],sol[7],err[8],err[9],sol[9],err[10],err[11],nodispchi,FORMAT='(I4,2F7.3,I5,22F9.3)'
;       printf,1,i,xoutbin[i],youtbin[i],npixels[i],measured_sn[i],sn[i],sol[6],sol[0],sol[1],sol[2],sol[3],sol[7],sol[9],nodispchi,FORMAT='(I4,2F7.3,I5,10F9.3)'
    binimage[ind]=i
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

ENDFOR
CLOSE,1

;WRITEFITS,'vor_out/voronoi_'+tag+'_mg.fits',mgimage
;WRITEFITS,'vor_out/voronoi_'+tag+'_brg.fits',brgimage
;WRITEFITS,'vor_out/voronoi_'+tag+'_na.fits',naimage
;WRITEFITS,'vor_out/voronoi_'+tag+'_ca.fits',caimage
;WRITEFITS,'vor_out/voronoi_'+tag+'_co.fits',coimage
;WRITEFITS,'vor_out/voronoi_'+tag+'_sn.fits',snimage
;WRITEFITS,'vor_out/voronoi_'+tag+'_sp.fits',spimage
;WRITEFITS,'vor_out/voronoi_'+tag+'_lum.fits',lumimage

binfile='vor_out/gnirs_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'_bin.fits'
WRITEFITS,binfile,binimage    
velfile='vor_out/vel_gnirs_sn'+STRTRIM(FIX(min_sn),2)+'.fits'
WRITEFITS,velfile,velimage
dispfile='vor_out/disp_gnirs_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'.fits'
WRITEFITS,dispfile,dispimage


vsys=velimage[FIX(xcen+0.5),FIX(ycen+0.5)]
velbin=velbin-vsys
velimage=velimage-vsys
loadct,13
plot_velfield,xbar,ybar,velbin,range=[-50,50]

cont=signal

DEVICE,/close
!P.MULTI=[0,1,1]

STOP

outpsfile='vor_out/voronoi_gnirs_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'.ps'
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
;velcontimage[bad]=!VALUES.F_NAN
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

