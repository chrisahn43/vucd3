@ppxf_gnirs_test.pro
PRO VORONOI_NOOUTPUT,GNIRS=gnirs,min_sn=min_sn,CHECK=check

GNIRS=1
IF NOT (KEYWORD_SET(min_sn)) THEN min_sn=45.

infile='ngc404_combine_near9.fits'
fwhmfile='../sky/ngc404/fullcube_disp_med.fits'
initvel=-74.32 ;measured central velocity GNIRS

ifu=READFITS(infile,head,ext=1)
var=READFITS(infile,ext=2)
fwhm=READFITS(fwhmfile)
tag='fwhmcor'
size=SIZE(ifu,/DIM)
lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
lambda=FINDGEN(size[2])*dlambda+lambda0
makex,ifu,x,y


;outplotfile='vor_out/gnirs_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'.ps'
;myplot,file=outplotfile,ysize=9,yoff=1,/inches
!P.MULTI=[0,1,2]

;calculate signal to noise
minfit=2.28e4  ;part of the target spectra to extract
maxfit=2.41e4
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
lowdispchibin=FLTARR(nbins)
ntemplate=28
tempsolbin=FLTARR(7,ntemplate,nbins)
IF KEYWORD_SET(GNIRS) THEN $
  outfile='vor_out/gnirs_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'.dat' ELSE $
  outfile='vor_out/nifs_'+tag+'_sn'+STRTRIM(FIX(min_sn),2)+'.dat'

!P.MULTI=[0,1,2]
FOR i=0,nbins-1 DO BEGIN
    ind=WHERE(binnum EQ i)
    IF (KEYWORD_SET(GNIRS)) THEN BEGIN
        ppxf_gnirs_test,lambda,REFORM(bin_spectra[i,*]),REFORM(bin_var[i,*]),bin_fwhm[i],outrv=vel,outdisp=disp,outsol=sol,binnum=i,initvel=initvel,lowdispchi=lowdispchi,tempsol=tempsol;,outerror=err,/doerror
        spimage[ind]=sol[7]
        spbin[i]=sol[7]
        spdispbin[i]=sol[8]
        lumimage[ind]=sol[9]
        lumbin[i]=sol[9]
        chiimage[ind]=sol[6]
        chibin[i]=sol[6]
        lowdispchibin[i]=lowdispchi
        tempsolbin[*,*,i]=tempsol
;        print,1,i,xoutbin[i],youtbin[i],npixels[i],measured_sn[i],sn[i],sol[6],sol[0],err[0],err[1],sol[1],err[2],err[3],sol[2],err[4],err[5],sol[3],err[6],err[7],sol[7],err[8],err[9],sol[9],err[10],err[11],FORMAT='(I4,2F7.3,I5,21F8.2)'
       print,i,xoutbin[i],youtbin[i],npixels[i],measured_sn[i],sn[i],sol[6],sol[0],sol[1],sol[2],sol[3],sol[7],sol[9],FORMAT='(I4,2F7.3,I5,7F8.2)'
    ENDIF ELSE BEGIN
        ppxf_nifs,lambda,REFORM(bin_spectra[i,*]),REFORM(bin_var[i,*]),outrv=vel,outdisp=disp,outsol=sol,outerror=err,/doerror
        chiimage[ind]=sol[6]
        chibin[i]=sol[6]
        printf,1,i,xoutbin[i],youtbin[i],npixels[i],measured_sn[i],sn[i],sol[6],sol[0],err[0],err[1],sol[1],err[2],err[3],sol[2],err[4],err[5],sol[3],err[6],err[7],FORMAT='(I4,2F7.3,I5,15F8.2)'
;         printf,1,i,xoutbin[i],youtbin[i],npixels[i],measured_sn[i],sn[i],sol[6],sol[0],sol[1],sol[2],sol[3],FORMAT='(I4,2F7.3,I5,5F8.2)'

    ENDELSE
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

tempbinnum=REBIN(INDGEN(ntemplate),ntemplate,nbins)
subvelarr=REBIN(REFORM(velbin,1,nbins),ntemplate,nbins)
tempvelbin=tempsolbin[0,*,*]
veldiff=REFORM(tempvelbin-subvelarr)
subdisparr=REBIN(REFORM(dispbin,1,nbins),ntemplate,nbins)
tempdispbin=tempsolbin[1,*,*]
dispdiff=REFORM(tempdispbin-subdisparr)
tempchibin=tempsolbin[6,*,*]
snarr=REBIN(REFORM(measured_sn,1,nbins),ntemplate,nbins)

meanveldiff=FLTARR(ntemplate)
sigmaveldiff=FLTARR(ntemplate)
meandispdiff=FLTARR(ntemplate)
sigmadispdiff=FLTARR(ntemplate)
chidiff=FLTARR(ntemplate)
FOR i=0,ntemplate-1 DO BEGIN
    ind=WHERE(tempbinnum EQ i AND snarr GT 30.)
    MEANCLIP,veldiff[ind],meanvel,sigmavel
    meanveldiff[i]=meanvel & sigmaveldiff[i]=sigmavel
    MEANCLIP,dispdiff[ind],meandisp,sigmadisp
    meandispdiff[i]=meandisp & sigmadispdiff[i]=sigmadisp
    chidiff[i]=MEDIAN(tempchibin[ind])
ENDFOR
plot,INDGEN(ntemplate),meanveldiff,psym=4,symsize=2,xrange=[-1,ntemplate],/xsty
oploterror,INDGEN(ntemplate),meanveldiff,sigmaveldiff,psym=4,symsize=2
plots,[-1,ntemplate],[0,0]  

plot,INDGEN(ntemplate),meandispdiff,psym=4,symsize=2,xrange=[-1,ntemplate],/xsty
oploterror,INDGEN(ntemplate),meandispdiff,sigmadispdiff,psym=4,symsize=2
plots,[-1,ntemplate],[0,0]  


plot,INDGEN(ntemplate),veldiff[*,0],xrange=[-1,ntemplate],/xsty,yrange=[-20,20]
ind=WHERE(measured_sn GT 40.,nind)
FOR i=1,nind-1 DO oplot,INDGEN(ntemplate),veldiff[*,[ind[i]]]

READCOL,'../winge/winge_rvcor_new.dat',wingeroot,rv,spclass,lumclass,FORMAT='A,F,F,I',/SILENT

OPENW,1,'winge_template_systematics.dat'
printf,1,'Root lumclass       dV  dVerr   dDisp  dDisperr  Chi'
FOR i=0,ntemplate-1 DO printf,1,wingeroot[i],lumclass[i],meanveldiff[i],sigmaveldiff[i],meandispdiff[i],sigmadispdiff[i],chidiff[i],FORMAT='(A12,I2,5F8.3)'
CLOSE,1
FREE_LUN,1    

!P.MULTI=[0,1,2]


STOP
END

