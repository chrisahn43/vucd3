;------------------------------------------------------------------------------
pro ppxf_gnirs_long,lambda,spec,var,fwhm,OUTRV=outrv,OUTDISP=outdisp,OUTSOL=sol,DOERROR=doerror,OUTERROR=outerror,BINNUM=BINNUM,INITVEL=initvel,nodispchi=nodispchi,BESTFIT=bestfit,loglamout=loglamout

loadct,12,/silent
minfit=2.285e4  ;part of the target spectra to extract
maxfit=2.42e4
mintemp=2.28e4 ;part of the template to use
maxtemp=2.425e4
velscale=10. ;for the data
velscaletemp=1. ;for the templates
tempbin=FIX(velscale/velscaletemp)


;prepare the galaxy data
ind=WHERE(lambda GT minfit AND lambda LT maxfit,nind)
gal_lin=spec[ind]
ini_error=SQRT(var[ind])
bad=WHERE(ini_error EQ 0.0,nbad)
IF (nbad GT 0) THEN ini_error[bad]=MEDIAN(ini_error)
lamRange = [lambda[ind[0]],lambda[ind[nind-1]]]
log_rebin, lamRange, gal_lin, galaxy, logLam1, VELSCALE=velScale
log_rebin, lamRange, ini_error^2, errorsquared,  VELSCALE=velScale
error=SQRT(errorsquared)

;now prepare the templates
wingesave='winge_templates_long.idl'
IF (FILE_TEST(wingesave) EQ 0) THEN BEGIN
;READCOL,'../winge/winge_table.dat',wingeroot,rv,FORMAT='A,X,F',/SILENT
;winge='../winge/'+STRLOWCASE(wingeroot)+'_v1r.fits'
READCOL,'../winge/winge_rvcor_final.dat',wingeroot,rv,spclass,lumclass,FORMAT='A,F,F,F'
winge='../winge/'+STRLOWCASE(wingeroot)
nfiles=N_ELEMENTS(wingeroot)
wingeflux=READFITS(winge[0],head)
wingelam=SXPAR(head,'CRVAL1')+FINDGEN(N_ELEMENTS(wingeflux))*SXPAR(head,'CD1_1')
ind=WHERE(wingelam GT mintemp AND wingelam LT maxtemp,nind)
ssp=wingeflux[ind]
lamRange = [wingelam[ind[0]],wingelam[ind[nind-1]]]
reflamrange=lamrange
log_rebin, lamRange, ssp, sspNew, logLam2, VELSCALE=velScaletemp
templates = dblarr(n_elements(sspNew),nfiles)
for j=0,nfiles-1 do begin
    wingeflux=READFITS(winge[j],head)
    wingelam=SXPAR(head,'CRVAL1')+FINDGEN(N_ELEMENTS(wingeflux))*SXPAR(head,'CD1_1')    
    ind=WHERE(wingelam GT mintemp AND wingelam LT maxtemp,nind)
    ssp=wingeflux[ind]
    lamRange = [wingelam[ind[0]],wingelam[ind[nind-1]]]
    log_rebin, lamRange, ssp, ssplogged, onelogLam, VELSCALE=velScaletemp
    sspnew=RESAMPLE_SPECTRA(ssplogged,oneloglam,loglam2,oversamp=2)
;    print,MIN(oneloglam),MIN(loglam2)
;    npix=N_ELEMENTS(sspNew)
;    inpix=FINDGEN(npix)
;    outpix=FINDGEN(npix)+rv[j]/velScale
;    sspvcor=RESAMPLE_SPECTRA(sspNew,inpix,outpix)
    sspvcor=SHIFT(sspnew,-(FIX(rv[j])+0.5))
    templates[*,j] = sspvcor/median(sspvcor)
;    templates[*,j] = convol(sspvcor,lsf)/median(sspvcor) ; Degrade template resolution
    IF (j EQ 0) THEN plot,sspvcor/MEDIAN(sspvcor),xrange=[2000,3000],yrange=[0.5,1.1] ELSE oplot,sspvcor/MEDIAN(sspvcor)
    print,wingeroot[j]
;    oplot,sspvcor,color=100
endfor
!P.MULTI=[0,1,2]
save,templates,logLam2,spclass,lumclass,file=wingesave
ENDIF ELSE RESTORE, file=wingesave

;now correct templates for difference in resolution
; Winge is 1.84? Angstrom pixels, FWHM of 1.82 pixels -- 3.35 Ang FWHM
;   at 2.3um, this is 43.5 km/sec
; FWHM of NIFS is 4.1 Angstroms, 53.3 km/sec
; HOWEVER, using templates the difference is ~13 km/sec.  
; DUH! this is because the FWHM=2.35*sigma.  The difference in sigmas
; is SQRT(22.7^2-18.5^2)=13.1 == the derived difference from the templates.
c=2.99792458d5
sigma2fwhm=2*SQRT(2*ALOG(2))
meanlambda=(maxfit+minfit)/2.
nifssigma=(fwhm/sigma2fwhm)/meanlambda*c ;in km/sec
; Winge is 1.84? Angstrom pixels, FWHM of 1.82 pixels -- 3.35 Ang FWHM
;   at 2.3um, this is 43.5 km/sec
gnirssigma=19.5;(1.84*1.82/sigma2fwhm)/meanlambda*c ;in km/sec
templatesize=SIZE(templates,/dim)
FOR i=0,templatesize[1]-1 DO BEGIN
    IF (nifssigma GT gnirssigma) THEN BEGIN
        sigma = SQRT(nifssigma^2-gnirssigma^2)/velScaletemp ; Quadratic sigma difference in pixels 
        lsf = psf_Gaussian(NPIXEL=(8*sigma > 5.), ST_DEV=sigma, /NORM, NDIM=1)
        convtemp = convol(templates[*,i],lsf) ; Degrade template resolution
    ENDIF ELSE convtemp=templates[*,i]
    outpix=N_ELEMENTS(convtemp) / tempbin
    inmaxpix=outpix*tempbin-1
    rebintemp=REBIN(convtemp[0:inmaxpix],outpix)
    loglamrebin=REBIN(loglam2[0:inmaxpix],outpix)
    IF (i EQ 0) THEN convtemplates=DBLARR(outpix,templatesize[1])
    convtemplates[*,i] = rebintemp
ENDFOR
 
; The galaxy and the template spectra do not have the same starting wavelength.
; For this reason an extra velocity shift DV has to be applied to the template
; to fit the galaxy spectrum. We remove this artificial shift by using the
; keyword VSYST in the call to PPXF below, so that all velocities are
; measured with respect to DV.
;
dv = (logLam2[0]-logLam1[0])*299792.458d ; km/s

; Here the actual fit starts. The best fit is plotted on the screen.
; Gas emission lines are excluded from the pPXF fit using the GOODPIXELS keyword.
;
;templates=templates[*,20]
start = [initvel, 3.*velscale] ; (km/s), starting guess for [V,sigma]
lambdalin=EXP(loglam1)
ind=WHERE((lambdalin GT 24050. AND lambdalin LT 24070.) OR (lambdalin GT 23205. AND lambdalin LT 23225.),comp=goodpix)
mdegree=10
ppxf, convtemplates, galaxy, error, velScale, start, sol, $
    /PLOT, MOMENTS=4, MDEGREE=mdegree, VSYST=dv,WEIGHTS=weights,goodpix=goodpix,bestfit=bestfit;,/clean;,/OVERSAMPLE
x = range(-1d,1d,n_elements(galaxy))
mpoly = 1d                      ; Multiplicative polynomial
for j=1,MDEGREE do mpoly += legendre(x,j)*sol[6+j]
oplot,mpoly*MEDIAN(bestfit),color=150


check=[sol[0],velscale/10.0]
ppxf, convtemplates, galaxy, error, velScale, check, checksol, MOMENTS=0, DEGREE=4, VSYST=dv,WEIGHTS=weights,goodpix=goodpix,/quiet;,/plot
nodispchi=checksol[6]
meanspec=TOTAL(weights*spclass)/TOTAL(weights)
sdevspec=SQRT(TOTAL(weights*(spclass-meanspec)^2)/TOTAL(weights))
meanlum=TOTAL(weights*lumclass)/TOTAL(weights)
sdevlum=SQRT(TOTAL(weights*(lumclass-meanlum)^2)/TOTAL(weights))
sol=[sol,meanspec,sdevspec,meanlum,sdevlum]
outrv=sol[0]
outdisp=sol[1]
loadct,0,/silent
loglamout=loglam1

IF KEYWORD_SET(DOERROR) THEN BEGIN
    nruns=100
    npixels=N_ELEMENTS(galaxy)
    allsol=FLTARR(11,nruns)
    FOR i=0,nruns-1 DO BEGIN
        newgal=galaxy+error*RANDOMN(seed,npixels)
        ppxf, convtemplates, newgal, error, velScale, start, onesol, $
          MOMENTS=4, DEGREE=10, VSYST=dv,bestfit=bestfit,bias=0.1,WEIGHTS=weights,goodpix=goodpix,/quiet
        meanspec=TOTAL(weights*spclass)/TOTAL(weights)
        sdevspec=SQRT(TOTAL(weights*(spclass-meanspec)^2)/TOTAL(weights))
        meanlum=TOTAL(weights*lumclass)/TOTAL(weights)
        sdevlum=SQRT(TOTAL(weights*(lumclass-meanlum)^2)/TOTAL(weights))
        onesol=[onesol,meanspec,sdevspec,meanlum,sdevlum]
        allsol[*,i]=onesol
    ENDFOR
    titstring='Bin '+STRTRIM(binnum,2)+' V='+STRING(sol[0],FORMAT='(F7.1)')+' !9s!x='+STRING(sol[1],FORMAT='(F7.1)')
    plothist,allsol[1,*],title=titstring,xtitle='!9s!x [km/sec]',ytitle='# of MC runs',bin=0.5
    print, MEAN(allsol[1,*]), STDDEV(allsol[1,*])
    outerror=[MEAN(allsol[0,*]),STDDEV(allsol[0,*]),MEAN(allsol[1,*]),STDDEV(allsol[1,*]),MEAN(allsol[2,*]),STDDEV(allsol[2,*]),MEAN(allsol[3,*]),STDDEV(allsol[3,*]),MEAN(allsol[7,*]),STDDEV(allsol[7,*]),MEAN(allsol[9,*]),STDDEV(allsol[9,*])]
ENDIF    


end
;------------------------------------------------------------------------------
