;------------------------------------------------------------------------------
pro PPXF_PHOENIX,lambda,spec,var,fwhm,OUTRV=outrv,OUTDISP=outdisp,OUTSOL=sol,DOERROR=doerror,OUTERROR=outerror,BINNUM=BINNUM,INITVEL=initvel,nodispchi=nodispchi,QUIET=QUIET,BESTFIT=BESTFIT,LOGLAMOUT=LOGLAMOUT,GALAXY=galaxy,RESIDSN=RESIDSN

loadct,12,/silent
minfit=2.295e4  ;part of the target spectra to extract
maxfit=2.395e4
mintemp=2.28e4 ;part of the template to use
maxtemp=2.40e4
velscale=25 ;for the data
velscaletemp=1. ;for the templates
tempbin=FIX(velscale/velscaletemp)

ind=WHERE(lambda GT minfit AND lambda LT maxfit,nind)
gal_lin=spec[ind]
ini_error=SQRT(var[ind])
bad=WHERE(ini_error EQ 0.0,nbad)
IF (nbad GT 0) THEN ini_error[bad]=MEDIAN(ini_error)
lamRange = [lambda[ind[0]],lambda[ind[nind-1]]]
log_rebin, lamRange, gal_lin, galaxy, logLam1, VELSCALE=velScale
log_rebin, lamRange, ini_error^2, errorsquared,  VELSCALE=velScale
error=SQRT(errorsquared)


phoenixsave='phoenix_templates.idl'
IF (FILE_TEST(phoenixsave) EQ 0) THEN BEGIN
READCOL,'./phoenix/templates.dat',phoenixfile,rv,spclass,lumclass,FORMAT='A,I,F,F',/SILENT
phoenix='./phoenix/'+STRLOWCASE(phoenixfile)
goodtemp=[0,1,2,3];phoenix templates to use
nfiles=N_ELEMENTS(goodtemp)
phoenixflux=READFITS(phoenix[0],head)
phoenixlam=READFITS('./phoenix/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
ind=WHERE(phoenixlam GT mintemp AND phoenixlam LT maxtemp,nind)
ssp=phoenixflux[ind]
lamRange = [phoenixlam[ind[0]],phoenixlam[ind[nind-1]]]
log_rebin, lamRange, ssp, sspNew, logLam2, VELSCALE=velScaletemp
print,lamrange,nind,N_ELEMENTS(sspnew)

templates = dblarr(n_elements(sspNew),nfiles)
for i=0,nfiles-1 do begin
    j=goodtemp[i]        
    phoenixflux=READFITS(phoenix[j],head)
    ind=WHERE(phoenixlam GT mintemp AND phoenixlam LT maxtemp,nind)
    ssp=phoenixflux[ind]
    lamRange = [phoenixlam[ind[0]],phoenixlam[ind[nind-1]]]
    log_rebin, lamRange, ssp, sspNew, logLam2, VELSCALE=velScaletemp

    print,lamrange,nind,N_ELEMENTS(sspnew)
    templates[*,i] = sspnew/median(sspnew)
;    templates[*,j] = convol(sspvcor,lsf)/median(sspvcor) ; Degrade template resolution
;    IF (i EQ 0) THEN plot,sspvcor,xrange=[1600,2800] ELSE oplot,sspvcor,color=i*20+40
    print,phoenixfile[j]
;    oplot,sspnew,color=100
;    a=GET_KBRD()
endfor
save,templates,logLam2,spclass,lumclass,goodtemp,file=phoenixsave
ENDIF ELSE RESTORE, file=phoenixsave
spclass=spclass[goodtemp]
lumclass=lumclass[goodtemp]

;now correct templates for difference in resolution
c=2.99792458d5
sigma2fwhm=2*SQRT(2*ALOG(2))
meanlambda=(maxfit+minfit)/2.
nifssigma=(fwhm/sigma2fwhm)/meanlambda*c ;in km/sec
; Phoenix FWHM resolution is R~100000, this is 3 km/sec
phoenixsigma=3./sigma2fwhm ;in km/sec
IF (nifssigma GT phoenixsigma) THEN BEGIN
    sigma = SQRT(nifssigma^2-phoenixsigma^2)/velScaletemp ; Quadratic sigma difference in pixels 
    lsf = psf_Gaussian(NPIXEL=(8*sigma > 5.), ST_DEV=sigma, /NORM, NDIM=1)
    templatesize=SIZE(templates,/dim)
;convtemplates=FLTARR(templatesize)
    FOR i=0,templatesize[1]-1 DO BEGIN
        convtemp=convol(templates[*,i],lsf) ; Degrade template resolution
        ;now rebin
        outpix=N_ELEMENTS(convtemp) / tempbin
        inmaxpix=outpix*tempbin-1
        rebintemp=REBIN(convtemp[0:inmaxpix],outpix)
        loglamrebin=REBIN(loglam2[0:inmaxpix],outpix)
        IF (i EQ 0) THEN BEGIN
            convtemplates=DBLARR(outpix,templatesize[1])
        ENDIF
        convtemplates[*,i] = rebintemp
    ENDFOR
ENDIF 
 
; The galaxy and the template spectra do not have the same starting wavelength.
; For this reason an extra velocity shift DV has to be applied to the template
; to fit the galaxy spectrum. We remove this artificial shift by using the
; keyword VSYST in the call to PPXF below, so that all velocities are
; measured with respect to DV.
;
dv = (logLamrebin[0]-logLam1[0])*299792.458d ; km/s

; Here the actual fit starts. The best fit is plotted on the screen.
; Gas emission lines are excluded from the pPXF fit using the GOODPIXELS keyword.
READCOL,'h2_wavnum.dat',wave,wn,strength,name,FORMAT='D,F,F,A',/SILENT
lamlin=EXP(loglam1)
lind=WHERE(strength GT 0.2,nlind)
badind=[0]
linewidth=18.
wave=wave*(1+initvel/2.9979d5)
FOR i=0,nlind-1 DO BEGIN
    inline=WHERE(lamlin GT wave[lind[i]]-linewidth/2. AND lamlin LT wave[lind[i]]+linewidth/2.)
    badind=[badind,inline]
    ENDFOR
ca8lambda=2.32204e4*(1+initvel/2.9979d5)
ca8ind=WHERE(ABS(lamlin-ca8lambda) LT 15.)
fe3lambda=2.3485e4*(1+initvel/2.9979d5)
fe3ind=WHERE(ABS(lamlin-fe3lambda) LT 25.)
badind=[badind,ca8ind,fe3ind]
badind=badind[1:*]
linearr=FLTARR(N_ELEMENTS(lamlin))
linearr[badind]=1


start = [initvel, 3.*velscale] ; (km/s), starting guess for [V,sigma]
;start = [initvel, 40.] ; (km/s), starting guess for [V,sigma]
ppxf, convtemplates, galaxy, error, velScale, start, sol,  MOMENTS=4, DEGREE=4, VSYST=dv,WEIGHTS=weights,QUIET=QUIET,BESTFIT=bestfit,/PLOT,/OVERSAMPLE
axis,xaxis=2,xrange=[EXP(MIN(loglam1))/1.e4,EXP(MAX(loglam1))/1.e4],chars=1.5,/xsty,xtitle='microns'
axis,xaxis=1,xrange=[EXP(MIN(loglam1))/1.e4,EXP(MAX(loglam1))/1.e4],chars=0.0001,/xsty
resid=galaxy-bestfit
MEANCLIP,resid,meanresid,4.0,subs=goodresid
;calculate error to use in MC tests
newerror=REPLICATE(STDDEV(resid[goodresid]),N_ELEMENTS(galaxy))
residsn=MEDIAN(galaxy/STDDEV((resid[goodresid])))
print,'S/N residual,median error',residsn,MEDIAN(galaxy/error)
error_ratio=residsn/MEDIAN(galaxy/error)


;now check what the chisq is if there is no dispersion
check=[sol[0],velscale/10.0]
ppxf, convtemplates, galaxy, error, velScale, check, checksol, MOMENTS=0, DEGREE=4, VSYST=dv,WEIGHTS=weights,goodpix=goodpix,/quiet;,/plot
nodispchi=checksol[6]
loglamout=loglam1

meanspec=TOTAL(weights*spclass)/TOTAL(weights)
sdevspec=SQRT(TOTAL(weights*(spclass-meanspec)^2)/TOTAL(weights))
meanlum=TOTAL(weights*lumclass)/TOTAL(weights)
sdevlum=SQRT(TOTAL(weights*(lumclass-meanlum)^2)/TOTAL(weights))
sol=[sol,meanspec,sdevspec,meanlum,sdevlum]
outrv=sol[0]
outdisp=sol[1]
loadct,0,/silent

IF KEYWORD_SET(DOERROR) THEN BEGIN
    nruns=50
    npixels=N_ELEMENTS(galaxy)
    allsol=FLTARR(11,nruns)
    FOR i=0,nruns-1 DO BEGIN
;        newgal=galaxy+newerror*RANDOMN(seed,npixels)
        newgal=galaxy+error*RANDOMN(seed,npixels)
        ppxf, convtemplates, newgal, error, velScale, start, onesol, $
          MOMENTS=4, DEGREE=4, VSYST=dv,bestfit=bestfit,bias=0.1,WEIGHTS=weights,goodpix=goodpix,/quiet
        meanspec=TOTAL(weights*spclass)/TOTAL(weights)
        sdevspec=SQRT(TOTAL(weights*(spclass-meanspec)^2)/TOTAL(weights))
        meanlum=TOTAL(weights*lumclass)/TOTAL(weights)
        sdevlum=SQRT(TOTAL(weights*(lumclass-meanlum)^2)/TOTAL(weights))
        onesol=[onesol,meanspec,sdevspec,meanlum,sdevlum]
        allsol[*,i]=onesol
    ENDFOR
    titstring='V='+STRING(sol[0],FORMAT='(F7.1)')+' !9s!x='+STRING(sol[1],FORMAT='(F7.1)')
    plothist,allsol[1,*],title=titstring,xtitle='!9s!x [km/sec]',ytitle='# of MC runs',bin=0.5,chars=1.5
    print, MEAN(allsol[1,*]), STDDEV(allsol[1,*])
    outerror=[MEAN(allsol[0,*]),STDDEV(allsol[0,*]),MEAN(allsol[1,*]),STDDEV(allsol[1,*]),MEAN(allsol[2,*]),STDDEV(allsol[2,*]),MEAN(allsol[3,*]),STDDEV(allsol[3,*]),MEAN(allsol[7,*]),STDDEV(allsol[7,*]),MEAN(allsol[9,*]),STDDEV(allsol[9,*])]
ENDIF    


end
;------------------------------------------------------------------------------
