;------------------------------------------------------------------------------
pro ppxf_wallace_full,lambda,spec,var,fwhm,OUTRV=outrv,OUTDISP=outdisp,OUTSOL=sol,DOERROR=doerror,OUTERROR=outerror,BINNUM=BINNUM,INITVEL=initvel,nodispchi=nodispchi,QUIET=QUIET,BESTFIT=bestfit,loglamout=loglamout

loadct,12,/silent
minfit=2.03e4  ;part of the target spectra to extract
maxfit=2.39e4
mintemp=2.02e4 ;part of the template to use
maxtemp=2.395e4
minfit=2.10e4  ;part of the target spectra to extract
maxfit=2.15e4
mintemp=2.08e4 ;part of the template to use
maxtemp=2.17e4
velscale=10. ;for the data
velscaletemp=1. ;for the templates
tempbin=FIX(velscale/velscaletemp)

;rv0=5 ;alpha boo's radial velocity (no V_LSR correction)
;currently no rv0 correction!
ind=WHERE(lambda GT minfit AND lambda LT maxfit,nind)
gal_lin=spec[ind]
ini_error=SQRT(var[ind])
bad=WHERE(ini_error EQ 0.0,nbad)
IF (nbad GT 0) THEN ini_error[bad]=MEDIAN(ini_error)
lamRange = [lambda[ind[0]],lambda[ind[nind-1]]]
log_rebin, lamRange, gal_lin, galaxy, logLam1, VELSCALE=velScale
log_rebin, lamRange, ini_error^2, errorsquared,  VELSCALE=velScale
error=SQRT(errorsquared)


wallacesave='wallace_templates_full.idl'
IF (FILE_TEST(wallacesave) EQ 0) THEN BEGIN
;READCOL,'../wallace/wallace_table.dat',wallaceroot,rv,FORMAT='A,X,F',/SILENT
;wallace='../wallace/'+STRLOWCASE(wallaceroot)+'_v1r.fits'
READCOL,'../wallace96/wallace_rvcor.dat',wallacefile,rv,spclass,lumclass,FORMAT='A,I,F,F',/SILENT
wallaceroot=STRMID(wallacefile,0,5)
wallace='../wallace96/'+STRLOWCASE(wallacefile)
wallacedqfiles='../wallace96/'+STRLOWCASE(wallaceroot)+'_err.fits'
goodtemp=[1,2,3,5,6,8,9,10];wallace templates to use
nfiles=N_ELEMENTS(goodtemp)
wallaceflux=READFITS(wallace[0],head)
wallacelam=SXPAR(head,'CRVAL1')+FINDGEN(N_ELEMENTS(wallaceflux))*SXPAR(head,'CD1_1')
ind=WHERE(wallacelam GT mintemp AND wallacelam LT maxtemp,nind)
ssp=wallaceflux[ind]
lamRange = [wallacelam[ind[0]],wallacelam[ind[nind-1]]]
log_rebin, lamRange, ssp, sspNew, logLam2, VELSCALE=velScaletemp
print,lamrange,nind,N_ELEMENTS(sspnew)

templates = dblarr(n_elements(sspNew),nfiles)
templatedq= FLTARR(n_elements(sspNew),nfiles)
for i=0,nfiles-1 do begin
    j=goodtemp[i]        
    wallaceflux=READFITS(wallace[j],head)
    wallacedq=READFITS(wallacedqfiles[j])
    wallacelam=SXPAR(head,'CRVAL1')+FINDGEN(N_ELEMENTS(wallaceflux))*SXPAR(head,'CD1_1')    
    ind=WHERE(wallacelam GT mintemp AND wallacelam LT maxtemp,nind)
    ssp=wallaceflux[ind]
    sspdq=wallacedq[ind]
    lamRange = [wallacelam[ind[0]],wallacelam[ind[nind-1]]]
    log_rebin, lamRange, ssp, sspNew, logLam2, VELSCALE=velScaletemp
    log_rebin, lamRange, sspdq, logsspdq, loglamfoo, VELSCALE=velScaletemp

    print,lamrange,nind,N_ELEMENTS(sspnew)
    sspvcor=SHIFT(sspnew,-rv[j])    
    templates[*,i] = sspvcor/median(sspvcor)
    templatedq[*,i] = logsspdq
;    templates[*,j] = convol(sspvcor,lsf)/median(sspvcor) ; Degrade template resolution
    IF (i EQ 0) THEN plot,sspvcor,psym=3 ELSE oplot,sspvcor,psym=3,color=i*20+40
    print,wallacefile[j]
;    oplot,sspnew,color=100
    a=GET_KBRD()
endfor
save,templates,templatedq,logLam2,spclass,lumclass,goodtemp,file=wallacesave
ENDIF ELSE RESTORE, file=wallacesave
spclass=spclass[goodtemp]
lumclass=lumclass[goodtemp]

;now correct templates for difference in resolution
c=2.99792458d5
sigma2fwhm=2*SQRT(2*ALOG(2))
meanlambda=(maxfit+minfit)/2.
nifssigma=(fwhm/sigma2fwhm)/meanlambda*c ;in km/sec
; Wallace FWHM resolution is R~45000, this is 6.67 km/sec
wallacesigma=6.67/sigma2fwhm ;in km/sec
IF (nifssigma GT wallacesigma) THEN BEGIN
    sigma = SQRT(nifssigma^2-wallacesigma^2)/velScaletemp ; Quadratic sigma difference in pixels 
    lsf = psf_Gaussian(NPIXEL=(8*sigma > 5.), ST_DEV=sigma, /NORM, NDIM=1)
    templatesize=SIZE(templates,/dim)
;convtemplates=FLTARR(templatesize)
    FOR i=0,templatesize[1]-1 DO BEGIN
        convtemp=convol(templates[*,i],lsf) ; Degrade template resolution
        convtempdq=convol(templatedq[*,i],lsf)
        ;now rebin
        outpix=N_ELEMENTS(convtemp) / tempbin
        inmaxpix=outpix*tempbin-1
        rebintemp=REBIN(convtemp[0:inmaxpix],outpix)
        loglamrebin=REBIN(loglam2[0:inmaxpix],outpix)
        rebindq=REBIN(convtempdq[0:inmaxpix],outpix)
        IF (i EQ 0) THEN BEGIN
            convtemplates=DBLARR(outpix,templatesize[1])
            convtemplatedq=FLTARR(outpix,templatesize[1])
        ENDIF
        convtemplates[*,i] = rebintemp
        convtemplatedq[*,i] = rebindq
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
;
;first mask out gas emission
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
humpind=WHERE(lamlin GT 20470. AND lamlin LT 20800.)
badind=[badind,humpind]
badind=badind[1:*]
linearr=FLTARR(N_ELEMENTS(lamlin))
linearr[badind]=1

;templates=templates[*,20]
start = [initvel, 3.*velscale] ; (km/s), starting guess for [V,sigma]
;start = [initvel, 40.] ; (km/s), starting guess for [V,sigma]
ppxf, convtemplates, galaxy, error, velScale, start, testsol, MOMENTS=4, DEGREE=4, VSYST=dv,WEIGHTS=testweights,/quiet;,/OVERSAMPLE
ind=WHERE(testweights GT 0.0,nind)
testweights=testweights/(TOTAL(testweights)*FLOAT(nind))
IF (nind GT 1) THEN $
  mask=SHIFT(TOTAL(convtemplatedq[*,ind],2)/FLOAT(nind),FIX(testsol[0]/velscale)) ELSE $
  mask=SHIFT(REFORM(convtemplatedq[*,0]),FIX(testsol[0]/velscale))
;plot,EXP(loglamrebin),mask
intmask=INTERPOL(mask,loglamrebin,loglam1)
;What should the clipping value be?
goodpix=WHERE(intmask LT 0.2 AND linearr EQ 0)
mdegree=15
ppxf, convtemplates, galaxy, error, velScale, start, sol, /PLOT, MOMENTS=4, MDEGREE=mdegree, VSYST=dv,WEIGHTS=weights,goodpix=goodpix,QUIET=QUIET,BESTFIT=bestfit;,/OVERSAMPLE
loglamout=loglam1
x = range(-1d,1d,n_elements(galaxy))
mpoly = 1d                      ; Multiplicative polynomial
for j=1,MDEGREE do mpoly += legendre(x,j)*sol[6+j]
oplot,mpoly*MEDIAN(bestfit)*0.7,color=150

axis,/xaxis,xrange=[EXP(MIN(loglam1)),EXP(MAX(loglam1))]
;now check what the chisq is if there is no dispersion
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

IF KEYWORD_SET(DOERROR) THEN BEGIN
    nruns=100
    npixels=N_ELEMENTS(galaxy)
    allsol=FLTARR(11,nruns)
    FOR i=0,nruns-1 DO BEGIN
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
    titstring='Bin '+STRTRIM(binnum,2)+' V='+STRING(sol[0],FORMAT='(F7.1)')+' !9s!x='+STRING(sol[1],FORMAT='(F7.1)')
    plothist,allsol[1,*],title=titstring,xtitle='!9s!x [km/sec]',ytitle='# of MC runs',bin=0.5
    print, MEAN(allsol[1,*]), STDDEV(allsol[1,*])
    outerror=[MEAN(allsol[0,*]),STDDEV(allsol[0,*]),MEAN(allsol[1,*]),STDDEV(allsol[1,*]),MEAN(allsol[2,*]),STDDEV(allsol[2,*]),MEAN(allsol[3,*]),STDDEV(allsol[3,*]),MEAN(allsol[7,*]),STDDEV(allsol[7,*]),MEAN(allsol[9,*]),STDDEV(allsol[9,*])]
ENDIF    


end
;------------------------------------------------------------------------------
