;------------------------------------------------------------------------------
PRO ppxf_nifs,lambda,spec,var,OUTRV=outrv,OUTDISP=outdisp,OUTSOL=sol,DOERROR=doerror,OUTERROR=outerror

loadct, 12,/SILENT ; Green=40, Blue=100, Magenta=130, Red=200
minfit=2.28e4  ;part of the target spectra to extract
maxfit=2.395e4
mintemp=2.26e4 ;part of the template to use
maxtemp=2.41e4
ind=WHERE(lambda GT minfit AND lambda LT maxfit,nind)
gal_lin=spec[ind]
IF (KEYWORD_SET(var) EQ 1) THEN BEGIN
    ini_error=SQRT(var[ind])
    bad=WHERE(ini_error EQ 0.0,nbad)
    IF (nbad GT 0) THEN ini_error[bad]=MEDIAN(ini_error)
ENDIF ELSE BEGIN
    snind=WHERE(lambda GT 22140. AND lambda LT 22866)
    MEANCLIP,spec[snind],meany,4.,subs=subs
    signalnoise=MEAN(spec[snind[subs]])/STDDEV(spec[snind[subs]])
    print,'S/N calculated',signalnoise
    ini_error = gal_lin*0.+MEDIAN(gal_lin)/signalnoise
ENDELSE

lamRange = [lambda[ind[0]],lambda[ind[nind-1]]]
log_rebin, lamRange, gal_lin, galaxy, logLam1, VELSCALE=velScale
log_rebin, lamRange, ini_error, error,  VELSCALE=velScale


;template = file_search('nifs_templates/*n.fits',COUNT=nfiles)
READCOL,'../kband/nifs_templates/template_table.dat',template,rv,FORMAT='A,F',/SILENT
template=template[3:5]
rv=rv[3:5]
nfiles=N_ELEMENTS(template)
template='../kband/nifs_templates/'+template

lambdascale=1.0;0.995
templateflux=READFITS(template[0],head,/SILENT);,ext=1)
templatelam=SXPAR(head,'CRVAL1')+FINDGEN(N_ELEMENTS(templateflux))*SXPAR(head,'CD1_1')*lambdascale
ind=WHERE(templatelam GT mintemp AND templatelam LT maxtemp,nind)
ssp=templateflux[ind]
lamRange = [templatelam[ind[0]],templatelam[ind[nind-1]]]
log_rebin, lamRange, ssp, sspNew, logLam2, VELSCALE=velScale
templates = dblarr(n_elements(sspNew),nfiles)
for j=0,nfiles-1 do begin
    templateflux=READFITS(template[j],head,/SILENT);,ext=1)
    templatelam=SXPAR(head,'CRVAL1')+FINDGEN(N_ELEMENTS(templateflux))*SXPAR(head,'CD1_1')*lambdascale    
    ind=WHERE(templatelam GT mintemp AND templatelam LT maxtemp,nind)
    ssp=templateflux[ind]
    lamRange = [templatelam[ind[0]],templatelam[ind[nind-1]]]
;    print,lamrange
    log_rebin, lamRange, ssp, sspNew, logLam2, VELSCALE=velScale
;    templates[*,j] = convol(sspNew,lsf)/median(sspNew) ; Degrade template resolution
    
    npix=N_ELEMENTS(sspNew)
    inpix=FINDGEN(npix)
    outpix=FINDGEN(npix)+rv[j]/velScale
    sspvcor=RESAMPLE_SPECTRA(sspNew,inpix,outpix)
    templates[*,j] = sspvcor/median(sspvcor)
endfor


;NOTE: need to add radial velocity to account for template offset
dv = (logLam2[0]-logLam1[0])*299792.458d0+rv[0] ; km/s

;templates=templates[*,0]
start = [-40., 30.] ; (km/s), starting guess for [V,sigma]
ppxf, templates, galaxy, error, velScale, start, sol, $
    /PLOT, MOMENTS=4, DEGREE=4, VSYST=dv,MDEGREE=4;,/OVERSAMPLE
outrv=sol[0]
outdisp=sol[1]
outsol=sol
loadct,0,/silent

IF KEYWORD_SET(DOERROR) THEN BEGIN
    nruns=250
    npixels=N_ELEMENTS(galaxy)
    allsol=FLTARR(7,nruns)
    FOR i=0,nruns-1 DO BEGIN
        newgal=galaxy+error*RANDOMN(seed,npixels)
        ppxf, templates, newgal, error, velScale, start, onesol, $
          MOMENTS=4, DEGREE=4, VSYST=dv,bestfit=bestfit,bias=0.1,WEIGHTS=weights,/quiet
        allsol[*,i]=onesol
    ENDFOR
    !P.MULTI=[0,1,2]
    plothist,allsol[0,*]
    plothist,allsol[1,*]
    !P.MULTI=[0,1,1]
    WAIT,5
    print, MEAN(allsol[1,*]), STDDEV(allsol[1,*])
    outerror=[MEAN(allsol[0,*]),STDDEV(allsol[0,*]),MEAN(allsol[1,*]),STDDEV(allsol[1,*]),MEAN(allsol[2,*]),STDDEV(allsol[2,*]),MEAN(allsol[3,*]),STDDEV(allsol[3,*])]
ENDIF    

end
;------------------------------------------------------------------------------
