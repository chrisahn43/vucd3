;------------------------------------------------------------------------------
pro rv_correct

mintemp=2.20e4
maxtemp=2.38e4
velscale=1.

READCOL,'wallacelist.dat',wallaceroot,sptype,lumclass,FORMAT='A,A,A',/SILENT
nfiles=N_ELEMENTS(wallaceroot)
wallace=STRLOWCASE(wallaceroot)+'_dop.fits'
wallaceflux=READFITS(wallace[5],head)
wallacelam=SXPAR(head,'CRVAL1')+FINDGEN(N_ELEMENTS(wallaceflux))*SXPAR(head,'CD1_1')
reflam=wallacelam

ind=WHERE(wallacelam GT mintemp AND wallacelam LT maxtemp,nind)
ssp=wallaceflux[ind]
lamRange = [wallacelam[ind[0]],wallacelam[ind[nind-1]]]
print,nind
log_rebin, lamRange, ssp, sspNew, logLam2, VELSCALE=velscale
templates = dblarr(n_elements(sspNew),nfiles)
templam = dblarr(nind,nfiles)
templin = dblarr(nind,nfiles)
for j=0,nfiles-1 do begin
    wallaceflux=READFITS(wallace[j],head,/SILENT)
    wallacelam=SXPAR(head,'CRVAL1')+FINDGEN(N_ELEMENTS(wallaceflux))*SXPAR(head,'CD1_1')    
    print,wallace[j],MIN(wallacelam),MAX(wallacelam),SXPAR(head,'CD1_1')
;    wallaceflux=RESAMPLE_SPECTRA(wallaceflux,wallacelam,reflam,oversamp=3)
    
    ind=WHERE(reflam GT mintemp AND reflam LT maxtemp,nind)
    ssp=wallaceflux[ind]
    templam[*,j]=reflam[ind]
    templin[*,j]=wallaceflux[ind]
    lamRange = [reflam[ind[0]],reflam[ind[nind-1]]]
;    print,nind, lamrange
;    print,SXPAR(head,'CD1_1')
    log_rebin, lamRange, ssp, sspNew, logLam2, VELSCALE=velscale
    
    templates[*,j]=sspNew/median(sspNew)
;    help,sspnew
endfor
;STOP
;dv = (logLam2[0]-logLam1[0])*299792.458d ; km/s
start = [10., 10.] ; (km/s), starting guess for [V,sigma]

refind=5;(WHERE(wallaceroot EQ 'HD4188'))[0]
reference=templates[*,refind]
error=reference*0.+1.
;OPENW,1,'wallace_rvcor.dat'
loadct,12,/sile
!P.MULTI=[0,1,2]
FOR i=0,nfiles-1 DO BEGIN
;    ppxf, reference, templates[*,i], error, velScale, start, sol, $
;      /PLOT, MOMENTS=2, DEGREE=4, MDEGREE=0 ;,/OVERSAMPLE
    velrange=RANGE(-150,150)
    result=c_correlate(reference,templates[*,i],velrange)
    plot,velrange,result,xtitle='dV [km/sec]',ytitle='XCor'
    plot,templam[*,i],templin[*,i],xrange=[2.292e4,2.3e4],/xsty;,xrange=[mintemp,maxtemp],/xsty
    tmp=max(result,maxpos)
    oplot,templam[*,i],SHIFT(templin[*,i],-velrange[maxpos]),color=100
    oplot,templam[*,refind],templin[*,5],color=200

    print,wallace[i],velrange[maxpos],' ',sptype[i],' ',lumclass[i]
;    printf,1,wallace[i],sol[0]

    a=get_kbrd()
ENDFOR
!P.MULTI=[0,1,1]
;CLOSE,1
;FREE_LUN,1
end
;------------------------------------------------------------------------------