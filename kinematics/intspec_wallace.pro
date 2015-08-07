@ppxf_wallace.pro
PRO INTSPEC_WALLACE,MIN_SN=min_sn

initvel=713.
READCOL,'intspec/spectra_fwhm.dat',rin,rout,avrad,fwhm,filename,FORMAT='F,F,F,F,A'

infiles='intspec/'+filename
test=READFITS(infiles[0],head)
lambda0=SXPAR(head,'CRVAL1') & dlambda=SXPAR(head,'CD1_1')
refpix=SXPAR(head,'CRPIX1')
lambda=(FINDGEN(N_ELEMENTS(test))*dlambda+lambda0-dlambda*refpix)
nlambda=N_ELEMENTS(lambda)
nfiles=N_ELEMENTS(infiles)
;read in files & calculate signal to noise
specarr=FLTARR(nfiles,nlambda)
snind=WHERE(lambda GT 22140. AND lambda LT 22600.)
snarr=FLTARR(nfiles)


outplotfile='vor_out/intspec_wallace_best8.ps'
myplot,file=outplotfile,ysize=9,yoff=1,/inches
!P.MULTI=[0,1,2]
outfile='vor_out/intspec_wallace_best8.dat'
OPENW,1,outfile,WIDTH=1000
FOR i=0,nfiles-1 DO BEGIN
   allspec=READFITS(infiles[i])
   flux=allspec[*,0]
   var=allspec[*,1]
   minfit=2.290e4               ;part of the target spectra to extract
   maxfit=2.395e4
   ind=WHERE(lambda GT minfit AND lambda LT maxfit)
   ppxf_wallace,lambda,flux,var,fwhm[i],outrv=vel,outdisp=disp,outsol=sol,binnum=i,initvel=initvel,nodispchi=nodispchi,outerror=err,/doerror,residsn=residsn
   printf,1,infiles[i],rin[i],rout[i],avrad[i],residsn,sol[6],sol[0],err[0],err[1],sol[1],err[2],err[3],sol[2],err[4],err[5],sol[3],err[6],err[7],sol[7],err[8],err[9],sol[9],err[10],err[11],nodispchi,FORMAT='(A30,24F9.3)'

;       printf,1,i,xoutbin[i],youtbin[i],npixels[i],measured_sn[i],sn[i],sol[6],sol[0],sol[1],sol[2],sol[3],sol[7],sol[9],nodispchi,FORMAT='(I4,2F7.3,I5,10F9.3)'

ENDFOR
CLOSE,1


DEVICE,/close
SET_PLOT,'x'
!P.MULTI=[0,1,1]


STOP
END

