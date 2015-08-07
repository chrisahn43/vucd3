@ppxf_wallace.pro
PRO INTSPEC_INTEGRATED

initvel=1290.

infiles='intspec/m60-ucd1_best9_int_fc0-2.fits'
allspec=READFITS(infiles[0],head)
lambda0=SXPAR(head,'CRVAL1') & dlambda=SXPAR(head,'CD1_1')
refpix=SXPAR(head,'CRPIX1')
lambda=(FINDGEN(N_ELEMENTS(allspec[*,0]))*dlambda+lambda0-dlambda*refpix)


outplotfile='vor_out/intspec_integrated_02.ps'
myplot,file=outplotfile,ysize=9,yoff=1,/inches
!P.MULTI=[0,1,2]
outfile='vor_out/intspec_integrated_02.dat'
OPENW,1,outfile,WIDTH=1000
flux=allspec[*,0]
var=allspec[*,1]
minfit=2.295e4                  ;part of the target spectra to extract
maxfit=2.395e4
ind=WHERE(lambda GT minfit AND lambda LT maxfit)
ppxf_wallace,lambda,flux,var,4.1,outrv=vel,outdisp=disp,outsol=sol,binnum=i,initvel=initvel,nodispchi=nodispchi,outerror=err,/doerror,residsn=residsn
printf,1,infiles,0,15,7.,residsn,sol[6],sol[0],err[0],err[1],sol[1],err[2],err[3],sol[2],err[4],err[5],sol[3],err[6],err[7],sol[7],err[8],err[9],sol[9],err[10],err[11],nodispchi,FORMAT='(A30,24F9.3)'

CLOSE,1


DEVICE,/close
SET_PLOT,'x'
!P.MULTI=[0,1,1]


STOP
END

