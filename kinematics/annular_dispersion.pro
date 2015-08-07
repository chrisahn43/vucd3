@ppxf_wallace.pro
@ppxf.pro
PRO ANNULAR_DISPERSION


a=READFITS('intspec/vcc1254_int_fc2-4.fits',head)
flux=DOUBLE(a[*,0])
var=DOUBLE(a[*,1])

lambda0=SXPAR(head,'CRVAL1') & dlambda=SXPAR(head,'CD1_1')
lambda=FINDGEN(N_ELEMENTS(flux))*dlambda+lambda0
;plot,lambda,flux/MEDIAN(flux),yrange=[0.7,1.3],/ysty
initvel=1278.
ind=WHERE(lambda GT 2.30e4 and lambda LT 2.4e4)
scale=MEDIAN(flux[ind])/2.7
flux=flux/scale
var=var/scale^2
print,MEDIAN(flux[ind])
;foovar=REPLICATE(5.e-4,N_ELEMENTS(flux))
myplot,file='annular_dispersion.ps'
ppxf_wallace,lambda,flux,var,4.3,outrv=vel,outdisp=disp,outsol=sol,binnum=i,initvel=initvel,nodispchi=nodispchi,outerror=err,/doerror
device,/close
set_plot,'x'    

STOP
END
