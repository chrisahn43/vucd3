PRO ASTROMETRY

gim=READFITS('HST_10137_03_ACS_HRC_F814W_drz.fits',ext=1,ghead)
zim=READFITS('HST_10137_03_ACS_HRC_F606W_drz.fits',ext=1,zhead)
kim=READFITS('vucd3_combine_best8_cont.fits',khead)

gxcen=628-1. & gycen=660.-1.
zxcen=628.-1. & zycen=660.-1.
kxcen=44.5-1. & kycen=43.-1.

gcntrd,gim,gxcen,gycen,gxcntrd,gycntrd,4
gcntrd,zim,zxcen,zycen,zxcntrd,zycntrd,4
gcntrd,kim,kxcen,kycen,kxcntrd,kycntrd,4

print,gxcntrd,gycntrd,zxcntrd,zycntrd
print,kxcntrd,kycntrd

xyad,ghead,gxcntrd,gycntrd,gra,gdec
xyad,khead,kxcntrd,kycntrd,kra,kdec
xyad,zhead,zxcntrd,zycntrd,zra,zdec
print,(gra-kra)*3600.,(gdec-kdec)*3600.
print,(gra-zra)*3600.,(gdec-zdec)*3600.

crval1=SXPAR(khead,'CRVAL1')+(gra-kra);+0.025d0/3600.
crval2=SXPAR(khead,'CRVAL2')+(gdec-kdec);-0.01d0/3600.
SXADDPAR,khead,'CRVAL1',crval1
SXADDPAR,khead,'CRVAL2',crval2
WRITEFITS,'vucd3_best8_cont_astcor.fits',kim,khead

crval1=SXPAR(zhead,'CRVAL1')+(gra-zra);+0.025d0/3600.
crval2=SXPAR(zhead,'CRVAL2')+(gdec-zdec);-0.01d0/3600.
SXADDPAR,zhead,'CRVAL1',crval1
SXADDPAR,zhead,'CRVAL2',crval2
WRITEFITS,'vucd3_hrc_f606w_astcor.fits',zim,zhead


STOP
END
