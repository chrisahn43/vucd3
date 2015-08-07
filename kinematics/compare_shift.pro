PRO COMPARE_SHIFT

finfile='CatfbrgnN20140220S0284.fits'
fifu=READFITS(finfile,head,ext=1)
size=SIZE(fifu,/DIM)
lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
flambda=FINDGEN(size[2])*dlambda+lambda0
fx=38 & fy=14

minfile='CatfbrgnN20140518S0055.fits'
mifu=READFITS(minfile,head,ext=1)
size=SIZE(mifu,/DIM)
lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
mlambda=FINDGEN(size[2])*dlambda+lambda0
mx=40 & my=15


mspec=TOTAL(TOTAL(mifu[mx-4:mx+4,my-4:my+4,*],1),1)
fspec=TOTAL(TOTAL(fifu[fx-4:fx+4,fy-4:fy+4,*],1),1)
plot,mlambda,mspec/MEDIAN(mspec),xrange=[2.28e4,2.35e4]
oplot,flambda*(1.+39./2.99d5),fspec/MEDIAN(fspec)/1.05,color=100
STOP
END
