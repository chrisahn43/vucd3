PRO MAKE_RESIM,resfile,im,model,weight,FITSKY=fitsky
;make IShape style residual image with four panels

resid=im-model
IF (KEYWORD_SET(fitsky)) THEN resid=resid+fitsky

xsize=(SIZE(im))[1]
ysize=(SIZE(im))[2]

outim=FLTARR(xsize*2,ysize*2)

;normalize the weight image so that it'll look ok.
weight=weight/MEDIAN(weight)*MEDIAN(im)

outim[0:xsize-1,0:ysize-1]=resid
outim[0:xsize-1,ysize:2*ysize-1]=model
outim[xsize:2*xsize-1,0:ysize-1]=weight
;outim[xsize:2*xsize-1,0:ysize-1]=resid/im
outim[xsize:2*xsize-1,ysize:2*ysize-1]=im

WRITEFITS,resfile,outim
END
