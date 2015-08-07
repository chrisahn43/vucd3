PRO COMPARE_NOISE

;filestem="kband_gemcombine"
infile='cube_combine.fits'

im=READFITS(infile,head,ext=1)
var=READFITS(infile,ext=2)
imsize=SIZE(im,/dim)

lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
lambda=FINDGEN(imsize[2])*dlambda+lambda0
snind=WHERE(lambda GT 22140. AND lambda LT 22600.)
mnoiseim=im
tnoiseim=im
contim=im
FOR i=0,imsize[0]-1 DO BEGIN
    FOR j=0,imsize[1]-1 DO BEGIN
        uspec=im[i,j,snind]
        uvar=var[i,j,snind]
        test=WHERE(FINITE(uspec) EQ 1,ntest)
        IF (ntest GT 2) THEN BEGIN
            MEANCLIP,uspec,meany,4.,subs=subs
            mnoiseim[i,j]=MEAN(uspec[subs])/STDDEV(uspec[subs])
            contim[i,j]=MEAN(uspec[subs])
;            tnoiseim[i,j]=MEAN(uspec[subs]/SQRT(uvar[subs]))        
        ENDIF ELSE mnoiseim[i,j]=-1.
        tnoiseim[i,j]=MEDIAN(uspec/SQRT(uvar))        
    ENDFOR
ENDFOR    

ratio=tnoiseim/mnoiseim
;tsn=contim/
makex,ratio,x,y
myplot,file='compare_noise.ps',ysize=9,yoff=1,/inches
!P.MULTI=[0,1,2]
plot,contim[10:50,10:50],ratio[10:50,10:50],psym=4,xtitle='Signal',ytitle='Variance/Measured S/N',/xlog,xrange=[10,2000],/xsty,charsize=1.8
loadct,12
plot,contim[10:50,10:50],tnoiseim[10:50,10:50],psym=4,/xlog,/ylog,yrange=[10,200],xrange=[10,2000],/ysty,/xsty,xtitle='Signal',ytitle='Signal/Noise',charsize=1.8,symsize=0.4
oplot,contim[10:50,10:50],mnoiseim[10:50,10:50],psym=1,color=200,symsize=0.4
legend,['Variance','Measured'],psym=[4,1],color=[0,200],charsize=1.8
loadct,0
device,/close
set_plot,'x'
STOP
END
