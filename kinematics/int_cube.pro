PRO INT_CUBE,filestem,CONT=cont,CALIB=calib,SN=sn,TSN=tsn

;filestem="kband_gemcombine"
infile=filestem+'.fits'
outfile='lum_model/'+filestem+'_int.fits'

junk=READFITS(infile,head0,ext=0,/SILENT)
im=READFITS(infile,head,ext=1,/SILENT)
var=READFITS(infile,ext=2,/SILENT)
imsize=SIZE(im,/dim)
outim=TOTAL(im,3)/FLOAT(imsize[2])

IF KEYWORD_SET(CONT) THEN BEGIN
    outfile='lum_model/'+filestem+'_cont.fits'
    lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
    lambda=FINDGEN(imsize[2])*dlambda+lambda0
    snind=WHERE(lambda GT 22140. AND lambda LT 22600.)
    FOR i=0,imsize[0]-1 DO BEGIN
        FOR j=0,imsize[1]-1 DO BEGIN
            uspec=im[i,j,snind]
            test=WHERE(FINITE(uspec) EQ 1,ntest)
            IF (ntest GT 2) THEN BEGIN
                MEANCLIP,uspec,meany,4.,subs=subs
                outim[i,j]=MEAN(uspec[subs])
            ENDIF ELSE outim[i,j]=-1.                
        ENDFOR
    ENDFOR    
ENDIF

IF KEYWORD_SET(SN) THEN BEGIN
    outfile='lum_model/'+filestem+'_sn.fits'
    lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
    lambda=FINDGEN(imsize[2])*dlambda+lambda0
    snind=WHERE(lambda GT 22140. AND lambda LT 22600.)
    FOR i=0,imsize[0]-1 DO BEGIN
        FOR j=0,imsize[1]-1 DO BEGIN
            uspec=im[i,j,snind]
            test=WHERE(FINITE(uspec) EQ 1,ntest)
            IF (ntest GT 2) THEN BEGIN
                MEANCLIP,uspec,meany,4.,subs=subs
                outim[i,j]=MEAN(uspec[subs])/STDDEV(uspec[subs])
            ENDIF ELSE outim[i,j]=-1.
        ENDFOR
    ENDFOR    
ENDIF

IF KEYWORD_SET(TSN) THEN BEGIN
    outfile='lum_model/'+filestem+'_tsn.fits'
    lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
    lambda=FINDGEN(imsize[2])*dlambda+lambda0
    snind=WHERE(lambda GT 22140. AND lambda LT 22600.)
    FOR i=0,imsize[0]-1 DO BEGIN
        FOR j=0,imsize[1]-1 DO BEGIN
            uspec=im[i,j,snind]
            uvar=var[i,j,snind]
            outim[i,j]=MEDIAN(uspec/SQRT(uvar))        
        ENDFOR
    ENDFOR
    print,filestem, ' Max S/N ',MAX(outim[10:50,10:50])
ENDIF

IF KEYWORD_SET(CALIB) THEN BEGIN
    outfile='lum_model/'+filestem+'_calib.fits'
    lambda0=SXPAR(head,'CRVAL3') & dlambda=SXPAR(head,'CD3_3')
    lambda=FINDGEN(imsize[2])*dlambda+lambda0
    READCOL,'../flux_cal/kshort_2mass_response.dat',tlambda,tresponse,FORMAT='X,F,F'
    foo=READPLAINFITS('../flux_cal/bb9500.fits',/noerr)
    bb=foo[1,*]
    response=INTERPOL(tresponse,tlambda*1.e4,lambda)
    scale=response*bb
    plot,lambda,scale/MEDIAN(scale)    
    scale=scale/MEAN(scale)
    correction=REBIN(REFORM(scale,1,1,imsize[2]),imsize[0],imsize[1],imsize[2])
    outim=TOTAL(im*correction,3)/FLOAT(imsize[2])
ENDIF

SXDELPAR,head,'EXTNAME'
SXDELPAR,head,'EXTVER'

pa=SXPAR(head0,'PA')
IF (ABS(pa) GT 0.01) THEN BEGIN
    rotat=-pa/!RADEG
    cdelt=0.043d0
    cd = (cdelt / 3600.)*[ [-cos(rotat),-sin(rotat)], [-sin(rotat), cos(rotat)] ]
    SXADDPAR,head,'CD1_1',cd[0,0]
    SXADDPAR,head,'CD1_2',cd[0,1]
    SXADDPAR,head,'CD2_1',cd[1,0]
    SXADDPAR,head,'CD2_2',cd[1,1]
ENDIF

writefits,outfile,outim,head


END
