PRO CUBE_COMBINE
;align images spatially and combine with cosmic ray rejection

thresh=3. ;threshold in sigma
outfile='arc_combine_best8.fits'

READCOL,'Catlist_best8',fitsfiles,FORMAT='A'
READCOL,'integer_offsets_best8.dat',xoff,yoff,FORMAT='I,I'
;fitsfiles=fitsfiles+'.fits'

zeroext=READFITS(fitsfiles[0],zerohead,ext=0)
im=READFITS(fitsfiles[0],head,ext=1)
outhead=head
imsize=SIZE(im,/DIM)
nim=N_ELEMENTS(fitsfiles)

;FOR i=0,nim-1 DO fits_info,fitsfiles[i]

outimsize=[imsize[0]+MAX(xoff),imsize[1]+MAX(yoff),imsize[2]]
imstack=FLTARR([nim,outimsize])
varstack=FLTARR([nim,outimsize])
dqstack=FLTARR([nim,outimsize])
pixelstack=INTARR([nim,outimsize])
FOR i=0,nim-1 DO BEGIN
    im=READFITS(fitsfiles[i],ext=1)
    var=READFITS(fitsfiles[i],ext=2)
    dq=READFITS(fitsfiles[i],ext=3)
    xmin=xoff[i] & xmax=xoff[i]+imsize[0]-1
    ymin=yoff[i] & ymax=yoff[i]+imsize[1]-1
    imstack[i,xmin:xmax,ymin:ymax,*]=im
    varstack[i,xmin:xmax,ymin:ymax,*]=var
    dqstack[i,xmin:xmax,ymin:ymax,*]=dq
    pixelstack[i,xmin:xmax,ymin:ymax,*]=1
ENDFOR

outim=FLTARR(outimsize)
varim=FLTARR(outimsize)
dqim=FLTARR(outimsize)
npixim=INTARR(outimsize)
sigmaim=FLTARR(outimsize)
flagim=INTARR(outimsize)
varim=FLTARR(outimsize)
FOR i=1,outimsize[0]-2 DO BEGIN
FOR j=1,outimsize[1]-2 DO BEGIN
    print,i,j
FOR k=0,outimsize[2]-1 DO BEGIN
    ;base cosmic ray rejection on neighboring x-pixels
;    allval=imstack[*,i-1:i+1,j,k]
;    allpix=pixelstack[*,i-1:i+1,j,k]
    allval=imstack[*,i-1:i+1,j-1:j+1,k]
    allpix=pixelstack[*,i-1:i+1,j-1:j+1,k]
    allind=WHERE(allpix EQ 1 AND FINITE(allval) EQ 1,nallind)    
    IF (nallind GT 1) THEN BEGIN
;        av=MEDIAN(allval[allind])
;        sigma=STDDEV(allval[allind]) > 0.1        
        meanclip,allval[allind],av,sigma,clipsig=thresh
        sigma=sigma > 0.1
        val=imstack[*,i,j,k]
        var=varstack[*,i,j,k]
        dq=dqstack[*,i,j,k]
        pix=pixelstack[*,i,j,k]
        IF (ABS(MEDIAN(val)-av) GT sigma) THEN flagim[i,j,k]=1
;        ind=WHERE(pix EQ 1 AND ABS(val-av) LT thresh*sigma,nind)
        ind=WHERE(pix EQ 1 AND ABS(val-av) LT thresh*sigma AND dq LT 0.2,nind)
        ;adding a variance clip causes problem on the red end in the center
        IF (nind GE 2) THEN BEGIN
            outim[i,j,k]=MEAN(val[ind])
            varim[i,j,k]=TOTAL(var[ind])/(FLOAT(nind))^2
            dqim[i,j,k]=MEAN(dq[ind])
            npixim[i,j,k]=nind
            sigmaim[i,j,k]=sigma
        ENDIF ELSE BEGIN
            IF (nind EQ 1) THEN BEGIN
                outim[i,j,k]=val[ind]
                varim[i,j,k]=var[ind]
                dqim[i,j,k]=dq[ind]
                npixim[i,j,k]=nind
                sigmaim[i,j,k]=sigma
            ENDIF
        ENDELSE
;        print,allval
;        STOP
    ENDIF

ENDFOR    
ENDFOR
ENDFOR

;WRITEFITS,'cube_combine.fits',outim,head
;WRITEFITS,'cc_npix.fits',npixim,head
;WRITEFITS,'cc_sigma.fits',sigmaim,head
;WRITEFITS,'cc_flag.fits',flagim,head

fits_write,outfile,zeroext,zerohead
fits_open,outfile,fcbout,/update
fits_write,fcbout,outim,outhead,extname='SCI'
fits_write,fcbout,varim,outhead,extname='VAR'
fits_write,fcbout,dqim,outhead,extname='DQ'
fits_write,fcbout,npixim,outhead,extname='NPIX'
fits_write,fcbout,sigmaim,outhead,extname='SIG'
fits_write,fcbout,flagim,outhead,extname='FLAG'

fits_close,fcbout


STOP
END
