PRO MEASURE_SKYVEL

COMMON skyblock, gauss_sep

infile='m60-ucd1_sky_combine_may18.fits'
cube=READFITS(infile,head,ext=1,/SILENT)
npix=READFITS(infile,head,ext=4,/SILENT)
cubesize=SIZE(cube,/dim)
lambda0=SXPAR(head,'CRVAL3')
dlambda=SXPAR(head,'CD3_3')
lambda=FINDGEN(cubesize[2])*dlambda+lambda0

sigma2fwhm=2*SQRT(2*ALOG(2))

READCOL,'../rousselot2000.dat',skylambda,skyint,FORMAT='F,F',skip=28
;these lines all have close doublets except 21955 line and are
;relatively isolated
uselambda=[20339.,20412.7,20563.5,20729.0,21176.,21802.3,21955.6,22052.,22125.5,22247.,22312.7,22460.0,22517.9]
nlines=N_ELEMENTS(uselambda)
reflambda=FLTARR(nlines,2)
FOR i=0,nlines-1 DO BEGIN
    ind=WHERE(ABS(skylambda-uselambda[i]) LT 2. AND skyint GT 1.,nind)
    IF (nind NE 2) THEN STOP
    reflambda[i,*]=skylambda[ind]
ENDFOR

c=2.99792458d5
intmap=FLTARR(cubesize[0],cubesize[1],nlines)
velmap=FLTARR(cubesize[0],cubesize[1],nlines)
dispmap=FLTARR(cubesize[0],cubesize[1],nlines)
intmap2=FLTARR(cubesize[0],cubesize[1],nlines)
velmap2=FLTARR(cubesize[0],cubesize[1],nlines)
dispmap2=FLTARR(cubesize[0],cubesize[1],nlines)

minlcont=-30.
maxlcont=-20.
minucont=20.
maxucont=30.
fitwidth=15.

loadct,5
FOR k=0,nlines-1 DO BEGIN
    centerwave=MEAN(reflambda[k,*])
    gauss_sep=reflambda[k,1]-reflambda[k,0]
    lind=WHERE(lambda GE centerwave+minlcont AND lambda LE centerwave+maxlcont,nlind)
    uind=WHERE(lambda GE centerwave+minucont AND lambda LE centerwave+maxucont,nuind)

    fitind=WHERE(lambda GE centerwave-fitwidth AND lambda LE centerwave+fitwidth,nfitind)

FOR i=0,cubesize[0]-1 DO BEGIN
FOR j=0,cubesize[1]-1 DO BEGIN
    spec=cube[i,j,*]
    IF (MAX(spec) EQ 0.) THEN GOTO, SKIPTHIS
    background=MEDIAN([spec[lind],spec[uind]])    
    estimates=[MAX(spec[fitind]),centerwave,4.2/sigma2fwhm,background]

;expresult=curvefit(rphysical[good],rgbdensity[good],weights[good],expa,sigma,function_name='exp',chisq=expchisq,/NODER)

    fit=GAUSSFIT(lambda[fitind],spec[fitind],par,nterms=4,estimates=estimates,chisq=chisq)
    par2=estimates
;    weights=1./(spec[fitind] > 1.)
;    weights=REPLICATE(1.,nfitind)
    fit2=CURVEFIT(lambda[fitind],spec[fitind],weights,par2,par2err,function_name='gauss_binary',chisq=chisq2,itmax=100)
;    print,par[0],par2[0],par[1],par2[1],par[2],par2[2],par[3],par2[3],chisq,chisq2
    intmap[i,j,k]=par[0]*par[2]
    velmap[i,j,k]=(par[1]-centerwave)/centerwave*c
    dispmap[i,j,k]=par[2]*sigma2fwhm
    intmap2[i,j,k]=par2[0]*par2[2]
    velmap2[i,j,k]=(par2[1]-centerwave)/centerwave*c
    dispmap2[i,j,k]=par2[2]*sigma2fwhm
    titstring=STRING(i)+'_'+STRING(j)+' '+STRING(chisq,FORMAT='(F6.1)')
    IF (j EQ 30) THEN BEGIN
    plot,lambda,spec,title=titstring,xrange=[centerwave+minlcont,centerwave+maxucont],/xsty
    oplot,lambda[fitind],fit,color=100
    oplot,lambda[fitind],fit2,color=200
    oplot,lambda[fitind],spec[fitind],psym=4
    WAIT,0.1
    ENDIF
SKIPTHIS:
ENDFOR
ENDFOR
intmap[*,*,k]=intmap[*,*,k]/MEDIAN(intmap[*,*,k])
intmap2[*,*,k]=intmap2[*,*,k]/MEDIAN(intmap2[*,*,k])
plot, dispmap[*,*,k],dispmap2[*,*,k],psym=4,xrange=[3,5],yrange=[3,5]
ENDFOR
FOR k=0,nlines-1 DO print,'Dispersion of',STRING(MEAN(reflambda[k,*]),'(F9.2)'),' is ',STRING(MEDIAN(dispmap[*,*,k]),FORMAT='(F6.3)'),STRING(MEDIAN(dispmap2[*,*,k]),FORMAT='(F6.3)')

ind=WHERE(ABS(dispmap) GT 10.)
dispmap[ind]=10.
ind=WHERE(ABS(dispmap2) GT 10.)
dispmap2[ind]=10.

plothist,dispmap2[*,*,1],xrange=[2,6],bin=0.1,yrange=[0,1000]
for i=1,nlines-1 do plothist,dispmap2[*,*,i],bin=0.1,/over,color=i*10


totalfrac=FLTARR(nlines)
lt4frac=FLTARR(nlines)
meanarr=FLTARR(nlines)
sigmaarr=FLTARR(nlines)
goodpix=WHERE(npix[*,*,1000] GE 1,totalpix)
dispsubmap=dispmap2
medianmap=MEDIAN(dispmap2,dim=3)
FOR i=0,nlines-1 DO BEGIN
    thisdisp=dispmap2[*,*,i]
    totalind=WHERE(dispmap2[*,*,i] GT 3. AND dispmap2[*,*,i] LT 5.5,ntotalind)
    lt4ind=WHERE(dispmap2[*,*,i] GT 3. AND dispmap2[*,*,i] LT 4.,nlt4ind)
    totalfrac[i]=ntotalind/FLOAT(totalpix)
    lt4frac[i]=nlt4ind/FLOAT(totalpix)    
    dispsubmap[*,*,i]=dispmap2[*,*,i]-medianmap
    sigmaarr[i]=STDDEV(thisdisp[totalind]-medianmap[totalind])
    meanarr[i]=MEAN(thisdisp[totalind])
ENDFOR

plot,uselambda,totalfrac,psym=4
oplot,uselambda,lt4frac,color=100,psym=4

plot,reflambda[*,1]-reflambda[*,0],meanarr,psym=4,ysty=16
plot,uselambda,sigmaarr
redind=WHERE(uselambda GT 2.15e4,comp=blueind)
WRITEFITS,'fullcube_int_may18.fits',intmap2
WRITEFITS,'fullcube_vel_may18.fits',velmap2
WRITEFITS,'fullcube_disp_may18.fits',dispmap2
WRITEFITS,'fullcube_dispsub_may18.fits',dispsubmap
WRITEFITS,'fullcube_int_med_may18.fits',MEDIAN(intmap2,dim=3)
WRITEFITS,'fullcube_vel_med_may18.fits',MEDIAN(velmap2,dim=3)
WRITEFITS,'fullcube_disp_med_may18.fits',MEDIAN(dispmap2,dim=3)
WRITEFITS,'redlines_disp_med_may18.fits',MEDIAN(dispmap2[*,*,redind],dim=3)
WRITEFITS,'bluelines_disp_med_may18.fits',MEDIAN(dispmap2[*,*,blueind],dim=3)

STOP
END


PRO doublegauss,x,a,f,pder
; a double gaussian where the separation between the components is
; known and the amplitude is the same.  
;   amplitude=a[0]
;   central wavelength=a[1]
;   sigma=a[2]
;   zeropoint=a[3]

COMMON skyblock, gauss_sep


;  f=a[3]+(a[0]/(SQRT(2*!PI)*a[2])*EXP(-(x-(a[1]-gauss_sep/2.))^2/(2.*a[2]^2)));+a[0]/(SQRT(2*!PI)*a[2])*EXP(-(x-(a[1]+gauss_sep/2.))^2/(2.*a[2]^2)))/2.

   f=a[3]+(a[0]*EXP(-(x-(a[1]-gauss_sep/2.))^2/(2.*a[2]^2))+a[0]*EXP(-(x-(a[1]+gauss_sep/2.))^2/(2.*a[2]^2)))/2.

   IF (N_PARAMS() GE 4) THEN $
     pder = [[(EXP(-(x-(a[1]-gauss_sep/2.))^2/(2.*a[2]^2))+EXP(-(x-(a[1]+gauss_sep/2.))^2/(2.*a[2]^2)))/2.], $
             [(a[0]*((x-a[1]+gauss_sep/2.)/a[2]^2)*EXP(-(x-(a[1]-gauss_sep/2.))^2/(2.*a[2]^2))+a[0]*((x-a[1]+gauss_sep/2.)/a[2]^2)*EXP(-(x-(a[1]+gauss_sep/2.))^2/(2.*a[2]^2)))/2.], $
             [(a[0]*(-(x-a[1]+gauss_sep/2.)^2/a[2]^3)*EXP(-(x-(a[1]-gauss_sep/2.))^2/(2.*a[2]^2))+a[0]*(-(x-a[1]-gauss_sep/2.)^2/a[2]^3)*EXP(-(x-(a[1]+gauss_sep/2.))^2/(2.*a[2]^2)))/2.], $
             [REPLICATE(1.,N_ELEMENTS(x))]]

   RETURN
END


; NAME:
;   GAUSS_BINARY
;   X = VALUES OF INDEPENDENT VARIABLE.
;   A = PARAMETERS OF EQUATION DESCRIBED BELOW.
; OUTPUTS:
;   F = VALUE OF FUNCTION AT EACH X(I).
;
; OPTIONAL OUTPUT PARAMETERS:
;   PDER = (N_ELEMENTS(X),6) ARRAY CONTAINING THE
;       PARTIAL DERIVATIVES.  P(I,J) = DERIVATIVE
;       AT ITH POINT W/RESPECT TO JTH PARAMETER.
; COMMON BLOCKS:
;   NONE.
; SIDE EFFECTS:
;   NONE.
; RESTRICTIONS:
;   NONE.
; PROCEDURE:
;   F = A(0)*EXP(-Z^2/2) + A(3)
;   Z = (X-A(1)+-gauss_sep)/A(2)

PRO GAUSS_BINARY,X,A,F,PDER

COMMON skyblock, gauss_sep

COMPILE_OPT idl2, hidden
ON_ERROR,2                      ;Return to caller if an error occurs
n = n_elements(a)
nx = N_ELEMENTS(x)

if a[2] ne 0.0 then begin
    Z1 = (X-A[1]-gauss_sep/2.)/A[2]           ;GET Z
    Z2 = (X-A[1]+gauss_sep/2.)/A[2]  
    EZ1 = EXP(-Z1^2/2.)           ;GAUSSIAN PART
    EZ2 = EXP(-Z2^2/2.)        
endif else begin
    z = REPLICATE(FIX(100, TYPE=SIZE(x,/TYPE)), nx)
    ez = z*0
endelse

 F = A[0]*(EZ1+EZ2)/2. + A[3]

 IF N_PARAMS(0) LE 3 THEN RETURN ;NEED PARTIAL?
;
 PDER = FLTARR(nx, n)           ;YES, MAKE ARRAY.
 PDER[*,0] = (EZ1+EZ2)/2.                 ;COMPUTE PARTIALS
 if a[2] ne 0. then PDER[*,1] = A[0] * (EZ1 * Z1/A[2] + EZ2 * Z2/A[2])/2.
 PDER[*,2] = A[0] * (EZ1 * Z1^2/A[2] + EZ2 * Z2^2/A[2])/2.
 PDER[*,3] = 1.
 RETURN
END


