PRO FIT_MULTIAGE_NGC404

;first get the model!
inpath='/Users/seth/programs/cb07/fits/' ;'/home/aseth/programs/cb07/fits/'
c=2.9979d5

testfile=inpath+'cb2007_hr_stelib_m52_chab_ssp.ised.fit'
;testfile=inpath+'bc2003_hr_m52_chab_ssp.ised.fit' 
edition='cb07'

test=MRDFITS(testfile,1)

outstem='siegel_m62young'
outfile='models/'+outstem+'.fit'

;useage=[5.,25,101,286,640,904,1600,2500,5000,10000]*1.e6
useage=[1.e6,1.e7,5.e7,1.e8,3.e8,6.e8,1.e9,2.5e9,5.e9,1.e10,1.3e10]
;useage=[5.e7,7.e7,1.e8]
nmodels=N_ELEMENTS(useage)
usemetal=['m62','m62','m62','m62','m62','m62','m62','m62','m52','m52','m32']
;usemetal=['m52','m52','m52','m52','m52','m52','m52','m62','m42','m42','m32']
;usemetal=REPLICATE('m42',nmodels)
ageind=LONARR(nmodels)

FOR i=0,nmodels-1 DO BEGIN
    min=MIN(ABS(useage[i]-test.age),minpos)
    ageind[i]=minpos
ENDFOR

;output wavelengths for Christy's code based on her bc_simple_models.pro
npix = round((alog10(8500) - alog10(3400)) / 1d-4)
logwl = dindgen(npix) * 1d-4 + alog10(3400.d0)
outc= {wave: 10.0^logwl, flux: fltarr(npix,nmodels), $
       age: test.age[ageind], id: strarr(nmodels), mlv:FLTARR(nmodels), $
       bv:FLTARR(nmodels), vi:FLTARR(nmodels), vj:FLTARR(nmodels), $
       vh:FLTARR(nmodels), vk:FLTARR(nmodels), norm:FLTARR(nmodels)}
FOR i=0,nmodels-1 DO BEGIN

    IF (edition EQ 'bc03') THEN infile=inpath+'bc2003_hr_'+usemetal[i]+'_chab_ssp.ised.fit' ELSE infile=inpath+'cb2007_hr_stelib_'+usemetal[i]+'_chab_ssp.ised.fit'    
    in=MRDFITS(infile,1)
    ;get data on colors and M/Ls
    IF edition EQ 'cb07' THEN BEGIN 
        suppath='/Users/seth/programs/cb07/stelib/'
        color2file=suppath+'cb2007_hr_stelib_'+usemetal[i]+'_chab_ssp.2color'
        color4file=suppath+'cb2007_hr_stelib_'+usemetal[i]+'_chab_ssp.4color'
    ENDIF ELSE BEGIN
        suppath='/Users/seth/programs/bc03/models/Padova1994/chabrier/'
        color2file=suppath+'bc2003_hr_'+usemetal[i]+'_chab_ssp.2color'
        color4file=suppath+'bc2003_hr_'+usemetal[i]+'_chab_ssp.4color'        
    ENDELSE
    READCOL,color2file,logage2,vi,vj,vk,jh,/SILENT,FORMAT='F,X,X,X,X,F,F,F,X,F'
    READCOL,color4file,logage4,bmag,vmag,mlb,mlv,/SILENT,FORMAT='F,X,F,F,F,F'
    bv=bmag-vmag
    vh=vj+jh

    ;apply empirical correction derived for NGC404
    READCOL,'empirical_dv.dat',dvlam,dv,/SILENT
    dv=dv+31. ;to correct for NGC404 velocity
    dvinterp=INTERPOL(dv,dvlam,in.wave,/spline)
    dlambda=dvinterp/c*in.wave
    dvind=WHERE(in.wave GT 3700. AND in.wave LT 6250.,ndvind,comp=nocor)
    in.wave[dvind]=in.wave[dvind]+dlambda[dvind]
;    pixsize=ABS(in.wave[0:*]-in.wave[1:*])
;    pixsize=[pixsize,pixsize[N_ELEMENTS(pixsize)-1]]
;    pixsize[nocor]=1.0 & pixsize[dvind[ndvind-1]]=1.0
 
    string=STRUPCASE(usemetal[i])+', Chabrier IMF'
    nwl = where(outc.wave gt 5450 and outc.wave lt 5550)
    inspeccor=in.flux[*,ageind[i]] ;/pixsize
    linterp, in.wave, inspeccor, outc.wave, spec
    spec=spec/MEDIAN(spec[nwl])
    outc.flux[*,i]=spec
    outc.mlv[i]=mlv[ageind[i]]
    outc.bv[i]=bv[ageind[i]]
    outc.vi[i]=vi[ageind[i]]
    outc.vj[i]=vj[ageind[i]]
    outc.vh[i]=vh[ageind[i]]
    outc.vk[i]=vk[ageind[i]]
ENDFOR
;STOP
IF (FILE_TEST(outfile)) THEN SPAWN,'rm '+outfile
MWRFITS,outc,outfile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;now moving on to the fitting!;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

instr_disp=90.
rv=-34.6;rvcor-cor
model = outc
;infile='NGC404_comb.fits'
infile='NGC404_bulge.fits'
minlambda=3740.
maxlambda=6200.


a=readplainfits(infile)
xwave=a[0,*]
calspec=a[1,*]/1.d-16
calspecerr=a[2,*]/1.d-16        

corwave=xwave / (1+rv/c)
ind=WHERE(corwave GE minlambda AND corwave LE maxlambda)
wave=corwave[ind]
calspec=calspec[ind]
calspecerr=calspecerr[ind]

;now try correcting for bad flux calibration
residarr=READFITS('resid/res_smooth.fits',/SILENT)
residlam=residarr[0,*]
resid=residarr[1,*]
intresid=INTERPOL(resid,residlam,wave)
calspec=calspec*(1.-intresid)

!P.MULTI=[0,1,2]
vdisp=instr_disp
;;;;;;;;;;;DO THE FITTING;;;;;;;
coef=my_bc_continuum(model, wave, calspec, calspecerr, vdisp, yfit=continuum,outssps=outssps)

;remove emission lines for calculating Chi^2
READCOL,'emlines.dat',emvac,emname,FORMAT='D,A',/SILENT
;http://www.sdss.org/dr7/products/spectra/vacwavelength.html
emlambda=emvac / (1.0 + 2.735182E-4 + 131.4182 / emvac^2 + 2.76249E8 / emvac^4)
nemlines=N_ELEMENTS(emlambda)
emlamcor=emlambda               ;no need to redshift data
badind=[0]
linewidth=10.
FOR j=0,nemlines-1 DO BEGIN
    inline=WHERE(wave GT emlamcor[j]-linewidth/2. AND wave LT emlamcor[j]+linewidth/2.)
    badind=[badind,inline]
ENDFOR
badind=badind[1:*]
linearr=FLTARR(N_ELEMENTS(wave))
linearr[badind]=1
calcind=WHERE(linearr EQ 0,ncalcind)

formatstring='(A6,'+STRTRIM(nmodels,2)+'I6,A4)'
print, 'AGE=', coef.model_age/1e6, ' Myr', format=formatstring
formatstring='(A6,'+STRTRIM(nmodels,2)+'F6.3,A4)'
print, 'LFRAC=',coef.light_frac/TOTAL(coef.light_frac), format=formatstring
print, 'TAUV=',coef.tauv, format='(A6, F6.3)'
chisq=TOTAL(((calspec-continuum)^2/calspecerr^2)[calcind])
dof=N_ELEMENTS(calspec[calcind])-(N_ELEMENTS(coef.model_age)+1)
redchisq=chisq/dof
print, 'RCHI^2=',redchisq, format='(A7, F7.3)'
lumage=TOTAL(coef.model_age*coef.light_frac)/TOTAL(coef.light_frac)
massfrac=coef.light_frac*model.mlv/TOTAL(coef.light_frac*model.mlv)
formatstring='(A6,'+STRTRIM(nmodels,2)+'F6.3,A4)'
print, 'MFRAC=',massfrac/TOTAL(massfrac), format=formatstring
massage=TOTAL(coef.model_age*massfrac)/TOTAL(massfrac)
print,'LUMAGE,MASSAGE',lumage/1.e9,massage/1.e9,FORMAT='(A,2F8.2)'
mlv=TOTAL(coef.light_frac*model.mlv)/TOTAL(coef.light_frac)
lightfrac_i=coef.light_frac*10.^(0.4*model.vi)
allmli=10.^(-0.4*(model.vi-(4.83-4.08)))*model.mlv
mli=TOTAL(lightfrac_i*allmli)/TOTAL(lightfrac_i)
lightfrac_h=coef.light_frac*10.^(0.4*model.vh)
allmlh=10.^(-0.4*(model.vh-(4.83-3.32)))*model.mlv
mlh=TOTAL(lightfrac_h*allmlh)/TOTAL(lightfrac_h)
print,'M/L_V,I,H=',mlv, mli, mlh,format='(A12, 3F7.2)'
vmag=-2.5*ALOG10(TOTAL(coef.light_frac))
imag=-2.5*ALOG10(TOTAL(lightfrac_i))
hmag=-2.5*ALOG10(TOTAL(lightfrac_h))
vi=TOTAL(model.vi*coef.light_frac)/TOTAL(coef.light_frac)
print,'V-I=',vmag-imag,format='(A7, 2F7.3)'
print,'I-H=',imag-hmag,format='(A7, 2F7.3)'
print,'RESID=',STDDEV(((calspec-continuum)/calspec)[calcind]), format='(A6, F6.3)'
print, ' '


formatstring='(A16,A16,A7,'+STRTRIM(nmodels,2)+'I7,2A7,A10)'
print,'File','Model','TauV',model.age/1.e6,'M/L_I','V-I','Chi2',FORMAT=formatstring
formatstring='(A16,A16,'+STRTRIM(nmodels+3,2)+'F7.3,F10.4)'
print,infile,outstem,coef.tauv,coef.light_frac/TOTAL(coef.light_frac),mli,vmag-imag,redchisq,FORMAT=formatstring

plot,wave[calcind],((calspec-continuum)/calspec)[calcind],/xsty,yrange=[-.1,.1],title=infile


outmodelfile='multiage/'+outstem+'.fits'
outlambda=wave
outmodel=continuum
outresid=(calspec-continuum)/calspec
outarr=[[outlambda],[outmodel],[outresid],[outssps]]
;WRITEFITS,outmodelfile,outarr


; measure absorption line indices
;icoef = absline_index(wave, calspec)
;mcoef = absline_index(wave, continuum, tag='_model') ; measure off model
;    IF (STRMATCH(filenames[i],"NGC404*") EQ 1) THEN BEGIN
;        outmodelfile='multiage/'+ztag+'_'+modeltag+'_model.fits'
;        outlambda=wave
;        outmodel=continuum
;        outresid=(calspec-continuum)/calspec
;        outarr=[[outlambda],[outmodel],[outresid],[outssps]]
;        WRITEFITS,outmodelfile,outarr
;    ENDIF
;CLOSE,1
;FREE_LUN,1
;device,/close
;set_plot,'x'
STOP
END
