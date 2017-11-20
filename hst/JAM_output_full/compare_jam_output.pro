PRO COMPARE_JAM_OUTPUT
  infile='../../kinematics/vor_out/intspec_wallace_best8.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  xbin=rav*0.05                 ;major axis*0.05 arcsec/pixel
  xbinerror=((rav-rin)/2.)*0.05
  ybin=fltarr(n_elements(rav)) ;minor axis
  distance=16.5
;  mbhs=[0.,10^(findgen(6)*0.2+6.)]
  mbhs=[0,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6]
  ml=findgen(40)*0.05+0.05
  inclinations=[60.,70.,80.,90.]
  betas=findgen(11)*0.1-0.2
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  ninclinations=n_elements(inclinations)
  nbetas=n_elements(betas)
  fixgauss_mass=mrdfits('fixgausspsf_mass.fits',1)
  fixnorm_mass=mrdfits('./temp/temp/fixnormpsf_mass.fits',1)
  fixold_mass=mrdfits('fixoldpsf_mass.fits',1)
  freenorm_mass=mrdfits('freenormpsf_mass.fits',1)
  freenorm=mrdfits('freenormpsf.fits',1)
  fixnorm=mrdfits('fixnormpsf.fits',1)


  minchi=dblarr(nmbhs,nbetas)
  data=reform(fixnorm_mass,nmbhs,nbetas,ninclinations,nmls)
  for i=0,nmbhs-1 do begin
     for j=0,nbetas-1 do begin
        minchi[i,j]=MIN(data[i,j,*,*].chi2)
     endfor
  endfor
  
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='betambh.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7;BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixnorm_mass.outmbh,fixnorm_mass.inbeta,YTITLE='!7b!3',XTITLE='M!IBH!N [M!D!9n!3!N]',XCHARS=1.2,/NODATA,YCHARS=1.3,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(betas)-0.1,max(betas)+0.1],CHARSIZE=1.25,XSTY = 1,XRANGE=[min(mbhs)-2d4,MAX(mbhs)+2d4],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.95,0.98] ;,/XLOG
;  minlev=MIN(minchi) & maxlev=MAX(minchi) & nlev=2d5
;  LOADCT,0,/SILENT
;  chi2im=(minchi-minlev)/(maxlev-minlev)
;  chi2im=chi2im/TOTAL(chi2im)*nlev
;  imgunder,chi2im
  
;  PLOTSYM,0,/FILL
  OPLOT,fixnorm_mass.outmbh,fixnorm_mass.inbeta,PSYM=6,THICK=10,SYMSIZE=0.2,color=100
  xyouts,[1.e5],[-0.25],['VUCD3'],charthick=3,charsize=1.5,/data
  min=MIN(fixnorm_mass.chi2,minpos)
  CONTOUR,minchi,mbhs,betas,/OVERPLOT,LEVELS=[fixnorm_mass[minpos].chi2+2.3,fixnorm_mass[minpos].chi2+6.18,fixnorm_mass[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]
;  CONTOUR,SMOOTH(minchi,4),mbhs,betas,/OVERPLOT,LEVELS=[fixnorm_mass[minpos].chi2+3.2],C_LINESTYLE=[0],C_THICK=[10]

  xx=REPLICATE(fixnorm_mass[minpos].outmbh,N_ELEMENTS(fixnorm_mass.outmbh))
  yy=REPLICATE(fixnorm_mass[minpos].inbeta,N_ELEMENTS(fixnorm_mass.outmbh))
  OPLOT,xx,yy,psym=6,thick=15,symsize=0.7,color=150
  set_plot,'x'
;  stop
  minchi=dblarr(nmbhs,nmls)
  for i=0,nmbhs-1 do begin
     for j=0,nmls-1 do begin
        minchi[i,j]=MIN(data[i,*,*,j].chi2)
     endfor
  endfor
    
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='mlmbh.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7;BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixnorm_mass.outmbh,fixnorm_mass.ml,YTITLE='!7C!3',XTITLE='M!IBH!N [M!D!9n!3!N]',XCHARS=1.2,/NODATA,YCHARS=1.3,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(ml),2.05],CHARSIZE=1.25,XSTY = 1,XRANGE=[min(mbhs)-2d4,MAX(mbhs)+2d4],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.951,0.98] 

  
;  PLOTSYM,0,/FILL
  OPLOT,fixnorm_mass.outmbh,fixnorm_mass.ml,PSYM=6,THICK=10,SYMSIZE=0.2,color=100
  xyouts,[1.e5],[.1],['VUCD3'],charthick=3,charsize=1.5,/data

  min=MIN(fixnorm_mass.chi2,minpos)
  CONTOUR,smooth(minchi,1),mbhs,ml,/OVERPLOT,LEVELS=[fixnorm_mass[minpos].chi2+2.3,fixnorm_mass[minpos].chi2+6.18,fixnorm_mass[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixnorm_mass[minpos].outmbh,N_ELEMENTS(fixnorm_mass.outmbh))
  yy=REPLICATE(fixnorm_mass[minpos].ml,N_ELEMENTS(fixnorm_mass.outmbh))
  OPLOT,xx,yy,psym=6,thick=15,symsize=0.7,color=150
  set_plot,'x'
  stop
  minchi=dblarr(nbetas,nmls)
  for i=0,nbetas-1 do begin
     for j=0,nmls-1 do begin
        minchi[i,j]=MIN(data[*,i,*,j].chi2)
     endfor
  endfor
    
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='betaml.ps',BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixnorm_mass.inbeta,fixnorm_mass.ml,YTITLE='M/L',XTITLE='Beta',XCHARS=1,/NODATA,YCHARS=1,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(ml)-0.1,1.2],CHARSIZE=1,XSTY = 1,XRANGE=[min(betas)-0.1,max(betas)+0.1],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] 

  
  PLOTSYM,0,/FILL
  OPLOT,fixnorm_mass.inbeta,fixnorm_mass.ml,PSYM=8,THICK=5,SYMSIZE=1,color=100

  min=MIN(fixnorm_mass.chi2,minpos)
  CONTOUR,minchi,betas,ml,/OVERPLOT,LEVELS=[fixnorm_mass[minpos].chi2+2.3,fixnorm_mass[minpos].chi2+6.18,fixnorm_mass[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixnorm_mass[minpos].inbeta,N_ELEMENTS(fixnorm_mass.inbeta))
  yy=REPLICATE(fixnorm_mass[minpos].ml,N_ELEMENTS(fixnorm_mass.inbeta))
  OPLOT,xx,yy,psym=8,thick=10,symsize=2,color=150


  ind=where(fixnorm_mass.inbeta eq 0.)
  fixnorm_mass=fixnorm_mass[ind]
  minchi=dblarr(nmbhs,ninclinations)
  for i=0,nmbhs-1 do begin
     for j=0,ninclinations-1 do begin
        mbh=mbhs[i] & inc=inclinations[j]
        ind=where(fixnorm_mass.ininc eq inc and fixnorm_mass.outmbh eq mbh)
        minchi[i,j]=MIN(fixnorm_mass[ind].chi2)
     endfor
  endfor
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='incmbh_nobeta.ps',BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixnorm_mass.outmbh,fixnorm_mass.ininc,YTITLE='Inclination',XTITLE='MBH',XCHARS=1,/NODATA,YCHARS=1,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(inclinations)-5,max(inclinations)+5],CHARSIZE=1,XSTY = 1,XRANGE=[min(mbhs)-2d4,MAX(mbhs)+2d4],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] 

  
  PLOTSYM,0,/FILL
  OPLOT,fixnorm_mass.outmbh,fixnorm_mass.ininc,PSYM=8,THICK=5,SYMSIZE=1,color=100

  min=MIN(fixnorm_mass.chi2,minpos)
  CONTOUR,minchi,mbhs,inclinations,/OVERPLOT,LEVELS=[fixnorm_mass[minpos].chi2+2.3,fixnorm_mass[minpos].chi2+6.18,fixnorm_mass[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixnorm_mass[minpos].outmbh,N_ELEMENTS(fixnorm_mass.outmbh))
  yy=REPLICATE(fixnorm_mass[minpos].ininc,N_ELEMENTS(fixnorm_mass.outmbh))
  OPLOT,xx,yy,psym=8,thick=10,symsize=2,color=150

  minchi=dblarr(nmls,ninclinations)
  for i=0,nmls-1 do begin
     for j=0,ninclinations-1 do begin
        mls=ml[i] & inc=inclinations[j]
        ind=where(fixnorm_mass.ininc eq inc and fixnorm_mass.ml eq mls)
        minchi[i,j]=MIN(fixnorm_mass[ind].chi2)
     endfor
  endfor
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='incml_nobeta.ps',BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixnorm_mass.ml,fixnorm_mass.ininc,YTITLE='Inclination',XTITLE='M/L',XCHARS=1,/NODATA,YCHARS=1,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(inclinations)-5,max(inclinations)+5],CHARSIZE=1,XSTY = 1,XRANGE=[min(ml)-0.1,1.2],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] 

  
  PLOTSYM,0,/FILL
  OPLOT,fixnorm_mass.ml,fixnorm_mass.ininc,PSYM=8,THICK=5,SYMSIZE=1,color=100

  min=MIN(fixnorm_mass.chi2,minpos)
  CONTOUR,minchi,ml,inclinations,/OVERPLOT,LEVELS=[fixnorm_mass[minpos].chi2+2.3,fixnorm_mass[minpos].chi2+6.18,fixnorm_mass[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixnorm_mass[minpos].ml,N_ELEMENTS(fixnorm_mass.ml))
  yy=REPLICATE(fixnorm_mass[minpos].ininc,N_ELEMENTS(fixnorm_mass.ml))
  OPLOT,xx,yy,psym=8,thick=10,symsize=2,color=150

  stop

  set_plot,'x'
  ind=where(fixgauss_mass.inbeta eq 0.)
  fixgauss_mass=fixgauss_mass[ind]
  ind=where(fixnorm_mass.inbeta eq 0.)
  fixnorm_mass=fixnorm_mass[ind]
  ind=where(fixold_mass.inbeta eq 0.)
  fixold_mass=fixold_mass[ind]
  ind=where(freenorm_mass.inbeta eq 0.)
  freenorm_mass=freenorm_mass[ind]
  ind=where(freenorm.inbeta eq 0.)
  freenorm=freenorm[ind]
  ind=where(fixnorm.inbeta eq 0.)
  fixnorm=fixnorm[ind]

  likelihood=dblarr(nmbhs)
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(fixnorm_mass.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixnorm_mass[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixnorm_mass[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixnorm_mass[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixnorm_mass=likelihood/max(likelihood)
  mbhsfixnorm_mass=interpol(mbhs,likefixnorm_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  sigma=[-3,-2,-1,0,1,2,3]
  ml=findgen(30)*0.05+0.05
  nmls=n_elements(ml)
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(fixold_mass.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixold_mass[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixold_mass[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixold_mass[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixold_mass=likelihood/max(likelihood)
  mbhsfixold_mass=interpol(mbhs,likefixold_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])

  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(fixgauss_mass.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixgauss_mass[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixgauss_mass[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixgauss_mass[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixgauss_mass=likelihood/max(likelihood)
  mbhsfixgauss_mass=interpol(mbhs,likefixgauss_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(freenorm_mass.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(freenorm_mass[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(freenorm_mass[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(freenorm_mass[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreenorm_mass=likelihood/max(likelihood)
  mbhsfreenorm_mass=interpol(mbhs,likefreenorm_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(freenorm.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(freenorm[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(freenorm[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(freenorm[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreenorm=likelihood/max(likelihood)
  mbhsfreenorm=interpol(mbhs,likefreenorm,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(fixnorm.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixnorm[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixnorm[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixnorm[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixnorm=likelihood/max(likelihood)
  mbhsfixnorm=interpol(mbhs,likefixnorm,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  stop
;  set_plot,'ps' & LOADCT,0,/SILENT
;  device,filename='./cummulike.ps',/color
  stop


  mbhs=[findgen(17)*1.e5+3.4e6]
  ml=findgen(40)*0.05+0.05
  inclinations=60.
  beta=0.
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  ninclinations=n_elements(inclinations)
  final=REPLICATE({inmbh:0.0,chi2:0.0,ml:0.0,rms:fltarr(n_elements(disp))},nmbhs,nmls)
  readcol,'~/research/code/gemini15/vucd3/kinematics/newkinematic_psf_moffat.dat',normpsf,sigmapsf,format='F,F'
  readcol,'../evstigneeva/vucd3_mge_outputsersic.dat',intensity,sigmaarc,q,pa,format='F,F,F,F'
  readcol,'../evstigneeva/vucd3_mge_outputsersic_mass.dat',mass,format='D'
  surf_lum = intensity
  sigma_lum = sigmaarc
  qobs_lum = q
  surf_pot = mass
  sigma_pot = sigmaarc
  qobs_pot = q
;  for i=0,nmbhs-1 do begin
;     for j=0,nmls-1 do begin
;        fitml=ml[j]
;        mbh=mbhs[i]/fitml
;        jam_axisymmetric_rms,surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, inclinations, mbh, distance, xbin, ybin, rmsModel, BETA=beta,RMS=disp,ERMS=disperr,SIGMAPSF=sigmapsf,NORMPSF=normpsf,PIXSIZE=0.05,ml=fitml,chi2=chi2,STEP=0.02
;        final[i,j].chi2=chi2*FLOAT(n_elements(xbin))
;        final[i,j].inmbh=mbh*fitml
;        final[i,j].ml=fitml
;        final[i,j].rms=rmsmodel
;     endfor
;  endfor
;  mwrfits,final,'finalgrid.fits'
  stop
  final=mrdfits('finalgrid.fits',1)
  final=reform(final,17,40)
  stop
  mbhs=[0,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,3.4e6,findgen(4)*1.e5+3.6e6,findgen(4)*1.e5+4.1e6,findgen(4)*1.e5+4.6e6]
  mbhs=mbhs[sort(mbhs)]
  nmbhs=n_elements(mbhs)
  ml=findgen(40)*0.05+0.05
  nmls=n_elements(ml)
  ind=where(fixnorm_mass.ininc eq 60.)
  fixnorm_mass=fixnorm_mass[ind]
  fixnorm_mass=reform(fixnorm_mass,16,40)
  likelihood=dblarr(nmbhs)
  temp=0D
  for i=0,nmbhs-1 do begin
     if (mbhs[i] lt 3.4e6 OR mbhs[i] gt 5.e6) then begin
        mbhind=where(fixnorm_mass.inmbh gt mbhs[i]-1 and fixnorm_mass.inmbh lt mbhs[i]+1)
        for j=0,nmls-1 do begin
           mlind=where(fixnorm_mass[mbhind].ml eq ml[j])
           temp+=exp(-0.5*(fixnorm_mass[mbhind[mlind]].chi2))
        endfor
        likelihood[i]=temp
     endif else begin
        mbhind=where(final.inmbh gt mbhs[i]-1 and final.inmbh lt mbhs[i]+1)
        for j=0,nmls-1 do begin
           mlind=where(final[mbhind].ml eq ml[j])
           temp+=exp(-0.5*(final[mbhind[mlind]].chi2))
        endfor
        likelihood[i]=temp
     endelse
  endfor
  likefixnorm_mass=likelihood/max(likelihood)
  mbhsfixnorm_mass=interpol(mbhs,likefixnorm_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  print,'Best fit= ',mbhsfixnorm_mass[3]
  print,'1-sigma errors = ',mbhsfixnorm_mass[4]-mbhsfixnorm_mass[3],mbhsfixnorm_mass[3]-mbhsfixnorm_mass[2]
  print,'3-sigma errors = ',mbhsfixnorm_mass[6]-mbhsfixnorm_mass[3],mbhsfixnorm_mass[3]-mbhsfixnorm_mass[0]
;  myplot,filename='./cummulike.ps'
  set_plot,'ps'
  device,filename='./cummulike.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7
  LOADCT,0,/SILENT
  djs_plot,mbhs,likefixnorm_mass,ytitle='Cumulative Likelihood',xtitle='M_{BH} [M_{\odot}]',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=12,/nodata
  arrow,mbhsfixnorm_mass[3],0.,mbhsfixnorm_mass[3],1.,/data,thick=8.,color=100,hsize=0.01
  arrow,mbhsfixnorm_mass[2],0.,mbhsfixnorm_mass[2],1.,/data,thick=4,color=100,hsize=0.01
  arrow,mbhsfixnorm_mass[4],0.,mbhsfixnorm_mass[4],1.,/data,thick=4,color=100,hsize=0.01
  arrow,mbhsfixnorm_mass[0],0.,mbhsfixnorm_mass[0],1.,/data,thick=2,color=100,hsize=0.01
  arrow,mbhsfixnorm_mass[6],0.,mbhsfixnorm_mass[6],1.,/data,thick=2,color=100,hsize=0.01
  djs_oplot,mbhs,likefixnorm_mass,thick=12
  mbhs=[0,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6]

  djs_oplot,mbhs,likefixold_mass,thick=4,color='red',linestyle=1
  djs_oplot,mbhs,likefixgauss_mass,color='red',thick=4
  djs_oplot,mbhs,likefreenorm_mass,thick=4,color='blue'
  djs_oplot,mbhs,likefixnorm,thick=4,color='blue',linestyle=1
  djs_oplot,mbhs,likefreenorm,thick=4,color='blue',linestyle=2
  
  items=['Best Fit Model','PSF Variations','Model Variations']
  lines=[0,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/left
  xyouts,[6.5e6],[0.05],['VUCD3'],charthick=3,charsize=1.5,/data
  device,/close
  set_plot,'x'
  mbhs=[0,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,3.4e6,findgen(4)*1.e5+3.6e6,findgen(4)*1.e5+4.1e6,findgen(4)*1.e5+4.6e6]
  mbhs=mbhs[sort(mbhs)]
  nmbhs=n_elements(mbhs)

  likelihood=dblarr(nmls)
  temp=0.D
  mlpop=2.51623
  for i=0,nmls-1 do begin
     for j=0,nmbhs-1 do begin
        if (mbhs[j] lt 3.4e6 OR mbhs[j] gt 5.e6) then begin
           mlind=where(fixnorm_mass.ml eq ml[i])
           mbhind=where(fixnorm_mass[mlind].inmbh gt mbhs[j]-1 and fixnorm_mass[mlind].inmbh lt mbhs[j]+1)
           temp+=exp(-0.5*(fixnorm_mass[mlind[mbhind]].chi2))
        endif else begin
           mlind=where(final.ml eq ml[i])
           mbhind=where(final[mlind].inmbh gt mbhs[j]-1 and final[mlind].inmbh lt mbhs[j]+1)
           temp+=exp(-0.5*(final[mlind[mbhind]].chi2))
        endelse
     endfor
     likelihood[i]=temp
  endfor
  
  likefixnorm_mass=likelihood/max(likelihood)
  mlfixnorm_mass=(interpol(ml,likefixnorm_mass,[0.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865]))*mlpop
  print,'Best fit M/L = ',mlfixnorm_mass[3]
  print,'1-sigma error M/L = ',mlfixnorm_mass[4]-mlfixnorm_mass[3],mlfixnorm_mass[3]-mlfixnorm_mass[2]
  print,'3-sigma error M/L = ',mlfixnorm_mass[6]-mlfixnorm_mass[3],mlfixnorm_mass[3]-mlfixnorm_mass[0]
  likelihood=dblarr(nmls)
  temp=0D
  a=where(fixnorm_mass.inmbh eq 0.)
  for i=0,nmls-1 do begin
     ind=where((fixnorm_mass[a].ml gt ml[i]-0.01 and fixnorm_mass[a].ml lt ml[i]+0.01) OR fixnorm_mass[a].ml eq ml[i],c)
     if (c gt 1 or c lt 1) then stop
     temp+=exp(-0.5*(fixnorm_mass[a[ind]].chi2))
     likelihood[i]=temp
  endfor
  likefixnorm_nomass=likelihood/max(likelihood)
  mlfixnorm_nomass=(interpol(ml,likefixnorm_nomass,[0.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])) ;*mlpop
  print,'Best fit M/L ratio= ', mlfixnorm_nomass[3]
  print,'3-sigma error M/L ratio= ', mlfixnorm_nomass[6]-mlfixnorm_nomass[3],mlfixnorm_nomass[3]-mlfixnorm_nomass[0]
  stop
  zeropoint=25.28697d
  scale=0.025d
  Msun=4.53d
  mliin=1.69459
  mliout=3.16480
  rad=100.;8.744
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/deconvolved_f814.fits',oldimg,head
  hrotate,oldimg,head,img,newhead,1
  find_galaxy,img,majoraxis,eps,ang,xci,yci

  inind=where(q lt 0.8)
  inintensity=intensity[inind]
  insigma=sigmaarc[inind]
  inq=q[inind]
  inpa=pa[inind]
  mge2image,img,xci,yci,inintensity,insigma,inq,inpa,inmodel,zeropoint=zeropoint,scale=scale,msun=msun

  outind=where(q gt 0.8)
  outintensity=intensity[outind]
  outsigma=sigmaarc[outind]
  outq=q[outind]
  outpa=pa[outind]
  mge2image,img,xci,yci,outintensity,outsigma,outq,outpa,outmodel,zeropoint=zeropoint,scale=scale,msun=msun

  aper,inmodel,xci,yci,influx,influxerr,0.,skyerr,1,rad,-1,[1,1],/silent,setskyval=0.,/flux,/exact
  aper,outmodel,xci,yci,outflux,outfluxerr,0.,skyerr,1,rad,-1,[1,1],/silent,setskyval=0.,/flux,/exact

  bhmlind=where(final.chi2 eq min(final.chi2))
  bhml=final[bhmlind].ml
  inmag=zeropoint-2.5*alog10(influx)
  inabsmag=inmag-5*(alog10(16.5e6)-1)
  inlum=10^((inabsmag-msun)/(-2.5))
  outmag=zeropoint-2.5*alog10(outflux)
  outabsmag=outmag-5*(alog10(16.5e6)-1)
  outlum=10^((outabsmag-msun)/(-2.5))
  bhmass=bhml*((mliin*inlum)+(mliout*outlum))
  temp=where(fixnorm_mass.inmbh eq 0.)
  nobhmlind=where(fixnorm_mass[temp].chi2 eq min(fixnorm_mass[temp].chi2))
  nobhml=fixnorm_mass[temp[nobhmlind]].ml
  nobhmass=nobhml*((mliin*inlum)+(mliout*outlum))
  aveml=((influx*mliin)+(outflux*mliout))/(influx+outflux)
  print,'Total mass with black hole = ', final[bhmlind].inmbh, ' = ',bhmass
  print,'Dynamical M/L with black hole = ',bhml*aveml
  print,'Total mass w/o black hole = ', nobhmass
  print,'Dynamical M/L w/o black hole = ',nobhml*aveml
  print,'Total Luminosity = ', inlum+outlum
  intdisp=39700.
  r=((16.5e6)*(3.086e16))*((rad*0.025)*(4.848e-6))
  g=6.67e-11
  dynmass=(((intdisp^2)*r)/(g))/(1.989e30)
  print,'Total mass from integrated dispersion = ', dynmass
  set_plot,'ps'
  device,filename='./onedbestfit_BEST.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7
;  myplot,filename='./onedbestfit_BEST.ps'
  plotsym,0,1/2.,/fill
  djs_plot,xbin,disp,psym=8,ytitle='Dispersion [km s^{-1}]',xtitle='Radius ["]',xran=[0.,0.5],yran=[25,55],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  oploterror,xbin,disp,xbinerror,disperr,/NOHAT,psym=3
  a=where(fixnorm_mass.inmbh eq 0.)
  c=where(fixnorm_mass[a].chi2 eq min(fixnorm_mass[a].chi2))
  djs_oplot,xbin,fixnorm_mass[a[c]].rms,color='red',thick=4
  b=where(final.chi2 eq min(final.chi2))
  djs_oplot,xbin,final[b].rms,color='blue',thick=4
  fixnorm_mass=mrdfits('./temp/temp/fixnormpsf_mass.fits',1)
  mbhs=[0,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6]
  ml=findgen(40)*0.05+0.05
  inclinations=[60.,70.,80.,90.]
  betas=findgen(11)*0.1-0.2
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  ninclinations=n_elements(inclinations)
  nbetas=n_elements(betas)
;  a=where(fixnorm_mass.ininc eq 60.)
;  fixnorm_mass=fixnorm_mass[a]
  fixnorm_mass=reform(fixnorm_mass,nmbhs,nbetas,ninclinations,nmls)
  a=where(fixnorm_mass.outmbh eq 0.)
  b=where(fixnorm_mass[a].chi2 eq min(fixnorm_mass[a].chi2))
  help,fixnorm_mass[a[b]]
  djs_oplot,xbin,fixnorm_mass[a[b]].rms,color='grey',thick=4
  xyouts,[0.01],[26],['VUCD3'],charthick=3,charsize=1.5,/data
  device,/close
  set_plot,'x'
  stop

END
