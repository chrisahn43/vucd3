PRO COMPARE_JAM_OUTPUT
  infile='../../../kinematics/vor_out/intspec_wallace_best8.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  xbin=rav*0.05                 ;major axis*0.05 arcsec/pixel
  xbinerror=((rav-rin)/2.)*0.05
  ybin=fltarr(n_elements(rav)) ;minor axis
  distance=16.5
;  mbhs=[0.,10^(findgen(6)*0.2+6.)]
  mbhs=[0,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6]
  ml=findgen(30)*0.1+0.1
  inclinations=[60.,70.,80.,90.]
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  ninclinations=n_elements(inclinations)
  fixgauss_mass=mrdfits('fixgausspsf_mass.fits',1)
  fixnorm_mass=mrdfits('fixnormpsf_mass.fits',1)
  fixold_mass=mrdfits('fixoldpsf_mass.fits',1)
  freenorm_mass=mrdfits('freenormpsf_mass.fits',1)
  freenorm=mrdfits('freenormpsf.fits',1)
  fixnorm=mrdfits('fixnormpsf.fits',1)

                                ;COMPUTE LIKELIHOOD FUNCTION FOR EACH
                                ;START WITH FIX MODEL WITH MASS AND NEW PSF
  likelihood=dblarr(nmbhs)
  temp=0D
  for i=0,nmbhs-1 do begin
 ;    temp=0D
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
;  stop
  set_plot,'ps'
;  device,filename='./onedbestfit.ps',/color
;  plotsym,0,1/2.,/fill
;  djs_plot,xbin,disp,psym=8,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.5],yran=[25,55],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
;  oploterror,xbin,disp,xbinerror,disperr,/NOHAT,psym=3
;  a=where(fixnorm_mass.outmbh eq 0.)
;  c=where(fixnorm_mass[a].chi2 eq min(fixnorm_mass[a].chi2))
;  djs_oplot,xbin,fixnorm_mass[a[c]].rms,color='red',thick=4
;  b=where(fixnorm_mass.chi2 eq min(fixnorm_mass.chi2))
;  djs_oplot,xbin,fixnorm_mass[b].rms,color='blue',thick=4
;  device,/close
  temp=0D
  for i=0,nmbhs-1 do begin
 ;    temp=0D
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


  ;stop
                                ;FIX MODEL WITH MASS GAUSS PSF
  
  temp=0D
  for i=0,nmbhs-1 do begin
 ;    temp=0D
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


  ;stop
                                ;FREE MODEL WITH MASS WITH NEW PSF
  temp=0D
  for i=0,nmbhs-1 do begin
 ;    temp=0D
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


                                ;stop
  ml=findgen(26)*0.1+0.5
  nmls=n_elements(ml)
                                ;FREE MODEL WITH NEW PSF NO MASS
  temp=0D
  for i=0,nmbhs-1 do begin
 ;    temp=0D
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


  ;stop
                                ;FIX MODEL WITH NEW PSF NO MASS
  temp=0D
  for i=0,nmbhs-1 do begin
 ;    temp=0D
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
  device,filename='./cummulike.ps',/color
                                ;FIX MODEL WITH MASS OLD PSF
  djs_plot,mbhs,likefixnorm_mass,ytitle='Cummulative Likelihood',xtitle='M_{BH}',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=12
  djs_oplot,mbhs,likefixold_mass,color='red',thick=4
  djs_oplot,mbhs,likefixgauss_mass,color='red',thick=4,linestyle=1
  djs_oplot,mbhs,likefreenorm_mass,color='blue',thick=4
  djs_oplot,mbhs,likefixnorm,color='blue',thick=4,linestyle=1
  djs_oplot,mbhs,likefreenorm,color='blue',thick=4,linestyle=2
  items=['Best Fit Model','PSF Variations','Model Variations']
  lines=[0,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/left
  device,/close

  
;  djs_oplot,mbhs,likefixold_mass,color='blue',thick=4
;  djs_oplot,mbhs,likefixgauss_mass,color='green',thick=4
;  djs_oplot,mbhs,likefreenorm_mass,color='purple',thick=4
;  djs_oplot,mbhs,likefreenorm,color='red',thick=4
;  djs_oplot,mbhs,likefixnorm,color='cyan',thick=4
;  items=['Fix Mass New PSF','Fix Mass Old PSF','Fix Mass Gauss PSF','Free Mass New PSF','Free New PSF','Fix New PSF']
;  lines=[0,0,0,0,0,0]
;  color=['black','blue','green','purple','red','cyan']
;  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/bottom,/right
;  device,/close
  device,filename='./sigmambhs_psf.ps',/color
  djs_plot,sigma,mbhsfixnorm_mass,ytitle='M_{BH}',xtitle='\sigma',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=12,xran=[-3.5,3.5],yran=[-150000,7.e6],/xsty,/ysty
  djs_oplot,sigma,mbhsfixold_mass,color='red',thick=4
  djs_oplot,sigma,mbhsfixgauss_mass,color='red',thick=4,linestyle=1
  djs_oplot,sigma,mbhsfreenorm_mass,color='blue',thick=4
  djs_oplot,sigma,mbhsfixnorm,color='blue',thick=4,linestyle=1
  djs_oplot,sigma,mbhsfreenorm,color='blue',thick=4,linestyle=2
  items=['Best Fit Model','PSF Variations','Model Variations']
  lines=[0,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/left
  device,/close
  
  
  stop

  mbhs=[0.,4254896.9]
  ml=findgen(30)*0.1+0.1
  inclinations=90.
  beta=0.
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  ninclinations=n_elements(inclinations)
  final=REPLICATE({inmbh:0.0,chi2:0.0,ml:0.0,rms:fltarr(n_elements(disp))},nmbhs,nmls)
  readcol,'~/research/code/gemini15/vucd3/kinematics/newkinematic_psf_moffat.dat',normpsf,sigmapsf,format='F,F'
  readcol,'../../evstigneeva/vucd3_mge_outputsersic.dat',intensity,sigmaarc,q,pa,format='F,F,F,F'
  readcol,'../../evstigneeva/vucd3_mge_outputsersic_mass.dat',mass,format='D'
  surf_lum = intensity
  sigma_lum = sigmaarc
  qobs_lum = q
  surf_pot = mass
  sigma_pot = sigmaarc
  qobs_pot = q
  for i=0,nmbhs-1 do begin
     for j=0,nmls-1 do begin
        fitml=ml[j]
        mbh=mbhs[i]/fitml
        jam_axisymmetric_rms,surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, inclinations, mbh, distance, xbin, ybin, rmsModel, BETA=beta,RMS=disp,ERMS=disperr,SIGMAPSF=sigmapsf,NORMPSF=normpsf,PIXSIZE=0.05,ml=fitml,chi2=chi2,STEP=0.02
        final[i,j].chi2=chi2*FLOAT(n_elements(xbin))
        final[i,j].inmbh=mbh*fitml
        final[i,j].ml=fitml
        final[i,j].rms=rmsmodel
     endfor
  endfor
;  stop
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
  temp=where(final.inmbh eq 0.)
  nobhmlind=where(final[temp].chi2 eq min(final[temp].chi2))
  nobhml=final[temp[nobhmlind]].ml
  nobhmass=nobhml*((mliin*inlum)+(mliout*outlum))
  aveml=((influx*mliin)+(outflux*mliout))/(influx+outflux)
  print,'Total mass with black hole = ', final[bhmlind].inmbh, ' = ',bhmass
  print,'Dynamical M/L with black hole = ',bhml*aveml
  print,'Total mass w/o black hole = ', nobhmass
  print,'Dynamical M/L w/o black hole = ',nobhml*aveml
  intdisp=39700.
  r=((16.5e6)*(3.086e16))*((rad*0.025)*(4.848e-6))
  g=6.67e-11
  dynmass=(((intdisp^2)*r)/(g))/(1.989e30)
  print,'Total mass from integrated dispersion = ', dynmass
;  stop
;  set_plot,'ps'
  device,filename='./onedbestfit_BEST.ps',/color
  plotsym,0,1/2.,/fill
  djs_plot,xbin,disp,psym=8,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.5],yran=[25,55],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  oploterror,xbin,disp,xbinerror,disperr,/NOHAT,psym=3
  a=where(final.inmbh eq 0.)
  c=where(final[a].chi2 eq min(final[a].chi2))
  djs_oplot,xbin,final[a[c]].rms,color='red',thick=4
  b=where(final.chi2 eq min(final.chi2))
  djs_oplot,xbin,final[b].rms,color='blue',thick=4
  device,/close
  set_plot,'x'
  stop

END
