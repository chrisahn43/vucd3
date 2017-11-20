PRO COMPARE_JAM_OUTPUT
  infile='../../kinematics/vor_out/intspec_wallace_best8.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  xbin=rav*0.05                ;major axis*0.05 arcsec/pixel
  ybin=fltarr(n_elements(rav)) ;minor axis
  distance=16.5
  mbhs=[0.,10^(findgen(6)*0.2+6.)]
  ml=findgen(26)*0.1+0.1
  inclinations=[60.,70.,80.,90.]
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  ninclinations=n_elements(inclinations)
  fixgauss_mass=mrdfits('first_fit/fixgausspsf_mass.fits',1)
  fixnorm_mass=mrdfits('first_fit/fixnormpsf_mass.fits',1)
  fixold_mass=mrdfits('first_fit/fixoldpsf_mass.fits',1)
  freenorm_mass=mrdfits('first_fit/freenormpsf_mass.fits',1)
  freenorm=mrdfits('first_fit/freenormpsf.fits',1)
  fixnorm=mrdfits('first_fit/fixnormpsf.fits',1)

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
           temp+=exp(-0.5*(fixnorm_mass[mbhind[mlind[incind]]].chi2)^2)
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixnorm_mass=likelihood/max(likelihood)
  mbhsfixnorm_mass=interpol(mbhs,likefixnorm_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  sigma=[-3,-2,-1,0,1,2,3]
  stop
  set_plot,'ps'
  device,filename='./first_fit/onedbestfit.ps',/color
  djs_plot,xbin,disp,psym=4,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.5],yran=[25,55],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  oploterr,xbin,disp,disperr
  a=where(fixnorm_mass.outmbh eq 0.)
  c=where(fixnorm_mass[a].chi2 eq min(fixnorm_mass[a].chi2))
  djs_oplot,xbin,fixnorm_mass[a[c]].rms,color='red',thick=4
  b=where(fixnorm_mass.chi2 eq min(fixnorm_mass.chi2))
  djs_oplot,xbin,fixnorm_mass[b].rms,color='blue',thick=4
  device,/close
  device,filename='./first_fit/cummulike.ps',/color
                                ;FIX MODEL WITH MASS OLD PSF
  djs_plot,mbhs,likefixnorm_mass,ytitle='Cummulative Likelihood',xtitle='M_{BH}',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=4
  temp=0D
  for i=0,nmbhs-1 do begin
 ;    temp=0D
     mbhind=where(fixold_mass.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixold_mass[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixold_mass[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixold_mass[mbhind[mlind[incind]]].chi2)^2)
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixold_mass=likelihood/max(likelihood)
  mbhsfixold_mass=interpol(mbhs,likefixold_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])

  djs_oplot,mbhs,likefixold_mass,color='blue',thick=4

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
           temp+=exp(-0.5*(fixgauss_mass[mbhind[mlind[incind]]].chi2)^2)
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixgauss_mass=likelihood/max(likelihood)
  mbhsfixgauss_mass=interpol(mbhs,likefixgauss_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])

  djs_oplot,mbhs,likefixgauss_mass,color='green',thick=4

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
           temp+=exp(-0.5*(freenorm_mass[mbhind[mlind[incind]]].chi2)^2)
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreenorm_mass=likelihood/max(likelihood)
  mbhsfreenorm_mass=interpol(mbhs,likefreenorm_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])

  djs_oplot,mbhs,likefreenorm_mass,color='purple',thick=4

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
           temp+=exp(-0.5*(freenorm[mbhind[mlind[incind]]].chi2)^2)
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreenorm=likelihood/max(likelihood)
  mbhsfreenorm=interpol(mbhs,likefreenorm,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])

  djs_oplot,mbhs,likefreenorm,color='red',thick=4

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
           temp+=exp(-0.5*(fixnorm[mbhind[mlind[incind]]].chi2)^2)
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixnorm=likelihood/max(likelihood)
  mbhsfixnorm=interpol(mbhs,likefixnorm,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])

  djs_oplot,mbhs,likefixnorm,color='cyan',thick=4
  items=['Fix Mass New PSF','Fix Mass Old PSF','Fix Mass Gauss PSF','Free Mass New PSF','Free New PSF','Fix New PSF']
  lines=[0,0,0,0,0,0]
  color=['black','blue','green','purple','red','cyan']
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/bottom,/right
  device,/close
  device,filename='./first_fit/sigmambhs.ps',/color
  djs_plot,sigma,mbhsfixnorm_mass,ytitle='M_{BH}',xtitle='\sigma',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=4,xran=[-3.5,3.5],yran=[-150000,7.e6],/xsty,/ysty
  djs_oplot,sigma,mbhsfixold_mass,color='blue',thick=4
  djs_oplot,sigma,mbhsfixgauss_mass,color='green',thick=4
  djs_oplot,sigma,mbhsfreenorm_mass,color='purple',thick=4
  djs_oplot,sigma,mbhsfreenorm,color='red',thick=4
  djs_oplot,sigma,mbhsfixnorm,color='cyan',thick=4
  sym=[0,0,0,0,0,0]
  al_legend,items,psym=sym,colors=color,background_color='white',charthick=4,thick=3,/top,/left
  device,/close
  
  set_plot,'x'
  
  stop
  
  
END
