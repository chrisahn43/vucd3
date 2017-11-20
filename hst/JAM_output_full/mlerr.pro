PRO MLERR
  infile='../../kinematics/vor_out/intspec_wallace_best8.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  xbin=rav*0.05                 ;major axis*0.05 arcsec/pixel
  xbinerror=((rav-rin)/2.)*0.05
  ybin=fltarr(n_elements(rav)) ;minor axis
  distance=16.5
;  mbhs=[0.,10^(findgen(6)*0.2+6.)]
  mbhs=[0,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6]
  ml=findgen(30)*0.05+0.05
  inclinations=[60.,70.,80.,90.]
  betas=findgen(11)*0.1-0.2
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  ninclinations=n_elements(inclinations)
  fixgauss_mass=mrdfits('fixgausspsf_mass.fits',1)
  fixnorm_mass=mrdfits('fixnormpsf_mass.fits',1)
  fixold_mass=mrdfits('fixoldpsf_mass.fits',1)
  freenorm_mass=mrdfits('freenormpsf_mass.fits',1)
  freenorm=mrdfits('freenormpsf.fits',1)
  fixnorm=mrdfits('fixnormpsf.fits',1)
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
  mlpop=2.51623
  likelihood=dblarr(nmls)
  temp=0D
  for i=0,nmls-1 do begin
     mlind=where(fixnorm_mass.ml eq ml[i])
     for j=0,nmbhs-1 do begin
        mbhind=where(fixnorm_mass[mlind].outmbh eq mbhs[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixnorm_mass[mlind[mbhind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixnorm_mass[mlind[mbhind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixnorm_mass=likelihood/max(likelihood)
  mlfixnorm_mass=(interpol(ml,likefixnorm_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865]))*mlpop
  sigma=[-3,-2,-1,0,1,2,3]
  temp=0D
  for i=0,nmls-1 do begin
     mlind=where(fixold_mass.ml eq ml[i])
     for j=0,nmbhs-1 do begin
        mbhind=where(fixold_mass[mlind].outmbh eq mbhs[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixold_mass[mlind[mbhind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixold_mass[mlind[mbhind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixold_mass=likelihood/max(likelihood)
  mlfixold_mass=(interpol(ml,likefixold_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865]))*mlpop
  temp=0D
  for i=0,nmls-1 do begin
     mlind=where(fixgauss_mass.ml eq ml[i])
     for j=0,nmbhs-1 do begin
        mbhind=where(fixgauss_mass[mlind].outmbh eq mbhs[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixgauss_mass[mlind[mbhind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixgauss_mass[mlind[mbhind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixgauss_mass=likelihood/max(likelihood)
  mlfixgauss_mass=(interpol(ml,likefixgauss_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865]))*mlpop
  temp=0D
  for i=0,nmls-1 do begin
     mlind=where(freenorm_mass.ml eq ml[i])
     for j=0,nmbhs-1 do begin
        mbhind=where(freenorm_mass[mlind].outmbh eq mbhs[j])
        for k=0,ninclinations-1 do begin
           incind=where(freenorm_mass[mlind[mbhind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(freenorm_mass[mlind[mbhind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreenorm_mass=likelihood/max(likelihood)
  mlfreenorm_mass=(interpol(ml,likefreenorm_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865]))*mlpop

;  ml=findgen(26)*0.1+0.5
;  nmls=n_elements(ml)
;  likelihood=dblarr(nmls)
  temp=0D
  for i=0,nmls-1 do begin
     mlind=where(freenorm.ml eq ml[i])
     for j=0,nmbhs-1 do begin
        mbhind=where(freenorm[mlind].outmbh eq mbhs[j])
        for k=0,ninclinations-1 do begin
           incind=where(freenorm[mlind[mbhind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(freenorm[mlind[mbhind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreenorm=likelihood/max(likelihood)
  mlfreenorm=interpol(ml,likefreenorm,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  temp=0D
  for i=0,nmls-1 do begin
     mlind=where(fixnorm.ml eq ml[i])
     for j=0,nmbhs-1 do begin
        mbhind=where(fixnorm[mlind].outmbh eq mbhs[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixnorm[mlind[mbhind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixnorm[mlind[mbhind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixnorm=likelihood/max(likelihood)
  mlfixnorm=interpol(ml,likefixnorm,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  stop
  set_plot,'ps'
  device,filename='./cummulike_ml.ps',/color
  
  djs_plot,ml*mlpop,likefixnorm_mass,ytitle='Cummulative Likelihood',xtitle='M/L',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=12,xran=[0.,3.],xstyle=8,/nodata
  mlratio=(ml/mlpop)
  xt='M/L_{dyn} / M/L_{pop}'
  xt=textoidl(xt)

  axis,xaxis=1,xtickv=mlratio,xcharsize=1.5,charthick=4,xthick=4,xtitle=xt,/xsty,xran=[0.0,1.19226]

  arrow,1.72,0.,1.72,1.,/data,thick=8,color=100,hsize=0.01
  arrow,1.25,0.,1.25,1.,/data,thick=4,color=100,hsize=0.01
  arrow,2.26,0.,2.26,1.,/data,thick=4.,color=100,hsize=0.01
  arrow,0.53,0.,0.53,1.,/data,thick=4.,color=100,hsize=0.01;0.84
  arrow,3.26,0.,3.26,1.,/data,thick=2.,color=100,hsize=0.01;0.00135

  djs_oplot,ml*mlpop,likefixnorm_mass,thick=12
  djs_oplot,ml*mlpop,likefixnorm_mass,thick=12
  djs_oplot,ml*mlpop,likefixold_mass,thick=4,color='red'
  djs_oplot,ml*mlpop,likefixgauss_mass,color='red',thick=4,linestyle=1
  djs_oplot,ml*mlpop,likefreenorm_mass,thick=4,color='blue'
  djs_oplot,ml,likefixnorm,thick=4,color='blue',linestyle=1
  djs_oplot,ml,likefreenorm,thick=4,color='blue',linestyle=2

  items=['Best Fit Model','PSF Variations','Model Variations']
  lines=[0,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/left
  device,/close
  stop
END
