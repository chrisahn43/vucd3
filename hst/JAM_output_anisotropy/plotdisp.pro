PRO PLOTDISP

  infile='../../kinematics/vor_out/intspec_wallace_best8.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  xbin=rav*0.05                 ;major axis*0.05 arcsec/pixel
  xbinerror=((rav-rin)/2.)*0.05

  data=mrdfits('fixnormpsf_mass.fits',1)
  set_plot,'ps'
  device,filename='disp_anisotropy.ps',/color
  plotsym,0,1/2.,/fill
  djs_plot,xbin,disp,psym=8,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.5],yran=[25,55],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  oploterror,xbin,disp,xbinerror,disperr,/NOHAT,psym=3
  a=where(data.chi2 eq min(data.chi2))
  rms_anis=data[a].rms
  b=where(data.inbeta eq 0.)
  c=where(data[b].chi2 eq min(data[b].chi2))
  rms_noanis=data[b[c]].rms

  djs_oplot,xbin,rms_noanis,color='red',thick=4
  djs_oplot,xbin,rms_anis,color='blue',thick=4
  device,/close
  set_plot,'x'
  
END
