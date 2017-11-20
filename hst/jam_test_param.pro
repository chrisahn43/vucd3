PRO JAM_TEST_PARAM
  infile='../kinematics/vor_out/intspec_wallace_best8.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  xbin=rav*0.05                ;major axis*0.05 arcsec/pixel
  ybin=fltarr(n_elements(rav)) ;minor axis
  distance=16.5
  ml=findgen(40)*0.05+.05
;  mbhs = [0.,10^(findgen(6)*0.2+6.)] ; Black hole mass in solar masses
  mbhs=[4.4e6];,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6]
  betas=[0.]
  inclinations=[60.];,70.,80.,90.]
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  nbetas=n_elements(betas)
  ninclinations=n_elements(inclinations)

  ;START WITH BEST PSF THAT WE DERIVED LUM ONLY
  readcol,'./newkinematic_psf_moffat.dat',normpsf,sigmapsf,format='F,F'

  fixnormpsf_mass=REPLICATE({inmbh:0.0,inbeta:0.0,ininc:0.0,chi2:0.0,ml:0.0,outmbh:0.0,rms:fltarr(n_elements(disp))},nmbhs,nbetas,ninclinations,nmls)

;  readcol,'./evstigneeva/vucd3_mge_outputsersic.dat',Intensity,sigmaarc,q,format='F,F,F'
  readcol,'./make_decon/vucd3_mge_outputsersic_k.dat',intensity,sigmaarc,q,format='F,F,F'
  readcol,'./evstigneeva/vucd3_mge_outputsersic_nmass.dat',mass,format='D'
  surf_lum = intensity;/10.
  sigma_lum = sigmaarc
  qobs_lum = q
  surf_pot = mass
  sigma_pot = sigmaarc
  qobs_pot = q

  for i=0,nmbhs-1 do begin
     for j=0,nbetas-1 do begin
        for k=0,ninclinations-1 do begin
           for l=0,nmls-1 do begin
              fitml=ml[l]
              mbh=mbhs[i]/fitml
              beta=betas[j]
              inc=inclinations[k]
              jam_axisymmetric_rms,surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, inc, mbh, distance, xbin, ybin, rmsModel, BETA=beta,RMS=disp,ERMS=disperr,SIGMAPSF=sigmapsf,NORMPSF=normpsf,PIXSIZE=0.05,ml=fitml,chi2=chi2,STEP=0.02
              fixnormpsf_mass[i,j,k,l].chi2=chi2*FLOAT(n_elements(xbin))
              fixnormpsf_mass[i,j,k,l].outmbh=mbh*fitml
              fixnormpsf_mass[i,j,k,l].ml=fitml
              fixnormpsf_mass[i,j,k,l].inmbh=mbhs[i]
              fixnormpsf_mass[i,j,k,l].inbeta=betas[j]
              fixnormpsf_mass[i,j,k,l].ininc=inclinations[k]
              fixnormpsf_mass[i,j,k,l].rms=rmsmodel
           endfor
        endfor
     endfor
  endfor
  stop
  mwrfits,fixnormpsf_mass,'./JAM_output_fixedml/fixnormpsf_mass.fits'
  a=where(fixnormpsf_mass.outmbh eq 0.)
  c=where(fixnormpsf_mass[a].chi2 eq min(fixnormpsf_mass[a].chi2))
  rms_nobh_freenorm=fixnormpsf_mass[a[c]].rms
  b=where(fixnormpsf_mass.chi2 eq min(fixnormpsf_mass.chi2))
  rms_freenorm=fixnormpsf_mass[b].rms
  djs_plot,xbin,disp,psym=4,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.52],yran=[25,55],/xsty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  oploterr,xbin,disp,disperr
  djs_oplot,xbin,rms_nobh_freenorm,color='red',thick=4
  djs_oplot,xbin,rms_freenorm,color='blue',thick=4
  help,fixnormpsf_mass[a[c]]
  help,fixnormpsf_mass[b]
  stop

END
