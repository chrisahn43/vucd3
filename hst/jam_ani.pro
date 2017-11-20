PRO JAM_ANI
;THIS CODE RUNS JAM CODE AND WITH VARYING ANISOTROPY TO SEE AT WHAT
;POINT WE WILL HAVE A ZERO MASS BLACK HOLE.

  infile='../kinematics/vor_out/intspec_wallace_best8.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  xbin=rav*0.05                 ;major axis*0.05 arcsec/pixel
  ybin=fltarr(n_elements(rav)) ;minor axis
  distance=16.5
  ml=findgen(30)*0.1+.1
  mbhs=[0.]
  betas=findgen(20)*0.1-1.
  inclinations=[60.,70.,80.,90.]
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  nbetas=n_elements(betas)
  ninclinations=n_elements(inclinations)
  readcol,'./../kinematics/newkinematic_psf_moffat.dat',normpsf,sigmapsf,format='F,F'
  fixnormpsf_mass=REPLICATE({inmbh:0.0,inbeta:0.0,ininc:0.0,chi2:0.0,ml:0.0,outmbh:0.0,rms:fltarr(n_elements(disp))},nmbhs,nbetas,ninclinations,nmls)

  readcol,'./evstigneeva/vucd3_mge_outputsersic.dat',Intensity,sigmaarc,q,format='F,F,F'
  readcol,'./evstigneeva/vucd3_mge_outputsersic_mass.dat',mass,format='D'
  surf_lum = intensity
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
  mwrfits,fixnormpsf_mass,'./JAM_output_anisotropy/fixnormpsf_mass.fits'
  a=where(fixnormpsf_mass.chi2 eq min(fixnormpsf_mass.chi2))
  rms_anis=fixnormpsf_mass[a].rms
  b=where(fixnormpsf_mass.inbeta eq 0.)
  c=where(fixnormpsf_mass[b].chi2 eq min(fixnormpsf_mass[b].chi2))
  rms_noanis=fixnormpsf_mass[b[c]].rms

  djs_plot,xbin,disp,psym=4,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.52],yran=[25,55],/xsty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  oploterr,xbin,disp,disperr
  djs_oplot,xbin,rms_noanis,color='red',thick=4
  djs_oplot,xbin,rms_anis,color='blue',thick=4
  help,fixnormpsf_mass[a]
  help,fixnormpsf_mass[b[c]]
  stop

END
