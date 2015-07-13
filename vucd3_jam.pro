pro vucd3_jam
  infile='vucd3_radial.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  readcol,'vucd3_mge_output.dat',Intensity,sigmaarc,q,format='F,F,F'

  xbin=rav*0.05                ;major axis*0.05 arcsec/pixel
  ybin=fltarr(n_elements(rav)) ;minor axis 
;Luminosity should be 10^7 L sun, mass 65 * 10^6
  surf_lum = intensity
  sigma_lum = sigmaarc
  qobs_lum = q
  surf_pot = intensity
  sigma_pot = sigmaarc
  qobs_pot = q
  inc=90d
  distance = 16.5              ; Assume Virgo distance in Mpc (Mei et al. 2007)
  mbhs = [0,10^(findgen(16)*0.2+5.0)] ; Black hole mass in solar masses
  betas=findgen(10)*0.2-1.
  inclinations=[50.,60.,70.,80.,90.]
  nmbhs=n_elements(mbhs)
  nbetas=n_elements(betas)
  ninclinations=n_elements(inclinations)
  out=REPLICATE({inmbh:0.0,inbeta:0.0,ininc:0.0,chi2:0.0,ml:0.0,outmbh:0.0},nmbhs,nbetas,ninclinations)
  inmbh=fltarr(nmbhs,nbetas,ninclinations)
  inbeta=FLTARR(nmbhs,nbetas,ninclinations)
  ininclination=FLTARR(nmbhs,nbetas,ninclinations)
  outchi2=FLTARR(nmbhs,nbetas,ninclinations)
  outmbh=FLTARR(nmbhs,nbetas,ninclinations)
  outml=FLTARR(nmbhs,nbetas,ninclinations)
  for i=0,nmbhs-1 do begin
     for j=0,nbetas-1 do begin
        for k=0,ninclinations-1 do begin
           fitml=-1
           mbh=mbhs[i]
           beta=betas[j]
           inc=inclinations[k]
           jam_axisymmetric_rms,surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, inc, mbh, distance, xbin, ybin, rmsModel, BETA=beta,RMS=disp,ERMS=disperr,SIGMAPSF=[0.118,0.562],NORMPSF=[0.63,0.37],PIXSIZE=0.05,ml=fitml,chi2=chi2
           out[i,j,k].chi2=chi2*FLOAT(n_elements(xbin))
           out[i,j,k].outmbh=mbh*fitml
           out[i,j,k].ml=fitml
           out[i,j,k].inmbh=mbhs[i]
           out[i,j,k].inbeta=betas[j]
           out[i,j,k].ininc=inclinations[k]
        endfor
     endfor
  endfor
  stop

  ;assume 5 degrees of freedom
  djs_plot,alog10(out.outmbh),out.inbeta,psym=3,xtitle='log(MBH)',ytitle='!9b!x',chars=2,xrange=[6.0,7.5],/xsty

  sig1ind=where(out.chi2 lt MPCHILIM(0.68,5))
  djs_oplot,ALOG10(out[sig1ind].outmbh),out[sig1ind].inbeta,psym=4,thick=5

  sig2ind=WHERE(out.chi2 LT MPCHILIM(0.95,5) AND out.chi2 GE MPCHILIM(0.68,5))
  djs_oplot,ALOG10(out[sig2ind].outmbh),out[sig2ind].inbeta,psym=5,thick=4

  sig3ind=WHERE(out.chi2 LT MPCHILIM(0.997,5) AND out.chi2 GE MPCHILIM(0.95,5))
  djs_oplot,ALOG10(out[sig3ind].outmbh),out[sig3ind].inbeta,psym=1,thick=2
  stop
  ;ML-ratio
  djs_plot,out.ml,alog10(out.outmbh),psym=3,ytitle='log(MBH)',xtitle='M/L_i',chars=2,ysty=1,yrange=[6.0,7.5]
  djs_oplot,out[sig1ind].ml,ALOG10(out[sig1ind].outmbh),psym=4;,thick=5
  djs_oplot,out[sig2ind].ml,ALOG10(out[sig2ind].outmbh),psym=5,thick=4
  djs_oplot,out[sig3ind].ml,ALOG10(out[sig3ind].outmbh),psym=1,thick=2
  stop
; The model is similar but not identical to the adopted kinematics!
;SigmapSF and NORMPSF are for kinematic data
;  jam_axisymmetric_rms,surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, inc, mbh, distance, xbin, ybin, rmsModel, BETA=0., RMS=disp, ERMS=disperr, SIGMAPSF=[0.118,0.562], NORMPSF=[0.63,0.37], PIXSIZE=0.05, FLUX=flux,ml=-1.0
;    jam_axisymmetric_rms,surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, inc, mbh, distance, xbin, ybin, rmsModel2, BETA=0., RMS=disp, ERMS=disperr, SIGMAPSF=[0.2], NORMPSF=[1.], PIXSIZE=0.05, FLUX=flux,ml=-1.0

;  print,rmsmodel,disp
;  djs_plot,xbin,disp,psym=2,xran=[-0.02,0.52],yran=[28,55],/xstyle,/ystyle
;  djs_oplot,xbin,rmsmodel,psym=2,color='blue'
;  djs_oplot,xbin,rmsmodel2,psym=2,color='green'
;  print, MAX(rmsmodel)/MIN(rmsmodel)
;  stop


END