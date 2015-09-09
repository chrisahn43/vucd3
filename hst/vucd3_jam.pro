pro vucd3_jam
  infile='../kinematics/vor_out/intspec_wallace_best8.dat'
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
  distance = 16.5              ; Assume Virgo distance in Mpc (Mei et al. 2007)
  ml=findgen(50)*0.1+0.1
  mbhs = [0,10^(findgen(11)*0.2+5.0)] ; Black hole mass in solar masses
  betas=[0.];findgen(10)*0.2-1.
  inclinations=[90.];[50.,60.,70.,80.,90.]
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  nbetas=n_elements(betas)
  ninclinations=n_elements(inclinations)
  out=REPLICATE({inmbh:0.0,inbeta:0.0,ininc:0.0,chi2:0.0,ml:0.0,outmbh:0.0,rms:fltarr(n_elements(disp))},nmbhs,nbetas,ninclinations,nmls)
  for i=0,nmbhs-1 do begin
     for j=0,nbetas-1 do begin
        for k=0,ninclinations-1 do begin
           for l=0,nmls-1 do begin
              fitml=ml[l]
              mbh=mbhs[i]/fitml
              beta=betas[j]
              inc=inclinations[k]
              jam_axisymmetric_rms,surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, inc, mbh, distance, xbin, ybin, rmsModel, BETA=beta,RMS=disp,ERMS=disperr,SIGMAPSF=[0.118,0.562],NORMPSF=[0.63,0.37],PIXSIZE=0.05,ml=fitml,chi2=chi2
              out[i,j,k,l].chi2=chi2*FLOAT(n_elements(xbin))
              out[i,j,k,l].outmbh=mbh*fitml
              out[i,j,k,l].ml=fitml
              out[i,j,k,l].inmbh=mbhs[i]
              out[i,j,k,l].inbeta=betas[j]
              out[i,j,k,l].ininc=inclinations[k]
              out[i,j,k,l].rms=rmsmodel
           endfor
        endfor
     endfor
  endfor
  a=where(out.outmbh eq 0.)
  c=where(out[a].chi2 eq min(out[a].chi2))
  rms_nobh=out[a[c]].rms
  b=where(out.chi2 eq min(out.chi2))
  set_plot,'ps'
  device,filename='oned_bestfit_rms.ps',/color
  djs_plot,xbin,disp,psym=4,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.52],yran=[25,55],/xsty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  oploterr,xbin,disp,disperr
  djs_oplot,xbin,rms_nobh,color='red',thick=4
  djs_oplot,xbin,out[b].rms,color='blue',thick=4
  help,out[a[c]]
  help,out[b]
  device,/close
  set_plot,'x'
  stop
  bestmbh=fltarr(n_elements(mbhs))
  bestml=fltarr(n_elements(mbhs))
  bestchi2=fltarr(n_elements(mbhs))
  for i=0,n_elements(mbhs)-1 do begin
     ind=where(out.outmbh eq mbhs[i])
     index=where(out[ind].chi2 eq min(out[ind].chi2))
     bestmbh[i]=out[ind[index]].outmbh
     bestml[i]=out[ind[index]].ml
     bestchi2[i]=out[ind[index]].chi2
     if (bestmbh[i] eq 0.) then bestmbh[i]=1.
   endfor

  set_plot,'ps'
  device,filename='best_blackholevsml.ps',/color
  djs_plot,bestml,alog10(bestmbh),psym=4,xran=[0.9,4.5],yran=[4.5,7.5],/xsty,/ysty,xtitle='M/L_i',ytitle='Log(MBH)',charsize=1.5,charthick=4,xthick=3,ythick=3
  device,/close
  device,filename='best_chi2vsml.ps',/color
  djs_plot,bestml,bestchi2,psym=4,xran=[0.9,4.6],yran=[0.,3.5],/xsty,/ysty,xtitle='M/L_I',ytitle='(\chi)^2',charsize=1.5,charthick=4,xthick=3,ythick=3
  device,/close
  set_plot,'x'
;  items=['chi2 ~1.132 Mbh = 6.3e6 M/L=2.1','chi2=0.636 Mbh = 6.3e6 M/L=2.3','Mbh = 0 M/L=2.1']
;  colors=['blue','green','red']
;  legend,items,color=colors
  stop
;  p=plot(xbin,disp,symbol='D',sym_size=1,linstyle=6)
;  err=plot(xbin,err,/overplot)
;  p1=plot(xbin,out[a].rms,color='blue')
;  p2=plot(xbin,out[b].rms,color='green')
;  p3=plot(xbin,rmsmodel_nobh,color='red')

                                ;assume 5 degrees of freedom
  ;set_plot,'ps'
  ;device,filename='mbhbeta.ps',/color
  djs_plot,alog10(out.outmbh),out.inbeta,psym=3,xtitle='log(MBH)',ytitle='\beta',chars=2,xrange=[6.0,7.5],/xsty,charsize=1.5,charthick=4,xthick=3,ythick=3

  sig1ind=where(out.chi2 lt MPCHILIM(0.68,5))
  djs_oplot,ALOG10(out[sig1ind].outmbh),out[sig1ind].inbeta,psym=4,thick=5

  sig2ind=WHERE(out.chi2 LT MPCHILIM(0.95,5) AND out.chi2 GE MPCHILIM(0.68,5))
  djs_oplot,ALOG10(out[sig2ind].outmbh),out[sig2ind].inbeta,psym=5,thick=4

  sig3ind=WHERE(out.chi2 LT MPCHILIM(0.997,5) AND out.chi2 GE MPCHILIM(0.95,5))
  djs_oplot,ALOG10(out[sig3ind].outmbh),out[sig3ind].inbeta,psym=1,thick=2
  ;device,/close
  stop
                                ;ML-ratio
  ;device,filename='mlmbh.ps',/color
  djs_plot,out.ml,alog10(out.outmbh),psym=3,ytitle='log(MBH)',xtitle='M/L_i',chars=2,ysty=1,yrange=[6.0,7.5],charsize=1.5,charthick=4,xthick=3,ythick=3
  djs_oplot,out[sig1ind].ml,ALOG10(out[sig1ind].outmbh),psym=4;,thick=5
  djs_oplot,out[sig2ind].ml,ALOG10(out[sig2ind].outmbh),psym=5,thick=4
  djs_oplot,out[sig3ind].ml,ALOG10(out[sig3ind].outmbh),psym=1,thick=2
  ;device,/close
  ;set_plot,'x'
  stop
;  inmbh=fltarr(nmbhs,nbetas,ninclinations,nmls)
;  inbeta=FLTARR(nmbhs,nbetas,ninclinations,nmls)
;  ininclination=FLTARR(nmbhs,nbetas,ninclinations,nmls)
;  outchi2=FLTARR(nmbhs,nbetas,ninclinations,nmls)
;  outmbh=FLTARR(nmbhs,nbetas,ninclinations,nmls)
;  outml=FLTARR(nmbhs,nbetas,ninclinations,nmls)
;  stop

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
