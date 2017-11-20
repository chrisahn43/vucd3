PRO jam_models
  ;SET UP KINEMATICS TO FIT
  disp=[51.1,49.8,47.3,39.6,32.0]
  disperr=[3.3,3.,2.7,2.1,3.1]
  ;SET UP KINEMATIC PSF 
  normpsf=[0.552,0.175,0.243,0.028];,0.002]
  sigmapsf=[0.047,0.056,0.175,0.422];,1.003]
  ;DEFINE PARAMETERS FOR JAM CODE
  xbin=[0.0326,0.069,0.118,0.217,0.437]
  ybin=fltarr(n_elements(xbin)) 
  inc=90.                  ;Since we're testing the q, want maximum inclination
  distance = 16.5          ; Assume Virgo distance in Mpc (Mei et al. 2007)
  fitml=2.
  mbh=1.e6

  ;READ IN THE MGEs
  readcol,'vucd3_mgesersic_output_orig.dat',intensity,sigmaarc,q,format='F,F,F'
  surf_lum = intensity
  sigma_lum = sigmaarc
  qobs_lum = q
  surf_pot = intensity   
  sigma_pot = sigmaarc
  qobs_pot = q

  ;THIS FIT FAILS
  jam_axisymmetric_rms,surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, inc, mbh, distance, xbin, ybin, rmsModel, BETA=0.,RMS=disp,ERMS=disperr,SIGMAPSF=sigmapsf,NORMPSF=normpsf,PIXSIZE=0.05,ml=fitml,chi2=chi2;,step=0.005
  
  ;CHANGE THE VALUE OF Q FOR INNER SERSIC COMPONENT
  ind=where(q eq 0.62)
  qobs_lum[ind]=0.65
  qobs_pot[ind]=0.65

  ;THIS FIT WORKS
  jam_axisymmetric_rms,surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, inc, mbh, distance, xbin, ybin, qrmsModel, BETA=0.,RMS=disp,ERMS=disperr,SIGMAPSF=sigmapsf,NORMPSF=normpsf,PIXSIZE=0.05,ml=fitml,chi2=chi2
  

  ;PLOT SHOWING THE FITS
  djs_plot,xbin,disp,psym=4,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.52],yran=[15,55],/xsty
  oploterr,xbin,disp,disperr
  djs_oplot,xbin,rmsModel,color='blue'
  djs_oplot,xbin,qrmsModel,linestyle=2,color='blue'
  stop
END
