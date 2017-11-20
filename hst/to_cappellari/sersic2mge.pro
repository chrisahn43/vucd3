; 
; 
; CALLING SEQUENCE: 
; mge=sersics2mge(sersics,Msun,A_B=A_B,fitrange=fitrange)
;
; Routine to convert multiple Sersic into an single MGE for dynamical 
; modelling of nearby galaxies.
;
; NOTES:
; Each Sersic is fitted with 1D MGE. The fitted ranges are adapted 
; automatically to capture 99.98% of the flux for each Sersic. The numbers of 
; Gaussians are chosen depeding on the Sersic n.
; The total flux of each Sersic-expansion is normalized to the input magnitude.
; There are typically ~9 Gaussians per Sersic.
;
; Todo: Write a 2 dimensional fitter. That will reduce the total number of Gaussians.
; ToDo: Replace fit_1d_mge with analytic expansion (Sersic -> MGE)
; 
;
; References for Solar luminosities http://www.ucolick.org/~cnaw/sun.html
;                              http://www.astro.umd.edu/~ssm/ASTR620/mags.html
;
;
; INPUTS: 
;        sersics: structure of type {mag:10,re:3,n:2.2,pa:90,q:0.5}
; 	     Msun: Solar luminosity
;
;         
; OPTIONAL INPUT:
;		A_B       - Forground extinction
;   fitrange  - two elements vector with min and max range that should at least be fitted.
;               (Usefull for ensuring the MGE extends to HST scales.)
;
;
; OUTPUT:  resulting MGE in solar luminosities per parsec^2
;          fltarr(NGAUSS,4): 1,2,3,4 : Lum, sigma,q,PA
;
;
; pro example
;  sersic1={mag:8,re:3,n:2.2,pa:90,q:0.5}
;  sersic2={mag:9,re:10,n:4.2,pa:90,q:0.8}
;  sersics=[sersic1,sersic2]
;  Msun=3.27 ; Solar luminsoity in Ks band
;  A_B=0.01 ; forground extinction (in units of K_S band)
;  	
;  mge=SERSICs2MGE(sersics,Msun,A_b=A_b)
; end
; 
; HISTORY:
; v0.1 Feb 2 2014, Remco van den Bosch, Heidelberg
; V0.2 May 22 2014, RvdB, HD
;      fixed 1D -> 2D flux conversion


function func_sersic_bn, bn_arr
; helper function for computeSersicBn
;; Implicit equation for bn (deV profile)
;  Copied from SYNMAG's mgsersic.pro by Kevin Bundy
  common func_bn, nser,bn,re
  n=nser
  return, gamma(2.0*n) - 2*igamma(2.0*n, bn_arr)*gamma(2.0*n)
end

function computeSersicBn, n
	; Compute Sersic b_n given Sersic n 
  common func_bn, nser,bn,re
  nser= n
  bnf = newton((2*n - 0.324)>0.1, 'func_sersic_bn')  ;; Exact solution, interpolated
  
  return,bnf
end

function func_sersic_accr,R
  ; helper function for sersic_findrange
	common func_sersic_acc, nser,bn,re,f
	x=bn*(R/Re)^(1d/nser)
	return,igamma(2.0*nser,x)- f 
end

function sersic_findrange,f_in,rein,n
	; find the radius that encloses fraction f of the mass, given Sersic parameters, rein and n
	common func_sersic_acc, nser,bn,re,f
	bn=computeSersicBn(n)
	nser=n
	re=rein
	f=f_in
  r=zbrent( re/1e5, re*1e5, FUNC = 'func_sersic_accr',TOLERANCE=f_in/10d)
	return,r
end

function fitmge2sersic,n_sersic, Re_sersic, NSTEP=NSTEP, fitbound=fitbound,NGAUSS=NGAUSS,frac=frac
	if not keyword_set(NGAUSS) then ngauss=long(poly(alog10(n_sersic),linfit(alog10([0.5,30]),[1,30]))+1)
	if not keyword_set(nstep) then nstep=300
	if not keyword_set(frac) then  frac=1d-4
	
	fitmin=sersic_findrange(frac,Re_sersic,n_sersic) 
	fitmax=sersic_findrange(1d -frac,Re_sersic,n_sersic)
	print,fitmin,fitmax
	
	if keyword_set(fitbound) then begin
		fitmin = fitmin < fitbound[0]
		fitmax = fitmax > fitbound[1]
	endif
	r=range(fitmin,fitmax,nstep,/log)
	bn=computeSersicBn(n_sersic)
	profile = exp( -bn*((R/Re_sersic)^(1.0/n_sersic)-1) )  ;; Makes I0 into Ie (surface brightness at Re)
  
	profile = (1+randomn(seed,nstep)*0.001)*profile
 ;       stop
  mge_fit_1d, R, profile, NGAUSS=NGAUSS, SOL=sol,outer_slope=n_sersic

  ;; Change first dimension output from 1D flux integral to (total) amplitudes
  lum = sol[0,*] / sqrt(2*!pi*sol[1,*]^2)
  sol[0,*]=lum
	return,sol
end


function sersic2mge,sersic,Msun,A_b=A_b,fitrange=fitrange,frac=frac
	 ;if not keyword_set(fitrange) then fitrange=[0.01,10]

	  mge=fitmge2sersic(sersic.n,sersic.re,fitbound=fitrange,frac=frac)
;          stop
		zero_pt=21.0  ; arbitrary zero_point
    
		totflux=total( mge[0,*] * (2*!pi) * mge[1,*]^2 * sersic.q)
		targetmag=sersic.mag
		if keyword_set(A_B) then targetmag -= a_b
    normfl= mag2flux( targetmag,zero_pt)
		mge[0,*] *= normfl/totflux
    mu=flux2mag(mge[0,*],zero_pt) 
		mag= flux2mag(mge[0,*] * (2*!pi) * mge[1,*]^2 * sersic.q,zero_pt)
		L_obs=(64800.0/!pi)^2*10^(0.4*(Msun-mu))
		;print,"mag,mu,L_obs,sigma"
		;forprint,mag,mu,L_obs,mge[1,*]
    
		mge_out=dblarr(N_ELEMENTS(mge[1,*]),4)
		mge_out[*,0]=L_obs	
		mge_out[*,1]=mge[1,*]	 ; sigma
		mge_out[*,2]=reform(sersic[0].q)
		mge_out[*,3]=reform(sersic[0].pa)

		return, MGE_out
end

function sersics2mge,sersics,Msun,A_B=A_B,fitrange=fitrange,frac=frac
	;if not keyword_set(fitrange) then fitrange=[0.01,10]
	
	nS=N_ELEMENTS(sersics.mag)
	
	mge=SERSIC2MGE(sersics[0],msun,A_B=A_B,fitrange=fitrange,frac=frac)
	for i=1l,ns-1 do $
		mge=[mge,SERSIC2MGE(sersics[i],msun,A_B=A_B,fitrange=fitrange)]

     
		forprint,mge[*,0],mge[*,1],mge[*,2],mge[*,3]
		
		!p.multi=[0,1,2]
		plot,mge[*,1],mge[*,0],/xlog,/ylog,psym=1
		plot,mge[*,1],mge[*,2],/xlog,psym=1
		!p.multi=[0,1,1]
	return, mge
end



