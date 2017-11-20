PRO KBANDMGE
  scale=0.025
  zeropoint_r_ab=25.9799
  zeropoint_i_ab=25.28697
  zeropoint_r_ve=25.90069
  zeropoint_i_ve=24.86147
  colorin=0.61-(0.061-0.034)
  colorout=0.74-(0.061-0.034)
  vsun=4.80
  isun=4.53
  readcol,'~/bc03/bc03/hst_add.1ABmag',logage,u,g,r,i,z,f606w,f814w,f475w,f850lp,format='F,F,F,F,F,F,F,F,F,F'
  hstcolor=f606w-f814w
  readcol,'~/bc03/bc03/hst_add.4color',logage4,mbol,bmag,vmag,kmag,mliv,mrem,mret,mgal,sfr,mtot,mlbtot,mlvtot,format='F,F,F,F,F,F,F,F,F,F,F,F,F'
  hstcolor=f606w-f814w
  stop
  inind=where(hstcolor gt colorin-0.001 and hstcolor lt colorin+0.001)
  outind=where(hstcolor gt colorout-0.05 and hstcolor lt colorout+0.05)
  ksun=3.28
  ik=f814w-kmag
  inscale=(10^(-0.4*(ik[inind]-(isun-ksun))))
  inscale=inscale[0]
  outscale=(10^(-0.4*(ik[outind]-(isun-ksun))))
  outscale=outscale[0]
  readcol,'../evstigneeva/vucd3_mge_outputsersic.dat',intensity,sigma,q,pa,format='D,D,D,D'
  a=where(q lt 0.8)
  intensity[a]=intensity[a]*inscale
  b=where(q gt 0.8)
  intensity[b]=intensity[b]*outscale
  forprint,intensity,sigma,q,pa,textout='./vucd3_mge_outputsersic_k.dat',format='D,D,D,D'
  stop
  
END
