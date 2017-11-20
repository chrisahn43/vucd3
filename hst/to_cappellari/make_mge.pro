@sersic2mge.pro
PRO make_mge
                                ;THIS CODE CREATES THE MGE'S
                                ;FOR USE IN THE JAM CODE. THE
                                ;ONLY DIFFERENCE BETWEEN THESE TWO MGE
                                ;OUTPUTS IS THE AXIS RATIO (q) of the
                                ;inner component Sersic. It appears
                                ;that any ratio less than 0.65
                                ;doesn't work in the JAM code. 

  scale=0.025                   ;pixel scale HST ACS/HRC
  Msun=4.53d      ;ACS F814W  AB mag from http://www.ucolick.org/~cnaw/sun.html
  A_B=0.034d      ;I band extinction from NED

  re1=2.74d*scale               ;effective radii
  re2=24.71d*scale
  sersic_orig1={mag:18.58d,re:re1,n:3.25d,pa:17.97d,q:0.62d}
  sersic_orig2={mag:17.99d,re:re2,n:1.74d,pa:20.65d,q:0.89d}

  sersics_orig=[sersic_orig1,sersic_orig2]
  mge=sersics2mge(sersics_orig,Msun,A_B=A_B)
  forprint,mge[*,0],mge[*,1],mge[*,2],mge[*,3],textout='vucd3_mgesersic_output_orig.dat',format='F,F,F,F'


  sersic_q1={mag:18.58d,re:re1,n:3.25d,pa:17.97d,q:0.65d}
  sersic_q2={mag:17.99d,re:re2,n:1.74d,pa:20.65d,q:0.89d}
  sersics_q=[sersic_q1,sersic_q2]
  mge=sersics2mge(sersics_q,Msun,A_B=A_B)
  forprint,mge[*,0],mge[*,1],mge[*,2],mge[*,3],textout='vucd3_mgesersic_output_qchange.dat',format='F,F,F,F'
  
  
  
  
END
