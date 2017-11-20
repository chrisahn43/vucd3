PRO AVEML

  zeropoint=25.28697d
  scale=0.025d
  Msun=4.53d
  mliin=1.38172d;AB;1.69459 VEGA OLD 2.053 VEGA NEW;
  mliout=2.72326d;AB;3.16480 VEGA OLD 4.045 VEGA NEW;
  rad=8.744
  readcol,'./evstigneeva/vucd3_mge_outputsersic.dat',intensity,sigmaarc,q,pa,format='F,F,F,F'
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/deconvolved_f814.fits',oldimg,head
  hrotate,oldimg,head,img,newhead,1
  find_galaxy,img,majoraxis,eps,ang,xci,yci

  inind=where(q lt 0.8)
  inintensity=intensity[inind]
  insigma=sigmaarc[inind]
  inq=q[inind]
  inpa=pa[inind]
  mge2image,img,xci,yci,inintensity,insigma,inq,inpa,inmodel,zeropoint=zeropoint,scale=scale,msun=msun

  outind=where(q gt 0.8)
  outintensity=intensity[outind]
  outsigma=sigmaarc[outind]
  outq=q[outind]
  outpa=pa[outind]
  mge2image,img,xci,yci,outintensity,outsigma,outq,outpa,outmodel,zeropoint=zeropoint,scale=scale,msun=msun

  aper,inmodel,xci,yci,influx,influxerr,0.,skyerr,1,100.,-1,[1,1],/silent,setskyval=0.,/flux,/exact
  aper,outmodel,xci,yci,outflux,outfluxerr,0.,skyerr,1,100.,-1,[1,1],/silent,setskyval=0.,/flux,/exact
  inmag=zeropoint-2.5*alog10(influx)
  inabsmag=inmag-5*(alog10(16.5e6)-1)
  inlum=10^((inabsmag-msun)/(-2.5))
  outmag=zeropoint-2.5*alog10(outflux)
  outabsmag=outmag-5*(alog10(16.5e6)-1)
  outlum=10^((outabsmag-msun)/(-2.5))
  aveml=((influx*mliin)+(outflux*mliout))/(influx+outflux)
  print,aveml
  stop
END
