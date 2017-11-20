PRO COMPARE_JAM_OUTPUT
  infile='../../kinematics/vor_out/intspec_wallace_best8.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  xbin=rav*0.05                 ;major axis*0.05 arcsec/pixel
  xbinerror=((rav-rin)/2.)*0.05
  ybin=fltarr(n_elements(rav)) ;minor axis
  distance=16.5
;  mbhs=[0.,10^(findgen(6)*0.2+6.)]
  mbhs=[0,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6]
  ml=findgen(40)*0.05+0.05
  inclinations=[60.,70.,80.,90.]
  betas=findgen(11)*0.1-0.2
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  ninclinations=n_elements(inclinations)
  nbetas=n_elements(betas)
  fixgauss_mass_k=mrdfits('fixgausspsf_mass_k.fits',1)
  fixnorm_mass_k=mrdfits('fixnormpsf_mass_k.fits',1)
  fixold_mass_k=mrdfits('fixoldpsf_mass_k.fits',1)
  freenorm_mass_k=mrdfits('freenormpsf_mass_k.fits',1)
  freenorm_k=mrdfits('freenormpsf_k.fits',1)
  fixnorm_k=mrdfits('fixnormpsf_k.fits',1)
  fixnorm_mass_i=mrdfits('fixnormpsf_mass_i.fits',1)

  minchi=dblarr(nmbhs,nbetas)
  data=reform(fixnorm_mass_k,nmbhs,nbetas,ninclinations,nmls)
  for i=0,nmbhs-1 do begin
     for j=0,nbetas-1 do begin
        minchi[i,j]=MIN(data[i,j,*,*].chi2)
     endfor
  endfor
  
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='betambh.ps',/color,/HELVETICA,bits=4,/cmyk,/encaps,/inches,ysize=7;BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixnorm_mass_k.outmbh,fixnorm_mass_k.inbeta,YTITLE='!7b!3',XTITLE='M!IBH!N [M!D!9n!3!N]',XCHARS=1.2,/NODATA,YCHARS=1.3,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(betas)-0.1,max(betas)+0.1],CHARSIZE=1.25,XSTY = 1,XRANGE=[min(mbhs)-2d4,MAX(mbhs)+2d4],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.95,0.98] ;,/XLOG
  OPLOT,fixnorm_mass_k.outmbh,fixnorm_mass_k.inbeta,PSYM=6,THICK=10,SYMSIZE=0.2,color=100
  min=MIN(fixnorm_mass_k.chi2,minpos)
  CONTOUR,minchi,mbhs,betas,/OVERPLOT,LEVELS=[fixnorm_mass_k[minpos].chi2+2.3,fixnorm_mass_k[minpos].chi2+6.18,fixnorm_mass_k[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]
  xx=REPLICATE(fixnorm_mass_k[minpos].outmbh,N_ELEMENTS(fixnorm_mass_k.outmbh))
  yy=REPLICATE(fixnorm_mass_k[minpos].inbeta,N_ELEMENTS(fixnorm_mass_k.outmbh))
  OPLOT,xx,yy,psym=6,thick=15,symsize=0.7,color=150
  a=where(fixnorm_mass_k.inmbh eq 5.e5)
  b=where(fixnorm_mass_k[a].chi2 eq min(fixnorm_mass_k[a].chi2))
  arrow,0.,fixnorm_mass_k[a[b]].inbeta,fixnorm_mass_k[a[b]].inmbh,fixnorm_mass_k[a[b]].inbeta,/data,thick=10,hsize=0.01,color=145
  arrow,fixnorm_mass_k[a[b]].inmbh,-.3,fixnorm_mass_k[a[b]].inmbh,fixnorm_mass_k[a[b]].inbeta,/data,thick=10,hsize=0.01,color=145
  onepercent=fixnorm_mass_k[a[b]].rms
;  stop
  a=where(fixnorm_mass_k.inmbh eq 1.5e6)
  b=where(fixnorm_mass_k[a].chi2 eq min(fixnorm_mass_k[a].chi2))
  arrow,0.,fixnorm_mass_k[a[b]].inbeta,fixnorm_mass_k[a[b]].inmbh,fixnorm_mass_k[a[b]].inbeta,/data,thick=10,hsize=0.01,color=220
  arrow,fixnorm_mass_k[a[b]].inmbh,-.3,fixnorm_mass_k[a[b]].inmbh,fixnorm_mass_k[a[b]].inbeta,/data,thick=10,hsize=0.01,color=220
  fivepercent=fixnorm_mass_k[a[b]].rms
;  stop
  a=where(fixnorm_mass_k.inmbh eq 3.e6)
  b=where(fixnorm_mass_k[a].chi2 eq min(fixnorm_mass_k[a].chi2))
  arrow,0.,fixnorm_mass_k[a[b]].inbeta,fixnorm_mass_k[a[b]].inmbh,fixnorm_mass_k[a[b]].inbeta,/data,thick=10,hsize=0.01,color=200
  arrow,fixnorm_mass_k[a[b]].inmbh,-.3,fixnorm_mass_k[a[b]].inmbh,fixnorm_mass_k[a[b]].inbeta,/data,thick=10,hsize=0.01,color=200
  tenpercent=fixnorm_mass_k[a[b]].rms
;  stop
  set_plot,'x'
  set_plot,'ps' & LOADCT,40,/silent
  xyouts,[6.e6],[0.83],['VUCD3'],charthick=8,charsize=1.5,/data

  device,filename='onedanisotropy.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7
  plotsym,0,1/2.,/fill
  djs_plot,xbin,disp,psym=8,ytitle='Dispersion [km s^{-1}]',xtitle='Radius ["]',xran=[0.,0.5],yran=[25,55],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  oploterror,xbin,disp,xbinerror,disperr,/NOHAT,psym=3
  djs_oplot,xbin,onepercent,thick=4,color=145
  djs_oplot,xbin,fivepercent,thick=4,color=220
  djs_oplot,xbin,tenpercent,thick=4,color=200
  items=['1% BH Mass','5% BH Mass', '10% BH Mass']
  lines=[0,0,0]
  color=[145,220,200]
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/right
  xyouts,[0.01],[26],['VUCD3'],charthick=3,charsize=1.5,/data

  set_plot,'x'
  minchi=dblarr(nmbhs,nmls)
  for i=0,nmbhs-1 do begin
     for j=0,nmls-1 do begin
        minchi[i,j]=MIN(data[i,*,*,j].chi2)
     endfor
  endfor
    
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='mlmbh.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7;BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixnorm_mass_k.outmbh,fixnorm_mass_k.ml,YTITLE='!7C!3',XTITLE='M!IBH!N [M!D!9n!3!N]',XCHARS=1.2,/NODATA,YCHARS=1.3,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(ml),2.05],CHARSIZE=1.25,XSTY = 1,XRANGE=[min(mbhs)-2d4,MAX(mbhs)+2d4],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.951,0.98] 
  OPLOT,fixnorm_mass_k.outmbh,fixnorm_mass_k.ml,PSYM=6,THICK=10,SYMSIZE=0.2,color=100

  min=MIN(fixnorm_mass_k.chi2,minpos)
  CONTOUR,smooth(minchi,1),mbhs,ml,/OVERPLOT,LEVELS=[fixnorm_mass_k[minpos].chi2+2.3,fixnorm_mass_k[minpos].chi2+6.18,fixnorm_mass_k[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixnorm_mass_k[minpos].outmbh,N_ELEMENTS(fixnorm_mass_k.outmbh))
  yy=REPLICATE(fixnorm_mass_k[minpos].ml,N_ELEMENTS(fixnorm_mass_k.outmbh))
  OPLOT,xx,yy,psym=6,thick=15,symsize=0.7,color=150
    a=where(fixnorm_mass_k.inmbh eq 5.e5)
  b=where(fixnorm_mass_k[a].chi2 eq min(fixnorm_mass_k[a].chi2))
  arrow,0.,fixnorm_mass_k[a[b]].ml,fixnorm_mass_k[a[b]].inmbh,fixnorm_mass_k[a[b]].ml,/data,thick=10,hsize=0.01,color=145
  arrow,fixnorm_mass_k[a[b]].inmbh,0.05,fixnorm_mass_k[a[b]].inmbh,fixnorm_mass_k[a[b]].ml,/data,thick=10,hsize=0.01,color=145
  onepercent=fixnorm_mass_k[a[b]].rms
;  stop
  a=where(fixnorm_mass_k.inmbh eq 1.5e6)
  b=where(fixnorm_mass_k[a].chi2 eq min(fixnorm_mass_k[a].chi2))
  arrow,0.,fixnorm_mass_k[a[b]].ml,fixnorm_mass_k[a[b]].inmbh,fixnorm_mass_k[a[b]].ml,/data,thick=10,hsize=0.01,color=220
  arrow,fixnorm_mass_k[a[b]].inmbh,0.05,fixnorm_mass_k[a[b]].inmbh,fixnorm_mass_k[a[b]].ml,/data,thick=10,hsize=0.01,color=220
  fivepercent=fixnorm_mass_k[a[b]].rms
;  stop
  a=where(fixnorm_mass_k.inmbh eq 3.e6)
  b=where(fixnorm_mass_k[a].chi2 eq min(fixnorm_mass_k[a].chi2))
  arrow,0.,fixnorm_mass_k[a[b]].ml,fixnorm_mass_k[a[b]].inmbh,fixnorm_mass_k[a[b]].ml,/data,thick=10,hsize=0.01,color=200
  arrow,fixnorm_mass_k[a[b]].inmbh,0.05,fixnorm_mass_k[a[b]].inmbh,fixnorm_mass_k[a[b]].ml,/data,thick=10,hsize=0.01,color=200
  tenpercent=fixnorm_mass_k[a[b]].rms
;  stop
  xyouts,[6.e6],[1.95],['VUCD3'],charthick=8,charsize=1.5,/data

  set_plot,'x'
  minchi=dblarr(nbetas,nmls)
  for i=0,nbetas-1 do begin
     for j=0,nmls-1 do begin
        minchi[i,j]=MIN(data[*,i,*,j].chi2)
     endfor
  endfor
    
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='betaml.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7 ;BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR

  plot,fixnorm_mass_k.inbeta,fixnorm_mass_k.ml,YTITLE='!7C!3',XTITLE='Beta',XCHARS=1,/NODATA,YCHARS=1,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(ml)-0.1,2.2],CHARSIZE=1,XSTY = 1,XRANGE=[min(betas)-0.1,max(betas)+0.1],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] 

  OPLOT,fixnorm_mass_k.inbeta,fixnorm_mass_k.ml,PSYM=6,THICK=10,SYMSIZE=0.2,color=100

  min=MIN(fixnorm_mass_k.chi2,minpos)
  CONTOUR,minchi,betas,ml,/OVERPLOT,LEVELS=[fixnorm_mass_k[minpos].chi2+2.3,fixnorm_mass_k[minpos].chi2+6.18,fixnorm_mass_k[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixnorm_mass_k[minpos].inbeta,N_ELEMENTS(fixnorm_mass_k.inbeta))
  yy=REPLICATE(fixnorm_mass_k[minpos].ml,N_ELEMENTS(fixnorm_mass_k.inbeta))
  OPLOT,xx,yy,psym=6,thick=15,symsize=0.7,color=150


  ind=where(fixnorm_mass_k.inbeta eq 0.)
  fixnorm_mass_k=fixnorm_mass_k[ind]
  minchi=dblarr(nmbhs,ninclinations)
  for i=0,nmbhs-1 do begin
     for j=0,ninclinations-1 do begin
        mbh=mbhs[i] & inc=inclinations[j]
        ind=where(fixnorm_mass_k.ininc eq inc and fixnorm_mass_k.outmbh eq mbh)
        minchi[i,j]=MIN(fixnorm_mass_k[ind].chi2)
     endfor
  endfor
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='incmbh_nobeta.ps',BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixnorm_mass_k.outmbh,fixnorm_mass_k.ininc,YTITLE='Inclination',XTITLE='MBH',XCHARS=1,/NODATA,YCHARS=1,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(inclinations)-5,max(inclinations)+5],CHARSIZE=1,XSTY = 1,XRANGE=[min(mbhs)-2d4,MAX(mbhs)+2d4],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] 

  
  PLOTSYM,0,/FILL
  OPLOT,fixnorm_mass_k.outmbh,fixnorm_mass_k.ininc,PSYM=8,THICK=5,SYMSIZE=1,color=100

  min=MIN(fixnorm_mass_k.chi2,minpos)
  CONTOUR,minchi,mbhs,inclinations,/OVERPLOT,LEVELS=[fixnorm_mass_k[minpos].chi2+2.3,fixnorm_mass_k[minpos].chi2+6.18,fixnorm_mass_k[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixnorm_mass_k[minpos].outmbh,N_ELEMENTS(fixnorm_mass_k.outmbh))
  yy=REPLICATE(fixnorm_mass_k[minpos].ininc,N_ELEMENTS(fixnorm_mass_k.outmbh))
  OPLOT,xx,yy,psym=8,thick=10,symsize=2,color=150

  minchi=dblarr(nmls,ninclinations)
  for i=0,nmls-1 do begin
     for j=0,ninclinations-1 do begin
        mls=ml[i] & inc=inclinations[j]
        ind=where(fixnorm_mass_k.ininc eq inc and fixnorm_mass_k.ml eq mls)
        minchi[i,j]=MIN(fixnorm_mass_k[ind].chi2)
     endfor
  endfor
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='incml_nobeta.ps',BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixnorm_mass_k.ml,fixnorm_mass_k.ininc,YTITLE='Inclination',XTITLE='M/L',XCHARS=1,/NODATA,YCHARS=1,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(inclinations)-5,max(inclinations)+5],CHARSIZE=1,XSTY = 1,XRANGE=[min(ml)-0.1,1.2],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] 

  
  PLOTSYM,0,/FILL
  OPLOT,fixnorm_mass_k.ml,fixnorm_mass_k.ininc,PSYM=8,THICK=5,SYMSIZE=1,color=100

  min=MIN(fixnorm_mass_k.chi2,minpos)
  CONTOUR,minchi,ml,inclinations,/OVERPLOT,LEVELS=[fixnorm_mass_k[minpos].chi2+2.3,fixnorm_mass_k[minpos].chi2+6.18,fixnorm_mass_k[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixnorm_mass_k[minpos].ml,N_ELEMENTS(fixnorm_mass_k.ml))
  yy=REPLICATE(fixnorm_mass_k[minpos].ininc,N_ELEMENTS(fixnorm_mass_k.ml))
  OPLOT,xx,yy,psym=8,thick=10,symsize=2,color=150

  stop

  set_plot,'x'
  ind=where(fixgauss_mass_k.inbeta eq 0.)
  fixgauss_mass_k=fixgauss_mass_k[ind]
  ind=where(fixnorm_mass_k.inbeta eq 0.)
  fixnorm_mass_k=fixnorm_mass_k[ind]
  ind=where(fixold_mass_k.inbeta eq 0.)
  fixold_mass_k=fixold_mass_k[ind]
  ind=where(freenorm_mass_k.inbeta eq 0.)
  freenorm_mass_k=freenorm_mass_k[ind]
  ind=where(freenorm_k.inbeta eq 0.)
  freenorm_k=freenorm_k[ind]
  ind=where(fixnorm_k.inbeta eq 0.)
  fixnorm_k=fixnorm_k[ind]
  ind=where(fixnorm_mass_i.inbeta eq 0.)
  fixnorm_mass_i=fixnorm_mass_i[ind]
  likelihood=dblarr(nmbhs)
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(fixnorm_mass_k.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixnorm_mass_k[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixnorm_mass_k[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixnorm_mass_k[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixnorm_mass_k=likelihood/max(likelihood)
  mbhsfixnorm_mass_k=interpol(mbhs,likefixnorm_mass_k,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  sigma=[-3,-2,-1,0,1,2,3]
  ml=findgen(30)*0.05+0.05
  nmls=n_elements(ml)
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(fixold_mass_k.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixold_mass_k[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixold_mass_k[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixold_mass_k[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixold_mass_k=likelihood/max(likelihood)
  mbhsfixold_mass_k=interpol(mbhs,likefixold_mass_k,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])

  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(fixgauss_mass_k.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixgauss_mass_k[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixgauss_mass_k[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixgauss_mass_k[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixgauss_mass_k=likelihood/max(likelihood)
  mbhsfixgauss_mass_k=interpol(mbhs,likefixgauss_mass_k,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(freenorm_mass_k.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(freenorm_mass_k[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(freenorm_mass_k[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(freenorm_mass_k[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreenorm_mass_k=likelihood/max(likelihood)
  mbhsfreenorm_mass_k=interpol(mbhs,likefreenorm_mass_k,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(freenorm_k.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(freenorm_k[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(freenorm_k[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(freenorm_k[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreenorm_k=likelihood/max(likelihood)
  mbhsfreenorm_k=interpol(mbhs,likefreenorm_k,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(fixnorm_k.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixnorm_k[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixnorm_k[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixnorm_k[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixnorm_k=likelihood/max(likelihood)
  mbhsfixnorm_k=interpol(mbhs,likefixnorm_k,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(fixnorm_mass_i.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixnorm_mass_i[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixnorm_mass_i[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixnorm_mass_i[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixnorm_mass_i=likelihood/max(likelihood)
  mbhsfixnorm_mass_i=interpol(mbhs,likefixnorm_mass_i,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  
  stop
;  mbhs=[findgen(17)*1.e5+3.4e6]
;  ml=findgen(40)*0.05+0.05
;  inclinations=60.
;  beta=0.
;  nmbhs=n_elements(mbhs)
;  nmls=n_elements(ml)
;  ninclinations=n_elements(inclinations)
;  final=REPLICATE({inmbh:0.0,chi2:0.0,ml:0.0,rms:fltarr(n_elements(disp))},nmbhs,nmls)
;  readcol,'~/research/code/gemini15/vucd3/kinematics/newkinematic_psf_moffat.dat',normpsf,sigmapsf,format='F,F'
  readcol,'../make_decon/vucd3_mge_outputsersic_k.dat',intensity,sigmaarc,q,pa,format='F,F,F,F'
  readcol,'../evstigneeva/vucd3_mge_outputsersic_nmass.dat',mass,format='D'
;  surf_lum = intensity
;  sigma_lum = sigmaarc
;  qobs_lum = q
;  surf_pot = mass
;  sigma_pot = sigmaarc
;  qobs_pot = q
;  for i=0,nmbhs-1 do begin
;     for j=0,nmls-1 do begin
;        fitml=ml[j]
;        mbh=mbhs[i]/fitml
;        jam_axisymmetric_rms,surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, inclinations, mbh, distance, xbin, ybin, rmsModel, BETA=beta,RMS=disp,ERMS=disperr,SIGMAPSF=sigmapsf,NORMPSF=normpsf,PIXSIZE=0.05,ml=fitml,chi2=chi2,STEP=0.02
;        final[i,j].chi2=chi2*FLOAT(n_elements(xbin))
;        final[i,j].inmbh=mbh*fitml
;        final[i,j].ml=fitml
;        final[i,j].rms=rmsmodel
;     endfor
;  endfor
;  mwrfits,final,'finalgrid.fits'
;  stop
  final=mrdfits('finalgrid.fits',1)
  final=reform(final,17,40)
  stop
  mbhs=[0,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,3.4e6,findgen(4)*1.e5+3.6e6,findgen(4)*1.e5+4.1e6,findgen(4)*1.e5+4.6e6]
  mbhs=mbhs[sort(mbhs)]
  nmbhs=n_elements(mbhs)
  ml=findgen(40)*0.05+0.05
  nmls=n_elements(ml)
  ind=where(fixnorm_mass_k.ininc eq 60.)
  fixnorm_mass_k=fixnorm_mass_k[ind]
  fixnorm_mass_k=reform(fixnorm_mass_k,16,40)
  likelihood=dblarr(nmbhs)
  temp=0D
  for i=0,nmbhs-1 do begin
     if (mbhs[i] lt 3.4e6 OR mbhs[i] gt 5.e6) then begin
        mbhind=where(fixnorm_mass_k.inmbh gt mbhs[i]-1 and fixnorm_mass_k.inmbh lt mbhs[i]+1)
        for j=0,nmls-1 do begin
           mlind=where(fixnorm_mass_k[mbhind].ml eq ml[j])
           temp+=exp(-0.5*(fixnorm_mass_k[mbhind[mlind]].chi2))
        endfor
        likelihood[i]=temp
     endif else begin
        mbhind=where(final.inmbh gt mbhs[i]-1 and final.inmbh lt mbhs[i]+1)
        for j=0,nmls-1 do begin
           mlind=where(final[mbhind].ml eq ml[j])
           temp+=exp(-0.5*(final[mbhind[mlind]].chi2))
        endfor
        likelihood[i]=temp
     endelse
  endfor
  likefixnorm_mass_k=likelihood/max(likelihood)
  mbhsfixnorm_mass_k=interpol(mbhs,likefixnorm_mass_k,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  print,'Best fit= ',mbhsfixnorm_mass_k[3]
  print,'1-sigma errors = ',mbhsfixnorm_mass_k[4]-mbhsfixnorm_mass_k[3],mbhsfixnorm_mass_k[3]-mbhsfixnorm_mass_k[2]
  print,'3-sigma errors = ',mbhsfixnorm_mass_k[6]-mbhsfixnorm_mass_k[3],mbhsfixnorm_mass_k[3]-mbhsfixnorm_mass_k[0]
;  myplot,filename='./cummulike.ps'
  set_plot,'ps'
  device,filename='./cummulike.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7
  LOADCT,0,/SILENT
  djs_plot,mbhs,likefixnorm_mass_k,ytitle='Cumulative Likelihood',xtitle='M_{BH} [M_{\odot}]',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=12,/nodata
  arrow,mbhsfixnorm_mass_k[3],0.,mbhsfixnorm_mass_k[3],1.,/data,thick=8.,color=100,hsize=0.01
  arrow,mbhsfixnorm_mass_k[2],0.,mbhsfixnorm_mass_k[2],1.,/data,thick=4,color=100,hsize=0.01
  arrow,mbhsfixnorm_mass_k[4],0.,mbhsfixnorm_mass_k[4],1.,/data,thick=4,color=100,hsize=0.01
  arrow,mbhsfixnorm_mass_k[0],0.,mbhsfixnorm_mass_k[0],1.,/data,thick=2,color=100,hsize=0.01
  arrow,mbhsfixnorm_mass_k[6],0.,mbhsfixnorm_mass_k[6],1.,/data,thick=2,color=100,hsize=0.01
  djs_oplot,mbhs,likefixnorm_mass_k,thick=12
  mbhs=[0,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6]

  djs_oplot,mbhs,likefixold_mass_k,thick=4,color='red',linestyle=1
  djs_oplot,mbhs,likefixgauss_mass_k,color='red',thick=4
  djs_oplot,mbhs,likefreenorm_mass_k,thick=4,color='blue'
  djs_oplot,mbhs,likefixnorm_k,thick=4,color='blue',linestyle=1
  djs_oplot,mbhs,likefreenorm_k,thick=4,color='blue',linestyle=2
  djs_oplot,mbhs,likefixnorm_mass_i,thick=4,color='cyan'
  items=['Best Fit Model','PSF Variations','Mass Model Variations', 'Luminosity Model Variations']
  lines=[0,0,0,0]
  color=['black','red','blue','cyan']
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/left
  xyouts,[6.5e6],[0.05],['VUCD3'],charthick=3,charsize=1.5,/data
  device,/close
  set_plot,'x'
  mbhs=[0,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,3.4e6,findgen(4)*1.e5+3.6e6,findgen(4)*1.e5+4.1e6,findgen(4)*1.e5+4.6e6]
  mbhs=mbhs[sort(mbhs)]
  nmbhs=n_elements(mbhs)

  likelihood=dblarr(nmls)
  temp=0.D
  mlpop=2.1314557
  for i=0,nmls-1 do begin
     for j=0,nmbhs-1 do begin
        if (mbhs[j] lt 3.4e6 OR mbhs[j] gt 5.e6) then begin
           mlind=where(fixnorm_mass_k.ml eq ml[i])
           mbhind=where(fixnorm_mass_k[mlind].inmbh gt mbhs[j]-1 and fixnorm_mass_k[mlind].inmbh lt mbhs[j]+1)
           temp+=exp(-0.5*(fixnorm_mass_k[mlind[mbhind]].chi2))
        endif else begin
           mlind=where(final.ml eq ml[i])
           mbhind=where(final[mlind].inmbh gt mbhs[j]-1 and final[mlind].inmbh lt mbhs[j]+1)
           temp+=exp(-0.5*(final[mlind[mbhind]].chi2))
        endelse
     endfor
     likelihood[i]=temp
  endfor
  
  likefixnorm_mass_k=likelihood/max(likelihood)
  mlfixnorm_mass_k=(interpol(ml,likefixnorm_mass_k,[0.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865]))*mlpop
  print,'Best fit M/L = ',mlfixnorm_mass_k[3]
  print,'1-sigma error M/L = ',mlfixnorm_mass_k[4]-mlfixnorm_mass_k[3],mlfixnorm_mass_k[3]-mlfixnorm_mass_k[2]
  print,'3-sigma error M/L = ',mlfixnorm_mass_k[6]-mlfixnorm_mass_k[3],mlfixnorm_mass_k[3]-mlfixnorm_mass_k[0]
  likelihood=dblarr(nmls)
  temp=0D
  a=where(fixnorm_mass_k.inmbh eq 0.)
  for i=0,nmls-1 do begin
     ind=where((fixnorm_mass_k[a].ml gt ml[i]-0.01 and fixnorm_mass_k[a].ml lt ml[i]+0.01) OR fixnorm_mass_k[a].ml eq ml[i],c)
     if (c gt 1 or c lt 1) then stop
     temp+=exp(-0.5*(fixnorm_mass_k[a[ind]].chi2))
     likelihood[i]=temp
  endfor
  likefixnorm_nomass=likelihood/max(likelihood)
  mlfixnorm_nomass=(interpol(ml,likefixnorm_nomass,[0.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865]))*mlpop
  print,'Best fit M/L ratio= ', mlfixnorm_nomass[3]
  print,'3-sigma error M/L ratio= ', mlfixnorm_nomass[6]-mlfixnorm_nomass[3],mlfixnorm_nomass[3]-mlfixnorm_nomass[0]
  stop
  zeropoint=25.28697d
  scale=0.025d
  Msun=4.53d
  mliin=1.382d
  mliout=2.723d
  rad=100.;8.744
  fits_read,'/Users/chrisahn/research/code/gemini15/vucd3/hst/evstigneeva/deconvolved_f814.fits',oldimg,head
  readcol,'../evstigneeva/vucd3_mge_outputsersic.dat',intensity,sigmaarc,q,pa,format='F,F,F,F'
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

  aper,inmodel,xci,yci,influx,influxerr,0.,skyerr,1,rad,-1,[1,1],/silent,setskyval=0.,/flux,/exact
  aper,outmodel,xci,yci,outflux,outfluxerr,0.,skyerr,1,rad,-1,[1,1],/silent,setskyval=0.,/flux,/exact

  bhmlind=where(final.chi2 eq min(final.chi2))
  bhml=final[bhmlind].ml
  inmag=zeropoint-2.5*alog10(influx)
  inabsmag=inmag-5*(alog10(16.5e6)-1)
  inlum=10^((inabsmag-msun)/(-2.5))
  outmag=zeropoint-2.5*alog10(outflux)
  outabsmag=outmag-5*(alog10(16.5e6)-1)
  outlum=10^((outabsmag-msun)/(-2.5))
  bhmass=bhml*((mliin*inlum)+(mliout*outlum))
  temp=where(fixnorm_mass_k.inmbh eq 0.)
  nobhmlind=where(fixnorm_mass_k[temp].chi2 eq min(fixnorm_mass_k[temp].chi2))
  nobhml=fixnorm_mass_k[temp[nobhmlind]].ml
  nobhmass=nobhml*((mliin*inlum)+(mliout*outlum))
  aveml=((influx*mliin)+(outflux*mliout))/(influx+outflux)
  print,'Total mass with black hole = ', final[bhmlind].inmbh, ' = ',bhmass
  print,'Dynamical M/L with black hole = ',bhml*aveml
  print,'Total mass w/o black hole = ', nobhmass
  print,'Dynamical M/L w/o black hole = ',nobhml*aveml
  print,'Total Luminosity = ', inlum+outlum
  intdisp=39700.
  r=((16.5e6)*(3.086e16))*((rad*0.025)*(4.848e-6))
  g=6.67e-11
  dynmass=(((intdisp^2)*r)/(g))/(1.989e30)
  print,'Total mass from integrated dispersion = ', dynmass
  set_plot,'ps'
  device,filename='./onedbestfit_BEST.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7
;  myplot,filename='./onedbestfit_BEST.ps'
  plotsym,0,1/2.,/fill
  djs_plot,xbin,disp,psym=8,ytitle='Dispersion [km s^{-1}]',xtitle='Radius ["]',xran=[0.,0.5],yran=[25,55],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  oploterror,xbin,disp,xbinerror,disperr,/NOHAT,psym=3
  a=where(fixnorm_mass_k.inmbh eq 0.)
  c=where(fixnorm_mass_k[a].chi2 eq min(fixnorm_mass_k[a].chi2))
  djs_oplot,xbin,fixnorm_mass_k[a[c]].rms,color='red',thick=4
  b=where(final.chi2 eq min(final.chi2))
  djs_oplot,xbin,final[b].rms,color='blue',thick=4
  fixnorm_mass_k=mrdfits('./fixnormpsf_mass_k.fits',1)
  mbhs=[0,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6]
  ml=findgen(40)*0.05+0.05
  inclinations=[60.,70.,80.,90.]
  betas=findgen(11)*0.1-0.2
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  ninclinations=n_elements(inclinations)
  nbetas=n_elements(betas)
;  a=where(fixnorm_mass_k.ininc eq 60.)
;  fixnorm_mass_k=fixnorm_mass_k[a]
  fixnorm_mass_k=reform(fixnorm_mass_k,nmbhs,nbetas,ninclinations,nmls)
  a=where(fixnorm_mass_k.outmbh eq 0.)
  b=where(fixnorm_mass_k[a].chi2 eq min(fixnorm_mass_k[a].chi2))
  help,fixnorm_mass_k[a[b]]
  djs_oplot,xbin,fixnorm_mass_k[a[b]].rms,color='grey',thick=4
  xyouts,[0.01],[26],['VUCD3'],charthick=3,charsize=1.5,/data
  device,/close
  set_plot,'x'
  stop


END
