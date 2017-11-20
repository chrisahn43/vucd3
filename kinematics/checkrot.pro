PRO CHECKROT

  infile='./vor_out/intspec_wallace_best8.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'

;  forprint,rav,vel,velerr,disp,disperr,format='F,F,F,F',textout='checkrot_down.dat'
  xbin=rav*0.05
  set_plot,'ps'
  device,filename='checkrot_down.ps',/color
  djs_plot,xbin,disp,psym=4,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.52],yran=[25,55],/xsty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  oploterr,xbin,disp,disperr
  djs_plot,xbin,vel,psym=4,ytitle='Velocity (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.52],yran=[min(vel)-5,max(vel)+5],/xsty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2,/ysty
  oploterr,xbin,vel,velerr
  device,/close
  set_plot,'x'

  stop

  readcol,'checkrot_right.dat',rav,rvel,rvelerr,rdisp,rdisperr,format='F,F,F,F,F'
  readcol,'checkrot_left.dat',rav,lvel,lvelerr,ldisp,ldisperr,format='F,F,F,F,F'

  readcol,'checkrot_up.dat',rav,uvel,uvelerr,udisp,udisperr,format='F,F,F,F,F'

  readcol,'checkrot_down.dat',rav,dvel,dvelerr,ddisp,ddisperr,format='F,F,F,F,F'

  xbin=rav*0.05
  set_plot,'ps'
  device,filename='checkrot_all.ps',/color
  djs_plot,xbin,rdisp,psym=2,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.52],yran=[25,55],/xsty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  djs_oplot,xbin,ldisp,psym=2,symsize=2,color='blue'
  djs_oplot,xbin,udisp,psym=2,symsize=2,color='red'
  djs_oplot,xbin,ddisp,psym=2,symsize=2,color='green'
  djs_plot,xbin,rvel,psym=2,ytitle='Velocity (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.52],yran=[min(rvel)-10,max(rvel)+10],/xsty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2,/ysty
  djs_oplot,xbin,lvel,psym=2,symsize=2,color='blue'
  djs_oplot,xbin,uvel,psym=2,symsize=2,color='red'
  djs_oplot,xbin,dvel,psym=2,symsize=2,color='green'
  items=['right','left','up','down']
  sym=[2,2,2,2]
  colors=['black','blue','red','green']
  al_legend,items,psym=sym,color=colors,/window,background_color='white',charthick=4,thick=3

  djs_plot,xbin,rvel-lvel,psym=2,ytitle='Velocity (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.52],yran=[-20,25],/xsty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2,/ysty
  djs_oplot,xbin,uvel-dvel,psym=2,symsize=2,color='blue'
  items=['right - left','up - down']
  sym=[2,2]
  colors=['black','blue']
  al_legend,items,psym=sym,color=colors,/window,background_color='white',charthick=4,thick=3

  
  device,/close
  set_plot,'x'
  stop

END
