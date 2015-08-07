PRO PLOT_INTSPEC

infile='intspec_wallace_best8.dat'
READCOLMORE,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'


;READCOL,'m59co_noimbh_seefw0.2',rad0,disp0,FORMAT='F,X,F'
;READCOL,'m59co_10imbh_seefw0.2',rad10,disp10,FORMAT='F,X,F'
;READCOL,'m59co_30imbh_seefw0.2',rad30,disp30,FORMAT='F,X,F'
;disp0=disp0*29./48.
;disp10=disp10*29./48.
;disp30=disp30*29./48.


;myplot,file='disp_profile_wvlt.ps',ysize=6,/inches
myplot,file='disp_profile_best8.ps',ysize=6,/inches
!P.MULTI=[0,1,2]
loadct,13
plot,rav*0.05,disp,psym=4,xtitle='Radius ["]',ytitle='Dispersion [km/s]',chars=2,yrange=[21,56],ymargin=[-3,3.5],xchars=0.0001,thick=3,/ysty,xrange=[0,0.6],xmargin=[7,1.5]
axis,/xaxis,xrange=!X.CRANGE*80.,/xsty,xtitle='Radius [pc]',chars=2
;legend,['0% BH','10% BH','30% BH'],thick=[4,4,4],color=[30,60,100],linestyle=[0,2,4],/bottom,chars=2
;oplot,rad0,disp0,thick=4,color=30
;oplot,rad10,disp10,thick=4,color=60,linestyle=2
;oplot,rad30,disp30,thick=4,color=100,linestyle=4
oploterror,rav*0.05,disp,disperr,thick=3,psym=4
;oplot,vrad*0.05,vdisp,psym=5,color=250,thick=3
;oploterror,vrad*0.05,vdisp,vdisperr,psym=5,color=250,thick=3
;plots,!X.CRANGE,[48,48],linestyle=2,thick=4
;legend,['NIFS','SINFONI'],psym=[4,5],color=[0,250],thick=[3,3],/bottom,chars=2

plot,rav*0.05,sn,xtitle='Radius ["]',ytitle='S/N',chars=2,ymargin=[3.5,3],psym=4,yrange=[0,39],/ysty,thick=5,xrange=[0,0.6],xmargin=[7,1.5]
;oplot,vrad*0.05,vsn,psym=5,thick=3,color=250

device,/close
set_plot,'x'

STOP
END
