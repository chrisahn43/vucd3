FUNCTION kband_indices,lambda,spec,VELOCITY=velocity,PLOT=plot,MARMOL=marmol,bb=bb

IF NOT KEYWORD_SET(velocity) THEN velocity=244.
redshift=1.+velocity/2.9979d5

;two sided indices first
indices=['MgI','BrG','NaI','CaI','MgIa','CO20','CaVIII']
;done for NGC404
linecenters=[21065.,21660.9,22078,22643.5,22813.5,22957.,23220.4]
linewidths=[30.,20.,50.,55.,20.,52.,14.]
;linecenters=[21065.,21657.5,22073,22638.5,22813.5,22957.]
;linewidths=[30.,20.,40.,51.,20.,52.]
;linecenters=[21075.,21662.5,22077,22636.5,22957.]
;linewidths=[70.,47.,48.,51.,52.]
;linecenters=[21075.,21662.5,22077,22636.5,22813.5,22945.]
;linewidths=[70.,47.,48.,51.,20.,25.]
cont1centers=[21020.,21610.,22010.,22555.,22750.,22885.,23200.]
cont1widths=[40.,50.,60.,60.,60.,70.,25.]
;cont1centers=[21020.,21610.,22010.,22555.,22750.,22885.]
;cont1widths=[40.,50.,60.,60.,60.,70.]
cont2centers=[21130,21735.,22170.,22750.,22885.,22885.,23240.]
cont2widths=[40.,40.,60.,60.,70.,70.,25.]

IF (KEYWORD_SET(MARMOL)) THEN BEGIN
    linecenters=[21065.,21657.5,22073,22638.5,22813.5,22945.]
    linewidths=[30.,20.,40.,51.,20.,130.]
    cont1centers=[21020.,21610.,22010.,22555.,22750.,22505.]
    cont1widths=[40.,50.,60.,60.,60.,90.]
    cont2centers=[21130,21735.,22170.,22750.,22885.,22740.]
    cont2widths=[40.,40.,60.,60.,70.,60.]
ENDIF

linecenters=linecenters*redshift
cont1centers=cont1centers;*redshift
cont2centers=cont2centers;*redshift

values=FLTARR(N_ELEMENTS(indices))
IF KEYWORD_SET(PLOT) THEN !P.MULTI=[0,2,4]
FOR i=0,4 DO BEGIN
    IF (linecenters[i]-linewidths[i]/2. GT MIN(lambda) AND cont1centers[i]-cont1widths[i]/2. GT MIN(lambda)) THEN BEGIN    
        incont1=WHERE(lambda GE cont1centers[i]-cont1widths[i]/2. AND lambda LE cont1centers[i]+cont1widths[i]/2.,nincont1)
        incont2=WHERE(lambda GE cont2centers[i]-cont2widths[i]/2. AND lambda LE cont2centers[i]+cont2widths[i]/2.,nincont2)
        inline=WHERE(lambda GE linecenters[i]-linewidths[i]/2. AND lambda LE linecenters[i]+linewidths[i]/2.,ninline)
        MEANCLIP,spec[incont1],cont1,3.
        MEANCLIP,spec[incont2],cont2,3.
;        cont2=MEAN(spec[incont2])
        cont=INTERPOL([cont1,cont2],[cont1centers[i],cont2centers[i]],lambda[inline])
        meaninline=MEAN(spec[inline])
        values[i]=-2.5*ALOG10(MEAN(spec[inline]/cont))
;        values[i]=1.-MEAN(spec[inline]/cont)
    IF KEYWORD_SET(PLOT) THEN BEGIN        
        IF (i EQ 1) THEN plot,lambda,spec/MEDIAN(spec),xrange=[linecenters[i]-150.,linecenters[i]+150],ysty=1,charsize=2,yrange=[0.85,1.3] $
        ELSE plot,lambda,spec/MEDIAN(spec),xrange=[linecenters[i]-150.,linecenters[i]+150],ysty=1,charsize=2,yrange=[0.7,1.15]
        
        plots,!X.CRANGE,[cont/MEDIAN(spec),cont/MEDIAN(spec)],linestyle=2

        plots,[cont1centers[i]-0.5*cont1widths[i],cont1centers[i]+0.5*cont1widths[i]],[cont1/MEDIAN(spec),cont1/MEDIAN(spec)],thick=2,color=100
        plots,[cont2centers[i]-0.5*cont1widths[i],cont2centers[i]+0.5*cont1widths[i]],[cont2/MEDIAN(spec),cont2/MEDIAN(spec)],thick=2,color=100
        plots,[linecenters[i]-0.5*linewidths[i],linecenters[i]+0.5*linewidths[i]],[meaninline/MEDIAN(spec),meaninline/MEDIAN(spec)],thick=2,color=100
        xyouts,linecenters[i],[1.08],indices[i]
    ENDIF
    ENDIF ELSE values[i]=-99.
ENDFOR


;now the CO indices (one-sided)
FOR i=5,5 DO BEGIN
    IF KEYWORD_SET(MARMOL) THEN BEGIN
        IF NOT KEYWORD_SET(bb) THEN $
          bb=READPLAINFITS('../flux_cal/bb9500.fits',/noerr)
        bbint=INTERPOL(bb[1,*],bb[0,*],lambda)
        newspec=spec*bbint
        newspec=newspec/MEDIAN(newspec)
        incont1=WHERE(lambda GE cont1centers[i]-cont1widths[i]/2. AND lambda LE cont1centers[i]+cont1widths[i]/2.,nincont1)
        incont2=WHERE(lambda GE cont2centers[i]-cont2widths[i]/2. AND lambda LE cont2centers[i]+cont2widths[i]/2.,nincont2)
        inline=WHERE(lambda GE linecenters[i]-linewidths[i]/2. AND lambda LE linecenters[i]+linewidths[i]/2.,ninline)
        MEANCLIP,newspec[incont1],cont1,3.
        MEANCLIP,newspec[incont2],cont2,3.
        values[i]=(cont1+cont2)/(2.*MEAN(newspec[inline]))
        IF KEYWORD_SET(PLOT) THEN BEGIN        
            meaninline=MEAN(newspec[inline])            
            plot,lambda,newspec,xrange=[linecenters[i]-500.,linecenters[i]+100],ysty=16,charsize=2
            plots,[cont1centers[i]-0.5*cont1widths[i],cont1centers[i]+0.5*cont1widths[i]],[cont1,cont1],thick=2,color=100
            plots,[cont2centers[i]-0.5*cont1widths[i],cont2centers[i]+0.5*cont1widths[i]],[cont2,cont2],thick=2,color=100
            plots,[linecenters[i]-0.5*linewidths[i],linecenters[i]+0.5*linewidths[i]],[meaninline,meaninline],thick=2,color=100
            xyouts,linecenters[i],[1.08],indices[i]
        ENDIF
    ENDIF ELSE BEGIN
        incont1=WHERE(lambda GE cont1centers[i]-cont1widths[i]/2. AND lambda LE cont1centers[i]+cont1widths[i]/2.,nincont1)
        inline=WHERE(lambda GE linecenters[i]-linewidths[i]/2. AND lambda LE linecenters[i]+linewidths[i]/2.,ninline)
        MEANCLIP,spec[incont1],cont1,3.
;    cont1=MEAN(spec[incont1])
        meaninline=MEAN(spec[inline])
        values[i]=-2.5*ALOG10(MEAN(spec[inline]/cont1))
;    values[i]=1.-(MEAN(spec[inline]/cont1))
        IF KEYWORD_SET(PLOT) THEN BEGIN        
            plot,lambda,spec/MEDIAN(spec),xrange=[linecenters[i]-100.,linecenters[i]+100],ysty=16,charsize=2
            plots,!X.CRANGE,[cont1/MEDIAN(spec),cont1/MEDIAN(spec)],linestyle=2
            plots,[cont1centers[i]-0.5*cont1widths[i],cont1centers[i]+0.5*cont1widths[i]],[cont1/MEDIAN(spec),cont1/MEDIAN(spec)],thick=2,color=100
            plots,[linecenters[i]-0.5*linewidths[i],linecenters[i]+0.5*linewidths[i]],[meaninline/MEDIAN(spec),meaninline/MEDIAN(spec)],thick=2,color=100
            xyouts,linecenters[i],[1.08],indices[i]
;            print, cont1centers[i]+0.5*cont1widths[i]
        ENDIF

    ENDELSE
ENDFOR


FOR i=6,6 DO BEGIN
    IF (linecenters[i]-linewidths[i]/2. GT MIN(lambda) AND cont1centers[i]-cont1widths[i]/2. GT MIN(lambda)) THEN BEGIN    
        incont1=WHERE(lambda GE cont1centers[i]-cont1widths[i]/2. AND lambda LE cont1centers[i]+cont1widths[i]/2.,nincont1)
        incont2=WHERE(lambda GE cont2centers[i]-cont2widths[i]/2. AND lambda LE cont2centers[i]+cont2widths[i]/2.,nincont2)
        inline=WHERE(lambda GE linecenters[i]-linewidths[i]/2. AND lambda LE linecenters[i]+linewidths[i]/2.,ninline)
        MEANCLIP,spec[incont1],cont1,3.
        MEANCLIP,spec[incont2],cont2,3.
;        cont2=MEAN(spec[incont2])
        cont=INTERPOL([cont1,cont2],[cont1centers[i],cont2centers[i]],lambda[inline])
        meaninline=MEAN(spec[inline])
        values[i]=-2.5*ALOG10(MEAN(spec[inline]/cont))
;        values[i]=1.-MEAN(spec[inline]/cont)
    IF KEYWORD_SET(PLOT) THEN BEGIN        
        IF (i EQ 1) THEN plot,lambda,spec/MEDIAN(spec),xrange=[linecenters[i]-150.,linecenters[i]+150],ysty=1,charsize=2,yrange=[0.85,1.3] $
        ELSE plot,lambda,spec/MEDIAN(spec),xrange=[linecenters[i]-150.,linecenters[i]+150],ysty=1,charsize=2,yrange=[0.7,1.15]
        
        plots,!X.CRANGE,[cont/MEDIAN(spec),cont/MEDIAN(spec)],linestyle=2

        plots,[cont1centers[i]-0.5*cont1widths[i],cont1centers[i]+0.5*cont1widths[i]],[cont1/MEDIAN(spec),cont1/MEDIAN(spec)],thick=2,color=100
        plots,[cont2centers[i]-0.5*cont1widths[i],cont2centers[i]+0.5*cont1widths[i]],[cont2/MEDIAN(spec),cont2/MEDIAN(spec)],thick=2,color=100
        plots,[linecenters[i]-0.5*linewidths[i],linecenters[i]+0.5*linewidths[i]],[meaninline/MEDIAN(spec),meaninline/MEDIAN(spec)],thick=2,color=100
        xyouts,linecenters[i],[1.08],indices[i]
    ENDIF

    ENDIF ELSE values[i]=-99.
ENDFOR




IF KEYWORD_SET(PLOT) THEN !P.MULTI=[0,1,1]
return,values
END
