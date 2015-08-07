PRO MAKE_NIFS15A,REMCO=REMCO, JEANS=JEANS

COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr  

velerrlim=100.
disperrlim=100.
velsnlim=1.
dispsnlim=1.

modelfile='remco_models/outputs/bh6.80q0.73a0.00b0.00r0.55dc1.00f-10.00/ml4.90/nn_p_kinem.out'
readcol,SKIPLIN=2, modelfile, id,fluxm,flux, velm ,vel,dvel,sigm,sig,dsig,h3m,h3,dh3,h4m, h4,dh4
id=id-1
id[0:32]=id[0:32]-1
vel=vel-1.

galaxies=['m60-ucd1']
vlsrcor=[-22.006]
vsys=[1300.]
velulim=[40.]
velllim=[-40.]
dispulim=[110.]
displlim=[40.]
distance=[16.5]
radii=[20.]
dms=DM(distance)
ngalaxies=N_ELEMENTS(galaxies)

FOR i=0,0 DO BEGIN
    print,i, galaxies[i]

;    datafile='vor_out/newbin_best9_sn25.dat'
    kimfile='lum_model/best9_subarr.fits'    
    binimfile='vor_out/newbin_best9_sn25_bin.fits'    
    ;read in data file
;    READCOL,datafile,binnum,x,y,npix,measured_sn,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err
;    vel=vel+vlsrcor[i]-vsys[i] ;correct velocities to systemic
    ;read in bin file
    binim=READFITS(binimfile)
    imsize=SIZE(binim,/dim)
    velim=FLTARR(imsize)
    dispim=FLTARR(imsize)
    h3im=FLTARR(imsize)
    h4im=FLTARR(imsize)
;    snim=FLTARR(imsize)
    mvelim=FLTARR(imsize)
    mdispim=FLTARR(imsize)
    mh3im=FLTARR(imsize)
    mh4im=FLTARR(imsize)
    nbins=N_ELEMENTS(id)
    FOR j=0,nbins-1 DO BEGIN
        ind=WHERE(binim EQ id[j])
        velim[ind]=vel[j]
        dispim[ind]=sig[j]
        h3im[ind]=h3[j]
        h4im[ind]=h4[j]
        mvelim[ind]=velm[j]
        mdispim[ind]=sigm[j]
        mh3im[ind]=h3m[j]
        mh4im[ind]=h4m[j]        
     ENDFOR   
;    STOP
    kim=READFITS(kimfile,khead)
    
    maxval=MAX(kim,maxpos)
    xinit=maxpos MOD imsize[0]
    yinit=maxpos/imsize[0]
    GCNTRD,kim,xinit,yinit,xcen,ycen,4.5
    print,xcen,ycen
    makex,kim,xarr,yarr
    rarr=SQRT((xarr-xcen)^2 + (yarr-ycen)^2)
    bad=WHERE(binim EQ 0 AND rarr GT 3.,nbad)
    IF (nbad GT 0) THEN BEGIN
        velim[bad]=-999.
        dispim[bad]=-999.
        h3im[bad]=-999.
        h4im[bad]=999.
    ENDIF
    rad=radii[i]
    minx=FIX(xcen-rad) & maxx=FIX(xcen+rad+0.5)
    miny=FIX(ycen-rad) & maxy=FIX(ycen+rad+0.5)
    subkim=kim[minx:maxx,miny:maxy]
    subvelim=velim[minx:maxx,miny:maxy]
;    snimage=snim[minx:maxx,miny:maxy]
    velimage=subvelim           ;-subvelim[rad+1,rad]
    h3image=h3im[minx:maxx,miny:maxy]
    h4image=h4im[minx:maxx,miny:maxy]
    dispimage=dispim[minx:maxx,miny:maxy]    
    mvelimage=mvelim[minx:maxx,miny:maxy]
    mdispimage=mdispim[minx:maxx,miny:maxy]
    mh3image=mh3im[minx:maxx,miny:maxy]
    mh4image=mh4im[minx:maxx,miny:maxy]
    cont=subkim



    psfile='m60-ucd1_nifs15a.ps'
;    myplot,file=psfile,xsize=12.4,ysize=9,yoff=1,/inches,bits=8,xoffset=0.25,yoffset=1.5,/cmyk

    !P.FONT=0
    !P.THICK=2
    !P.CHARTHICK=2
    !X.THICK=2
    !Y.THICK=2
;    set_plot,'ps'
;    device,filename=psfile,/color,xsize=6.5,ysize=2.5,bits=8,/inches
    myplot, filename=psfile,/color,xsize=6.5,ysize=2.5,bits=8,/inches

    !P.MULTI=[0,2,1]
    loadct,0,/silent
    size=SIZE(subkim,/dim)
    xcoord=(FINDGEN(size[0])-rad)*0.05
    ycoord=(FINDGEN(size[1])-rad)*0.05

    
    klevels=[0.016,0.04,0.1,0.25,0.625]*MAX(subkim)
;    klevels=[0.15,0.2,0.3,0.4,0.6,0.8]*MAX(subkim)
;    IF (i EQ 4) THEN klevels=[0.101,0.127,0.167,0.25,0.506,0.80]*MAX(subkim)


;first do velocities   
    minlev=velllim[i] & maxlev=velulim[i]
    nlev=200.
    dellev=(maxlev-minlev)/nlev
    vellevels=findgen(nlev)*dellev+minlev
    velcontimage=velimage
    low=WHERE(velimage LE minlev,nlow)
    IF (nlow GT 0) THEN velcontimage[low]=!VALUES.F_NAN
    high=WHERE(velimage GT maxlev,nhigh)
    IF (nhigh GT 0) THEN velcontimage[high]=!VALUES.F_NAN
    CONTOUR, cont,xcoord,ycoord, levels=klevels,/fill,/xs,/ys,/nodata,charsize=1.0,xtitle='RA Offset ["]',ytitle='Dec Offset ["]',xminor=2,/iso,xmargin=[9,4],ymargin=[3.5,3.5]
    pixelscale=SIN(1.0D0/3600./180.*!DPI)*10.^((dms[i]+5.)/5.)
    print,galaxies[i],pixelscale
    xmin=MIN(xcoord)*pixelscale
    xmax=MAX(xcoord)*pixelscale

;make the 255th color gray
    LOADCT, 33,/silent ;used to use 13
    R_CURR[255]=150 & G_CURR[255]=150 & B_CURR[255]=150
    R_CURR[0]=0 & G_CURR[0]=0 & B_CURR[0]=0
    tvlct,r_curr,g_curr,b_curr
    velcolimage=(velcontimage-minlev)/(maxlev-minlev)*255

    replace=WHERE(FINITE(velcolimage) EQ 0,nreplace)
    IF (nreplace GT 0) THEN velcolimage[replace]=255
    useimage=REBIN(velcolimage,size[0]*5,size[1]*5,/sample)
    imgunder,useimage

    axis,yaxis=0,ychars=0.00000001,/ysty
    axis,xaxis=0,xchars=0.00000001,/xsty
    axis,/xaxis,xrange=[xmax,xmin],/xsty,xtitle='RA Offset [pc]',charsize=1.0

    contour,cont,xcoord,ycoord,/overplot,levels=klevels,thick=4
    COLORBAR,range=[min(vellevels),max(vellevels)+dellev],format='(I4)',/vertical,right=1,divisions=5,charsize=1.0,title='!C Relative Velocity [km/s]',POSITION=[!X.WINDOW[1],!Y.WINDOW[0],!X.WINDOW[1]+0.03,!Y.WINDOW[1]]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;now do dispersions
    minlev=displlim[i] & maxlev=dispulim[i];displim[i]
    IF (galaxies[i] EQ 'ngc205') THEN minlev=0.0
    nlev=100.
    dellev=(maxlev-minlev)/nlev
    vellevels=findgen(nlev)*dellev+minlev
    dispcontimage=dispimage
    ind=WHERE(dispcontimage GT maxlev,nhigh)
    IF (nhigh GT 0) THEN dispcontimage[ind]=maxlev-0.2
    ind=WHERE(dispcontimage LT minlev,nlow)
    IF (nlow GT 0) THEN dispcontimage[ind]=!VALUES.F_NAN

    CONTOUR, cont,xcoord,ycoord, levels=levels,/fill,/xs,/ys,/nodata,charsize=1.0,xtitle='RA Offset ["]',ytitle='Dec Offset ["]',xminor=2,/iso,xmargin=[9,4],ymargin=[3.5,3.5]
;make the 255th color gray
    LOADCT,33,/silent
    R_CURR[255]=150 & G_CURR[255]=150 & B_CURR[255]=150
    R_CURR[0]=0 & G_CURR[0]=0 & B_CURR[0]=0
    tvlct,r_curr,g_curr,b_curr
    dispcolimage=(dispcontimage-minlev)/(maxlev-minlev)*255
    clipped=WHERE(FINITE(dispcolimage) EQ 0,nclipped)
    IF (nclipped GT 0) THEN dispcolimage[clipped]=255
;    useimage=dispcolimage
    useimage=REBIN(dispcolimage,size[0]*5,size[1]*5,/sample)
    imgunder,useimage
    axis,yaxis=0,ychars=0.00000001,/ysty
    axis,xaxis=0,xchars=0.00000001,/xsty
    axis,/xaxis,xrange=[xmax,xmin],/xsty,xtitle='RA Offset [pc]',charsize=1.0

;levels at 1-mag intervals
    contour,cont,xcoord,ycoord,/overplot,levels=klevels,thick=4
    COLORBAR,range=[min(vellevels),max(vellevels)+dellev],format='(I4)',/vertical,right=1,divisions=5,charsize=1.0,title='Dispersion [km/s]',POSITION=[!X.WINDOW[1],!Y.WINDOW[0],!X.WINDOW[1]+0.03,!Y.WINDOW[1]]



    loadct,0,/silent
    device,/close
    set_plot,'x'
ENDFOR
STOP
END
