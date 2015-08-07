PRO CALC_OFFSETS
;to be run after analyze_cube and going through with the IRAF data.

READCOl,'offsets.dat',xidl,yidl,xiraf,yiraf,FORMAT='X,F,F,F,F'
x=(xidl+xiraf)/2.
y=(yidl+yiraf)/2.
nfiles=N_ELEMENTS(xidl)
FOR i=0,nfiles-1 DO print,x[i],y[i],MAX(x)-x[i],MAX(y)-y[i],FORMAT='(4F7.2)'
STOP
END
