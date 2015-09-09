@read_wallace.pro
PRO RUN_WALLACE
;create all fits files for the wallace spectra.
READCOL,'wallacenames.txt',names,FORMAT='A'
nspec=N_ELEMENTS(names)

FOR i=0,nspec-1 DO read_wallace,names[i]
STOP
END
