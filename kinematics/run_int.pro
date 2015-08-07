@int_cube.pro
PRO RUN_INT

READCOL, 'Catlist',files,FORMAT='A'
files=STRMID(files,0,22)
nfiles=N_ELEMENTS(files)
FOR i=0,nfiles-1 DO int_cube,files[i],/tsn
print,files
STOP
END
