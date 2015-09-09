PRO VLSR
;http://www.dark-cosmology.dk/~justyn/iraf/obsdb.dat
;KPNO
longitude=111.+36.d0/60.
latitude=31.+57.8d0/60.
altitude = 2120.
;gzcep
ra=332.7136529/15.
dec=58.2012608
;http://www.csgnetwork.com/julianmodifdateconv.html
jd=44656.5 ;gzcep Feb 21 1981
helcorr,longitude,latitude,altitude,ra,dec,jd,corr,hjd;,/debug
print,'gzcep ',corr


;gasge
ra=295.0241325/15.
dec=18.0138906
;http://www.csgnetwork.com/julianmodifdateconv.html
jd=49914.5 ;gzcep Jul 16 1995
helcorr,longitude,latitude,altitude,ra,dec,jd,corr,hjd;,/debug
print,'gasge ',corr

STOP
END
