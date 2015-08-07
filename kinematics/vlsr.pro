PRO VLSR
;http://www.dark-cosmology.dk/~justyn/iraf/obsdb.dat
;Gemini north
longitude=155.+28.142805d0/60.
latitude=19.+49.42809d0/60.
altitude = 4213.4
ra=187.73916667/15.
dec=12.42911111
;Feb 1
jd=57054.611263
helcorr,longitude,latitude,altitude,ra,dec,jd,corr,hjd
print,'Feb 1',corr
;May 3
jd=57145.398618
helcorr,longitude,latitude,altitude,ra,dec,jd,corr,hjd
print,'May 3',corr
;within 0.03 km/sec of the rvcorrect value. Yippee!
;based on rvcorrect, it seems like this value should be ADDED to the
;radial velocities.  
;Feb 1       22.726699
;May 3      -17.869727

STOP
END
