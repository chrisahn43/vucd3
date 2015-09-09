pro plotiso
  readcol,'output73983041177.dat',Z,logage,m_init,m_act,L,teff,logg,mbol,umag,bmag,vmag,rmag,imag,jmag,hmag,kmag,imfstage,blank,format='F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,I'
  age=findgen(23)*0.05+9.
  rimag=rmag - imag
;  color=['black','blue','green','yellow','purple','teal','pink','red','orange','brown','grey',
  for i=0,n_elements(age)-1 do begin
     a=where(logage eq age[i])
     djs_plot,rimag[a],imag[a],psym=2,yran=[6,-3],xran=[0.3,1.3]
     
     print,age[i]
     if age[i] eq 10.1 then begin
        ;1st RGB
        b=where((rimag[a] gt 0.45) and (rimag[a] lt 0.55))
        c=where(imag[a[b]] lt 2. and imag[a[b]] gt 1.)
        djs_oplot,rimag[a[b[c]]],imag[a[b[c]]],psym=2,color='red'
        print,teff[a[b[c]]]
        print,logg[a[b[c]]]
        stop
        ;2nd RGB
        b=where(rimag[a] gt 0.8 and rimag[a] lt 0.9)
        c=where(imag[a[b]] lt -1.5 and imag[a[b]] gt -1.8)
        djs_oplot,rimag[a[b[c]]],imag[a[b[c]]],psym=2,color='red'
        print,teff[a[b[c]]]
        print,logg[a[b[c]]]
        stop
        ;main sequence
        b=where(rimag[a] gt 0.3 and rimag[a] lt 0.37)
        c=where(imag[a[b]] lt 4.)
        djs_oplot,rimag[a[b[c]]],imag[a[b[c]]],psym=2,color='red'
        print,teff[a[b[c]]]
        print,logg[a[b[c]]]
        stop
        ;AGB
        b=where(rimag[a] gt 0.55 and rimag[a] lt 0.6)
        c=where(imag[a[b]] gt -1. and imag[a[b]] lt -0.5)
        djs_oplot,rimag[a[b[c]]],imag[a[b[c]]],psym=2,color='red'
        print,teff[a[b[c]]]
        print,logg[a[b[c]]]
        stop
     endif
 
  endfor
  
  stop
END
