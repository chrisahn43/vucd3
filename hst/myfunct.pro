function myfunct, x,y,p
  model=p[0]+p[1]*x+p[2]*y;+p[3]*x*y+p[4]*x^2+p[5]*y^2+p[6]*x^2*y+p[7]*x*y^2+p[8]*x^3+p[9]*y^3
  return,model
END

  

  
