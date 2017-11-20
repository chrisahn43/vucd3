pro myfunct, x, A, model,PDER
  model=A[0]*SIN(A[1]*x + A[2])
  IF N_PARAMS() GE 4 THEN $
    pder = [[SIN(A[1]*x + A[2])], [A[0] * X * COS(A[1]*X+A[2])], [A[0]*COS(A[1]*X+A[2])]]

end
