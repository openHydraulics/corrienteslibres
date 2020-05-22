%Cálculo del calado en régimen uniforme con la Ec. de Manning
function y=yManning(Q,n,b,z,I)
  yant=0.1;
  y=0;
  while abs(y-yant)>1e-4
    yant=y;
    y=(Q*n/I^(1/2)*(b+2*yant*(1+z^2)^(1/2))^(2/3)/(b+z*yant)^(5/3))^(3/5);
    abs(y-yant);
  endwhile
endfunction
