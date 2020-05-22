%Cálculo del calado en régimen uniforme con la Ec. de Manning
function yc=ycritico(Q,b,z)
  yant=0.1;
  yc=0;
  while abs(yc-yant)>1e-4
    yant=yc;
    yc=(Q^2/9.80665*(b+2*z*yant)/(b+z*yant)^3)^(1/3);
  endwhile
endfunction
