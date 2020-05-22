clear;clc
close all
#Cálculo del régimen gradualmente variado en canales de sección uniforme

#Las condiciones de contorno y de signo deltay deben usarse conociendo previamente el tipo de curva (S1, S2, S3, F1, F2, F3, C1 ó C2) que se está integrando.

addpath('./src');

x=[];
y=[];
H0=[];
I=[];

Q=55;%Caudal
I0=0.01;%Pendiente topográfica rasante del tramo del canal
n=0.012;%Aspereza de Manning
b=5;%Ancho de la solera del canal
z=0;%Inclinación talud quijeros
deltay=-0.001;%Importante introducir el signo (+ ó -) según la curva a integrar

yc=ycritico(Q,b,z);
Ic=(Q*n*(b+2*yc*sqrt(1+z^2))^(2/3)/(b*yc+z*yc^2)^(5/3))^2;
y0=yManning(Q,n,b,z,I0);

%Condición de contorno.
x(1)=0;%Valor x de la condición de contorno.
y(1)=yc;
H0(1)=y(1)+(Q/(b*y(1)+z*y(1)^2))^2/2/9.80665;
I(1)=(Q*n*(b+2*y(1)*(1+z^2)^(1/2))^(2/3)/(b*y(1)+z*y(1)^2)^(5/3))^2;

numpasos=floor((y0-yc)/deltay)

if I0>Ic
  display('Régimen rápido')
else
  display('Régimen lento')
endif

for i=1:numpasos
    y(i+1)=y(i)+deltay;
    H0(i+1)=y(i+1)+(Q/(b*y(i+1)+z*y(i+1)^2))^2/2/9.80665;
    I(i+1)=(Q*n*(b+2*y(i+1)*(1+z^2)^(1/2))^(2/3)/(b*y(i+1)+z*y(i+1)^2)^(5/3))^2;
    deltaH0=H0(i+1)-H0(i);
    Imed=(I(i)+I(i+1))/2;
    deltax=deltaH0/(I0-Imed);
    x(i+1)=x(i)+deltax; 
endfor

plot(x,y,x,H0)
