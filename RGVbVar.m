clear;clc
close all
#Cálculo del régimen gradualmente variado en canales de sección no uniforme

x=[];
y=[];
H0=[];
I=[];
yc=[];
Q=[];
b=[];


I0=0.003;
n=0.012;
z=0;%Sección trapecial, con z=0 es rectangular
L=40;%Longitud de tramo de canal
Laliv=23;%Longitud del aliviadero
%Naliv=23;%Número de secciones del aliviadero
N=40;%Número de secciones a calcular
deltax=-L/N;%Distancia entre secciones, negativo hacia aguas arriba de la sección de control
q0=55;q1=(q0+q0/Laliv)/Laliv;%Acumulación de caudal con x, en la forma Q=q0+q1*x
tol=1e-6;%Parámetro para considerar alcanzada la convergencia del cálculo iterativo

x(1)=0;%Posción x de la sección de control
b0=5;b1=(b0-2+(b0-2)/Laliv)/Laliv;%Cambio de la sección con x, en la forma b=b0+b1*x
b(1)=b0;
y(1)=((q0/b0)^2/9.80665)^(1/3);
Q(1)=q0;
yc(1)=y(1);

for i=1:N
  %Se calculan N secciones separadas una distancia deltax
  deltay=1e-2;deltayant=0;
  while abs(deltay-deltayant)>tol
    deltayant=deltay;
    
    x(i+1)=x(i)+deltax;
    if x(i+1)>=(Laliv-L)
      Q(i+1)=q0;
      b(i+1)=b0;
    else
      Q(i+1)=q0+q1*(x(i)+(L-Laliv));
      b(i+1)=b0+b1*(x(i)+(L-Laliv));
    endif
    H0(i)=y(i)+(Q(i)/b(i)/y(i))^2/2/9.80665;
    I(i)=(Q(i)*n*(b(i)+2*y(i)*(1+z^2)^(1/2))^(2/3)/(b(i)*y(i)+z*y(i)^2)^(5/3))^2;
    y(i+1)=y(i)+deltay;
    H0(i+1)=y(i+1)+(Q(i+1)/b(i+1)/y(i+1))^2/2/9.80665;
    I(i+1)=(Q(i+1)*n*(b(i+1)+2*y(i+1)*(1+z^2)^(1/2))^(2/3)/(b(i+1)*y(i+1)+z*y(i+1)^2)^(5/3))^2;
    deltaH0=H0(i+1)-H0(i);
    Imed=(I(i)+I(i+1))/2;
    dH0dbmed=-1/3*((Q(i)+Q(i+1))/2)^2/9.80665*(1/(b(i)^3+y(i)^2)+1/(b(i+1)^3+y(i+1)^2))/2;
    if x(i+1)>=(Laliv-L)
      deltaxdif=deltaH0/(I0-Imed);
    else
      deltaxdif=deltaH0/(I0-Imed-b1*dH0dbmed);
    endif
    deltay=deltay*(abs(deltax/deltaxdif))^(1/2);
  endwhile
  yc(i+1)=((Q(i)/b(i+1))^2/9.80665)^(1/3);
endfor

subplot(1,2,1)
plot(x,y,'+-',x,H0,x,yc)
xlabel('x(m)')
ylabel('yc(m), y(m), Ho(m)')

subplot(1,2,2)
plot(x,b,x,Q)
xlabel('x(m)')
ylabel('b(m), Q(m3/s)')