clear;clc
#Cálculo del régimen gradualmente variado en canales de sección no uniforme

x=[];
y=[];
H0=[];
I=[];
yc=[];
Q=[];


I0=0.003;
n=0.012;
z=0;%Sección trapecial, con z=0 es rectangular
L=23;%Longitud de tramo de canal
N=23;%Número de secciones
deltax=-L/N;%Distancia entre secciones, negativo hacia aguas arriba de la sección de control
q0=55;q1=(0-(q0))/(N*deltax);%Acumulación de caudal con x, en la forma Q=q0+q1*x
tol=1e-6;%Parámetro para considerar alcanzada la convergencia del cálculo iterativo

x(1)=0;%Posción x de la sección de control
b0=5;b1=(b0-3)/N;%Cambio de la sección con x, en la forma b=b0+b1*x
y(1)=(((q0+q1*x(1))/(b0+b1*x(1)))^2/9.80665)^(1/3);
yc(1)=y(1);

for i=1:N
  %Se calculan N secciones separadas una distancia deltax
  i  
  deltay=1e-2;deltayant=0;
  while abs(deltay-deltayant)>tol
    deltayant=deltay;
    Q(i)=q0+q1*x(i);
    H0(i)=y(i)+(Q(i)/(b0+b1*x(i))/y(i))^2/2/9.80665;
    I(i)=(Q(i)*n*((b0+b1*x(i))+2*y(i)*(1+z^2)^(1/2))^(2/3)/((b0+b1*x(i))*y(i)+z*y(i)^2)^(5/3))^2;
    x(i+1)=x(i)+deltax;
    y(i+1)=y(i)+deltay;
    Q(i+1)=q0+q1*x(i+1);
    H0(i+1)=y(i+1)+(Q(i+1)/(b0+b1*x(i+1))/y(i+1))^2/2/9.80665;
    I(i+1)=(Q(i+1)*n*((b0+b1*x(i+1))+2*y(i+1)*(1+z^2)^(1/2))^(2/3)/((b0+b1*x(i+1))*y(i+1)+z*y(i+1)^2)^(5/3))^2;
    deltaH0=H0(i+1)-H0(i);
    Imed=(I(i)+I(i+1))/2;
    dH0dbmed=-1/3*((Q(i)+Q(i+1))/2)^2/9.80665*(1/((b0+b1*x(i))^3+y(i)^2)+1/((b0+b1*x(i+1))^3+y(i+1)^2))/2;
    deltaxdif=deltaH0/(I0-Imed-b1*dH0dbmed);
    deltay=deltay*(abs(deltax/deltaxdif))^(1/2)
  endwhile
  yc(i+1)=(((q0+q1*x(i+1))/(b0+b1*x(i+1)))^2/9.80665)^(1/3);
endfor

plot(x,y,'+-',x,H0,x,yc)
xlabel('x(m)')
ylabel('yc(m), y(m), Ho(m)')
