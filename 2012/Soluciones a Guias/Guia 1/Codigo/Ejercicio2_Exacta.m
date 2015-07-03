function [ T ] = Ejercicio2_Exacta( CantX )
%Solucion exacta a la ODE T'' + T = 4(x-1/2)^2 -1
%Con T'(0) = -10 y T'(1) = 0

dx = 1/(CantX-1);
x = 0:dx:1;

k1 = -6;
k2 = -(6*cos(1)-4)/(sin(1));

T = k1*sin(x) + k2*cos(x) + 4*x.^2 - 4*x -8;

end

