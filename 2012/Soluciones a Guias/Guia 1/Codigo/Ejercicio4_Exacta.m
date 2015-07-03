function [y] = Ejercicio4_Exacta(x)
%Solucion exacta a phi'' + phi - 10 x (x+1) = 0
%con condiciones de borde phi(0) = 1 y phi(1) = 0

k1 = -21*cos(1)/sin(1);
k2 = 21;

y = k1*sin(x) + k2*cos(x) + 10*x.^2 + 10*x - 20;

end
