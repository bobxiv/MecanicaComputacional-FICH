% Ejercicio 4 Guia MDF 2012

% Definiendo la ecuacion diferencial a resolver
%----------------------------------------------

%   d(T)/dt + k*Lap(T) + c*(T-Tamb) + Q = 0

EqDiffParam.Dimension = '1D';%1D o 2D Si es 1D se consideran solo CC
                             %de Izquierda y Derecha

EqDiffParam.Q         = sym('-10*x*(x+1)');
EqDiffParam.k         = sym('1');
EqDiffParam.c         = sym('1');
EqDiffParam.Tamb      = sym('0');

%ACA ESTA LA CLAVE DEL EJERCICIO
%Dominio del problema
X = [0, 0.05, 0.15, 0.35, 0.5, 0.75, 1];
Y = [1];

EqDiffParam.EqEspaciado = false;%nuestros nodos estaran equiespaciados

% Definiendo las CC
%------------------

% Nomenclatura:
%           CC                  P           Izquierda
%Condicion de Contorno      Parametros      Izquierda

CCPIzquierda.Tipo = 'Dirichlet';
CCPDerecha.Tipo   = 'Dirichlet';
CCPInferior.Tipo  = 'Dirichlet';
CCPSuperior.Tipo  = 'Dirichlet';

CCPIzquierda.Timpuesta = sym(1);
CCPDerecha.Timpuesta   = sym(0);

% Llamada a resolucion
%---------------------

figure(1);
%           vecX vecY 
T = SolveFD( X  , Y  , EqDiffParam, CCPIzquierda, CCPDerecha, CCPInferior, CCPSuperior, false );
plot(X,T,'r');
xlabel('x');
ylabel('T');
title ('Temperaturas(x)');
hold on;
plot(X,T,'r*');


minX = 0; maxX = 1; CantX = length(X);
X = linspace(minX, maxX, CantX);
EqDiffParam.EqEspaciado = true;%nuestros nodos estaran equiespaciados

hold on;
%           vecX vecY 
T = SolveFD( X  , Y  , EqDiffParam, CCPIzquierda, CCPDerecha, CCPInferior, CCPSuperior, false );
plot(X,T,'b');
xlabel('x');
ylabel('T');
title ('Temperaturas(x)');
hold on;
plot(X,T,'b*');

T = Ejercicio4_Exacta(linspace(0,1,100));
hold on;
plot(linspace(0,1,100),T,'k');
xlabel('x');
ylabel('T');
title ('Temperaturas(x)');

legend('Espaciado heterogeneo','','Espaciado homogeneo','','Solucion Analitica');