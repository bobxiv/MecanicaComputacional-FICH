%Este script resuelve la ecuacion de calor que se muestra abajo realizando
%un mapeo de un dominio en forma de circulo o corona (posiblemente abierto)
%a un rectangulo estructurado y cartesiano en las variables (r,theta),
%de forma tal que se puede usar el metodo de diferencias finitas sin mas
%consideraciones

% Ploteara en un grafico la solucion en el espacio (r,theta)
% Y en otro grafico la solucion en (x,t), osea en el circulo/corona


% Definiendo la ecuacion diferencial a resolver
%----------------------------------------------
% No es la ecuacion de calor completa y le faltan coeficientes y terminos
% pero solo de esta hice la transformacion... es dificil

%   d2(T)/dx2 + d2(T)/dy2 + Q = 0

EqDiffParam.Q         = sym('0');%en funcion de x e y

%Dominio del problema
minR     = 1; maxR     = 2;        CantR     = 10;
minTheta = 0; maxTheta = (3/2)*3.1415; CantTheta = 10;
R     = linspace(minR, maxR, CantR);
Theta = linspace(minTheta, maxTheta, CantTheta);

% Definiendo las CC
%------------------

% Nomenclatura:
%           CC                  P           Izquierda
%Condicion de Contorno      Parametros      Izquierda

CCPIzquierda.Tipo = 'Dirichlet';%'Neumann';
CCPDerecha.Tipo   = 'Dirichlet';%'Neumann';
CCPInferior.Tipo  = 'Dirichlet';%'Neumann';
CCPSuperior.Tipo  = 'Dirichlet';%'Neumann';

CCPIzquierda.Timpuesta = sym('100');%Radio minimo
CCPDerecha.Timpuesta   = sym('50');%Radio maximo
CCPInferior.Timpuesta  = sym('120');%Apertura de rosquilla
CCPSuperior.Timpuesta  = sym('80');%Apertura de rosquilla

CCPIzquierda.q = sym('100');
CCPDerecha.q   = sym('5');
CCPInferior.q  = sym('-100');
CCPSuperior.q  = sym('0');

% Llamada a resolucion
%---------------------

figure(1);
%           vecX vecY 
T = SolveFD( R  , Theta  , EqDiffParam, CCPIzquierda, CCPDerecha, CCPInferior, CCPSuperior, true );