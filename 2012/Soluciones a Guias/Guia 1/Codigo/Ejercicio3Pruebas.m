% Ejercicio 3 Guia MDF 2012

% Definiendo la ecuacion diferencial a resolver
%----------------------------------------------

%   d(T)/dt + k*Lap(T) + c*(T-Tamb) + Q = 0

EqDiffParam.Dimension = '2D';%1D o 2D Si es 1D se consideran solo CC
                             %de Izquierda y Derecha

EqDiffParam.Q         = sym('-2*(x^2+y^2)');
EqDiffParam.k         = sym('1');
EqDiffParam.c         = sym('0');
EqDiffParam.Tamb      = sym('0');

%Dominio del problema
minX = 0; maxX = 3; CantX = 7;
minY = 0; maxY = 1; CantY = 5;
X = linspace(minX, maxX, CantX);
Y = linspace(minY, maxY, CantY);

EqDiffParam.EqEspaciado = true;%nuestros nodos estaran equiespaciados

% Definiendo las CC
%------------------

% Nomenclatura:
%           CC                  P           Izquierda
%Condicion de Contorno      Parametros      Izquierda

CCPIzquierda.Tipo = 'Dirichlet';
CCPDerecha.Tipo   = 'Dirichlet';
CCPInferior.Tipo  = 'Dirichlet';
CCPSuperior.Tipo  = 'Dirichlet';

CCPIzquierda.Timpuesta = sym(0);
CCPDerecha.Timpuesta   = sym(0);
CCPInferior.Timpuesta  = sym(0);
CCPSuperior.Timpuesta  = sym(0);

% Llamada a resolucion (a)
%-------------------------

figure(1);
Tmdf      = SolveFD( X, Y, EqDiffParam, CCPIzquierda, CCPDerecha, CCPInferior, CCPSuperior, true );
title('Dirichlet T=0 en CC')


% Llamada a resolucion (b)
%-------------------------

CCPIzquierda.Timpuesta = sym(100);
CCPDerecha.Timpuesta   = sym(100);
CCPInferior.Timpuesta  = sym(100);
CCPSuperior.Timpuesta  = sym(100);

figure(2);
Tmdf      = SolveFD( X, Y, EqDiffParam, CCPIzquierda, CCPDerecha, CCPInferior, CCPSuperior, true );
title('Dirichlet T=100 en CC')

% Llamada a resolucion (c)
%-------------------------

CCPIzquierda.Tipo = 'Dirichlet';
CCPDerecha.Tipo   = 'Dirichlet';
CCPInferior.Tipo  = 'Neumann';
CCPSuperior.Tipo  = 'Dirichlet';

CCPIzquierda.Timpuesta = sym(100);
CCPDerecha.Timpuesta   = sym(100);
CCPInferior.q          = sym(0);
CCPSuperior.Timpuesta  = sym(100);

figure(3);
Tmdf      = SolveFD( X, Y, EqDiffParam, CCPIzquierda, CCPDerecha, CCPInferior, CCPSuperior, true );
title('CC Mixtas')
