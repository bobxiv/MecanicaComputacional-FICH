% Definiendo la ecuacion diferencial a resolver
%----------------------------------------------

%   d(T)/dt + k*Lap(T) + c*(T-Tamb) + Q = 0

EqDiffParam.Dimension = '1D';%1D o 2D Si es 1D se consideran solo CC
                             %de Izquierda y Derecha

EqDiffParam.Q         = sym('x');
EqDiffParam.k         = sym('10');
EqDiffParam.c         = sym('0');
EqDiffParam.Tamb      = sym('0');

%Dominio del problema
minX = 0; maxX = 1; CantX = 4;
minY = 0; maxY = 1; CantY = 1;
X = linspace(minX, maxX, CantX);
Y = linspace(minY, maxY, CantY);

EqDiffParam.EqEspaciado = true;%nuestros nodos estaran equiespaciados

% Definiendo las CC
%------------------

% Nomenclatura:
%           CC                  P           Izquierda
%Condicion de Contorno      Parametros      Izquierda

CCPIzquierda.Tipo = 'Dirichlet';%'Neumann';
CCPDerecha.Tipo   = 'Neumann';%'Neumann';
CCPInferior.Tipo  = 'Neumann';%'Neumann';
CCPSuperior.Tipo  = 'Neumann';%'Neumann';

CCPIzquierda.Timpuesta = sym('1');
CCPDerecha.Timpuesta   = sym('0');
CCPInferior.Timpuesta  = sym(10);
CCPSuperior.Timpuesta  = sym(10);

CCPIzquierda.q = sym('0');
CCPDerecha.q   = sym('1');
CCPInferior.q  = sym(0);
CCPSuperior.q  = sym(0);

% Llamada a resolucion
%---------------------

figure(1);
%           vecX vecY 
T = SolveFD( X  , Y  , EqDiffParam, CCPIzquierda, CCPDerecha, CCPInferior, CCPSuperior, true );