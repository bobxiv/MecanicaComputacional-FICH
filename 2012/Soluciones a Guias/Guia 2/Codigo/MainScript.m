% Definiendo la ecuacion diferencial a resolver
%----------------------------------------------

%   k*Lap(T) + c*(T-Tamb) + Q = 0

EqDiffParam.Dimension = '2D';%1D o 2D Si es 1D se consideran solo CC
                             %de Izquierda y Derecha

EqDiffParam.Q         = sym('0');
EqDiffParam.k         = sym('1');
EqDiffParam.c         = sym('0');
EqDiffParam.Tamb      = sym('0');
EqDiffParam.h         = sym('0');

%Dominio del problema
EqDiffParam.minX = -1; EqDiffParam.maxX = 1;
EqDiffParam.minY = -1; EqDiffParam.maxY = 1;

EqDiffParam.EqEspaciado = true;%nuestros nodos estaran equiespaciados

% Definiendo las CC
%------------------

% Nomenclatura:
%           CC                  P           Izquierda
%Condicion de Contorno      Parametros      Izquierda

CCPIzquierda.Tipo = 'Dirichlet';%'Neumann';
CCPDerecha.Tipo   = 'Dirichlet';%'Neumann';
CCPInferior.Tipo  = 'Dirichlet';%'Dirichlet';
CCPSuperior.Tipo  = 'Dirichlet';%'Neumann';

CCPIzquierda.Timpuesta = sym('1-y^2');
CCPDerecha.Timpuesta   = sym('1-y^2');
CCPInferior.Timpuesta  = sym('1-x^2');
CCPSuperior.Timpuesta  = sym('1-x^2');

CCPIzquierda.q = sym('0');
CCPDerecha.q   = sym('5');
CCPInferior.q  = sym('0');
CCPSuperior.q  = sym('0');

% Llamada a resolucion
%---------------------
syms symm syml x y;
funN = sym('x^(symm-1)*y^(symm-1)');%es fundamental que  se tenga a la
funW = sym('x^(syml-1)*y^(syml-1)');%cte en la funcion de forma!!!
%funN = sym('x^(symm-1)');
%funW = sym('x^(syml-1)');

Cant = 5;

Plot.hacer = true;
Plot.titulo = [char(EqDiffParam.k), '*Lap(T) + ', char(EqDiffParam.c), '*(T-', char(EqDiffParam.Tamb), ') + ', char(EqDiffParam.Q), ' = 0'];
Plot.X = linspace(EqDiffParam.minX, EqDiffParam.maxX, 12);
Plot.Y = linspace(EqDiffParam.minY, EqDiffParam.maxY, 12);

figure(1);
%           forma    peso 
T = SolveRP( funN  , funW , Cant, EqDiffParam, CCPIzquierda, CCPDerecha, CCPInferior, CCPSuperior, Plot );