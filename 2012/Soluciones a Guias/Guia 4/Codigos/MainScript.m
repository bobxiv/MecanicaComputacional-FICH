% Definiendo la ecuacion diferencial a resolver
%----------------------------------------------

%   k*Lap(T) + c*(T-Tamb) + Q = 0

ProbElasticidadViga.F = [0, 0];%Fuerza de volumen

% Definiendo las CC
%------------------

ProbElasticidadViga.FlexDerecha = [0, -0.1];

%Inicializamos la malla que vamos a usar para el calculo
%-------------------------------------------------------

%Malla de triangulos
NombreArchivoMalla = 'mallaVigaTriangulos';

%Malla de cuadrangulo
%NombreArchivoMalla = 'mallaVigaCuadrilateros';

% Llamada a resolucion
%---------------------

Plot.hacer = true;
Plot.titulo = 'Elasticidad';

figure(1);
T = SolveFEM_Elasticidad( ProbElasticidadViga, NombreArchivoMalla, Plot );