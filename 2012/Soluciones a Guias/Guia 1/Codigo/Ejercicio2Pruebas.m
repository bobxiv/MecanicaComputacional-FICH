% Ejercicio 2 Guia MDF 2012

% Definiendo la ecuacion diferencial a resolver
%----------------------------------------------

EqDiffParam.Dimension = '1D';%1D o 2D Si es 1D se consideran solo CC
                             %de Izquierda y Derecha

EqDiffParam.Q         = sym('1-4*(x-1/2)^2');
EqDiffParam.k         = sym('1');
EqDiffParam.c         = sym('1');
EqDiffParam.Tamb      = sym('0');

EqDiffParam.EqEspaciado = true;%nuestros nodos estaran equiespaciados

% Definiendo las CC
%------------------

% Nomenclatura:
%           CC                  P           Izquierda
%Condicion de Contorno      Parametros      Izquierda

CCPIzquierda.Tipo = 'Neumann';
CCPDerecha.Tipo   = 'Neumann';
CCPInferior.Tipo  = 'Dirichlet';%'Neumann';
CCPSuperior.Tipo  = 'Dirichlet';%'Neumann';

CCPIzquierda.q = sym(10);
CCPDerecha.q   = sym(0);


% Llamada a resolucion
%---------------------

minX = 0; maxX = 1;

NCantX= 10:10:100;
Errors = zeros(length(NCantX),1);
for i = 1:length(NCantX)
    %Dominio del problema
    
    X = linspace(minX, maxX, NCantX(i));
    Y = [1];

    Tmdf      = SolveFD( X, Y, EqDiffParam, CCPIzquierda, CCPDerecha, CCPInferior, CCPSuperior, false );
    Tanalitic = Ejercicio2_Exacta( NCantX(i) );
    e         = ErrorNorma2(Tanalitic, Tmdf);
    Errors(i) = e;
end


% Ploteo
%-------

%Errores
figure(1)
plot(NCantX, Errors,'r');
title('Desarrollo del error');
xlabel('Cantidad Nodos');
ylabel('RMS al cuadrado');

%Errores escala logaritmica
figure(2);
plot(log(NCantX), log(Errors),'r');
title('Desarrollo del error - Escala Logaritmica');
xlabel 'Cantidad Nodos'
ylabel 'RMS al cuadrado'

%Curvas de temperatura
figure(3)
plot(0:1/(NCantX(length(NCantX))-1):1, Tmdf,'r');
hold on;
Tanalitic = Ejercicio2_Exacta( length(0:1/10:1) );
plot(0:1/10:1, Tanalitic,'g*');
title('Curva de temperaturas');
xlabel('X');
ylabel('Temperatura');
legend('MDF', 'Analitica');