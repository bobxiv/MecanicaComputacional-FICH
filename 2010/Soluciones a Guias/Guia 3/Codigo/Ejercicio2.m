function [] = Ejercicio2(  )

    %Analitic Solution: T = - (Q*t^2)/2 + (Q/2 - 1)*t + 1
    TReal = ResolverProblema(8000, false);
    dif = zeros(length(10:2:76),1);
    it = 1;
    for i= 10:2:76
        TApprox = ResolverProblema(i, false);
        dif(it) = ErrorCuadraticoMedio(TReal , TApprox);
        it=it+1;
    end
    
    plot(10:2:76,dif);
    title('Error');
    xlabel('Nodos');
    ylabel('Error');

    %ResolverProblema(5000,true);

end


%Nn es la cantidad de nodos con que se resulve
function [a] = ResolverProblema( Nn , DoPlot)

malla = GeneradorMallaUnidimensional(Nn, 1);

%Seteo de Ctes
malla.k = 1.0;

for m=1:ceil(malla.Nn/2.0)
    malla.Q(m) = 1.0;
end
for m=(ceil(malla.Nn/2.0)+1):malla.Nn
    malla.Q(m) = 0.0;
end

%Seteo de Condiciones Iniciales

malla.IzquierdaCI = 'Dirichlet';
malla.IzquierdaPhi  = 1.0;
malla.DerechaCI = 'Dirichlet';
malla.DerechaPhi = 0.0;
%malla.DerechaCI = 'Newmann';%Para correr con condiciones Newmann
%descomentar aca
malla.Derechaq = 0.0;

%Comienzo de resolucion

[elemCount,nodesIndexs] = size(malla.e);

%Matrices globales
K = zeros(malla.Nn,malla.Nn);
K = sparse(K);
f = zeros(malla.Nn,1);
a = zeros(malla.Nn,1);

%Calculo de la matriz elemental
%Y Ensamblado de la matriz global
for m=1:elemCount
    
    K_e = [  1/malla.h^2  , -1/malla.h^2; %i o 1
            -1/malla.h^2  ,  1/malla.h^2; %j o 2
            ];
        
    K_e = K_e * malla.k * malla.h;%malla.area(m);
        
    %Se agrega la contribucion elemental a la matriz global
    K( malla.e(m,1) , malla.e(m,1) ) = K( malla.e(m,1) , malla.e(m,1) ) + K_e(1,1);
    K( malla.e(m,1) , malla.e(m,2) ) = K( malla.e(m,1) , malla.e(m,2) ) + K_e(1,2);
    
    K( malla.e(m,2) , malla.e(m,1)) = K( malla.e(m,2) , malla.e(m,1) ) + K_e(2,1);
    K( malla.e(m,2) , malla.e(m,2)) = K( malla.e(m,2) , malla.e(m,2) ) + K_e(2,2);
    
    %Contribucion de fe con termino fuente
    Qprom = (malla.Q(malla.e(m,1)) + malla.Q(malla.e(m,2))) / 2.0;
    f(malla.e(m,1)) = f(malla.e(m,1)) + 1/3.0 * malla.h * Qprom;%malla.area(m) * Qprom;
    f(malla.e(m,2)) = f(malla.e(m,2)) + 1/3.0 * malla.h * Qprom;%malla.area(m) * Qprom;

end

%Contribucion de Condiciones de contorno a f

%Izquierda CI
if strcmp(malla.IzquierdaCI , 'Dirichlet' ) == 1
    %Condicion Dirichlet izquierda
    K( 1 , : ) = zeros(1,malla.Nn);
    K( 1 , 1 ) = 1;
    f( 1 ) = malla.IzquierdaPhi;
end

%Derecha CI
if strcmp(malla.DerechaCI , 'Dirichlet' ) == 1
    %Condicion Dirichlet derecha
    K( malla.Nn , :       ) = zeros(1,malla.Nn);
    K( malla.Nn , malla.Nn ) = 1;
    f( malla.Nn ) = malla.DerechaPhi;
else
    %Condicion Newmann derecha
    f( malla.Nn ) = f( malla.Nn ) - 1/2 * malla.Derechaq * malla.h;
end


%Solucion del problema
a = K\f;

%Graficando
if DoPlot
    plot(malla.n,a);

    xlabel('x');
    zlabel('T');
    title('k*d^2T/dx^2+Q=0')
end

%NOTA: NO TIENE SENTIDO REFINAR DENTRO DE LOS ELEMENTOS, OSEA HACER UNA
%MALLA MAS CHICA Y EVALUAR LA FUNCION EN ESA MALLA MAS CHICA PORQUE POR
%USAR UNA FUNCIONES DE FORMA LINEALES DENTRO DE 1 ELEMENTO HACE UNA
%INTERPOLACION LINEAL ENTRE LOS 2 NODOS... osea que el grafico no se veria
%mas refinado


end

