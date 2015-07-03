%Resuelve la ecuacion de calor
%   k*Lap(T) + c*(T-Tamb) + Q = 0
%                                          EN 1D o 2D
%Donde k, c, Tamb y Q pueden ser ctes o funciones
%Sujeta a cualquier condicion de borde Neumann, Dirichlet o mixta
%El dominio de solucion es rectangular y solo se puede aplicar una CC a
%un borde (Izquierdo, Derecho, Inferior y Superior)
%Utiliza el metodo de Residuos Ponderados
%Requiere de un estado inicial. Los parametros son:
%
%                 - funN        funciones de forma
%                 - funW        funciones de peso
%                 - Cant        cantidad de ecuaciones
%                 - EqDiffParam parametros de la ecuacion diferencial
%                 - CCPXXX      parametros de la CC XXX
%                 - Plot        estructura con informacion de ploteo
%
%               Para ayuda ver el script de ejemplo MainScript.m
function [ a ] = SolveRP( funN, funW, Cant, EqDiffParam, CCPIzquierda, CCPDerecha, CCPInferior, CCPSuperior, Plot )

    %Creando la matriz del problema y variables basicas
    %--------------------------------------------------
    global K;
    %global Nx;
    %global Ny;
    %Nx = length(X);
    %Ny = length(Y);
    K = sparse(Cant,Cant);
    
%     if( strcmp(EqDiffParam.Dimension, '1D') && Ny>1 )
%         error('Intenta realizar una simulacon 1D pero el vector e y posee mas de 1 elemento! Si el problema es 1D se consideran solo los nodos en x')
%     end
%         
%     if( strcmp(EqDiffParam.Dimension, '1D') )
%         dy = 1;
%         Ny = 1;
%     end
    
%     global vecX;
%     global vecY;
%     vecX = X;
%     vecY = Y;
    
%     global EqEspaciado;%si los nodos estan equiespaciados
%     EqEspaciado = EqDiffParam.EqEspaciado;
    
    %Parametros de Ecuacion de Calor
    %-------------------------------
    global Q;%Fuente de Calor
    global k;%Difusividad
    global c;%Perdida de calor al ambiente
    global Tamb;%Temperatura ambiente
    global h;
    
    Q     = EqDiffParam.Q;
    k     = EqDiffParam.k;
    c     = EqDiffParam.c;
    Tamb  = EqDiffParam.Tamb;
    h     = EqDiffParam.h;
    
    global Dimension;%Dimension del problema
    Dimension = EqDiffParam.Dimension;
    
    global minX;
    global maxX;
    global minY;
    global maxY;
    minX = EqDiffParam.minX;
    maxX = EqDiffParam.maxX;
    minY = EqDiffParam.minY;
    maxY = EqDiffParam.maxY;
    
    %Parametros de CC
    %----------------
    global Tipocci;
    global Tipoccd;
    global Tipoccin;
    global Tipoccs;
    Tipocci  = CCPIzquierda.Tipo;%'Dirichlet' o 'Neumann'
    Tipoccd  = CCPDerecha.Tipo;  %'Dirichlet' o 'Neumann'
    Tipoccin = CCPInferior.Tipo; %'Dirichlet' o 'Neumann'
    Tipoccs  = CCPSuperior.Tipo; %'Dirichlet' o 'Neumann'
    
    global qcci; %conductividad en los borde izquierdo
    global qccd; %conductividad en los borde derecho
    global qccin;%conductividad del borde inferior
    global qccs; %conductividad del borde superior
    global Tcci; %Temperatura impuesta en el borde izquierdo
    global Tccd; %Temperatura impuesta en el borde derecho
    global Tccin;%Temperatura impuesta en el borde inferior
    global Tccs; %Temperatura impuesta en el borde superior
    qcci  = CCPIzquierda.q;
    qccd  = CCPDerecha.q;
    qccin = CCPInferior.q;
    qccs  = CCPSuperior.q;
    Tcci  = CCPIzquierda.Timpuesta;
    Tccd  = CCPDerecha.Timpuesta;
    Tccin = CCPInferior.Timpuesta;
    Tccs  = CCPSuperior.Timpuesta;
    
    %Parametros de MRP
    %-----------------
    
    global N;%Funciones de forma
    global W;%Funciones de peso
    N = funN;
    W = funW;
    
    global M;
    M = Cant;
  
    %Ensablado de la matriz global para nodos internos
    %-------------------------------------------------
    
    %Recorre por cada ecuacion y termino
    %if( strcmp(Dimension, '2D') )
        for l= 1:M
            for m= 1:M
                ContribucionNodal(l,m);%Agrega la contribucion de (l,m) a K
            end
        end
    %else
        %for j= 2:Nx-1
        %    ContribucionNodal(1,j);%Agrega la contribucion de (1,j) a K
        %end
    %end
    
    %Ensablado de las CC y vector independiente
    %------------------------------------------
    
    %Creacion del vector independiente
    global b;
    b = zeros(M,1);        %Tramo 2: d^2Fi/dx^2     = 0
    
    %Calculo de las condiciones de borde. Sus contribuciones a la matriz
    %y al vector independiente
    
%     if( strcmp(Dimension, '2D') )
%         %Recorre CC inferior
%         for i= 1:Nx
%             CCInferior(i);
%         end
%         %Recorre CC superior
%         for i= 1:Nx
%             CCSuperior(i);
%         end
%     end
%     %Recorre CC izquierda
%     for j= 1:Ny
%         CCIzquierda(j);
%     end
%     %Recorre CC derecha
%     for j= 1:Ny
%         CCDerecha(j);
%     end

    %CondicionesDirichlet();%Agrega (si hay) condiciones Dirichlet
    
    %Calculo del vector independiente de los elementos internos
    %if( strcmp(Dimension, '2D') )
        for l= 1:M
           %for m= 1:Cant
                ContribucionIndependiente(l);%Agrega la contribucion de (l,m) a b
           %end
        end
%     else
%          for j= 2:Nx-1
%             ContribucionIndependiente(1,j);%Agrega la contribucion de (1,j) a b
%         end
    %end
    
    %Resolucion de la matriz global
    %------------------------------
    
    %spy(K);
    %p = det(K);
    
    %Soluciona el sistema 
    %   K*a=b -> a = inv(K)*b
    a = K\b;
    
    %q = det(K);
    
    %Plotear si se lo pide
    %---------------------
    
    if( Plot.hacer )
        %dibujar:
        if( strcmp(Dimension, '2D') )
            
            Tval = zeros(length(Plot.Y)*length(Plot.X),1);

            for y=1:length(Plot.Y)
                for x=1:length(Plot.X)
                    Tval(length(Plot.X)*(y-1)+x) = T2d(Plot.X(x),Plot.Y(y),a,N,M);
                end
            end
            Tval = reshape(Tval, length(Plot.X), length(Plot.Y));
            Tval = Tval';%Porque reshape devuelve column-wise y yo uso row-wise
            
            surf(Plot.X,Plot.Y,Tval);
            xlabel('x');
            ylabel('y');
            zlabel('T');
            
        else
           
            Tval = zeros(length(Plot.X),1);

            for x=1:length(Plot.X)
                Tval(x) = T1d(Plot.X(x),a,N,M);
            end
            
            plot(Plot.X,Tval);
            xlabel('x');
            ylabel('T');
            
        end
        title (Plot.titulo);
    end

end

function val = T1d(valX,a,N,Cant) 
    syms symm x; 
    val = 0; 
    for k=1:Cant
        val = val + a(k)*subs(N,{symm,x},{k,valX}); 
    end
end

function val = T2d(valX,valY,a,N,Cant) 
    syms symm x y; 
    val = 0; 
    for k=1:Cant
        val = val + a(k)*subs(N,{symm,x,y},{k,valX,valY}); 
    end
end