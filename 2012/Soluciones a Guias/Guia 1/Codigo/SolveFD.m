%Resuelve la ecuacion de calor
%   d(T)/dt + k*Lap(T) + c*(T-Tamb) + Q = 0
%                                          EN 1D o 2D
%Donde k, c, Tamb y Q pueden ser ctes o funciones
%Sujeta a cualquier condicion de borde Neumann, Dirichlet o mixta
%El dominio de solucion es rectangular y solo se puede aplicar una CC a
%un borde (Izquierdo, Derecho, Inferior y Superior)
%Utiliza el metodo de Diferencias Finitas
%Requiere de un estado inicial. Los parametros son:
%
%                 - X           vector de las posiciones x
%                 - Y           vector de las posiciones y
%                 - EqDiffParam parametros de la ecuacion diferencial
%                 - CCPXXX      parametros de la CC XXX
%                 - IfPlot      verdadero si se desea que plotee
%
%               Para ayuda ver el script de ejemplo MainScript.m
function [ a ] = SolveFD( X, Y, EqDiffParam, CCPIzquierda, CCPDerecha, CCPInferior, CCPSuperior, IfPlot )

    %Creando la matriz del problema y variables basicas
    %--------------------------------------------------
    global K;
    global Nx;
    global Ny;
    Nx = length(X);
    Ny = length(Y);
    K = sparse(Nx*Ny,Nx*Ny);
    
    if( strcmp(EqDiffParam.Dimension, '1D') && Ny>1 )
        error('Intenta realizar una simulacon 1D pero el vector e y posee mas de 1 elemento! Si el problema es 1D se consideran solo los nodos en x')
    end
        
    if( strcmp(EqDiffParam.Dimension, '1D') )
        dy = 1;
        Ny = 1;
    end
    
    global vecX;
    global vecY;
    vecX = X;
    vecY = Y;
    
    global EqEspaciado;%si los nodos estan equiespaciados
    EqEspaciado = EqDiffParam.EqEspaciado;
    
    %Parametros de Ecuacion de Calor
    %-------------------------------
    global Q;%Fuente de Calor
    global k;%Difusividad
    global c;%Perdida de calor al ambiente
    global Tamb;%Temperatura ambiente
    
    Q     = EqDiffParam.Q;
    k     = EqDiffParam.k;
    c     = EqDiffParam.c;
    Tamb  = EqDiffParam.Tamb;
    
    global Dimension;%Dimension del problema
    Dimension = EqDiffParam.Dimension;
    
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
  
    %Ensablado de la matriz global para nodos internos
    %-------------------------------------------------
    
    %Recorre por nodos internos y acumula 
    if( strcmp(Dimension, '2D') )
        for i= 2:Ny-1
            for j= 2:Nx-1
                ContribucionNodal(i,j);%Agrega la contribucion de (i,j) a K
            end
        end
    else
        for j= 2:Nx-1
            ContribucionNodal(1,j);%Agrega la contribucion de (1,j) a K
        end
    end
    
    %Ensablado de las CC y vector independiente
    %------------------------------------------
    
    %Creacion del vector independiente
    global b;
    b = zeros(Nx*Ny,1);        %Tramo 2: d^2Fi/dx^2     = 0
    
    %Calculo de las condiciones de borde. Sus contribuciones a la matriz
    %y al vector independiente
    
    if( strcmp(Dimension, '2D') )
        %Recorre CC inferior
        for i= 1:Nx
            CCInferior(i);
        end
        %Recorre CC superior
        for i= 1:Nx
            CCSuperior(i);
        end
    end
    %Recorre CC izquierda
    for j= 1:Ny
        CCIzquierda(j);
    end
    %Recorre CC derecha
    for j= 1:Ny
        CCDerecha(j);
    end
    
    %Calculo del vector independiente de los elementos internos
    if( strcmp(Dimension, '2D') )
        for i= 2:Ny-1
           for j= 2:Nx-1
                ContribucionIndependiente(i,j);%Agrega la contribucion de (i,j) a b
           end
        end
    else
         for j= 2:Nx-1
            ContribucionIndependiente(1,j);%Agrega la contribucion de (1,j) a b
        end
    end
    
    %Resolucion de la matriz global
    %------------------------------
    
    %Soluciona el sistema 
    %   K*a=b -> a = inv(K)*b
    a = K\b;
    
    %Plotear si se lo pide
    %---------------------
    
    if( IfPlot )
        %dibujar:
        if( strcmp(Dimension, '2D') )
            a = reshape(a, Nx, Ny);
            a = a';%Porque reshape devuelve column-wise y yo uso row-wise
            
            surf(X,Y,a);
            xlabel('x');
            ylabel('y');
            zlabel('T');
            title ('Temperaturas(x, y)');
        else
            plot(X,a);
            xlabel('x');
            ylabel('T');
            title ('Temperaturas(x)');
        end
    end

end