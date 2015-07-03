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
function [ a ] = SolveFD( R, Theta, EqDiffParam, CCPIzquierda, CCPDerecha, CCPInferior, CCPSuperior, IfPlot )

    %Creando la matriz del problema y variables basicas
    %--------------------------------------------------
    global K;
    global Nx;
    global Ny;
    Nx = length(R);
    Ny = length(Theta);
    K = sparse(Nx*Ny,Nx*Ny);
       
    global vecR;
    global vecTheta;
    vecR     = R;
    vecTheta = Theta;
    
    %Parametros de Ecuacion de Calor
    %-------------------------------
    global Q;%Fuente de Calor
    
    Q     = EqDiffParam.Q;
    
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
    for i= 2:Ny-1
        for j= 2:Nx-1
            ContribucionNodal(i,j);%Pone la contribucion de (i,j) a K
        end
    end
    
    %Agregacion de las CC y vector independiente
    %-------------------------------------------
    
    %Creacion del vector independiente
    global b;
    b = zeros(Nx*Ny,1);        %Tramo 2: d^2Fi/dx^2     = 0
    
    %Calculo de las condiciones de borde. Sus contribuciones a la matriz
    %y al vector independiente
    
    %Recorre CC inferior
    for i= 1:Nx
        CCInferior(i);
    end
    %Recorre CC superior
    for i= 1:Nx
        CCSuperior(i);
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
    for i= 2:Ny-1
       for j= 2:Nx-1
            ContribucionIndependiente(i,j);%Agrega la contribucion de (i,j) a b
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
        a = reshape(a, Nx, Ny);
        a = a';%Porque reshape devuelve column-wise y yo uso row-wise

        surf(R,Theta,a);
        xlabel('r');
        ylabel('theta');
        zlabel('T');
        title ('Temperaturas(r, theta)');
        
        
        figure
        [r,t]=meshgrid(R,Theta);
        surf(r.*cos(t),r.*sin(t),a);
        xlabel('x');
        ylabel('y');
        zlabel('T');
        title ('Temperaturas(x, y)');
    end

end