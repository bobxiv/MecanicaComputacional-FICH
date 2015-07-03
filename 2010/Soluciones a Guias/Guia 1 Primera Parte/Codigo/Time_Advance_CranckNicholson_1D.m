%La ecuacion que resulve es:
%          -dT/dt + d^2T/dx^2 - T + Q = 0
%                                          EN 1D
%Requiere de un estado inicial. Los parametros son:
%
%                      - X  el vector de puntos en el espacio
%                      - T0 como estado inicial
%                      - tf el tiempo final
%                      - dx el intervalo espacial
%                      - dt el intervalo de tiempo
%                      - k    cte de conductividad
%                      - c    cte de perdida de calor al medio ambiente
%                      - Tamb cte de temperatura ambiente 
%                      - CC1  la condicion de contorno 1
%                      - CC2  la condicion de contorno 2
%                      - calcHistory si es verdadero devuelve el historial
%                        de los estados en TH
%
%                      Q debe ser una funcion de x definida 
%                      en este entorno!!!
%                      CC1 y CC2 deben estar pasadas en esta estructura:
%                      struct('Tipo', [0,1] , 'Valor', [numerico])
%                      Con Tipo 0 -> Dirichlet y Tipo 1 -> Newmann
%
%                       E = O(dt^2)
function [Tn, TH] = Time_Advance_CranckNicholson_1D(X, T0, tf, dx, dt, k, c, T_amb, CC1, CC2, caclcHistory)

    %Armado de Matriz de solucion temporal de Crank-Nicholson

    N = length(X);
    A = eye(N);
    %Stencil de la matriz
    for i= 2:N-1
        A(i,i)   = -(dx^2/dt)-k-(dx^2*c)/2;
        A(i,i+1) = k/2;
        A(i,i-1) = k/2;
    end

    
    %Condicion de Contorno 1
    if( CC1.Tipo == 0 )%Es Dirichlet?
        A(1,1) = CC1.Valor;
    else if( CC1.Tipo == 1 )%Es Newmann?
            A(1,1) = -1.0/dx;   %Derivada 
            A(1,2) =  1.0/dx;    %hacia adelante
        end
    end
    
    %Condicion de Contorno 2
    if( CC2.Tipo == 0 )%Es Dirichlet?
        A(N) = CC2.Valor;
    else if( CC2.Tipo == 1 )%Es Newmann?
            A(N,N-1) = -1.0/dx;   %Derivada 
            A(N,N)   =  1.0/dx;    %hacia atras 
        end
    end
    
   
    sparse(A);
    
    
    
    Tn   = T0;
    Tn_1 = zeros(N,1);
    if( caclcHistory )
        TH = zeros(length(0:dt:tf),N);%T History
        TH(1,:) = Tn';
    else
        TH(1,1) = 0;
    end
    
    for timeIt=2:length(0:dt:tf)
    
        %b = MakeBforCranckNicholson(Tn,dx,dt*(k-1),N);% Con esto y sin
        %actualizar Tn me da una hoya que se achica, con la linea de abajo
        %y actualizando Tn me da lo mismo que el explicito... ver que es lo
        %correcto!
        b = MakebforCranckNicholson(X,Tn,dx,dt,k,c,T_amb);
        
        Tn_1 = A\b;
        
        if( caclcHistory )
            TH(timeIt,:) = Tn_1';
        end
        
        Tn = Tn_1;%updatear Temperatura
    
    end

end

%Evalue el vector independiente del metodo Cranck Nicholson
%Este se debe re calcular cuando se itera. La ecuacion para la cual sirve
%es la que se muestra arriba.
function b = MakebforCranckNicholson(X,Tn,dx,dt,k, c, T_amb)

    N = length(Tn);
    
    b = zeros(N,1);
    b(1) = 1;% CC:    T(x=0,t) = 1   
    %for i=2:ceil(N/2)
    %    b(i) = -Tn(i-1)/(2*dx^2)+(+1.0/2-1.0/dt+1.0/(dx^2))*Tn(i)-Tn(i+1)/(2*dx^2) - 1;
    %end
    %for i= ceil(N/2)+1:N-1
    %    b(i) = -Tn(i-1)/(2*dx^2)+(+1.0/2-1.0/dt+1.0/(dx^2))*Tn(i)-Tn(i+1)/(2*dx^2);
    %end
    for i= 2:N-1
        b(i) = (-k/2)*Tn(i-1) + ((-dx^2/dt)+k+(c*dx^2)/2)*Tn(i) + (-k/2)*Tn(i+1) - Q(X(i))*dx^2 - c*dx^2*T_amb;
    end
    
    b(N,1) = 0;% CC:    q      = 0
    
end

