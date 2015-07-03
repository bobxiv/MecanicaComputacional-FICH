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
%                       E = O(dt)
 function [Tn, TH] = Time_Advance_Explicit_1D(X, T0, tf, dx, dt, k, c, T_amb, CC1, CC2, caclcHistory)
 
    N    = length(T0);
    Tn   = T0;
    Tn_1 = zeros(N,1);
    
    if( caclcHistory )
        TH = zeros(length(0:dt:tf), N);%T History
        TH(1,:) = Tn';
    else
        TH(1,1) = 0;
    end
    
    r = dt/dx^2;%mi aux
    maxIter = length(0:dt:tf);
    %Iteramos en el tiempo
    for timeIt=2:maxIter
        
        %Stencil de iteracion en el cuerpo
        for i=2:N-1
            Tn_1(i) = Tn(i+1)*k*r +(-c*dt-2*r+1)*Tn(i)+Tn(i-1)*k*r+Q(X(i))*dt+T_amb*dt;
        end
        
        %Condicion de Contorno 1
        if( CC1.Tipo == 0 )%Es Dirichlet?
            Tn_1(1) = CC1.Valor;%T_i = T_0(x=0) se mantiene Dirichlet
        else if( CC1.Tipo == 1 )%Es Newmann?
                Tn_1(1) = Tn_1(2) + (dx/k) * CC1.Valor;%T_i = I_{i+1} + (dx/k)*q
            end
        end
        
        %Condicion de Contorno 2
        if( CC2.Tipo == 0 )%Es Dirichlet?
            Tn_1(N) = CC2.Valor;%T_i = T_0(x=0) se mantiene Dirichlet
        else if( CC2.Tipo == 1 )%Es Newmann?
                Tn_1(N) = Tn_1(N-1) - (dx/k) * CC2.Valor;%T_i = I_{i-1} - (dx/k)*q
            end
        end
       
        if( caclcHistory )
            TH(timeIt,:) = Tn_1';
        end
        
        Tn = Tn_1;%updatear Temperatura
        
    end
 
 end
