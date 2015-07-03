function [ ] = Ejericio2(  )

    disp('Ecuacion del calor con termino temporal. Solucion por metodo Explicito:')
    disp('Calculando solucion pseudo-analitica')

    %                           N ,tf ,      dt ,  F ,feedback? ,useFourier
    Preparar_Resolver(21,  1, 0.001, 0.3,     true, true );
    Tanalitic=Preparar_Resolver(21,  1, 0.000001, 0.1,     false, false );
    %Tanalitic=Preparar_Resolver(101,  1, 0.001, 0.2,     true, true );
     
    disp('Solucion pseudo-analitica calculada...')
    disp('Calculando el error contra esta:')
    
    di = 1/50000;
    df = 1/10000;
    dt = di  :  (df-di)/(100.0-1)  :  df;
    e = zeros( length(dt), 1 );
    for i= 1:length(e)
        T = Preparar_Resolver( 21, 1, dt(i), 1, false, false );
        
        e(i) = ErrorNorma2(Tanalitic, T);
        if( mod(i,(length(e)/10.0)) == 0 )
            disp( [num2str((i/length(e))* 100) , '%'] );
        end
    end
    
    disp('Fin del calculo')
    
    figure(4);
    plot(dt, e);
    title('dt vs Error');
    xlabel 'dt';
    ylabel 'e';
    
    figure(5);
    plot(log(dt), log(e),'y');
    title 'Desarrollo del error - Escala Logaritmica'
    xlabel 'dt'
    ylabel 'e'
    
end

%La ecuacion que resulve es:
% -dT/dt + d^2T/dx^2 - T + Q = 0
%Con dT/dt negativo es decir que tiende a estacionarse
%
%Se usa un esquema explicito para avanzar en el tiempo
function [ Tn ] = Preparar_Resolver( N, tf, dt, F, feedback, useFourier )
   
    X = 0    :1.0/(N-1):    1;%produce vector de N elementos
    dx = 1.0/(N-1);         
 
    if( feedback )%Graficas de curva a la que tendera la ecuacion
        %                                                   N, dx, k, c, Tamb, Q
        Ttiende = Heat_Convection_Diffusion_Reaction_Source(N, dx, 1, 1,    0, 1);
        
        figure(1);
        plot(X,Ttiende,'y');
        title('T(x,t) en t = +inf');
        xlabel 'x';
        ylabel 'T';
        legend('A lo que tiende T');
    end
    
    %Analisis de estabilidad con numero de Fourier
    %El maximo fourier esta en el orden de F = 0.5, para este problema
    %dt = (F*dx^2)/k;   Pero k = 1 aca ENTONCES
    if( useFourier )
        dt = F*dx^2;
    end
        
    %Calcular pasos con Esquema Explicito
    %Tn   = Ttiende;%A ESTO VA A TENDER, PERO PARA QUE SE MUEVA EMPIEZO EN ALGO DISTINTO
    T0 = zeros(N,1); T0(1) = 1;
    
    CC1 = struct('Tipo', 0 , 'Valor', 1);%Tipo 0 es Dirichlet
    CC2 = struct('Tipo', 1 , 'Valor', 0);%Tipo 1 es Newmann
    
    %                                   X,T0, tf, dx, dt, k, c, T_amb, CC1, CC2, caclcHistory
    [Tn, TH] = Time_Advance_Explicit_1D(X,T0, tf, dx, dt, 1, 1, 0    , CC1, CC2, feedback);
      
    if( feedback )%Graficas de resultado final
        %Dibujar superficie de historial
        figure(2);
        t = 0:dt:tf;
        surf(X',t',TH);
        title(['T(x,t) en t dentro [0,',num2str(tf),']']);
        xlabel 'x';
        ylabel 't';
        zlabel 'T';

        %Dibujar curva final e inicial superpuestas
        figure(3);
        plot(X,TH(length(TH(:,1)),:),'b');
        hold on;
        plot(X,TH(1,:),'r');
        title(['Curvas en T(x,0) y T(x,',num2str(tf),')']);
        xlabel 'x';
        ylabel 'T';
        legend('Solucion Tf','Solucion T0');
    end

end

 %La ecuacion que resulve es:
 %          k * d^2T/dx^2 - c * (T - T_amb) + Q = 0
 %
 %                      - N    Cantidad de puntos a calcular
 %                      - dx   el intervalo en el espacio
 %                      - k    cte de conductividad
 %                      - c    cte de perdida de calor al medio ambiente
 %                      - Tamb cte de temperatura ambiente 
%                       - Q    vector de fuente de calor
 function T = Heat_Convection_Diffusion_Reaction_Source(N,dx, k, c, Tamb, Q)

    %Armado de Matriz
    r = k/dx^2;%my aux
    A = eye(N);
    %Condicion de Dirichlet
    A(1,1) = 1;
    for i= 2:N-1
        A(i,i-1) = 1.0*r;   %stencil para i-1
        A(i,i)   = (-2.0)*r-c*1;%stencil para i=j
        A(i,i+1) = 1.0*r;   %stencil para i+1
    end
    %Condicion de Von Neumann
    A(N,N-1) = -1.0/dx;   %Derivada 
    A(N,N)   = 1.0/dx;      %hacia atras 
   
    sparse(A);
    
    b = zeros(N,1);
    b(1,1) = 1;% CC:    T(x=0) = 1   
    for i=2:ceil(N/2)
       b(i) = -Q-c*Tamb;
    end
    b(N,1) = 0;% CC:    q      = 0

    T = A\b;

 end
