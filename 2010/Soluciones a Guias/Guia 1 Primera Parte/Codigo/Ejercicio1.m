function [  ] = Ejercicio1( )


    %Preparar_Resolver(100, true); %Con esto se puede graficar los
                                  %resultados
    Preparar_Resolver(110, true);

    Npar= 10:10:1000;
    Nimpar = 11:10:1000;
    ErrorsPar = zeros(length(Npar),1);    %NOTAR puntos num pares pero intervalos impares
    ErrorsImpar = zeros(length(Nimpar),1);%NOTAR puntos num impares pero intervalos pares
    for i = 1:length(Npar)
       ErrorsPar(i) = Preparar_Resolver(Npar(i), false);
    end
    for i = 1:length(Nimpar)
        ErrorsImpar(i) = Preparar_Resolver(Nimpar(i), false);
    end
    
    figure(1)
    plot(Npar, ErrorsPar,'g');
    hold on;
    Cuadratica = Npar.^(-2);%-0.000001*N.^2;
    plot(Npar, Cuadratica,'y');
    plot(Nimpar, ErrorsImpar,'r');
    title 'Desarrollo del error'
    xlabel 'Puntos'
    ylabel 'Errores'
    %legend('Error calculado', 'Puntos^2');
    legend('Intervalos impares', 'Intervalos pares');
    %Notar que son intervalos pares o impares NO PUNTOS PARES O IMPARES

    %figure(2);
    %plot(log(N), log(Errors),'r');
    %title 'Desarrollo del error - Escala Logaritmica'
    %xlabel 'Puntos'
    %ylabel 'Errores'

end

function [ e ] = Preparar_Resolver( N, feedback )

    close all;
    
    X = 0    :1.0/(N-1):    1;%produce vector de N
    dx = 1.0/(N-1);           

    T1 = -X(1:ceil(N/2)).^2./2 + (-5./8).*X(1:ceil(N/2)) +1;%Solucion Analitica, tramo 1
    T2 = 9/8.*(-X(ceil(N/2)+1:N)+1);                        %Solucion Analitica, tramo 2
    
    Tanalitic = [T1' ; T2'];
    
    if( feedback )
        plot(X,Tanalitic,'b');
    end
    
    %Creando la matriz del problema
    A = eye(N);
    A = A.*(-2);
    A(1,1) = 1;
    for i= 2:N-1
        A(i,i+1) = 1;
        A(i,i-1) = 1;
    end
    A(N,N) = 1;

    sparse(A);
    
    %Creacion del vector independiente
    b = zeros(N,1);        %Tramo 2: d^2Fi/dx^2     = 0
    b(1) = 1;
    for i=2:ceil(N/2)
       b(i) = -1*(dx^2);%Tramo 1: d^2Fi/dx^2 + 1 = 0
    end
    b(N) = 0;
    
    T = A\b;%Solucionando

    if( feedback )
        hold on;
        plot(X,T,'r');
        xlabel 'x';
        ylabel 'Y';
        legend('Analitica','MDF');
    end
    
    e = CalcularError(Tanalitic, T);

end
