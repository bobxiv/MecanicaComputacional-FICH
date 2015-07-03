% Ensambla la contribucion nodal de (i,j) al vector independiente
% solo para elementos (i,j) internos!!! Para CC usar las funciones de CC
function ContribucionIndependiente(i,j)

    global b;
    
    %Parametros de Ecuacion de Calor
    global Q;%Fuente de Calor
    global k;%Difusividad
    global c;%Perdida de calor al ambiente
    global Tamb;%Temperatura ambiente
    
    global Dimension;%Dimension del problema
    
    global vecX;
    global vecY;
    
    global EqEspaciado;%si los nodos estan equiespaciados
    
    % Matriz K se accede ( indice de ecuacion a modificar, indice de termino en la ecuacion )
    % Un nodo (i,j) en el mundo se corresponde con un (solo) indice Tij
    
    Tij   = WorldToMat(i  ,j);
    
    dx_m = vecX(j  ) - vecX(j-1);%atras
    dx_p = vecX(j+1) - vecX(j  );%adelante
    
    if( strcmp(Dimension, '2D') )
        dy_m = vecY(i  ) - vecY(i-1);%atras
        dy_p = vecY(i+1) - vecY(i  );%adelante
        
        syms x y;
        b(Tij) = -subs(Q, [x y], [vecX(j), vecY(i)] )*(dx_m*dx_p)*(dy_m*dy_p)  +  subs(c)*subs(Tamb)*(dx_m*dx_p)*(dy_m*dy_p);
    else
        syms x;
        b(Tij) = -subs(Q, x, vecX(j) )  +  subs(c)*subs(Tamb);
    end
    
end

