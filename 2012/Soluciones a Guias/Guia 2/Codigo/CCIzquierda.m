% Ensambla la contribucion nodal de (i,1) a la matriz global y al vector 
% independiente. Solo para nodos de CC Izquierda!!!
function CCIzquierda(i)

    global K;
    global Nx;
    global Ny;
    global b;
    
    %Parametros de Ecuacion de Calor
    global Q;%Fuente de Calor
    global k;%Difusividad
    global c;%Perdida de calor al ambiente
    global Tamb;%Temperatura ambiente
    global Tipocci;%tipo de cc
    global qcci;%conductividad del borde izquierdo
    global Tcci;%Temperatura impuesta en el borde derecho
    
    global Xrang;
    global Yrang;
    
    % Matriz K se accede ( indice de ecuacion a modificar, indice de termino en la ecuacion )
    % Un nodo (i,j) en el mundo se corresponde con un (solo) indice Tij
    
    Tij   = WorldToMat(i  ,1  );
    
    if( 
    
end