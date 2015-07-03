% Ensambla la contribucion nodal de (i,j) al vector independiente
% solo para elementos (i,j) internos!!! Para CC usar las funciones de CC
function ContribucionIndependiente(i,j)

    global b;
    
    %Parametros de Ecuacion de Calor
    global Q;%Fuente de Calor
    
    global vecR;
    global vecTheta;
    
    % Matriz K se accede ( indice de ecuacion a modificar, indice de termino en la ecuacion )
    % Un nodo (i,j) en el mundo se corresponde con un (solo) indice Tij
    
    Tij   = WorldToMat(i  ,j);

    syms x y;
    b(Tij) = -subs(Q, [x y], [vecR(j)*cos(vecTheta(i)), vecR(j)*sin(vecTheta(i))] );
    
end

