% Obtiene el indice de la matriz para el nodo (i,j)
% En 2D siempre ordena por filas
function [Tij] = WorldToMat(l,m)
    global M;
    %global Ny;
    
    Tij = M*(l-1)+m;%contando por fila
end