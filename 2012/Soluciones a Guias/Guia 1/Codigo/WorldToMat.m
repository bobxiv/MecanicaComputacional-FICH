% Obtiene el indice de la matriz para el nodo (i,j)
% En 2D siempre ordena por filas
function [Tij] = WorldToMat(i,j)
    global Nx;
    %global Ny;
    
    Tij = Nx*(i-1)+j;%contando por fila
end