% Ensambla la contribucion nodal de (i,j) a la matriz global
% solo para elementos (i,j) internos!!! Para CC usar las funciones de CC
function ContribucionNodal(i,j)

    global K;
    
    %Parametros de Ecuacion de Calor
    global Q;%Fuente de Calor
    
    global vecR;
    global vecTheta;
    
    % Matriz K se accede ( indice de ecuacion a modificar, indice de termino en la ecuacion )
    % Un nodo (i,j) en el mundo se corresponde con un (solo) indice Tij
    
    Tij   = WorldToMat(i  ,j  );
    Tijp1 = WorldToMat(i  ,j+1);
    Tijm1 = WorldToMat(i  ,j-1);
    Tip1j = WorldToMat(i+1,j  );
    Tim1j = WorldToMat(i-1,j  );
    
    dr     = vecR(j)-vecR(j-1);
    dtheta = vecTheta(i)-vecTheta(i-1);
    
    r = vecR(j);

    K(Tij,Tij  ) = - 2/dr^2 - 2/(dtheta^2*r^2);
    K(Tij,Tim1j) = 1/dr^2 - 1/(r*2*dr);
    K(Tij,Tip1j) = 1/dr^2 + 1/(r*2*dr);
    K(Tij,Tijm1) = 1/(dtheta^2*r^2);
    K(Tij,Tijp1) = 1/(dtheta^2*r^2);
    
end