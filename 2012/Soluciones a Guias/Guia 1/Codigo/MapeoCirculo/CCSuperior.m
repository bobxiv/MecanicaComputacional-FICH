% Ensambla la contribucion nodal de (Ny,j) a la matriz global y al vector 
% independiente. Solo para nodos de CC Superior!!!
function CCSuperior(j)

    global K;
    global Nx;
    global Ny;
    global b;
    
    %Parametros de Ecuacion de Calor
    global Q;%Fuente de Calor

    global Tipoccs;%tipo de cc
    global qccs;%conductividad del borde superior
    global Tccs;%Temperatura impuesta en el borde derecho
    
    global vecR;
    global vecTheta;
    
    % Matriz K se accede ( indice de ecuacion a modificar, indice de termino en la ecuacion )
    % Un nodo (i,j) en el mundo se corresponde con un (solo) indice Tij
    
    Tij   = WorldToMat(Ny  ,j  );
    Tim1j = WorldToMat(Ny-1,j  );
    Tim2j = WorldToMat(Ny-2,j  );

    if( strcmp(Tipoccs, 'Neumann') )
        dtheta = vecTheta(Ny  ) - vecTheta(Ny-1);%adelante
    
        K(Tij,Tij  ) = 3/(2*dtheta);
        K(Tij,Tim1j) = -4/(2*dtheta);
        K(Tij,Tim2j) = 1/(2*dtheta);
        b(Tij      ) = -subs(qccs);
    end
    
    if( strcmp(Tipoccs, 'Dirichlet') )
        K(Tij,Tij  ) = 1;
        b(Tij      ) = subs(Tccs);
    end
    
end