% Ensambla la contribucion nodal de (i,Nx) a la matriz global y al vector 
% independiente. Solo para nodos de CC Derecha!!!
function CCDerecha(i)

    global K;
    global Nx;
    global Ny;
    global b;
    
    %Parametros de Ecuacion de Calor
    global Q;%Fuente de Calor

    global Tipoccd;%tipo de cc
    global qccd;%conductividad del borde derecho
    global Tccd;%Temperatura impuesta en el borde derecho
    
    global vecR;
    global vecTheta;
    
    % Matriz K se accede ( indice de ecuacion a modificar, indice de termino en la ecuacion )
    % Un nodo (i,j) en el mundo se corresponde con un (solo) indice Tij
    
    Tij   = WorldToMat(i  ,Nx  );
    Tijm1 = WorldToMat(i  ,Nx-1);
    Tijm2 = WorldToMat(i  ,Nx-2);
    
    if( strcmp(Tipoccd, 'Neumann') )
        dr = vecR(Nx)-vecR(Nx-1);

        K(Tij,Tij  ) = 3/(2*dr);
        K(Tij,Tijm1) = - 4/(2*dr);
        K(Tij,Tijm2) = 1/(2*dr);
        b(Tij      ) = -subs(qccd);
    end
    
    if( strcmp(Tipoccd, 'Dirichlet') )
        K(Tij,Tij  ) = 1;
        b(Tij      ) = subs(Tccd);
    end
    
end