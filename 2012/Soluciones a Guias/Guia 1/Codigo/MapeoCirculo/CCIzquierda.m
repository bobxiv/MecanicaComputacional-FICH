% Ensambla la contribucion nodal de (i,1) a la matriz global y al vector 
% independiente. Solo para nodos de CC Izquierda!!!
function CCIzquierda(i)

    global K;
    global Nx;
    global Ny;
    global b;
    
    %Parametros de Ecuacion de Calor
    global Q;%Fuente de Calor

    global Tipocci;%tipo de cc
    global qcci;%conductividad del borde izquierdo
    global Tcci;%Temperatura impuesta en el borde derecho
    
    global vecR;
    global vecTheta;
    
    % Matriz K se accede ( indice de ecuacion a modificar, indice de termino en la ecuacion )
    % Un nodo (i,j) en el mundo se corresponde con un (solo) indice Tij
    
    Tij   = WorldToMat(i  ,1  );
    Tijp1 = WorldToMat(i  ,1+1  );
    Tijp2 = WorldToMat(i  ,1+2  );
    
    if( strcmp(Tipocci, 'Neumann') )
        dr = vecR(2)-vecR(1);
        
        K(Tij,Tij  ) = - 3/(2*dr);
        K(Tij,Tijp1) = 4/(2*dr);
        K(Tij,Tijp2) = - 1/(2*dr);
        b(Tij      ) = -subs(qcci);
    end
    
    if( strcmp(Tipocci, 'Dirichlet') )
        K(Tij,Tij  ) = 1;
        b(Tij      ) = subs(Tcci);
    end
    
end